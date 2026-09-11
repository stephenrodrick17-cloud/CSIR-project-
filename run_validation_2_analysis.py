#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pan-Fibrotic Gene Signature Project - Validation 2 Pipeline
-----------------------------------------------------------
Runs end-to-end Validation 2 pipeline across 4 organs (kidney, liver, lung, skin):

Tasks:
1. Differential Expression (Welch's t-test, BH adjustment, adj.P.Val < 0.05, |logFC| > 0.585)
2. Narrow to Validation 2 Signature (Intersect DEGs with organ's Validation 1 signature)
3. Spearman Correlation with Severity Numeric (where available)
4. Functional Enrichment Analysis (GO Biological Process + KEGG via gseapy)
5. Functional mRNA Interaction Network (STRING REST API)
6. Cross-organ Summary & All-4-Organ Overlap Analysis
7. Signature Heatmaps (Z-scored expression)
"""

import os
import sys
import time
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

try:
    import gseapy as gp
except ImportError:
    gp = None

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(BASE_DIR, "organ_validation2_data")
VAL1_DIR = os.path.join(BASE_DIR, "validation_1_all_organs")
OUT_DIR = os.path.join(BASE_DIR, "validation_2")

ORGANS = ["kidney", "liver", "lung", "skin"]

def run_task1_de(expr_df, sample_groups):
    """
    TASK 1: Differential Expression Analysis (Welch's t-test + BH-adjusted p-values)
    """
    ctrl_samples = sample_groups[sample_groups["group"] == "control"]["sample_id"].tolist()
    fib_samples = sample_groups[sample_groups["group"] == "fibrotic"]["sample_id"].tolist()

    ctrl_expr = expr_df[ctrl_samples].values
    fib_expr = expr_df[fib_samples].values

    ctrl_mean = np.mean(ctrl_expr, axis=1)
    fib_mean = np.mean(fib_expr, axis=1)
    logfc = fib_mean - ctrl_mean

    t_stats, p_vals = stats.ttest_ind(fib_expr, ctrl_expr, axis=1, equal_var=False, nan_policy='omit')
    p_vals = np.nan_to_num(p_vals, nan=1.0)
    
    _, adj_p_vals, _, _ = multipletests(p_vals, method='fdr_bh')

    top_table = pd.DataFrame({
        "gene": expr_df["gene"],
        "control_mean": ctrl_mean,
        "fibrotic_mean": fib_mean,
        "logFC": logfc,
        "P.Value": p_vals,
        "adj.P.Val": adj_p_vals
    })

    sig_degs = top_table[(top_table["adj.P.Val"] < 0.05) & (top_table["logFC"].abs() > 0.585)]
    return top_table, sig_degs

def run_task2_narrow(sig_degs, val1_path):
    """
    TASK 2: Intersect Task 1 DEGs with Validation 1 signature for this organ
    """
    if os.path.exists(val1_path):
        val1_df = pd.read_csv(val1_path)
        val1_genes = set(val1_df["gene"].astype(str))
    else:
        val1_genes = set()

    sig_deg_genes = set(sig_degs["gene"].astype(str))
    val2_genes = sorted(list(val1_genes.intersection(sig_deg_genes)))
    return val2_genes, len(val1_genes)

def run_task3_correlation(expr_df, sample_groups, val2_genes, organ_out_dir, organ):
    """
    TASK 3: Spearman Correlation Analysis with Severity Numeric
    """
    corr_dir = os.path.join(organ_out_dir, "correlation")
    os.makedirs(corr_dir, exist_ok=True)

    if "severity_numeric" not in sample_groups.columns:
        print(f"[{organ.upper()}] No severity_numeric column in sample_groups.csv. Skipping Spearman correlation.")
        return [], None

    valid_samples = sample_groups.dropna(subset=["severity_numeric"])
    if len(valid_samples) == 0:
        return [], None

    samples = valid_samples["sample_id"].tolist()
    severities = valid_samples["severity_numeric"].values

    expr_sub = expr_df[expr_df["gene"].isin(val2_genes)].set_index("gene")[samples]

    corr_list = []
    sig_corr_genes = []
    for g in expr_sub.index:
        g_vals = expr_sub.loc[g].values
        rho, pval = stats.spearmanr(g_vals, severities)
        corr_list.append({"gene": g, "spearman_rho": rho, "p_value": pval})
        if pval < 0.05:
            sig_corr_genes.append(g)

    corr_df = pd.DataFrame(corr_list)
    corr_df.to_csv(os.path.join(corr_dir, f"{organ}_severity_correlation.csv"), index=False)

    # Plot boxplot / scatter plot per gene
    if len(corr_df) > 0:
        top_genes = corr_df.sort_values("p_value").head(6)["gene"].tolist()
        fig, axes = plt.subplots(nrows=(len(top_genes)+1)//2, ncols=min(2, len(top_genes)), figsize=(10, 4 * ((len(top_genes)+1)//2)))
        axes = np.array(axes).flatten() if len(top_genes) > 1 else [axes]
        
        for idx, g in enumerate(top_genes):
            g_vals = expr_sub.loc[g].values
            rho_val = corr_df.loc[corr_df["gene"] == g, "spearman_rho"].values[0]
            p_val = corr_df.loc[corr_df["gene"] == g, "p_value"].values[0]
            
            sns.regplot(x=severities, y=g_vals, ax=axes[idx], color="#2b5c8f", scatter_kws={'s': 50})
            axes[idx].set_title(f"{g} (rho={rho_val:.2f}, p={p_val:.3e})", fontsize=11)
            axes[idx].set_xlabel("Severity Numeric")
            axes[idx].set_ylabel("Expression Level")
            
        plt.tight_layout()
        plt.savefig(os.path.join(corr_dir, f"{organ}_severity_correlation.png"), dpi=300)
        plt.close()

    return sig_corr_genes, corr_df

def run_task4_enrichment(gene_list, enrichment_out_dir, prefix):
    """
    TASK 4: Functional Enrichment Analysis (GO Biological Process + KEGG)
    """
    os.makedirs(enrichment_out_dir, exist_ok=True)
    if len(gene_list) == 0:
        print(f"[{prefix}] No genes to run enrichment.")
        return

    if gp is None:
        print(f"[{prefix}] gseapy module not available. Skipping enrichment.")
        return

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
            organism='human',
            outdir=enrichment_out_dir,
            cutoff=0.05
        )
        res = enr.results
        res.to_csv(os.path.join(enrichment_out_dir, f"{prefix}_enrichment_results.csv"), index=False)

        if len(res) > 0:
            top_terms = res.sort_values("Adjusted P-value").head(10).copy()
            top_terms["-log10(P-value)"] = -np.log10(top_terms["P-value"] + 1e-15)
            
            plt.figure(figsize=(9, 5))
            sns.barplot(x="-log10(P-value)", y="Term", data=top_terms, palette="crest")
            plt.title(f"Top Enriched Functional Terms - {prefix.capitalize()}", fontsize=12)
            plt.xlabel("-log10(P-value)")
            plt.ylabel("Term")
            plt.tight_layout()
            plt.savefig(os.path.join(enrichment_out_dir, f"{prefix}_enrichment_barplot.png"), dpi=300)
            plt.close()
    except Exception as e:
        print(f"[{prefix}] Enrichment error: {e}")

def run_task5_network(gene_list, network_out_dir, prefix):
    """
    TASK 5: STRING REST API mRNA Interaction Network Query
    """
    os.makedirs(network_out_dir, exist_ok=True)
    if len(gene_list) == 0:
        print(f"[{prefix}] No genes for network query.")
        return

    string_url = "https://string-db.org/api"
    ids_param = "%0d".join(gene_list)

    # a & b: Fetch network PNG image
    try:
        img_url = f"{string_url}/png/network?identifiers={ids_param}&species=9606&required_score=400"
        img_resp = requests.get(img_url)
        if img_resp.status_code == 200:
            with open(os.path.join(network_out_dir, f"{prefix}_network.png"), "wb") as f:
                f.write(img_resp.content)
            print(f"[{prefix}] Saved STRING network image.")
        else:
            print(f"[{prefix}] STRING image query failed with status code {img_resp.status_code}")
    except Exception as e:
        print(f"[{prefix}] STRING image error: {e}")
        
    time.sleep(1)

    # c: Fetch raw edge list TSV
    try:
        tsv_url = f"{string_url}/tsv/network?identifiers={ids_param}&species=9606&required_score=400"
        tsv_resp = requests.get(tsv_url)
        if tsv_resp.status_code == 200:
            with open(os.path.join(network_out_dir, f"{prefix}_network_edges.tsv"), "w", encoding="utf-8") as f:
                f.write(tsv_resp.text)
            print(f"[{prefix}] Saved STRING network edge list TSV.")
        else:
            print(f"[{prefix}] STRING edge list query failed with status code {tsv_resp.status_code}")
    except Exception as e:
        print(f"[{prefix}] STRING TSV error: {e}")

    time.sleep(1)

def run_task7_heatmap(expr_df, sample_groups, val2_genes, organ_out_dir, organ):
    """
    TASK 7 / TASK 6 Heatmap: Z-scored expression heatmap ordered by group/severity
    """
    if len(val2_genes) == 0:
        return

    sg = sample_groups.copy()
    if "severity_numeric" in sg.columns:
        sg = sg.sort_values(by=["group", "severity_numeric"])
    else:
        sg = sg.sort_values(by=["group"])

    sorted_samples = sg["sample_id"].tolist()
    expr_sub = expr_df[expr_df["gene"].isin(val2_genes)].set_index("gene")[sorted_samples]

    z_scores = expr_sub.apply(lambda x: (x - x.mean()) / (x.std() + 1e-8), axis=1)

    plt.figure(figsize=(10, max(4, len(val2_genes) * 0.4)))
    sns.heatmap(z_scores, cmap="coolwarm", center=0, cbar_kws={'label': 'Z-score'})
    plt.title(f"Validation 2 Signature Heatmap - {organ.capitalize()}", fontsize=12)
    plt.xlabel("Samples (Control -> Fibrotic / Severity)")
    plt.ylabel("Gene Symbol")
    plt.tight_layout()
    plt.savefig(os.path.join(organ_out_dir, f"{organ}_heatmap.png"), dpi=300)
    plt.close()

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    summary_data = []
    organ_val2_dict = {}

    print("=================================================================")
    print("STARTING VALIDATION 2 PIPELINE FOR ALL 4 ORGANS")
    print("=================================================================\n")

    for organ in ORGANS:
        print(f">>> Processing Organ: {organ.upper()} <<<")
        organ_out_dir = os.path.join(OUT_DIR, organ)
        os.makedirs(organ_out_dir, exist_ok=True)

        expr_file = os.path.join(INPUT_DIR, organ, "expr_matrix.csv")
        group_file = os.path.join(INPUT_DIR, organ, "sample_groups.csv")
        val1_file = os.path.join(VAL1_DIR, organ, f"{organ}_validated_signature.csv")

        expr_df = pd.read_csv(expr_file)
        sample_groups = pd.read_csv(group_file)

        n_input_genes = len(expr_df)

        # TASK 1: Differential Expression
        top_table, sig_degs = run_task1_de(expr_df, sample_groups)
        top_table.to_csv(os.path.join(organ_out_dir, f"{organ}_deg_top_table.csv"), index=False)
        sig_degs.to_csv(os.path.join(organ_out_dir, f"{organ}_sig_degs.csv"), index=False)

        # TASK 2: Intersect with Validation 1 Signature
        val2_genes, n_val1_genes = run_task2_narrow(sig_degs, val1_file)
        val2_df = pd.DataFrame({"gene": val2_genes})
        val2_df.to_csv(os.path.join(organ_out_dir, f"{organ}_validation2_signature.csv"), index=False)
        organ_val2_dict[organ] = val2_genes

        # TASK 3: Spearman Correlation
        sig_corr_genes, corr_df = run_task3_correlation(expr_df, sample_groups, val2_genes, organ_out_dir, organ)

        # TASK 4: Functional Enrichment
        enrich_dir = os.path.join(OUT_DIR, "enrichment", organ)
        run_task4_enrichment(val2_genes, enrich_dir, organ)

        # TASK 5: STRING mRNA Interaction Network
        net_dir = os.path.join(OUT_DIR, "network", organ)
        run_task5_network(val2_genes, net_dir, organ)

        # TASK 7 / TASK 6.1: Heatmap
        run_task7_heatmap(expr_df, sample_groups, val2_genes, organ_out_dir, organ)

        summary_data.append({
            "organ": organ,
            "n_input_genes_task1": n_input_genes,
            "n_Validation1_genes": n_val1_genes,
            "n_Validation2_genes": len(val2_genes),
            "sig_correlation_genes": ", ".join(sig_corr_genes) if sig_corr_genes else "None"
        })

        # TASK 6.2: Print per-organ summary
        print(f"\n--- {organ.upper()} SUMMARY ---")
        print(f"Total genes evaluated in Task 1: {n_input_genes}")
        print(f"Validation 1 signature genes: {n_val1_genes}")
        print(f"Validation 2 signature genes (Task 2): {len(val2_genes)} ({', '.join(val2_genes)})")
        print(f"Significant severity correlation genes (Task 3): {', '.join(sig_corr_genes) if sig_corr_genes else 'None'}\n")

    # CROSS-ORGAN SUMMARY (TASK 6)
    print("=================================================================")
    print("CROSS-ORGAN COMBINED ANALYSIS & OVERLAP")
    print("=================================================================")

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(OUT_DIR, "validation2_summary.csv"), index=False)

    # Combined Funnel Bar Chart
    plt.figure(figsize=(9, 5))
    x = np.arange(len(ORGANS))
    width = 0.35
    plt.bar(x - width/2, summary_df["n_Validation1_genes"], width, label="Validation 1", color="#3470a3")
    plt.bar(x + width/2, summary_df["n_Validation2_genes"], width, label="Validation 2", color="#3ea066")
    plt.xticks(x, [o.capitalize() for o in ORGANS])
    plt.ylabel("Gene Count")
    plt.title("Pan-Fibrotic Gene Signature Funnel (Validation 1 vs Validation 2)", fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "validation2_funnel_chart.png"), dpi=300)
    plt.close()

    # Compute All-4-Organ Overlap
    all4_genes = set.intersection(*[set(organ_val2_dict[o]) for o in ORGANS])
    all4_list = sorted(list(all4_genes))

    all4_df = pd.DataFrame({"gene": all4_list})
    all4_df.to_csv(os.path.join(OUT_DIR, "validation2_all4organs_overlap.csv"), index=False)

    print(f"\nAll 4 Organs Validation 2 Overlapping Genes ({len(all4_list)}): {', '.join(all4_list)}")

    if len(all4_list) > 0:
        run_task4_enrichment(all4_list, os.path.join(OUT_DIR, "enrichment", "all4organs"), "all4organs")
        run_task5_network(all4_list, os.path.join(OUT_DIR, "network", "all4organs"), "all4organs")

    print("\n=================================================================")
    print("VALIDATION 2 PIPELINE SUCCESSFULLY COMPLETED FOR ALL ORGANS")
    print("=================================================================")

if __name__ == "__main__":
    main()
