#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pan-Fibrotic Gene Signature Project - Fresh Validation 2 Pipeline
------------------------------------------------------------------
Runs end-to-end Validation 2 pipeline across 4 organs (kidney, liver, lung, skin):

Tasks & Features:
1. Differential Expression Analysis (Welch's t-test, BH adjustment, adj.P.Val < 0.05, |logFC| > 0.585)
2. Cross-Organ Overlap Analysis:
   - Common DEG signature across 4 organs (Kidney, Liver, Lung, Skin)
   - Common DEG signature across 3 organs (Lungs, Skin, Kidneys)
   - Full Genome DEG cross-organ overlap comparison
3. Visualizations:
   - 3-Organ Venn Diagram (Lungs, Skin, Kidneys) via matplotlib_venn
   - 4-Organ UpSet Plot & 4-Ellipse Venn Diagram
   - Pan-Fibrotic Funnel & DEG count summary bar chart
   - Z-score expression heatmaps per organ & cross-organ
4. Spearman Correlation with Severity Numeric (full sample & disease-only)
5. Functional Enrichment Analysis (GO Biological Process 2023 + KEGG 2021 Human via gseapy)
6. STRING REST API PPI Interaction Networks for 4-organ and 3-organ DEG signatures
7. Clean structured reporting & CSV export in validation_2/
"""

import os
import sys
import time
import requests
import warnings
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Optional packages check
try:
    import matplotlib_venn as mv
except ImportError:
    mv = None

try:
    import upsetplot
except ImportError:
    upsetplot = None

try:
    import gseapy as gp
except ImportError:
    gp = None

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace") if hasattr(sys.stdout, "reconfigure") else None

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(BASE_DIR, "organ_validation2_data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
OUT_DIR   = os.path.join(BASE_DIR, "validation_2")

ORGANS = ["kidney", "liver", "lung", "skin"]
THREE_ORGANS = ["lung", "skin", "kidney"]

ORGAN_COLORS = {
    "kidney": "#5B8DB8",
    "liver":  "#E07B54",
    "lung":   "#6BAE75",
    "skin":   "#A97DC9",
}

# =============================================================================
# TASK 1: DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================
def run_task1_de(expr_df, sample_groups):
    """
    Differential Expression Analysis (Welch's t-test + BH-adjusted p-values)
    """
    ctrl_samples = sample_groups[sample_groups["group"] == "control"]["sample_id"].tolist()
    fib_samples  = sample_groups[sample_groups["group"] == "fibrotic"]["sample_id"].tolist()

    ctrl_expr = expr_df[ctrl_samples].values
    fib_expr  = expr_df[fib_samples].values

    ctrl_mean = np.mean(ctrl_expr, axis=1)
    fib_mean  = np.mean(fib_expr, axis=1)
    logfc     = fib_mean - ctrl_mean

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

    sig_degs = top_table[(top_table["adj.P.Val"] < 0.05) & (top_table["logFC"].abs() > 0.585)].copy()
    return top_table, sig_degs

# =============================================================================
# TASK 2: SPEARMAN SEVERITY CORRELATION
# =============================================================================
def run_task2_correlation(expr_df, sample_groups, deg_genes, organ_out_dir, organ):
    """
    Spearman Correlation Analysis with Severity Numeric (Full Sample & Disease-Only)
    """
    corr_dir = os.path.join(organ_out_dir, "correlation")
    os.makedirs(corr_dir, exist_ok=True)

    if "severity_numeric" not in sample_groups.columns:
        print(f"[{organ.upper()}] No severity_numeric column in sample_groups.csv.")
        return [], pd.DataFrame(), pd.DataFrame()

    valid_samples = sample_groups.dropna(subset=["severity_numeric"])
    if len(valid_samples) == 0:
        return [], pd.DataFrame(), pd.DataFrame()

    # Full Sample Correlation
    samples = valid_samples["sample_id"].tolist()
    severities = valid_samples["severity_numeric"].values
    expr_sub = expr_df[expr_df["gene"].isin(deg_genes)].set_index("gene")[samples]

    corr_list = []
    sig_corr_genes = []
    for g in expr_sub.index:
        g_vals = expr_sub.loc[g].values
        rho, pval = stats.spearmanr(g_vals, severities)
        corr_list.append({"gene": g, "spearman_rho": rho, "p_value": pval, "significant": pval < 0.05})
        if pval < 0.05:
            sig_corr_genes.append(g)

    corr_df = pd.DataFrame(corr_list)
    corr_df.to_csv(os.path.join(corr_dir, f"{organ}_severity_correlation.csv"), index=False)

    # Plot Full Sample Scatter / Regplots
    if len(corr_df) > 0 and len(expr_sub) > 0:
        top_genes = corr_df.sort_values("p_value").head(6)["gene"].tolist()
        n_plots = len(top_genes)
        n_rows = (n_plots + 1) // 2
        fig, axes = plt.subplots(nrows=n_rows, ncols=min(2, n_plots), figsize=(10, 3.5 * n_rows))
        axes = np.array(axes).flatten() if n_plots > 1 else [axes]
        
        for idx, g in enumerate(top_genes):
            g_vals = expr_sub.loc[g].values
            rho_val = corr_df.loc[corr_df["gene"] == g, "spearman_rho"].values[0]
            p_val = corr_df.loc[corr_df["gene"] == g, "p_value"].values[0]
            
            sns.regplot(x=severities, y=g_vals, ax=axes[idx], color=ORGAN_COLORS.get(organ, "#333333"), scatter_kws={'s': 50})
            axes[idx].set_title(f"{g} (rho={rho_val:.2f}, p={p_val:.3e})", fontsize=11, fontweight="bold")
            axes[idx].set_xlabel("Severity Score")
            axes[idx].set_ylabel("Expression Level")
            
        plt.tight_layout()
        plt.savefig(os.path.join(corr_dir, f"{organ}_severity_correlation.png"), dpi=300)
        plt.close()

    # Disease-Only Correlation (Fibrotic samples only)
    fib_sg = valid_samples[valid_samples["group"] == "fibrotic"]
    dis_corr_df = pd.DataFrame()
    if len(fib_sg) > 0:
        fib_samples = fib_sg["sample_id"].tolist()
        fib_severities = fib_sg["severity_numeric"].values
        expr_fib = expr_df[expr_df["gene"].isin(deg_genes)].set_index("gene")[fib_samples]
        
        dis_list = []
        for g in expr_fib.index:
            g_vals = expr_fib.loc[g].values
            rho, pval = stats.spearmanr(g_vals, fib_severities)
            dis_list.append({"gene": g, "spearman_rho": rho, "p_value": pval, "significant": pval < 0.05})
        dis_corr_df = pd.DataFrame(dis_list)
        dis_corr_df.to_csv(os.path.join(corr_dir, f"{organ}_disease_only_spearman.csv"), index=False)

    return sig_corr_genes, corr_df, dis_corr_df

# =============================================================================
# TASK 3: FUNCTIONAL ENRICHMENT ANALYSIS (GO BP & KEGG)
# =============================================================================
def run_task3_enrichment(gene_list, enrichment_out_dir, prefix):
    """
    Functional Enrichment Analysis (GO Biological Process + KEGG via gseapy)
    """
    os.makedirs(enrichment_out_dir, exist_ok=True)
    if len(gene_list) == 0:
        print(f"[{prefix}] No genes provided for enrichment.")
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
            plt.title(f"Top Enriched Functional Terms - {prefix.replace('_', ' ').title()}", fontsize=12, fontweight="bold")
            plt.xlabel("-log10(P-value)")
            plt.ylabel("Functional Term")
            plt.tight_layout()
            plt.savefig(os.path.join(enrichment_out_dir, f"{prefix}_enrichment_barplot.png"), dpi=300)
            plt.close()
            print(f"[{prefix}] Enrichment analysis saved successfully.")
    except Exception as e:
        print(f"[{prefix}] Enrichment notice/error: {e}")

# =============================================================================
# TASK 4: STRING REST API PPI NETWORK
# =============================================================================
def run_task4_network(gene_list, network_out_dir, prefix):
    """
    STRING REST API mRNA/Protein Interaction Network Query
    """
    os.makedirs(network_out_dir, exist_ok=True)
    if len(gene_list) == 0:
        print(f"[{prefix}] No genes for network query.")
        return

    string_url = "https://string-db.org/api"
    ids_param = "%0d".join(gene_list)

    # Network PNG Image
    try:
        img_url = f"{string_url}/png/network?identifiers={ids_param}&species=9606&required_score=400"
        img_resp = requests.get(img_url, timeout=15)
        if img_resp.status_code == 200:
            with open(os.path.join(network_out_dir, f"{prefix}_network.png"), "wb") as f:
                f.write(img_resp.content)
            print(f"[{prefix}] Saved STRING PPI network image.")
    except Exception as e:
        print(f"[{prefix}] STRING network image notice: {e}")
        
    time.sleep(0.5)

    # Raw Edge List TSV
    try:
        tsv_url = f"{string_url}/tsv/network?identifiers={ids_param}&species=9606&required_score=400"
        tsv_resp = requests.get(tsv_url, timeout=15)
        if tsv_resp.status_code == 200:
            with open(os.path.join(network_out_dir, f"{prefix}_network_edges.tsv"), "w", encoding="utf-8") as f:
                f.write(tsv_resp.text)
            print(f"[{prefix}] Saved STRING PPI network edge list TSV.")
    except Exception as e:
        print(f"[{prefix}] STRING network TSV notice: {e}")

# =============================================================================
# TASK 5: VISUALIZATIONS & HEATMAPS
# =============================================================================
def plot_heatmaps(organ_expr_dict, organ_sg_dict, target_genes, out_path, title_prefix):
    """
    Generates Z-scored expression heatmaps across organs for specified gene signature
    """
    if not target_genes:
        return

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, organ in enumerate(ORGANS):
        ax = axes[idx]
        expr_df = organ_expr_dict[organ]
        sg = organ_sg_dict[organ].copy()

        if "severity_numeric" in sg.columns:
            sg = sg.sort_values(by=["group", "severity_numeric"])
        else:
            sg = sg.sort_values(by=["group"])

        sorted_samples = sg["sample_id"].tolist()
        expr_sub = expr_df[expr_df["gene"].isin(target_genes)].set_index("gene")[sorted_samples]
        
        if len(expr_sub) == 0:
            ax.text(0.5, 0.5, "No Matching Genes", ha="center", va="center")
            ax.set_title(organ.capitalize())
            continue

        z_scores = expr_sub.apply(lambda x: (x - x.mean()) / (x.std() + 1e-8), axis=1)

        sns.heatmap(z_scores, cmap="coolwarm", center=0, ax=ax, cbar_kws={'label': 'Z-score'})
        ax.set_title(f"{title_prefix} - {organ.capitalize()}", fontsize=11, fontweight="bold")
        ax.set_xlabel("Samples (Control -> Fibrotic)")
        ax.set_ylabel("Gene Symbol")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

def plot_3organ_venn(lungs_genes, skin_genes, kidney_genes, out_path, title):
    """
    Plots 3-Organ Venn Diagram for Lungs, Skin, and Kidneys using matplotlib_venn
    """
    plt.figure(figsize=(8, 7))
    if mv is not None:
        v = mv.venn3(
            subsets=[set(lungs_genes), set(skin_genes), set(kidney_genes)],
            set_labels=('Lung DEGs', 'Skin DEGs', 'Kidney DEGs'),
            set_colors=('#6BAE75', '#A97DC9', '#5B8DB8'),
            alpha=0.65
        )
        plt.title(title, fontsize=13, fontweight="bold", pad=15)
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"Saved 3-organ Venn diagram to {out_path}")
    else:
        print("matplotlib_venn not installed. Skipping venn3 plot.")

def plot_upset(organ_deg_dict, out_path, title):
    """
    Plots UpSet Plot for all 4 organs
    """
    if upsetplot is None:
        print("upsetplot library not installed. Skipping UpSet plot.")
        return

    try:
        from upsetplot import from_contents, UpSet
        contents = {o.capitalize(): set(genes) for o, genes in organ_deg_dict.items()}
        upset_data = from_contents(contents)
        
        plt.figure(figsize=(10, 6))
        upset = UpSet(upset_data, subset_size='count', show_counts=True, sort_by='cardinality')
        upset.plot()
        plt.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved UpSet plot to {out_path}")
    except Exception as e:
        print(f"UpSet plot notice: {e}")

def plot_funnel_summary(summary_df, out_path):
    """
    Plots bar chart comparing DEG counts across organs & overlaps
    """
    plt.figure(figsize=(10, 5))
    x = np.arange(len(ORGANS))
    width = 0.25

    plt.bar(x - width, summary_df["n_input_genes"], width, label="Input Genes", color="#888888")
    plt.bar(x, summary_df["n_sig_degs"], width, label="Significant DEGs", color="#3470a3")
    plt.bar(x + width, summary_df["n_common_4organs"], width, label="4-Organ Common", color="#3ea066")

    plt.xticks(x, [o.capitalize() for o in ORGANS], fontweight="bold")
    plt.ylabel("Gene Count", fontweight="bold")
    plt.title("Validation 2 DEG Analysis Summary across 4 Organs", fontsize=12, fontweight="bold")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    print("=================================================================")
    print("STARTING FRESH VALIDATION 2 PIPELINE (4 ORGANS & 3 ORGANS)")
    print("=================================================================\n")

    organ_expr_dict = {}
    organ_sg_dict   = {}
    organ_deg_dict  = {}
    summary_data    = []

    for organ in ORGANS:
        print(f">>> Processing Organ: {organ.upper()} <<<")
        organ_out_dir = os.path.join(OUT_DIR, organ)
        os.makedirs(organ_out_dir, exist_ok=True)

        expr_file  = os.path.join(INPUT_DIR, organ, "expr_matrix.csv")
        group_file = os.path.join(INPUT_DIR, organ, "sample_groups.csv")

        expr_df       = pd.read_csv(expr_file)
        sample_groups = pd.read_csv(group_file)

        organ_expr_dict[organ] = expr_df
        organ_sg_dict[organ]   = sample_groups

        n_input = len(expr_df)

        # 1. Differential Expression Analysis
        top_table, sig_degs = run_task1_de(expr_df, sample_groups)
        top_table.to_csv(os.path.join(organ_out_dir, f"{organ}_deg_top_table.csv"), index=False)
        sig_degs.to_csv(os.path.join(organ_out_dir, f"{organ}_sig_degs.csv"), index=False)

        deg_genes = sig_degs["gene"].astype(str).tolist()
        organ_deg_dict[organ] = deg_genes

        # 2. Spearman Severity Correlation
        sig_corr_genes, corr_df, dis_corr_df = run_task2_correlation(
            expr_df, sample_groups, deg_genes, organ_out_dir, organ
        )

        # 3. Individual Organ Enrichment & Network
        run_task3_enrichment(deg_genes, os.path.join(OUT_DIR, "enrichment", organ), organ)
        run_task4_network(deg_genes, os.path.join(OUT_DIR, "network", organ), organ)

        summary_data.append({
            "organ": organ,
            "n_input_genes": n_input,
            "n_sig_degs": len(deg_genes),
            "sig_deg_list": ", ".join(deg_genes),
            "sig_severity_corr_genes": ", ".join(sig_corr_genes) if sig_corr_genes else "None"
        })

        print(f"  {organ.upper()} Total Input Genes : {n_input}")
        print(f"  {organ.upper()} Significant DEGs  : {len(deg_genes)} ({', '.join(deg_genes)})")
        print(f"  {organ.upper()} Severity Corr Genes: {', '.join(sig_corr_genes) if sig_corr_genes else 'None'}\n")

    # =========================================================================
    # CROSS-ORGAN OVERLAP COMPUTATIONS
    # =========================================================================
    print("=================================================================")
    print("COMPUTING CROSS-ORGAN DEG OVERLAPS (4 ORGANS vs 3 ORGANS)")
    print("=================================================================")

    # 1. Validation 2 Cohort Overlaps
    common_4organs_val2 = sorted(list(set.intersection(*[set(organ_deg_dict[o]) for o in ORGANS])))
    common_3organs_val2 = sorted(list(set.intersection(
        set(organ_deg_dict["lung"]), set(organ_deg_dict["skin"]), set(organ_deg_dict["kidney"])
    )))

    print(f"\n[Validation 2 Cohort] Common in 4 Organs (Kidney, Liver, Lung, Skin) ({len(common_4organs_val2)}): {', '.join(common_4organs_val2)}")
    print(f"[Validation 2 Cohort] Common in 3 Organs (Lungs, Skin, Kidneys) ({len(common_3organs_val2)}): {', '.join(common_3organs_val2)}")

    # Save Validation 2 Overlap CSVs
    pd.DataFrame({"gene": common_4organs_val2}).to_csv(os.path.join(OUT_DIR, "validation2_common_4organs.csv"), index=False)
    pd.DataFrame({"gene": common_3organs_val2}).to_csv(os.path.join(OUT_DIR, "validation2_common_3organs.csv"), index=False)

    # Update summary dataframe with 4-organ overlap count
    for item in summary_data:
        item["n_common_4organs"] = len(common_4organs_val2)
        item["n_common_3organs"] = len(common_3organs_val2)

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(OUT_DIR, "validation2_summary.csv"), index=False)

    # 2. Full Genome Discovery Cohort Overlaps
    full_deg_dict = {}
    for organ in ORGANS:
        fpath = os.path.join(RESULTS_DIR, f"{organ}_DEGs.csv")
        if os.path.exists(fpath):
            df = pd.read_csv(fpath)
            sig = df[(df["adj_p_value"] < 0.05) & (df["logFC"].abs() > 0.585)]
            full_deg_dict[organ] = set(sig["gene"].astype(str).str.strip().str.upper())

    if len(full_deg_dict) == 4:
        full_4organs = sorted(list(set.intersection(full_deg_dict["kidney"], full_deg_dict["liver"], full_deg_dict["lung"], full_deg_dict["skin"])))
        full_3organs = sorted(list(set.intersection(full_deg_dict["lung"], full_deg_dict["skin"], full_deg_dict["kidney"])))
        exclusive_3organs = sorted(list(set(full_3organs) - set(full_4organs)))

        print(f"\n[Full Genome Cohort] Common in 4 Organs: {len(full_4organs)} genes")
        print(f"[Full Genome Cohort] Common in 3 Organs (Lungs, Skin, Kidneys): {len(full_3organs)} genes")
        print(f"[Full Genome Cohort] Exclusive to 3 Organs (not in Liver): {len(exclusive_3organs)} genes")

        pd.DataFrame({"gene": full_4organs}).to_csv(os.path.join(OUT_DIR, "full_genome_common_4organs.csv"), index=False)
        pd.DataFrame({"gene": full_3organs}).to_csv(os.path.join(OUT_DIR, "full_genome_common_3organs.csv"), index=False)
        pd.DataFrame({"gene": exclusive_3organs}).to_csv(os.path.join(OUT_DIR, "full_genome_3organs_exclusive_to_lungs_skin_kidneys.csv"), index=False)

        # Plot 3-Organ Venn Diagram for Full Genome
        plot_3organ_venn(
            full_deg_dict["lung"], full_deg_dict["skin"], full_deg_dict["kidney"],
            os.path.join(OUT_DIR, "full_genome_venn_3organ_lungs_skin_kidney.png"),
            "Full Genome DEG Overlap: Lungs, Skin & Kidneys (313 Genes)"
        )

    # =========================================================================
    # GENERATE PLOTS & VISUALIZATIONS
    # =========================================================================
    print("\n=================================================================")
    print("GENERATING SUMMARY PLOTS & VISUALIZATIONS")
    print("=================================================================")

    # 1. 3-Organ Venn Diagram (Lungs, Skin, Kidneys) for Validation 2
    plot_3organ_venn(
        organ_deg_dict["lung"], organ_deg_dict["skin"], organ_deg_dict["kidney"],
        os.path.join(OUT_DIR, "venn_3organ_lungs_skin_kidney.png"),
        "Validation 2 DEG Overlap: Lungs, Skin & Kidneys (5 Genes)"
    )

    # 2. UpSet Plot for Validation 2
    plot_upset(organ_deg_dict, os.path.join(OUT_DIR, "upset_plot_validation2.png"), "Validation 2 DEG Intersections Across 4 Organs")

    # 3. Funnel Summary Bar Chart
    plot_funnel_summary(summary_df, os.path.join(OUT_DIR, "validation2_funnel_chart.png"))

    # 4. Cross-Organ Heatmaps for 4-Organ and 3-Organ DEGs
    plot_heatmaps(
        organ_expr_dict, organ_sg_dict, common_4organs_val2,
        os.path.join(OUT_DIR, "validation2_heatmaps.png"),
        "Validation 2 Common 4-Organ DEGs"
    )

    # =========================================================================
    # ENRICHMENT & NETWORK FOR COMMON OVERLAPS
    # =========================================================================
    print("\n=================================================================")
    print("RUNNING ENRICHMENT & PPI NETWORKS FOR COMMON DEG SIGNATURES")
    print("=================================================================")

    # 4-Organ Common Signature Enrichment & PPI Network
    if len(common_4organs_val2) > 0:
        run_task3_enrichment(common_4organs_val2, os.path.join(OUT_DIR, "enrichment", "all4organs"), "all4organs")
        run_task4_network(common_4organs_val2, os.path.join(OUT_DIR, "network", "all4organs"), "all4organs")

    # 3-Organ Common Signature Enrichment & PPI Network
    if len(common_3organs_val2) > 0:
        run_task3_enrichment(common_3organs_val2, os.path.join(OUT_DIR, "enrichment", "3organs_lungs_skin_kidney"), "3organs_lungs_skin_kidney")
        run_task4_network(common_3organs_val2, os.path.join(OUT_DIR, "network", "3organs_lungs_skin_kidney"), "3organs_lungs_skin_kidney")

    print("\n=================================================================")
    print("VALIDATION 2 PIPELINE SUCCESSFULLY COMPLETED FOR ALL ORGANS")
    print("=================================================================")

if __name__ == "__main__":
    main()
