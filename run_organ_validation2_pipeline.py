#!/usr/bin/env python3
"""
Multi-Organ Validation 2 Pipeline
=================================
Automated scanning, inspection, differential expression validation (Mann-Whitney U),
and severity/dose-response validation (Spearman correlation & Kruskal-Wallis) across 
Liver, Kidney, Lung, and Skin cohorts.

Authoritative Assumptions (Configurable below):
- Organ directories: ['kidney', 'liver', 'lung', 'skin']
- Validation2 search paths: '{organ}/validation2', '{organ}/validate 2', 'organ_validation2_data/{organ}'
- File extensions supported: .tsv, .csv, .txt
- Group labels: 'fibrotic' vs 'healthy' (or 'control', 'normal')
- Severity columns: 'severity', 'severity_numeric', 'stage', 'fibrosis_stage', 'score'
"""

import os
import sys
import glob
import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# Force UTF-8 output stdout encoding if possible
if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

# =============================================================================
# USER CONFIGURATION & ASSUMPTIONS
# =============================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "validation_2")
PLOTS_DIR = os.path.join(OUT_DIR, "plots")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

ORGANS = ["kidney", "liver", "lung", "skin"]

# Column name candidates for metadata auto-detection
GROUP_COLS = ["group", "condition", "status", "sample_group", "disease_state", "class"]
SAMPLE_ID_COLS = ["sample_id", "sample", "sample_name", "id", "gsm", "accession"]
SEVERITY_COLS = ["severity_numeric", "severity", "stage", "fibrosis_stage", "score", "tif_percent"]

# Values mapping to control vs fibrotic
CONTROL_VALS = ["control", "healthy", "ctrl", "normal", "f0", "stage0", "0"]
FIBROTIC_VALS = ["fibrotic", "fibrosis", "disease", "ckd", "cirrhosis", "ipf", "ssc"]

# Gene Panel & Discovery Direction file fallbacks
GENE_PANEL_CANDIDATES = [
    os.path.join(BASE_DIR, "final_gene_panel.csv"),
    os.path.join(RESULTS_DIR, "final_validated_pan_fibrotic_genes.csv"),
    os.path.join(RESULTS_DIR, "common_all_4_tissues_ecm_genes.csv"),
    os.path.join(RESULTS_DIR, "common_all_4_tissues_genes.csv"),
]

DISCOVERY_DIR_CANDIDATES = [
    os.path.join(BASE_DIR, "discovery_directions.csv"),
    os.path.join(RESULTS_DIR, "discovery_directions.csv"),
]

# Gene alias map for common symbol mismatches
GENE_ALIASES = {
    "CCN2": ["CTGF"],
    "CTGF": ["CCN2"],
    "COL1A1": ["COL1A1"],
    "COL1A2": ["COL1A2"],
    "COL3A1": ["COL3A1"],
    "AEBP1": ["ACL2"],
    "VWF": ["F8VWF"],
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_gene_panel():
    """Load final consensus gene panel."""
    for path in GENE_PANEL_CANDIDATES:
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                col = next((c for c in ["gene", "Gene", "symbol", "Gene Symbol"] if c in df.columns), df.columns[0])
                genes = df[col].dropna().astype(str).str.strip().str.upper().unique().tolist()
                print(f"[PANEL] Loaded {len(genes)} panel genes from: {os.path.basename(path)}")
                return genes, df
            except Exception as e:
                print(f"[PANEL WARNING] Could not read {path}: {e}")
    
    fallback = ["AEBP1", "COL1A1", "COL1A2", "COL3A1", "VWF"]
    print(f"[PANEL] Using default 5-gene ECM panel: {fallback}")
    return fallback, pd.DataFrame({"gene": fallback})


def load_discovery_directions(organ):
    """Load expected discovery directions (up/down) per gene."""
    for path in DISCOVERY_DIR_CANDIDATES:
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                if "organ" in df.columns:
                    df = df[df["organ"].astype(str).str.lower() == organ.lower()]
                gene_col = next((c for c in ["gene", "Gene", "symbol"] if c in df.columns), None)
                dir_col = next((c for c in ["direction", "discovery_direction", "logFC", "log2FoldChange"] if c in df.columns), None)
                
                if gene_col and dir_col:
                    res = {}
                    for g, val in zip(df[gene_col], df[dir_col]):
                        g_str = str(g).strip().upper()
                        if isinstance(val, (int, float)):
                            res[g_str] = "up" if val > 0 else "down"
                        else:
                            val_s = str(val).strip().lower()
                            res[g_str] = "up" if "up" in val_s or "+" in val_s else "down"
                    return res
            except Exception as e:
                pass

    return {"AEBP1": "up", "COL1A1": "up", "COL1A2": "up", "COL3A1": "up", "VWF": "up"}


def find_organ_validation_files(organ):
    """Scan and find validation2 files for a given organ."""
    search_patterns = [
        os.path.join(BASE_DIR, organ, "validation2", "*.*"),
        os.path.join(BASE_DIR, organ, "validate 2", "*.*"),
        os.path.join(BASE_DIR, organ, "Validation 2", "*.*"),
        os.path.join(BASE_DIR, organ.capitalize(), "validation2", "*.*"),
        os.path.join(BASE_DIR, organ.capitalize(), "validate 2", "*.*"),
        os.path.join(BASE_DIR, organ.capitalize(), "Validation 2", "*.*"),
        os.path.join(BASE_DIR, "organ_validation2_data", organ, "*.*"),
        os.path.join(BASE_DIR, "organ_validation2_data", organ.capitalize(), "*.*"),
    ]
    
    files = []
    for pat in search_patterns:
        matches = glob.glob(pat)
        for m in matches:
            if os.path.isfile(m) and m.endswith((".tsv", ".csv", ".txt")) and not m.endswith(".png"):
                files.append(m)
                
    return sorted(list(set(files)))


def inspect_and_parse_files(file_list):
    """
    Inspect files, print shape/columns/head, auto-detect expression matrix vs metadata.
    Returns: (expr_df, meta_df)
    """
    expr_df = None
    meta_df = None
    
    print(f"\n  Found {len(file_list)} candidate file(s):")
    for fpath in file_list:
        fname = os.path.basename(fpath)
        sep = "\t" if fpath.endswith((".tsv", ".txt")) else ","
        try:
            df = pd.read_csv(fpath, sep=sep)
        except Exception as e:
            print(f"    [READ ERROR] Could not parse {fname}: {e}")
            continue
            
        print(f"\n  -------------------------------------------------------------")
        print(f"  FILE INSPECTION: {fname}")
        print(f"    Path:  {fpath}")
        print(f"    Shape: {df.shape}")
        print(f"    Cols:  {list(df.columns[:8])}{'...' if len(df.columns)>8 else ''}")
        print(f"    First 5 rows:")
        print(df.head(5).to_string())
        print(f"  -------------------------------------------------------------")
        
        # Auto-detect Metadata vs Expression Matrix
        col_lower = [str(c).lower() for c in df.columns]
        is_metadata = any(gc in col_lower for gc in GROUP_COLS) or any(sc in col_lower for sc in SEVERITY_COLS)
        
        if is_metadata:
            meta_df = df
            print(f"    => AUTO-DETECTED as Sample Metadata file [OK]")
        else:
            # If numeric matrix or top.table with gene column
            expr_df = df
            print(f"    => AUTO-DETECTED as Gene Expression Data [OK]")
            
    return expr_df, meta_df


def prepare_expression_matrix(expr_df):
    """
    Auto-detect orientation (genes in rows vs genes in columns) and format dataframe.
    Returns: DataFrame with genes as index, samples as columns.
    """
    if expr_df is None or len(expr_df) == 0:
        return None
        
    df = expr_df.copy()
    
    # Check for top.table style dataframe (ID, logFC, Gene.symbol)
    sym_col = next((c for c in df.columns if str(c).lower() in ["gene", "symbol", "gene.symbol", "gene_symbol"]), None)
    
    if sym_col and len(df.columns) <= 15:
        # It is a top table / DE summary table. Convert to synthetic matrix or extraction
        df[sym_col] = df[sym_col].astype(str).str.strip().str.upper()
        df = df.drop_duplicates(subset=sym_col, keep="first").set_index(sym_col)
        return df
        
    first_col = df.columns[0]
    if df[first_col].dtype == object or isinstance(df[first_col].iloc[0], str):
        df[first_col] = df[first_col].astype(str).str.strip().str.upper()
        df = df.drop_duplicates(subset=first_col, keep="first").set_index(first_col)
    else:
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) < len(df.columns) / 2:
            df = df.transpose()
            df.index = df.index.astype(str).str.strip().str.upper()
            
    return df


def prepare_sample_metadata(meta_df, sample_names):
    """
    Extract sample_id, group (healthy/fibrotic), and severity_numeric if available.
    """
    if meta_df is None or len(meta_df) == 0:
        records = []
        for s in sample_names:
            s_low = str(s).lower()
            group = "fibrotic" if any(f in s_low for f in ["fib", "dis", "case", "stage", "patient"]) else "control"
            records.append({"sample_id": s, "group": group, "severity_numeric": np.nan})
        return pd.DataFrame(records)
        
    df = meta_df.copy()
    cols_low = {str(c).lower(): c for c in df.columns}
    
    id_col = next((cols_low[c] for c in SAMPLE_ID_COLS if c in cols_low), df.columns[0])
    group_col = next((cols_low[c] for c in GROUP_COLS if c in cols_low), None)
    sev_col = next((cols_low[c] for c in SEVERITY_COLS if c in cols_low), None)
    
    res = pd.DataFrame()
    res["sample_id"] = df[id_col].astype(str).str.strip()
    
    if group_col:
        raw_grp = df[group_col].astype(str).str.strip().str.lower()
        res["group"] = raw_grp.apply(lambda x: "control" if any(c in x for c in CONTROL_VALS) else "fibrotic")
    else:
        res["group"] = res["sample_id"].apply(lambda x: "control" if any(c in str(x).lower() for c in CONTROL_VALS) else "fibrotic")
        
    if sev_col:
        res["severity_numeric"] = pd.to_numeric(df[sev_col], errors="coerce")
    else:
        res["severity_numeric"] = np.nan
        
    return res


def match_gene_symbols(panel_genes, expr_genes):
    """
    Match panel genes against expression matrix genes using case-insensitive & alias matching.
    Returns: dict {panel_gene -> expr_gene_name_or_None}
    """
    expr_upper_map = {str(g).strip().upper(): g for g in expr_genes}
    matched = {}
    missing = []
    
    for p_gene in panel_genes:
        p_upper = str(p_gene).strip().upper()
        if p_upper in expr_upper_map:
            matched[p_gene] = expr_upper_map[p_upper]
        else:
            found = False
            aliases = GENE_ALIASES.get(p_upper, [])
            for alias in aliases:
                alias_up = alias.upper()
                if alias_up in expr_upper_map:
                    matched[p_gene] = expr_upper_map[alias_up]
                    found = True
                    print(f"    [ALIAS MATCH] '{p_gene}' matched to '{expr_upper_map[alias_up]}'")
                    break
            if not found:
                matched[p_gene] = None
                missing.append(p_gene)
                
    if missing:
        print(f"    [MISSING GENES] {len(missing)} panel gene(s) missing in expression data: {missing}")
    else:
        print(f"    [GENE MATCH] 100% of panel genes ({len(panel_genes)}) found in expression data [OK]")
        
    return matched


# =============================================================================
# MAIN ANALYSIS PIPELINE PER ORGAN
# =============================================================================

def process_single_organ(organ, panel_genes, disc_directions):
    print(f"\n" + "="*80)
    print(f"PROCESSING ORGAN: {organ.upper()}")
    print("="*80)
    
    files = find_organ_validation_files(organ)
    if not files:
        print(f"  [ERROR] No validation2 files found for organ '{organ}'! Skipping.")
        return None, None
        
    expr_df_raw, meta_df_raw = inspect_and_parse_files(files)
    
    if expr_df_raw is None:
        print(f"  [ERROR] Expression matrix could not be parsed for '{organ}'! Skipping.")
        return None, None
        
    expr_mat = prepare_expression_matrix(expr_df_raw)
    sample_meta = prepare_sample_metadata(meta_df_raw, expr_mat.columns)
    
    # Align samples between expression matrix and metadata
    valid_samples = [s for s in sample_meta["sample_id"] if s in expr_mat.columns]
    if len(valid_samples) == 0:
        valid_samples = list(expr_mat.columns)
        sample_meta = prepare_sample_metadata(None, valid_samples)
        
    # Check if expr_mat has expression numeric columns
    expr_numeric_cols = [c for c in valid_samples if c in expr_mat.columns]
    
    if len(expr_numeric_cols) > 0:
        expr_mat = expr_mat[expr_numeric_cols]
        sample_meta = sample_meta[sample_meta["sample_id"].isin(expr_numeric_cols)].reset_index(drop=True)
    
    print(f"\n  Final Aligned Data for {organ.upper()}:")
    print(f"    Samples: {len(sample_meta)} (Fibrotic: {(sample_meta['group']=='fibrotic').sum()}, Healthy: {(sample_meta['group']=='control').sum()})")
    print(f"    Genes in dataset: {len(expr_mat):,}")
    
    # Gene Matching
    gene_map = match_gene_symbols(panel_genes, expr_mat.index)
    
    # Check severity availability
    has_severity = sample_meta["severity_numeric"].notna().sum() >= 5
    sev_status = "AVAILABLE" if has_severity else "NOT AVAILABLE for this organ's validation2 data"
    print(f"    Severity metadata status: {sev_status}")
    
    results = []
    ctrl_samples = sample_meta[sample_meta["group"] == "control"]["sample_id"].tolist()
    fib_samples  = sample_meta[sample_meta["group"] == "fibrotic"]["sample_id"].tolist()
    
    for p_gene in panel_genes:
        mapped_symbol = gene_map[p_gene]
        expected_dir = disc_directions.get(p_gene, "up")
        
        if mapped_symbol is None or mapped_symbol not in expr_mat.index:
            results.append({
                "gene": p_gene,
                "discovery_direction": expected_dir,
                "validation2_direction": "N/A",
                "direction_match": False,
                "de_pvalue_raw": np.nan,
                "de_pvalue_adj": np.nan,
                "de_significant": False,
                "severity_rho_full": np.nan,
                "severity_p_full_raw": np.nan,
                "severity_p_full_adj": np.nan,
                "severity_rho_disease_only": np.nan,
                "severity_p_disease_only_raw": np.nan,
                "severity_p_disease_only_adj": np.nan,
                "kruskal_pvalue": np.nan,
                "severity_test_available": has_severity,
                "present_in_data": False,
            })
            continue
            
        gene_row = expr_mat.loc[mapped_symbol]
        
        # Check if row is expression series or DE top table row
        if len(expr_numeric_cols) > 0 and len(ctrl_samples) > 0 and len(fib_samples) > 0:
            ctrl_vals = pd.to_numeric(gene_row[ctrl_samples], errors="coerce").dropna().values
            fib_vals  = pd.to_numeric(gene_row[fib_samples], errors="coerce").dropna().values
            
            # --- TEST 1: Differential Expression (Mann-Whitney U Test) ---
            if len(ctrl_vals) > 0 and len(fib_vals) > 0:
                mwu_stat, mwu_p = stats.mannwhitneyu(fib_vals, ctrl_vals, alternative="two-sided")
                lfc = float(np.mean(fib_vals) - np.mean(ctrl_vals))
                val_dir = "up" if lfc > 0 else "down"
                dir_match = (val_dir == expected_dir)
            else:
                mwu_p, val_dir, dir_match = np.nan, "N/A", False
        else:
            # Fallback for top-table DE output
            mwu_p_col = next((c for c in ["P.Value", "pvalue", "PValue", "adj.P.Val", "padj"] if c in gene_row), None)
            lfc_col = next((c for c in ["logFC", "log2FoldChange", "LFC"] if c in gene_row), None)
            
            mwu_p = float(gene_row[mwu_p_col]) if mwu_p_col else np.nan
            lfc_val = float(gene_row[lfc_col]) if lfc_col else 0.0
            val_dir = "up" if lfc_val > 0 else "down"
            dir_match = (val_dir == expected_dir)
            
        # --- TEST 2: Severity / Dose-Response Validation ---
        rho_full, p_full = np.nan, np.nan
        rho_dis, p_dis = np.nan, np.nan
        kw_p = np.nan
        
        if has_severity and len(expr_numeric_cols) > 0:
            df_sev = sample_meta[["sample_id", "group", "severity_numeric"]].copy()
            df_sev["expr"] = pd.to_numeric(gene_row[df_sev["sample_id"]].values, errors="coerce")
            df_sev = df_sev.dropna(subset=["expr", "severity_numeric"])
            
            # Full-sample Spearman
            if len(df_sev) >= 5:
                rho_full, p_full = stats.spearmanr(df_sev["severity_numeric"], df_sev["expr"])
                
            # Disease-only Spearman (excluding controls)
            df_dis = df_sev[df_sev["group"] == "fibrotic"]
            if len(df_dis) >= 4 and df_dis["severity_numeric"].nunique() > 1:
                rho_dis, p_dis = stats.spearmanr(df_dis["severity_numeric"], df_dis["expr"])
                
            # Kruskal-Wallis across severity stages
            sev_groups = [g["expr"].values for _, g in df_sev.groupby("severity_numeric") if len(g) >= 2]
            if len(sev_groups) >= 2:
                kw_stat, kw_p = stats.kruskal(*sev_groups)
                
        results.append({
            "gene": p_gene,
            "discovery_direction": expected_dir,
            "validation2_direction": val_dir,
            "direction_match": dir_match,
            "de_pvalue_raw": mwu_p,
            "de_pvalue_adj": np.nan,
            "de_significant": False,
            "severity_rho_full": rho_full,
            "severity_p_full_raw": p_full,
            "severity_p_full_adj": np.nan,
            "severity_rho_disease_only": rho_dis,
            "severity_p_disease_only_raw": p_dis,
            "severity_p_disease_only_adj": np.nan,
            "kruskal_pvalue": kw_p,
            "severity_test_available": has_severity,
            "present_in_data": True,
        })
        
    res_df = pd.DataFrame(results)
    
    # --- Benjamini-Hochberg FDR Corrections ---
    valid_de_mask = res_df["de_pvalue_raw"].notna()
    if valid_de_mask.sum() > 0:
        _, res_df.loc[valid_de_mask, "de_pvalue_adj"], _, _ = multipletests(
            res_df.loc[valid_de_mask, "de_pvalue_raw"], method="fdr_bh"
        )
        res_df["de_significant"] = (res_df["de_pvalue_adj"] < 0.05) & (res_df["direction_match"] == True)
        
    valid_full_mask = res_df["severity_p_full_raw"].notna()
    if valid_full_mask.sum() > 0:
        _, res_df.loc[valid_full_mask, "severity_p_full_adj"], _, _ = multipletests(
            res_df.loc[valid_full_mask, "severity_p_full_raw"], method="fdr_bh"
        )
        
    valid_dis_mask = res_df["severity_p_disease_only_raw"].notna()
    if valid_dis_mask.sum() > 0:
        _, res_df.loc[valid_dis_mask, "severity_p_disease_only_adj"], _, _ = multipletests(
            res_df.loc[valid_dis_mask, "severity_p_disease_only_raw"], method="fdr_bh"
        )
        
    out_cols = [
        "gene", "discovery_direction", "validation2_direction", "direction_match",
        "de_pvalue_adj", "de_significant", "severity_rho_full", "severity_p_full_adj",
        "severity_rho_disease_only", "severity_p_disease_only_adj", "severity_test_available"
    ]
    
    out_df = res_df[out_cols].copy()
    out_csv = os.path.join(OUT_DIR, f"{organ}_validation2_results.csv")
    out_df.to_csv(out_csv, index=False)
    print(f"\n  [SAVED] {organ.upper()} summary CSV: {out_csv}")
    
    # --- PLOTTING FOR TOP 3 GENES ---
    if len(expr_numeric_cols) > 0:
        top_genes = res_df.sort_values("de_pvalue_raw", ascending=True).head(3)["gene"].tolist()
        print(f"  Generating Box Plots & Scatter Plots for top 3 genes: {top_genes}...")
        
        for g in top_genes:
            mapped_s = gene_map[g]
            if mapped_s is None or mapped_s not in expr_mat.index:
                continue
                
            g_series = expr_mat.loc[mapped_s]
            df_plot = sample_meta.copy()
            df_plot["Expression"] = pd.to_numeric(g_series[df_plot["sample_id"]].values, errors="coerce")
            df_plot = df_plot.dropna(subset=["Expression"])
            
            if len(df_plot) == 0:
                continue
                
            # Box Plot
            fig, ax = plt.subplots(figsize=(6, 5))
            sns.boxplot(data=df_plot, x="group", y="Expression", palette={"control": "#5B8DB8", "fibrotic": "#E07B54"}, ax=ax, width=0.4)
            sns.stripplot(data=df_plot, x="group", y="Expression", color="black", alpha=0.6, jitter=0.2, ax=ax)
            p_adj_val = res_df.loc[res_df['gene']==g, 'de_pvalue_adj'].values[0]
            p_str = f"{p_adj_val:.4e}" if pd.notna(p_adj_val) else "N/A"
            ax.set_title(f"{organ.upper()} - {g} Expression by Group\n(Mann-Whitney p_adj = {p_str})", fontweight="bold")
            ax.set_xlabel("Group", fontweight="bold")
            ax.set_ylabel("Gene Expression (log scale)", fontweight="bold")
            box_png = os.path.join(PLOTS_DIR, f"{organ}_{g}_expression_boxplot.png")
            fig.savefig(box_png, dpi=200, bbox_inches="tight")
            plt.close(fig)
            
            # Scatter Plot with Regression (Severity)
            if has_severity and df_plot["severity_numeric"].notna().sum() >= 5:
                df_scat = df_plot.dropna(subset=["severity_numeric"])
                fig, ax = plt.subplots(figsize=(6.5, 5))
                sns.regplot(data=df_scat, x="severity_numeric", y="Expression", color="#2c3e50",
                            scatter_kws={"s": 50, "alpha": 0.7}, line_kws={"color": "#e74c3c", "lw": 2}, ax=ax)
                rho_val = res_df.loc[res_df['gene']==g, 'severity_rho_full'].values[0]
                p_val = res_df.loc[res_df['gene']==g, 'severity_p_full_adj'].values[0]
                rho_str = f"{rho_val:.3f}" if pd.notna(rho_val) else "N/A"
                p_str2 = f"{p_val:.4e}" if pd.notna(p_val) else "N/A"
                ax.set_title(f"{organ.upper()} - {g} Expression vs Clinical Severity\n(Spearman rho = {rho_str}, p_adj = {p_str2})", fontweight="bold")
                ax.set_xlabel("Clinical Severity Score / Stage", fontweight="bold")
                ax.set_ylabel("Gene Expression", fontweight="bold")
                scat_png = os.path.join(PLOTS_DIR, f"{organ}_{g}_severity_scatterplot.png")
                fig.savefig(scat_png, dpi=200, bbox_inches="tight")
                plt.close(fig)
                
    return res_df, out_df


# =============================================================================
# PIPELINE EXECUTION & FUNNEL SUMMARY
# =============================================================================

def main():
    print("="*80)
    print("RUNNING MULTI-ORGAN VALIDATION 2 PIPELINE")
    print("="*80)
    
    panel_genes, _ = load_gene_panel()
    
    all_organ_results = {}
    all_summary_rows = []
    funnel_stats = {}
    
    for organ in ORGANS:
        disc_dirs = load_discovery_directions(organ)
        full_res, summary_df = process_single_organ(organ, panel_genes, disc_dirs)
        
        if summary_df is not None:
            summary_df["organ"] = organ
            all_organ_results[organ] = full_res
            all_summary_rows.append(summary_df)
            
            n_panel = len(panel_genes)
            n_present = full_res["present_in_data"].sum()
            n_de_sig = ((full_res["de_pvalue_adj"] < 0.05) & (full_res["direction_match"] == True)).sum()
            n_sev_full = (full_res["severity_p_full_adj"] < 0.05).sum() if full_res["severity_p_full_adj"].notna().sum() > 0 else 0
            n_sev_dis = (full_res["severity_p_disease_only_adj"] < 0.05).sum() if full_res["severity_p_disease_only_adj"].notna().sum() > 0 else 0
            
            funnel_stats[organ] = {
                "panel": n_panel,
                "present": n_present,
                "de_sig": n_de_sig,
                "sev_full": n_sev_full,
                "sev_dis": n_sev_dis,
            }

    # Combined Cross-Organ Summary CSV
    if all_summary_rows:
        combined_df = pd.concat(all_summary_rows, ignore_index=True)
        cols = ["organ"] + [c for c in combined_df.columns if c != "organ"]
        combined_df = combined_df[cols]
        combined_csv = os.path.join(OUT_DIR, "all_organs_validation2_summary.csv")
        combined_df.to_csv(combined_csv, index=False)
        print(f"\n[SAVED] Combined Cross-Organ Summary CSV: {combined_csv}")

    # --- FINAL FUNNEL SUMMARY TO CONSOLE ---
    print("\n" + "="*80)
    print("FINAL VALIDATION 2 FUNNEL SUMMARY REPORT")
    print("="*80)
    print(f"{'Organ':<10} {'(a) Panel Genes':>16} {'(b) Present in Val2':>20} {'(c) DE Sig + Direction':>22} {'(d) Full Severity Sig':>22} {'(e) Disease Severity Sig':>24}")
    print("-" * 115)
    for organ in ORGANS:
        if organ in funnel_stats:
            st = funnel_stats[organ]
            print(f"{organ.capitalize():<10} {st['panel']:>16,} {st['present']:>20,} {st['de_sig']:>22,} {st['sev_full']:>22,} {st['sev_dis']:>24,}")
        else:
            print(f"{organ.capitalize():<10} {'N/A':>16} {'N/A':>20} {'N/A':>22} {'N/A':>22} {'N/A':>24}")
    print("="*115)
    print("\nVALIDATION 2 PIPELINE COMPLETED SUCCESSFULLY.")


if __name__ == "__main__":
    main()
