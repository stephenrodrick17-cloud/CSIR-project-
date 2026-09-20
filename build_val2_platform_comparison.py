# -*- coding: utf-8 -*-
"""
Validation 2 Platform-Stratified Comparison (Microarray vs RNA-seq)
===================================================================
Produces:
1. results/val2_platform_stratified_comparison.csv
2. plots/val2_platform_stratified_comparison.png
"""
import os
import glob
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = r"D:\CSIR"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
PLOTS_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

# 1. Load 24 Clean Core ECM Genes
ecm_df = pd.read_csv(os.path.join(RESULTS_DIR, "ecm_clean_genes.csv"))
genes = sorted(ecm_df["gene"].tolist())

# 2. Build Probe Dictionary
global_bio_map = {}
for f in glob.glob(os.path.join(BASE_DIR, "**", "bioDBnet*.txt"), recursive=True):
    with open(f, "r", encoding="utf-8", errors="ignore") as fp:
        for line in fp:
            parts = line.strip().split("\t")
            if len(parts) >= 2 and parts[1].strip().upper() not in ("-", "NONE", "NAN", "NA"):
                global_bio_map[parts[0].strip()] = parts[1].strip().upper()

acc_map_file = os.path.join(BASE_DIR, "geo_cache", "gb_acc_to_symbol_map.json")
if os.path.exists(acc_map_file):
    with open(acc_map_file, "r") as f:
        global_bio_map.update(json.load(f))

def prep_table(fpath):
    df = pd.read_csv(fpath, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    cols_l = {c.lower(): c for c in df.columns}
    lfc_col = next((cols_l[c] for c in ["logfc", "log2foldchange", "lfc"] if c in cols_l), None)
    p_col = next((cols_l[c] for c in ["adj.p.val", "padj", "adj_p_value", "fdr"] if c in cols_l), None)
    sym_col = next((cols_l[c] for c in ["symbol", "gene.symbol", "gene_symbol", "genesymbol"] if c in cols_l), None)
    id_col = cols_l.get("id", df.columns[0])
    if sym_col:
        df["gene"] = df[sym_col].astype(str).str.strip().str.upper()
    else:
        df["gene"] = df[id_col].astype(str).str.strip().map(global_bio_map)
        df["gene"] = df["gene"].fillna(df[id_col].astype(str).str.strip().str.upper())
    df["lfc_clean"] = pd.to_numeric(df[lfc_col], errors="coerce")
    df["padj_clean"] = pd.to_numeric(df[p_col], errors="coerce")
    return df.dropna(subset=["gene", "lfc_clean", "padj_clean"]).sort_values("padj_clean").drop_duplicates("gene", keep="first").set_index("gene")

# Load Val 2 Datasets
val2_k = prep_table(r"D:\CSIR\Kidney\Validation 2\GSE30529.top.table.tsv")
val2_l = prep_table(r"D:\CSIR\Liver\validate 2\GSE14323.top.table.tsv")
val2_lu = prep_table(r"D:\CSIR\Lungs\validate 2\GSE83717.top.table.tsv")
val2_s = prep_table(r"D:\CSIR\Skin\validation 2\GSE125362.top.table.tsv")

rows = []
for _, r in ecm_df.iterrows():
    g = r["gene"]
    cat = r["matrisome_category"]
    
    # Microarray: Kidney (GSE30529)
    k_lfc = val2_k.loc[g, "lfc_clean"] if g in val2_k.index else np.nan
    k_p = val2_k.loc[g, "padj_clean"] if g in val2_k.index else np.nan
    k_sig = bool(k_p < 0.05) if pd.notna(k_p) else False
    
    # Microarray: Liver (GSE14323)
    l_lfc = val2_l.loc[g, "lfc_clean"] if g in val2_l.index else np.nan
    l_p = val2_l.loc[g, "padj_clean"] if g in val2_l.index else np.nan
    l_sig = bool(l_p < 0.05) if pd.notna(l_p) else False
    
    # Microarray Concordance
    if pd.notna(k_lfc) and pd.notna(l_lfc):
        m_conc = bool(np.sign(k_lfc) == np.sign(l_lfc))
        m_dir = np.sign(k_lfc) if m_conc else 0
    else:
        m_conc = False
        m_dir = 0
        
    # RNA-seq: Lung (GSE83717)
    lu_lfc = val2_lu.loc[g, "lfc_clean"] if g in val2_lu.index else np.nan
    lu_p = val2_lu.loc[g, "padj_clean"] if g in val2_lu.index else np.nan
    lu_sig = bool(lu_p < 0.05) if pd.notna(lu_p) else False
    
    # RNA-seq: Skin (GSE125362)
    s_lfc = val2_s.loc[g, "lfc_clean"] if g in val2_s.index else np.nan
    s_p = val2_s.loc[g, "padj_clean"] if g in val2_s.index else np.nan
    s_sig = bool(s_p < 0.05) if pd.notna(s_p) else False
    
    # RNA-seq Concordance
    if pd.notna(lu_lfc) and pd.notna(s_lfc):
        r_conc = bool(np.sign(lu_lfc) == np.sign(s_lfc))
        r_dir = np.sign(lu_lfc) if r_conc else 0
    elif pd.notna(lu_lfc):
        # Single observation available
        r_conc = True
        r_dir = np.sign(lu_lfc)
    elif pd.notna(s_lfc):
        r_conc = True
        r_dir = np.sign(s_lfc)
    else:
        r_conc = False
        r_dir = 0
        
    # Cross-Platform Concordance (Microarray matches RNA-seq)
    if m_conc and r_conc and m_dir != 0 and r_dir != 0:
        cross_conc = bool(m_dir == r_dir)
    else:
        cross_conc = False
        
    rows.append({
        "gene": g,
        "matrisome_category": cat,
        "kidney_microarray_logFC": round(k_lfc, 3) if pd.notna(k_lfc) else np.nan,
        "kidney_microarray_adjp": f"{k_p:.3e}" if pd.notna(k_p) else np.nan,
        "kidney_microarray_sig": k_sig,
        "liver_microarray_logFC": round(l_lfc, 3) if pd.notna(l_lfc) else np.nan,
        "liver_microarray_adjp": f"{l_p:.3e}" if pd.notna(l_p) else np.nan,
        "liver_microarray_sig": l_sig,
        "microarray_direction_concordant": m_conc,
        "lung_rnaseq_logFC": round(lu_lfc, 3) if pd.notna(lu_lfc) else np.nan,
        "lung_rnaseq_adjp": f"{lu_p:.3e}" if pd.notna(lu_p) else np.nan,
        "lung_rnaseq_sig": lu_sig,
        "skin_rnaseq_logFC": round(s_lfc, 3) if pd.notna(s_lfc) else np.nan,
        "skin_rnaseq_adjp": f"{s_p:.3e}" if pd.notna(s_p) else np.nan,
        "skin_rnaseq_sig": s_sig,
        "rnaseq_direction_concordant": r_conc,
        "cross_platform_concordant": cross_conc
    })

df_table = pd.DataFrame(rows)

# Save to results
out_csv = os.path.join(RESULTS_DIR, "val2_platform_stratified_comparison.csv")
df_table.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
n_m_conc = df_table["microarray_direction_concordant"].sum()
n_r_conc = df_table["rnaseq_direction_concordant"].sum()
n_cross_conc = df_table["cross_platform_concordant"].sum()

# Identify interesting cases: concordant in Microarray but discordant across platforms
# E.g. m_conc is True, but cross_platform_concordant is False
disagree_genes = df_table[df_table["microarray_direction_concordant"] & (~df_table["cross_platform_concordant"])]["gene"].tolist()

print("\n" + "=" * 80)
print("PLATFORM STRATIFICATION SUMMARY STATISTICS")
print("=" * 80)
print(f"Total Clean Core ECM Genes Evaluated        : {len(df_table)}")
print(f"Concordant WITHIN Microarray (Kidney vs Liver): {n_m_conc} / {len(df_table)} ({n_m_conc/len(df_table)*100:.1f}%)")
print(f"Concordant WITHIN RNA-seq (Lung vs Skin)     : {n_r_conc} / {len(df_table)} ({n_r_conc/len(df_table)*100:.1f}%)")
print(f"Concordant ACROSS Platforms (Array == RNA-seq): {n_cross_conc} / {len(df_table)} ({n_cross_conc/len(df_table)*100:.1f}%)")
print(f"\nGenes Concordant in Microarray but DISAGREE with RNA-seq ({len(disagree_genes)} genes):")
print(", ".join(disagree_genes))

# =============================================================================
# VISUALIZATION: Grouped Heatmap with Microarray vs RNA-seq Divider
# =============================================================================
# Prepare matrix for heatmap
mat_data = []
for _, r in df_table.iterrows():
    mat_data.append([
        r["kidney_microarray_logFC"],
        r["liver_microarray_logFC"],
        r["lung_rnaseq_logFC"],
        r["skin_rnaseq_logFC"]
    ])

heatmap_df = pd.DataFrame(mat_data, index=df_table["gene"], 
                          columns=["Kidney\n(GSE30529)", "Liver\n(GSE14323)", 
                                   "Lung\n(GSE83717)", "Skin\n(GSE125362)"])

# Add significance asterisk annotations
annot_matrix = []
for _, r in df_table.iterrows():
    k_ast = "*" if r["kidney_microarray_sig"] else ""
    l_ast = "*" if r["liver_microarray_sig"] else ""
    lu_ast = "*" if r["lung_rnaseq_sig"] else ""
    s_ast = "*" if r["skin_rnaseq_sig"] else ""
    
    k_txt = f"{r['kidney_microarray_logFC']:.2f}{k_ast}" if pd.notna(r['kidney_microarray_logFC']) else "NA"
    l_txt = f"{r['liver_microarray_logFC']:.2f}{l_ast}" if pd.notna(r['liver_microarray_logFC']) else "NA"
    lu_txt = f"{r['lung_rnaseq_logFC']:.2f}{lu_ast}" if pd.notna(r['lung_rnaseq_logFC']) else "NA"
    s_txt = f"{r['skin_rnaseq_logFC']:.2f}{s_ast}" if pd.notna(r['skin_rnaseq_logFC']) else "NA"
    annot_matrix.append([k_txt, l_txt, lu_txt, s_txt])

annot_df = pd.DataFrame(annot_matrix, index=df_table["gene"], columns=heatmap_df.columns)

fig, ax = plt.subplots(figsize=(11, 14))

# Diverging color map centered at 0
cmap = sns.diverging_palette(220, 20, as_cmap=True)

sns.heatmap(heatmap_df, annot=annot_df, fmt="", cmap="vlag", center=0,
            vmin=-3.0, vmax=3.0, cbar_kws={"label": "log2 Fold Change (Fibrosis vs Control)", "shrink": 0.8},
            linewidths=1.0, linecolor="#f1f5f9", ax=ax, annot_kws={"fontsize": 9, "fontweight": "bold"})

# Add Platform Headers & Shading / Vertical Divider
ax.axvline(2, color="#0f172a", linewidth=3.5, linestyle="-")

# Header banner texts
ax.text(1.0, -0.6, "MICROARRAY PLATFORMS", ha="center", va="center", fontsize=12, fontweight="bold", 
        color="#1e3a8a", bbox=dict(boxstyle="round,pad=0.4", fc="#dbeafe", ec="#1e3a8a", lw=1.5))

ax.text(3.0, -0.6, "RNA-SEQ PLATFORMS", ha="center", va="center", fontsize=12, fontweight="bold", 
        color="#065f46", bbox=dict(boxstyle="round,pad=0.4", fc="#d1fae5", ec="#065f46", lw=1.5))

plt.title("Validation 2 (Held-Out) Cross-Platform Stratification: Microarray vs RNA-seq\n(* indicates FDR adjusted p < 0.05)", 
          fontsize=14, fontweight="bold", pad=35)
plt.xlabel("Validation 2 Cohorts Stratified by Sequencing Technology", fontsize=12, fontweight="bold", labelpad=12)
plt.ylabel("24 Clean Core Extracellular Matrix (ECM) Genes", fontsize=12, fontweight="bold", labelpad=12)

plt.tight_layout()
out_png = os.path.join(PLOTS_DIR, "val2_platform_stratified_comparison.png")
plt.savefig(out_png, dpi=300)
plt.close()
print(f"Generated plot: {out_png}")
