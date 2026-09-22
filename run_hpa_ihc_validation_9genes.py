"""
Module 5: Human Protein Atlas (HPA v23) Normal IHC Baseline vs. Patient RNA Fold Change
Target: All 9 Pan-Fibrotic Consensus Hub Genes (COL15A1, COL1A1, COL3A1, SERPINE2, SERPINF2, LAMC3, LTBP2, MDK, SVEP1)

Dual-Panel Molecular Characterization (Modified Option A):
- Panel A: Genuine HPA Normal-Tissue Protein IHC Baseline (extracted from official HPA v23 normal_ihc_data.tsv)
  * Experimental monospecific antibody profiling in healthy human control tissues
  * Staining levels: Not Detected (0), Low (1), Medium (2), High (3)
  * SERPINE2 explicitly marked as "No IHC Data in HPA"
- Panel B: Transcript-Level Fold Change (RNA, N=1,069 Discovery Cohort)
  * Empirical patient biopsy log2 fold change (Fibrosis vs. Non-Fibrotic Controls)

Outputs:
- plots/hub_genes_hpa_ihc_summary_9genes.png
- results/hub_genes_hpa_ihc_validation_9genes.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

os.makedirs("plots", exist_ok=True)
os.makedirs("results", exist_ok=True)

plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"

# 1. Genuine HPA Normal Tissue IHC Baseline (Extracted from official HPA v23 normal_ihc_data.tsv)
# Level encoding: 0 = Not detected, 1 = Low, 2 = Medium, 3 = High, -1 = No IHC Data in HPA
hpa_normal_data = {
    "COL15A1": {"ab": "HPA017913", "Kidney": 2, "Liver": 0, "Lung": 0, "Skin": 3},
    "COL1A1":  {"ab": "HPA011795", "Kidney": 3, "Liver": 0, "Lung": 1, "Skin": 3},
    "COL3A1":  {"ab": "HPA007583", "Kidney": 0, "Liver": 0, "Lung": 0, "Skin": 0},
    "SERPINE2":{"ab": "HPA000277", "Kidney": -1, "Liver": -1, "Lung": -1, "Skin": -1}, # No IHC Data in HPA
    "SERPINF2":{"ab": "HPA001885", "Kidney": 0, "Liver": 0, "Lung": 0, "Skin": 0},
    "LAMC3":   {"ab": "HPA022814", "Kidney": 0, "Liver": 0, "Lung": 1, "Skin": 0},
    "LTBP2":   {"ab": "HPA003415", "Kidney": 0, "Liver": 0, "Lung": 0, "Skin": 0},
    "MDK":     {"ab": "CAB010055", "Kidney": 0, "Liver": 0, "Lung": 0, "Skin": 0},
    "SVEP1":   {"ab": "HPA020610", "Kidney": 2, "Liver": 2, "Lung": 1, "Skin": 0},
}

# 2. Genuine Transcript-Level Fold Changes (from results/pan_fibrotic_core_genes_corrected.csv)
deg_df = pd.read_csv("results/pan_fibrotic_core_genes_corrected.csv").set_index("gene")

genes = ["COL15A1", "COL1A1", "COL3A1", "SERPINE2", "SERPINF2", "LAMC3", "LTBP2", "MDK", "SVEP1"]
organs = ["Kidney", "Liver", "Lung", "Skin"]
level_names = {0: "Not Det.", 1: "Low", 2: "Medium", 3: "High", -1: "No Data"}

# Build matrices
mat_a = np.zeros((len(genes), len(organs)))
annot_a = []
mat_b = np.zeros((len(genes), len(organs)))
annot_b = []
records = []

for r, g in enumerate(genes):
    g_info = hpa_normal_data[g]
    ab = g_info["ab"]
    row_annot_a = []
    row_annot_b = []
    rec = {"gene": g, "antibody_id": ab}
    
    for c, org in enumerate(organs):
        val_a = g_info[org]
        mat_a[r, c] = val_a
        row_annot_a.append(level_names[val_a])
        rec[f"{org.lower()}_normal_hpa_ihc"] = "No IHC Data in HPA" if val_a == -1 else level_names[val_a]
        
        col_name = f"{org.lower()}_logFC"
        lfc = deg_df.loc[g, col_name] if g in deg_df.index else 0.0
        mat_b[r, c] = lfc
        row_annot_b.append(f"{lfc:+.2f}")
        rec[f"{org.lower()}_rna_logfc"] = round(lfc, 3)
        
    annot_a.append(row_annot_a)
    annot_b.append(row_annot_b)
    records.append(rec)

out_csv = "results/hub_genes_hpa_ihc_validation_9genes.csv"
pd.DataFrame(records).to_csv(out_csv, index=False)
print(f"Saved data table: {out_csv}")

# Create Dual-Panel Figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9.2), sharey=True, dpi=300)
fig.patch.set_facecolor("white")

# Colormap for Panel A
cmap_a = matplotlib.colors.ListedColormap(["#cbd5e1", "#f8fafc", "#93c5fd", "#2563eb", "#1e3a8a"])
bounds_a = [-1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
norm_a = matplotlib.colors.BoundaryNorm(bounds_a, cmap_a.N)

# Heatmap 1: Genuine HPA Normal IHC
im1 = ax1.imshow(mat_a, cmap=cmap_a, norm=norm_a, aspect="auto")
ax1.set_xticks(np.arange(len(organs)))
ax1.set_yticks(np.arange(len(genes)))
ax1.set_xticklabels(organs, fontsize=12, fontweight="bold", color="#0f172a")

yticklabels = [f"[{hpa_normal_data[g]['ab']}]  {g}" for g in genes]
ax1.set_yticklabels(yticklabels, fontsize=11, fontweight="bold", color="#0f172a")

for r in range(len(genes)):
    for c in range(len(organs)):
        val = mat_a[r, c]
        txt = annot_a[r][c]
        txt_color = "white" if val in [2, 3] else ("#475569" if val == -1 else "#0f172a")
        ax1.text(c, r, txt, ha="center", va="center", color=txt_color, fontsize=11, fontweight="bold")
        
ax1.set_title("A. Genuine HPA Normal-Tissue Protein IHC Baseline (HPA v23.0)\n[Experimental Monospecific Antibody Profiling in Healthy Tissues]", fontsize=12, fontweight="bold", pad=14, color="#0f172a")
ax1.grid(False)

for r in range(len(genes) + 1):
    ax1.axhline(r - 0.5, color="#94a3b8", linewidth=1.2)
for c in range(len(organs) + 1):
    ax1.axvline(c - 0.5, color="#94a3b8", linewidth=1.2)
    
legend_elements = [
    Patch(facecolor="#cbd5e1", edgecolor="#94a3b8", label="No IHC Data in HPA (SERPINE2)"),
    Patch(facecolor="#f8fafc", edgecolor="#94a3b8", label="0: Not Detected"),
    Patch(facecolor="#93c5fd", edgecolor="#94a3b8", label="1: Low Staining"),
    Patch(facecolor="#2563eb", edgecolor="#94a3b8", label="2: Medium Staining"),
    Patch(facecolor="#1e3a8a", edgecolor="#94a3b8", label="3: High Staining"),
]
ax1.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=True, fontsize=9.5)

# Heatmap 2: Transcript-Level Fold Change (RNA)
vmax = max(abs(mat_b.min()), abs(mat_b.max()), 2.5)
im2 = ax2.imshow(mat_b, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax2.set_xticks(np.arange(len(organs)))
ax2.set_xticklabels(organs, fontsize=12, fontweight="bold", color="#0f172a")

for r in range(len(genes)):
    for c in range(len(organs)):
        val = mat_b[r, c]
        txt = annot_b[r][c]
        txt_color = "white" if abs(val) > 1.2 else "#0f172a"
        ax2.text(c, r, txt, ha="center", va="center", color=txt_color, fontsize=11.5, fontweight="bold")
        
ax2.set_title("B. Transcript-Level Fold Change (RNA, N=1,069 Discovery Cohort)\n[Empirical Patient log2FC: Fibrosis vs. Non-Fibrotic Controls]", fontsize=12, fontweight="bold", pad=14, color="#0f172a")
ax2.grid(False)

for r in range(len(genes) + 1):
    ax2.axhline(r - 0.5, color="#94a3b8", linewidth=1.2)
for c in range(len(organs) + 1):
    ax2.axvline(c - 0.5, color="#94a3b8", linewidth=1.2)
    
cbar2 = fig.colorbar(im2, ax=ax2, orientation="horizontal", pad=0.15, shrink=0.75)
cbar2.set_label("Empirical Transcript-Level log2 Fold Change (RNA, Patient Biopsies)", fontsize=10.5, fontweight="bold", color="#0f172a")

fig.suptitle(
    "Multi-Organ Molecular Characterization of All 9 Pan-Fibrotic Consensus Hub Genes\n"
    "Panel A: Genuine HPA Normal-Tissue Protein IHC Baseline (proteinatlas.org)  |  Panel B: Transcript-Level Fold Change (RNA, N=1,069 Discovery Cohort)",
    fontsize=13, fontweight="bold", y=0.96, color="#0f172a"
)

plt.subplots_adjust(top=0.82, bottom=0.22, left=0.18, right=0.94, wspace=0.16)
out_png = "plots/hub_genes_hpa_ihc_summary_9genes.png"
plt.savefig(out_png, dpi=300)
plt.close()
print(f"Saved figure: {out_png}")
