# -*- coding: utf-8 -*-
"""
Generate Final High-Resolution Publication Figures for CSIR Pan-Fibrotic Study
- plots/study_design_funnel_corrected.png
- plots/ml_4model_consensus_hub_biomarkers.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

BASE_DIR = r"D:\CSIR"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
PLOTS_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# Load Master Evidence Table
master_csv = os.path.join(RESULTS_DIR, "final_master_evidence_table.csv")
df_master = pd.read_csv(master_csv)

# =============================================================================
# FIGURE 1: ML 4-Model Consensus Hub Biomarkers & Master Evidence Ranking
# =============================================================================
plt.figure(figsize=(15, 10))

# For plotting, sort by ML diagnostic AUC (fill NaN for TNXB as 0.50 with special annotation)
plot_df = df_master.copy()
plot_df["auc_plot"] = pd.to_numeric(plot_df["ml_diagnostic_auc"], errors="coerce").fillna(0.50)
plot_df = plot_df.sort_values(by="auc_plot", ascending=True).reset_index(drop=True)

tier_colors = {
    "Tier 1 (Full Spectrum)": "#1b4965",      # Deep Navy
    "Tier 2 (Validation-Only)": "#2a9d8f",    # Sea Green / Teal
    "Not Supported": "#adb5bd"                # Slate Grey
}

colors = [tier_colors.get(t, "#6c757d") for t in plot_df["final_evidence_tier"]]
y_pos = np.arange(len(plot_df))
bars = plt.barh(y_pos, plot_df["auc_plot"], color=colors, edgecolor="black", linewidth=0.8, height=0.68)

plt.yticks(y_pos, plot_df["gene"], fontsize=11, fontweight="bold")
plt.xlabel("Mean Diagnostic ROC AUC (5-Seed Average on Pure Discovery ComBat-Corrected Human Cohorts, N=1,069)", 
           fontsize=12, fontweight="bold", labelpad=10)
plt.title("Master Evidence Ranking: 24 Clean Core ECM Genes + TNXB Across All Validation Layers", 
          fontsize=15, fontweight="bold", pad=20)
plt.xlim(0.48, 0.92)
plt.axvline(x=0.70, color="#e63946", linestyle="--", linewidth=1.5, label="High Diagnostic Performance Threshold (AUC >= 0.70)")

# Value annotations on bars
for idx, (bar, row) in enumerate(zip(bars, plot_df.iterrows())):
    r = row[1]
    w = bar.get_width()
    gene = r["gene"]
    tier = r["final_evidence_tier"]
    stab = r["ml_stability_score"]
    
    if gene == "TNXB":
        label_text = "Unassessed in ML Matrix (4/4 Val2 Concordant; Reverse-MR Causal p=5.4e-7)"
        plt.text(0.505, bar.get_y() + bar.get_height()/2.0, label_text, 
                 va="center", ha="left", fontsize=9, fontweight="bold", color="#264653", style="italic")
    else:
        auc_val = r["ml_diagnostic_auc"]
        label_text = f"AUC = {auc_val:.3f} | Stability: {stab} | {tier}"
        plt.text(w + 0.005, bar.get_y() + bar.get_height()/2.0, label_text, 
                 va="center", ha="left", fontsize=9, fontweight="bold", color="#1d3557")

# Custom Legend
legend_elements = [
    patches.Patch(facecolor=tier_colors["Tier 1 (Full Spectrum)"], edgecolor="black", label="Tier 1 (Full Spectrum: Val2 >= 2/4 & ML Stability >= 4/5) [n=4]"),
    patches.Patch(facecolor=tier_colors["Tier 2 (Validation-Only)"], edgecolor="black", label="Tier 2 (Validation-Only: Val2 >= 2/4 Concordant Multi-Organ) [n=15]"),
    patches.Patch(facecolor=tier_colors["Not Supported"], edgecolor="black", label="Not Supported (Failed Multi-Organ Val2 Replication, <= 1/4) [n=6]"),
    plt.Line2D([0], [0], color="#e63946", linestyle="--", linewidth=1.5, label="Diagnostic Benchmark (AUC = 0.70)")
]
plt.legend(handles=legend_elements, loc="lower right", frameon=True, facecolor="white", edgecolor="#ced4da", fontsize=10)
plt.grid(axis="x", linestyle=":", alpha=0.6)
plt.tight_layout()

out_ml_png = os.path.join(PLOTS_DIR, "ml_4model_consensus_hub_biomarkers.png")
plt.savefig(out_ml_png, dpi=300)
plt.close()
print(f"Generated: {out_ml_png}")

# =============================================================================
# FIGURE 2: Study Design Funnel Diagram (Locked 3-Tier Multi-Organ Discovery)
# =============================================================================
fig, ax = plt.subplots(figsize=(14, 10))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# Title
ax.text(50, 96, "Locked Multi-Organ Discovery & Multi-Layer Validation Study Funnel", 
        fontsize=17, fontweight="bold", ha="center", va="center", color="#0f172a")
ax.text(50, 92.5, "Kidney, Liver, Lung, and Skin Fibrosis | Zero Data Leakage Enforced Across All Tiers", 
        fontsize=11, ha="center", va="center", color="#475569", style="italic")

# Funnel levels data
levels = [
    {
        "y": 80, "width": 88, "height": 8.5, "color": "#1e3a8a", "text_color": "white",
        "title": "TIER 1: 4-ORGAN DISCOVERY COHORTS (14 Datasets, N=1,069 Pure Human Biopsies)",
        "details": "Kidney (GSE66494, GSE104066, GSE104948, GSE104954) | Liver (GSE164760, GSE89377, GSE77627)\nLungs (GSE10667, GSE110147, GSE32537, GSE53845) | Skin (GSE130955, GSE181549, GSE95065)\nPer-Organ Significant DEGs: Kidney (12,442) | Liver (13,096) | Lungs (10,778) | Skin (3,079)"
    },
    {
        "y": 67, "width": 74, "height": 7.5, "color": "#2563eb", "text_color": "white",
        "title": "CONSERVED PAN-FIBROTIC CORE DEGs (N = 86 Genes)",
        "details": "Conserved across 4/4 Organs at |log2FC| >= 1.0 (>= 2-fold change) and Benjamini-Hochberg FDR p < 0.05"
    },
    {
        "y": 55, "width": 62, "height": 7.5, "color": "#0284c7", "text_color": "white",
        "title": "MATRISOME / ECM FILTERING & TARGET INCLUSION (N = 25 Candidate Genes)",
        "details": "24 Clean Core Human Extracellular Matrix (Matrisome) Genes + 1 High-Interest Target (TNXB)\nCategories: Collagens, ECM Glycoproteins, ECM Regulators, Secreted Factors, ECM-affiliated"
    },
    {
        "y": 43, "width": 52, "height": 7.5, "color": "#0d9488", "text_color": "white",
        "title": "TIER 2: VALIDATION LAYER 1 INDEPENDENT REPLICATION (4 Cohorts)",
        "details": "GSE200818 (Kidney), GSE162694 (Liver), GSE24206 (Lungs), GSE58095 (Skin)\n100% Direction Concordance (24/24 Genes) | 12/24 Statistically Significant in >= 2/4 Organs"
    },
    {
        "y": 31, "width": 44, "height": 7.5, "color": "#059669", "text_color": "white",
        "title": "TIER 3: VALIDATION LAYER 2 HELD-OUT REPLICATION (4 Multi-Platform Cohorts)",
        "details": "GSE30529 (Kidney), GSE14323 (Liver), GSE83717 (Lungs), GSE125362 (Skin)\n19 Genes Validated in >= 2/4 Organs (18 ECM + TNXB) | 6 Dropped as Not Supported (<= 1/4 Organs)"
    },
    {
        "y": 16, "width": 38, "height": 11.5, "color": "#0f766e", "text_color": "white",
        "title": "FINAL MASTER EVIDENCE STRATIFICATION (4-Model ML & Causal Multi-Omics)",
        "details": "• TIER 1 (Full Spectrum, n=4): COL15A1 (AUC 0.843), COL1A1 (0.790), SERPINE2 (0.786), SERPINF2 (0.764)\n• TIER 2 (Validation-Only, n=15): COL3A1, COL1A2, LTBP2, LAMC3, CLEC2D, PDGFD, CCL2,\n   SERPINH1, AEBP1, VWF, CCL5, CCL19, SVEP1, MFAP4, TNXB (Causal SSc MR p=5.4e-7)\n• NOT SUPPORTED (n=6): MDK, FGF14, SPARCL1, BMP1, CCL21, COLEC11 (Failed Val2 Replication)"
    }
]

for idx, lvl in enumerate(levels):
    w, h, y = lvl["width"], lvl["height"], lvl["y"]
    x = (100 - w) / 2
    
    # Rounded box
    rect = patches.FancyBboxPatch((x, y - h/2), w, h, boxstyle="round,pad=0.8,rounding_size=1.5",
                                  facecolor=lvl["color"], edgecolor="#0f172a", linewidth=1.5, zorder=2)
    ax.add_patch(rect)
    
    # Text
    ax.text(50, y + h*0.22, lvl["title"], fontsize=10.5, fontweight="bold", ha="center", va="center", color=lvl["text_color"], zorder=3)
    ax.text(50, y - h*0.20, lvl["details"], fontsize=8.5, ha="center", va="center", color="#f8fafc", zorder=3)
    
    # Connective arrow down
    if idx < len(levels) - 1:
        next_lvl = levels[idx + 1]
        ax.annotate("", xy=(50, next_lvl["y"] + next_lvl["height"]/2), xytext=(50, y - h/2),
                    arrowprops=dict(arrowstyle="-|>", color="#334155", lw=2, mutation_scale=15), zorder=1)

plt.tight_layout()
out_funnel_png = os.path.join(PLOTS_DIR, "study_design_funnel_corrected.png")
plt.savefig(out_funnel_png, dpi=300)
plt.close()
print(f"Generated: {out_funnel_png}")
