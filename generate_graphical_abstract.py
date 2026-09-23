# -*- coding: utf-8 -*-
"""
Publication-Grade Graphical Abstract Generator for CSIR Pan-Fibrotic Study.
Renders high-resolution (300 DPI) multi-panel summary of the end-to-end study architecture,
exact gene numbers, 5-method ML consensus, dual-platform validation tiers, and downstream modules.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec

os.makedirs("graphical_abstract", exist_ok=True)
os.makedirs("plots", exist_ok=True)

# Set high-end font styling
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"

fig = plt.figure(figsize=(20, 11.25), dpi=300)
fig.patch.set_facecolor("#f8fafc")

gs = GridSpec(2, 3, figure=fig, hspace=0.25, wspace=0.22, 
              left=0.04, right=0.96, top=0.91, bottom=0.06)

# Palette
C_HEADER = "#0f172a"
C_BLUE = "#1e40af"
C_GREEN = "#059669"
C_AMBER = "#d97706"
C_RED = "#dc2626"
C_PURPLE = "#7c3aed"
C_SLATE = "#475569"
C_LIGHT_BG = "#ffffff"
C_BORDER = "#cbd5e1"

# Title Header
fig.text(0.5, 0.958, "Conserved Pan-Fibrotic Core Transcriptomic Program & Multi-Model Consensus ML Biomarkers",
         fontsize=17, fontweight="bold", ha="center", va="center", color=C_HEADER)
fig.text(0.5, 0.932, "Universal Extracellular Matrix Architecture Across Human Kidney, Liver, Lung, and Skin Biopsies (N = 1,069)",
         fontsize=11.5, fontstyle="italic", ha="center", va="center", color="#334155")

# -------------------------------------------------------------------------
# PANEL A: STAGE 1 DISCOVERY & MATRISOME FILTERING
# -------------------------------------------------------------------------
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_facecolor(C_LIGHT_BG)
for spine in ax_a.spines.values():
    spine.set_color(C_BORDER)
    spine.set_linewidth(1.2)
ax_a.set_xticks([])
ax_a.set_yticks([])

ax_a.text(0.04, 0.92, "A. Cross-Organ Discovery & Core Filtering", fontsize=12, fontweight="bold", color=C_HEADER)

# Flow boxes
boxes_a = [
    (0.06, 0.70, 0.88, 0.16, "#eff6ff", "#3b82f6", "14 Human Discovery Cohorts (N = 1,069 Biopsies)", 
     "Kidney (4) | Liver (3) | Lung (4) | Skin (3)\nModerated eBayes |log2FC| >= 0.585, FDR p < 0.05"),
    (0.12, 0.42, 0.76, 0.16, "#fef2f2", "#ef4444", "Conserved 4-Organ Core: 86 DEGs",
     "Simultaneous statistical significance across all 4\nanatomically distinct organ systems"),
    (0.18, 0.14, 0.64, 0.16, "#ecfdf5", "#10b981", "Matrisome Clean Core: 24 ECM Genes",
     "Collagens, Regulators, Secreted Factors, Glycoproteins\n(Human Matrisome Masterlist Overlap)")
]

for x, y, w, h, bg, bc, title, desc in boxes_a:
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015,rounding_size=0.03",
                                  facecolor=bg, edgecolor=bc, linewidth=1.4)
    ax_a.add_patch(rect)
    ax_a.text(x + w/2, y + h*0.68, title, fontsize=9.2, fontweight="bold", ha="center", va="center", color=C_HEADER)
    ax_a.text(x + w/2, y + h*0.28, desc, fontsize=7.6, ha="center", va="center", color="#334155", linespacing=1.2)

# Arrows
ax_a.annotate("", xy=(0.5, 0.58), xytext=(0.5, 0.70), arrowprops=dict(arrowstyle="-|>", color="#64748b", lw=2.0))
ax_a.annotate("", xy=(0.5, 0.30), xytext=(0.5, 0.42), arrowprops=dict(arrowstyle="-|>", color="#64748b", lw=2.0))


# -------------------------------------------------------------------------
# PANEL B: INDEPENDENT PRE-ML REPLICATION
# -------------------------------------------------------------------------
ax_b = fig.add_subplot(gs[0, 1])
ax_b.set_facecolor(C_LIGHT_BG)
for spine in ax_b.spines.values():
    spine.set_color(C_BORDER)
    spine.set_linewidth(1.2)
ax_b.set_xticks([])
ax_b.set_yticks([])

ax_b.text(0.04, 0.92, "B. Pre-ML Multi-Center Replication", fontsize=12, fontweight="bold", color=C_HEADER)

boxes_b = [
    (0.05, 0.55, 0.90, 0.30, "#faf5ff", "#a855f7", "Validation 1: Internal Replication (N = 4, n = 384)",
     "• 24 / 24 Genes (100%) Concordant Direction in All 4 Organs\n• GSE200818 (Kidney) · GSE162694 (Liver)\n• GSE24206 (Lung) · GSE58095 (Skin)"),
    (0.05, 0.14, 0.90, 0.30, "#fff7ed", "#f97316", "Validation 2: Held-Out Replication (N = 4, n = 105)",
     "• 18 / 24 Genes Replicated in >= 2/4 Organs (Mann-Whitney U)\n• GSE30529 (Kidney) · GSE14323 (Liver)\n• GSE83717 (Lung RNA-seq) · GSE125362 (Skin Microarray)")
]

for x, y, w, h, bg, bc, title, desc in boxes_b:
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015,rounding_size=0.03",
                                  facecolor=bg, edgecolor=bc, linewidth=1.4)
    ax_b.add_patch(rect)
    ax_b.text(x + 0.03, y + h*0.78, title, fontsize=9.2, fontweight="bold", ha="left", va="center", color=C_HEADER)
    ax_b.text(x + 0.03, y + h*0.36, desc, fontsize=7.8, ha="left", va="center", color="#334155", linespacing=1.3)

ax_b.annotate("", xy=(0.5, 0.44), xytext=(0.5, 0.55), arrowprops=dict(arrowstyle="-|>", color="#64748b", lw=2.0))


# -------------------------------------------------------------------------
# PANEL C: 5-METHOD CONSENSUS ML & WGCNA
# -------------------------------------------------------------------------
ax_c = fig.add_subplot(gs[0, 2])
ax_c.set_facecolor(C_LIGHT_BG)
for spine in ax_c.spines.values():
    spine.set_color(C_BORDER)
    spine.set_linewidth(1.2)
ax_c.set_xticks([])
ax_c.set_yticks([])

ax_c.text(0.04, 0.92, "C. Consensus ML & WGCNA Discovery", fontsize=12, fontweight="bold", color=C_HEADER)

# ML methods subgrid
ml_box = patches.FancyBboxPatch((0.05, 0.52), 0.90, 0.33, boxstyle="round,pad=0.015,rounding_size=0.03",
                                facecolor="#f8fafc", edgecolor="#94a3b8", linewidth=1.4)
ax_c.add_patch(ml_box)
ax_c.text(0.5, 0.78, "5 Complementary Selection Architectures (5 Seeds)", fontsize=8.6, fontweight="bold", ha="center", va="center", color=C_HEADER)
ax_c.text(0.5, 0.63, "1. LASSO (L1 Regularization)     2. SVM-RFE (Recursive Elimination)\n3. Random Forest (Gini/Perm)   4. XGBoost (Gradient Boosting)\n5. TOM Co-expression Clustering / WGCNA (Soft-thresholding)",
          fontsize=7.6, ha="center", va="center", color="#334155", linespacing=1.35)

# Consensus output box
out_ml = patches.FancyBboxPatch((0.05, 0.12), 0.90, 0.33, boxstyle="round,pad=0.015,rounding_size=0.03",
                                facecolor="#eff6ff", edgecolor="#2563eb", linewidth=1.6)
ax_c.add_patch(out_ml)
ax_c.text(0.5, 0.35, "9 Consensus Hub Genes (>= 3/5 Votes)", fontsize=9.5, fontweight="bold", ha="center", va="center", color="#1e40af")
ax_c.text(0.5, 0.21, "5 Unanimous Hubs (5/5 Votes):\nCOL15A1 · COL1A1 · COL3A1 · SERPINE2 · SERPINF2\n+ 4 Extended Hubs: LAMC3 · LTBP2 · MDK · SVEP1",
          fontsize=7.8, ha="center", va="center", color="#0f172a", linespacing=1.2)


# -------------------------------------------------------------------------
# PANEL D: LOCKED 3-TIER BIOMARKER HIERARCHY
# -------------------------------------------------------------------------
ax_d = fig.add_subplot(gs[1, :2])
ax_d.set_facecolor(C_LIGHT_BG)
for spine in ax_d.spines.values():
    spine.set_color(C_BORDER)
    spine.set_linewidth(1.2)
ax_d.set_xticks([])
ax_d.set_yticks([])

ax_d.text(0.02, 0.93, "D. Stage 3 Dual-Platform Validation & Locked 3-Tier Hub Hierarchy", fontsize=12, fontweight="bold", color=C_HEADER)

# Tier 1 Card
t1 = patches.FancyBboxPatch((0.03, 0.12), 0.29, 0.74, boxstyle="round,pad=0.015,rounding_size=0.03",
                            facecolor="#eff6ff", edgecolor="#1d4ed8", linewidth=1.8)
ax_d.add_patch(t1)
ax_d.text(0.175, 0.78, "3 Strict Validated Hubs", fontsize=10.2, fontweight="bold", ha="center", va="center", color="#1e40af")
ax_d.text(0.175, 0.70, "(p < 0.05 on Microarray & RNA-seq)", fontsize=7.4, fontstyle="italic", ha="center", va="center", color="#3b82f6")
ax_d.text(0.175, 0.44, "• COL15A1 (AUC: 0.843)\n  Basement membrane organizer\n• COL3A1 (AUC: 0.852)\n  Early fibrillar repair collagen\n• SERPINE2 (AUC: 0.786)\n  Antiprotease/matrix protection",
          fontsize=8.2, ha="center", va="center", color="#0f172a", linespacing=1.35)

# Tier 2 Card
t2 = patches.FancyBboxPatch((0.355, 0.12), 0.29, 0.74, boxstyle="round,pad=0.015,rounding_size=0.03",
                            facecolor="#ecfdf5", edgecolor="#059669", linewidth=1.8)
ax_d.add_patch(t2)
ax_d.text(0.50, 0.78, "2 Direction-Concordant Hubs", fontsize=10.2, fontweight="bold", ha="center", va="center", color="#047857")
ax_d.text(0.50, 0.70, "(5/5 ML Votes; RNA-seq p >= 0.05)", fontsize=7.4, fontstyle="italic", ha="center", va="center", color="#059669")
ax_d.text(0.50, 0.44, "• COL1A1 (AUC: 0.790)\n  Canonical fibrillar collagen\n  RNA-seq p = 0.101 (Concordant Up)\n• SERPINF2 (AUC: 0.764)\n  Alpha-2-antiplasmin inhibitor\n  RNA-seq p = 0.788 (Concordant Down)",
          fontsize=8.2, ha="center", va="center", color="#0f172a", linespacing=1.35)

# Tier 3 Card
t3 = patches.FancyBboxPatch((0.68, 0.12), 0.29, 0.74, boxstyle="round,pad=0.015,rounding_size=0.03",
                            facecolor="#f8fafc", edgecolor="#64748b", linewidth=1.4)
ax_d.add_patch(t3)
ax_d.text(0.825, 0.78, "4 Stage-2 Candidates", fontsize=10.2, fontweight="bold", ha="center", va="center", color="#475569")
ax_d.text(0.825, 0.70, "(Failed Stage 3 Dual-Platform)", fontsize=7.4, fontstyle="italic", ha="center", va="center", color="#64748b")
ax_d.text(0.825, 0.44, "• LAMC3 (Discordant RNA-seq)\n• LTBP2 (Discordant RNA-seq)\n• SVEP1 (Discordant RNA-seq)\n• MDK (RNA-seq p = 0.760)\n\nRetained as pilot candidates",
          fontsize=8.2, ha="center", va="center", color="#334155", linespacing=1.35)


# -------------------------------------------------------------------------
# PANEL E: CLINICAL TRANSLATION & VALIDATION MODULES
# -------------------------------------------------------------------------
ax_e = fig.add_subplot(gs[1, 2])
ax_e.set_facecolor(C_LIGHT_BG)
for spine in ax_e.spines.values():
    spine.set_color(C_BORDER)
    spine.set_linewidth(1.2)
ax_e.set_xticks([])
ax_e.set_yticks([])

ax_e.text(0.04, 0.93, "E. Downstream Validation Modules", fontsize=12, fontweight="bold", color=C_HEADER)

modules_text = (
    "Mod 1: PPI Network (STRING v12)\n"
    "  • High-confidence Collagen Triad (COL1A1-COL3A1-COL15A1)\n\n"
    "Mod 2: Early-Detection Subgroup ROC (GSE84044)\n"
    "  • S0 vs S1-S2 Subclinical AUC = 0.716 (5-Hub) / 0.719 (9-Hub)\n\n"
    "Mod 3: Diagnostic Nomogram & DCA\n"
    "  • Multivariable AUC = 0.9002 (5-Hub) / 0.9238 (9-Hub)\n\n"
    "Mod 4: Immune & Stromal Infiltration\n"
    "  • ACTA2 (Myofibroblast), POSTN (Matrix), PECAM1 (Endothelial)\n\n"
    "Mod 5: HPA Baseline IHC & Biopsy RNA\n"
    "  • Normal protein localization vs. patient transcript log2FC\n\n"
    "Mod 6: GO / KEGG / DO Enrichment (Enrichr)\n"
    "  • Extracellular matrix & collagen fibril organization"
)

mod_box = patches.FancyBboxPatch((0.04, 0.06), 0.92, 0.82, boxstyle="round,pad=0.015,rounding_size=0.03",
                                 facecolor="#faf5ff", edgecolor="#c084fc", linewidth=1.4)
ax_e.add_patch(mod_box)
ax_e.text(0.07, 0.46, modules_text, fontsize=7.6, ha="left", va="center", color="#1e1b4b", linespacing=1.2)

# Save Outputs
out_fig = "graphical_abstract/graphical_abstract_publication_plate.png"
out_pdf = "graphical_abstract/graphical_abstract_publication_plate.pdf"
plt.savefig(out_fig, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()
print(f"Successfully generated: {out_fig} and {out_pdf}")
