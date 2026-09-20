#!/usr/bin/env python3
"""
Venn Diagram Visualizations
========================
Creates a multi-panel figure showing overlap between:
  A) Each organ's significant genes  vs  ECM reference genes (2-way Venns, 2x2 grid)
  B) ECM-subset Venn: Kidney-ECM ∩ Liver-ECM ∩ Lung-ECM  (3-way Venn, with Skin-ECM annotated)
  C) Pan-Fibrotic Core Genes (175) vs ECM Reference (1027)
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

plt.rcParams.update({"font.family": "Arial"})

# ---------- Paths ----------
BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "results")
REF = os.path.join(BASE, "reference")
ECM_XLSX = os.path.join(BASE, "ECM genes all.xlsx")
if not os.path.exists(ECM_XLSX):
    ECM_XLSX = os.path.join(REF, "ECM_genes_all.xlsx")

ADJ_P_THR = 0.05
LFC_THR = 0.585

# ---------- Load organ DEGs ----------
organs = {}
for org, fname in [("Kidney", "kidney_DEGs.csv"),
                  ("Liver",  "liver_DEGs.csv"),
                  ("Lung",   "lung_DEGs.csv"),
                  ("Skin",   "skin_DEGs.csv")]:
    df = pd.read_csv(os.path.join(RESULTS, fname))
    sig = df[(df["adj_p_value"] < ADJ_P_THR) & (df["logFC"].abs() > LFC_THR)].copy()
    sig["gene"] = sig["gene"].astype(str).str.strip().str.upper()
    organs[org] = set(sig["gene"].values)

# ---------- Load ECM masterlist ----------
ecm_df = pd.read_excel(ECM_XLSX, sheet_name="Hs_ECM_Masterlist", header=1, engine="openpyxl")
ecm_df["Gene Symbol"] = ecm_df["Gene Symbol"].astype(str).str.strip().str.upper()
ecm_set = set(ecm_df["Gene Symbol"].dropna().values)
ecm_set.discard("")
ecm_set.discard("NAN")

# For each organ, compute the organ ∩ ECM subset
org_ecm = {o: (organs[o] & ecm_set) for o in organs}

# Pan-fibrotic core (4-way organ intersection)
core_set = organs["Kidney"] & organs["Liver"] & organs["Lung"] & organs["Skin"]
core_ecm = core_set & ecm_set

# =============================================================================
# FIGURE – Multi-panel Venns
# =============================================================================
fig = plt.figure(figsize=(20, 14))
gs = fig.add_gridspec(3, 4, hspace=0.55, wspace=0.35)

# ------- ROW 1+2: 2x2 grid of 2-way Venns (Organ sig genes vs ECM) --------
venn_colors_organ = {"Kidney": "#1f77b4", "Liver": "#2ca02c", "Lung": "#ff7f0e", "Skin": "#d62728"}
positions = [(0, 0), (0, 2), (1, 0), (1, 2)]
grid_spec_cells = [(slice(0, 1), slice(0, 2)),
                  (slice(0, 1), slice(2, 4)),
                  (slice(1, 2), slice(0, 2)),
                  (slice(1, 2), slice(2, 4))]

for org_name, (row_s, col_s) in zip(organs.keys(), grid_spec_cells):
    ax = fig.add_subplot(gs[row_s, col_s])
    organ_genes = organs[org_name]
    o_and_e = organ_genes & ecm_set
    o_only = len(organ_genes - ecm_set)
    e_only = len(ecm_set - organ_genes)
    both = len(o_and_e)

    v = venn2(subsets=(o_only, e_only, both),
               set_labels=(f"{org_name}\nSig Genes\n(n={len(organ_genes):,})",
                           f"ECM\nMasterlist\n(n={len(ecm_set):,})"),
               set_colors=(venn_colors_organ[org_name], "#9467bd"),
               ax=ax, alpha=0.75)

    if v.get_label_by_id("11") is not None:
        v.get_label_by_id("11").set_text(f"{both}\n({both/len(organ_genes)*100:.1f}%\nof {org_name}")
        v.get_label_by_id("11").set_fontweight("bold")
    if v.get_label_by_id("10") is not None:
        v.get_label_by_id("10").set_text(f"{o_only}\n({o_only/len(organ_genes)*100:.1f}%)")
    if v.get_label_by_id("01") is not None:
        v.get_label_by_id("01").set_text(f"{e_only}\nECM-only")

    for text in v.set_labels:
        if text is not None:
            text.set_fontsize(10)
            text.set_fontweight("bold")

    # Highlight ECM-core genes (shared across all 4 organs that are also in ECM)
    # within each organ's overlap
    core_in_slice = len(core_ecm & organ_genes)
    ax.set_title(f"{org_name} Significant Genes  ×  ECM Masterlist\n"
                 f"  [of {both} ECM-overlapping, {core_in_slice} are in 4-org core]",
                 fontsize=11, fontweight="bold", pad=12)

# ------- ROW 3 LEFT: 3-way Venn of ECM-subset genes (Kidney ∩ Liver ∩ Lung) -------
ax3 = fig.add_subplot(gs[2, 0:2])
k_e = org_ecm["Kidney"]
li_e = org_ecm["Liver"]
lu_e = org_ecm["Lung"]
s_e = org_ecm["Skin"]

v3 = venn3([k_e, li_e, lu_e],
           set_labels=(f"Kidney-ECM\n(n={len(k_e)})",
                       f"Liver-ECM\n(n={len(li_e)})",
                       f"Lung-ECM\n(n={len(lu_e)})"),
           set_colors=("#1f77b4", "#2ca02c", "#ff7f0e"),
           ax=ax3, alpha=0.75)

for text in v3.set_labels:
    if text is not None:
        text.set_fontsize(9)
        text.set_fontweight("bold")
for txt_id in ["111", "110", "101", "011", "100", "010", "001"]:
    lab = v3.get_label_by_id(txt_id)
    if lab is not None:
        lab.set_fontsize(9)

# The 3-way intersection PLUS skin == 4-way ECM-core (core_ecm set)
three_way_intersect = k_e & li_e & lu_e
in_3_plus_skin = three_way_intersect & s_e
only_3_not_4 = three_way_intersect - s_e

ax3.set_title("ECM-Overlapping Genes: Kidney ∩ Liver ∩ Lung\n"
              f"(+ Skin-ECM n={len(s_e)}; 4-Org ECM-Core = {len(core_ecm)})",
              fontsize=10, fontweight="bold", pad=10)
# Annotate Skin count manually.
# Add Skin info callout
ax3.annotate(f"  Skin-ECM (not shown):\n"
            f"  • total = {len(s_e)} ECM sig genes\n"
            f"  • 3-way (K∩Li∩Lu) = {len(three_way_intersect)}\n"
            f"  • +Skin → 4-org core = {len(core_ecm)}",
            xy=(0.98, 0.02), xycoords="axes fraction",
            fontsize=8.5, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f7e1e4", ec="#d62728", lw=1.2),
            color="#7a0e0e", fontweight="bold")

# ------- ROW 3 RIGHT: 2-way Venn Pan-Fibrotic Core Genes vs ECM -------
ax_core = fig.add_subplot(gs[2, 2:4])
core_only = len(core_set - ecm_set)
ecm_only_vs_core = len(ecm_set - core_set)
both_core_ecm = len(core_ecm)

vc = venn2(subsets=(core_only, ecm_only_vs_core, both_core_ecm),
           set_labels=(f"Pan-Fibrotic Core\n(4-Org Shared)\n(n={len(core_set)})",
                       f"ECM\nMasterlist\n(n={len(ecm_set):,})"),
           set_colors=("#e377c2", "#9467bd"),
           ax=ax_core, alpha=0.80)
if vc.get_label_by_id("11") is not None:
    vc.get_label_by_id("11").set_text(f"{both_core_ecm}\n({both_core_ecm/len(core_set)*100:.1f}%\nof core)")
    vc.get_label_by_id("11").set_fontweight("bold")
    vc.get_label_by_id("11").set_color("#4a0e4e")
if vc.get_label_by_id("10") is not None:
    vc.get_label_by_id("10").set_text(f"{core_only}\nnon-ECM\ncore genes")
if vc.get_label_by_id("01") is not None:
    vc.get_label_by_id("01").set_text(f"{ecm_only_vs_core:,}\nECM-only\n(no sig in all 4)")

for text in vc.set_labels:
    if text is not None:
        text.set_fontsize(10)
        text.set_fontweight("bold")

# List top ECM core gene examples inside plot area below inlaid category breakdown:
cat_counts = {}
for g in core_ecm:
    row_match = ecm_df[ecm_df["Gene Symbol"] == g]["Matrisome Category"]
    if not row_match.empty:
        cat = str(row_match.values[0])
        cat_counts[cat] = cat_counts.get(cat, 0) + 1

cat_str = "\n".join([f"  • {c}: {n}" for c, n in sorted(cat_counts.items(), key=lambda x:-x[1])])
ax_core.set_title("Pan-Fibrotic Core Genes  ×  ECM Masterlist\n"
                 f"({both_core_ecm}/{len(core_set)} = {both_core_ecm/len(core_set)*100:.1f}% ECM-enriched)",
                 fontsize=11, fontweight="bold", pad=12)
ax_core.annotate(f"Matrisome Categories\nof {both_core_ecm} core-ECM genes:\n{cat_str}",
               xy=(0.02, 0.02), xycoords="axes fraction",
               fontsize=8, ha="left", va="bottom",
               bbox=dict(boxstyle="round,pad=0.4", fc="#f0f5ff", ec="#5e4e8c", lw=1),
               color="#2c1e5e", fontweight="bold")

# =============================================================================
# MAIN FIGURE TITLE
# =============================================================================
fig.suptitle("Cross-Organ Fibrosis DEGs × ECM Matrisome Overlap – Venn Compendium",
             fontsize=17, fontweight="bold", y=0.992,
             color="#1a1a4a")
fig.text(0.5, 0.965,
         f"Significance filter: adj_p_value < 0.05 AND |logFC| > 0.585"
         f"   |   ECM reference: Hs_ECM_Masterlist (N={len(ecm_set):,})"
         f"   |   4-org pan-fibrotic core: {len(core_set)} genes ({len(core_ecm)} ECM)",
         ha="center", fontsize=10.5, style="italic", color="#444444")

out_png = os.path.join(BASE, "venn_organ_vs_ecm_compendium.png")
plt.savefig(out_png, dpi=220, bbox_inches="tight", facecolor="white")
print(f"[SAVED] Multi-panel Venn figure: {out_png}")

plots_dir = os.path.join(BASE, "plots")
os.makedirs(plots_dir, exist_ok=True)
out_png_plots = os.path.join(plots_dir, "venn_organ_vs_ecm_compendium.png")
plt.savefig(out_png_plots, dpi=220, bbox_inches="tight", facecolor="white")
print(f"[SAVED] Multi-panel Venn figure: {out_png_plots}")
plt.close(fig)

# =============================================================================
# PRINT SUMMARY TO CONSOLE
# =============================================================================
print("\n" + "="*72)
print("VENN OVERLAP SUMMARY")
print("="*72)
print(f"\n{'Organ':<10} {'Sig Genes':>10} {'int ECM':>8} {'% of Organ':>11} {'In 4-Core int ECM':>20}")
print("-"*72)
for o in organs:
    n_sig = len(organs[o])
    n_oe = len(org_ecm[o])
    pct = n_oe / n_sig * 100 if n_sig > 0 else 0
    n_core_here = len(org_ecm[o] & core_set)
    print(f"{o:<10} {n_sig:>10,} {n_oe:>8,} {pct:>10.1f}% {n_core_here:>20}")

print(f"\nPan-Fibrotic 4-org core: {len(core_set)} genes")
print(f"  of which ECM:         {len(core_ecm)} ({len(core_ecm)/len(core_set)*100:.1f}%)")

three = len(org_ecm["Kidney"] & org_ecm["Liver"] & org_ecm["Lung"])
four  = len(org_ecm["Kidney"] & org_ecm["Liver"] & org_ecm["Lung"] & org_ecm["Skin"])
print(f"\nECM-subset overlaps:")
print(f"  Kidney int Liver int Lung (3-way ECM): {three} genes")
print(f"  +Skin -> 4-way ECM core:              {four} genes (= {four/max(1,three)*100:.1f}% of 3-way)")
print("="*72)
