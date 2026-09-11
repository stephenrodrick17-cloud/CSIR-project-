#!/usr/bin/env python3
"""
4-Organ Venn Diagram
====================
Draws a 4-set Venn diagram with four ellipses (the standard way to show
all 16 pairwise intersections between 4 gene sets, since true circles
cannot show all 4-set overlaps geometrically).

Shows the overlap of SIGNIFICANT genes between:
    Kidney  •  Liver  •  Lung  •  Skin

Genes are considered significant:  adj_p_value < 0.05  AND  |logFC| > 0.585

Uses the "venn" library (pip install venn) which supports venn4.
If that package is unavailable, falls back to matplotlib_venn-style venn3
panels for pairwise/three-way overlaps so the script still works.
"""

import os
import itertools
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

# ---------- Paths & thresholds ----------
BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "results")
ADJ_P_THR = 0.05
LFC_THR = 0.585

# ---------- Load organ significant genes ----------
print("Loading organ DEG results and filtering for significant genes ...")
organs_list = ["Kidney", "Liver", "Lung", "Skin"]
organ_colors = {
    "Kidney": "#1f77b4",
    "Liver":  "#2ca02c",
    "Lung":   "#ff7f0e",
    "Skin":   "#d62728",
}
organ_csv = {
    "Kidney": "kidney_DEGs.csv",
    "Liver":  "liver_DEGs.csv",
    "Lung":   "lung_DEGs.csv",
    "Skin":   "skin_DEGs.csv",
}

gene_sets = {}
for org in organs_list:
    df = pd.read_csv(os.path.join(RESULTS, organ_csv[org]))
    sig = df[(df["adj_p_value"] < ADJ_P_THR) & (df["logFC"].abs() > LFC_THR)].copy()
    sig["gene"] = sig["gene"].astype(str).str.strip().str.upper()
    gene_sets[org] = set(sig["gene"].values)
    print(f"  {org:<7}  {len(gene_sets[org]):>5,} significant genes")

# ---------- Compute ALL 15 intersection regions ----------
def compute_venn_counts(sets_dict, order):
    """
    For the four sets in `order`, compute the exact count of each
    Venn region (exclusive membership in that combination). Returns a
    dict keyed by a 4-bit mask string 'ABCD' (0 or 1 per set).
    """
    masks = {}
    all_genes = set().union(*[sets_dict[o] for o in order])
    for g in all_genes:
        bit = "".join(["1" if g in sets_dict[o] else "0" for o in order])
        masks[bit] = masks.get(bit, 0) + 1
    return masks

order = organs_list  # Kidney, Liver, Lung, Skin
counts = compute_venn_counts(gene_sets, order)

# Pretty names for non-empty subsets for printing
print("\nAll Venn regions (4-org):")
regions_report = []
for bits, cnt in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
    present = [order[i] for i, b in enumerate(bits) if b == "1"]
    absent = [order[i] for i, b in enumerate(bits) if b == "0"]
    label = "+".join(present)
    if absent:
        label += f"  (not in {','.join(absent)})"
    regions_report.append((bits, cnt, len(present), label))
    if cnt > 0:
        print(f"  {bits}  n={cnt:>5}   {label}")

# ---------- VENN DIAGRAM 1: Use `venn` library venn4 (4 ellipses) ----------
venn4_ok = False
try:
    from venn import venn
    # `venn` library venn4 function takes labels in the set order
    # We pass a dict of sets -> it auto-computes everything
    print("\nUsing `venn` library for 4-way ellipse Venn...")

    fig, ax = plt.subplots(figsize=(14, 12))
    # Prepare sets dict in desired order
    venn_sets = {org: gene_sets[org] for org in order}
    venn(venn_sets, ax=ax,
         colors=[organ_colors[o] for o in order],
         alpha=0.45,
         fmt="{size:,}",
         legend_loc="best",
         fontsize=10,
         figsize=(14,12))
    ax.set_title("4-Organ Venn Diagram — Significant DEG Overlap\n"
                 "(adj_p < 0.05  AND  |logFC| > 0.585)",
                 fontsize=17, fontweight="bold", pad=22)

    # Annotate the 4-way intersection
    four_way = counts.get("1111", 0)
    ax.annotate(f"ALL 4 ORGANS\nSHARED\nn = {four_way:,}",
                xy=(0.5, 0.5), xycoords="axes fraction",
                ha="center", va="center", fontsize=13,
                fontweight="bold", color="#1a1a1a",
                bbox=dict(boxstyle="round,pad=0.7", fc="#fff8dc",
                          ec="#b8860b", lw=2.5))

    out1 = os.path.join(BASE, "venn_4organ_ellipses.png")
    fig.savefig(out1, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[SAVED] 4-way ellipse Venn: {out1}")
    venn4_ok = True
except Exception as e:
    print(f"`venn` library not usable: {e}. Switching to manual 4-ellipse Venn.")

# ---------- FALLBACK / ADDITIONAL: Manual 4-ellipse Venn with labels ----------
print("\nBuilding manual 4-ellipse Venn (alternative layout with full legend) ...")

def draw_4ellipse_venn(sets, order, colors, ax):
    """
    Custom 4-ellipse Venn layout — 4 rotated ellipses centered roughly around
    a common center, with counts annotated in each of the 16 regions.
    """
    # Draw 4 overlapping ellipses at different rotations
    cx, cy = 0.5, 0.5
    w, h = 0.55, 0.30
    rotations = [45, -45, 135, -135]   # 4 angles -> all 16 regions formed
    offsets_x = [-0.08, 0.08, -0.08, 0.08]
    offsets_y = [-0.08, -0.08, 0.08, 0.08]

    for i, org in enumerate(order):
        ell = Ellipse(xy=(cx + offsets_x[i], cy + offsets_y[i]),
                      width=w, height=h, angle=rotations[i],
                      facecolor=colors[org], alpha=0.35,
                      edgecolor=colors[org], linewidth=2.0, linestyle="-",
                      label=f"{org}  (n={len(sets[org]):,})")
        ax.add_patch(ell)

    # Label each ellipse outside
    label_xy = [
        (0.09, 0.10),
        (0.91, 0.10),
        (0.09, 0.90),
        (0.91, 0.90),
    ]
    for i, org in enumerate(order):
        ax.annotate(org, xy=label_xy[i], fontsize=13, fontweight="bold",
                    color=colors[org], ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.4", fc="white",
                              ec=colors[org], lw=1.8))

    # Sample the 16 regions by grid test (simple Monte Carlo approach)
    # to label each region's count with an (x,y) at the region's centroid
    n_sample = 80000
    rng = np.random.default_rng(42)
    xs = rng.uniform(0, 1, n_sample)
    ys = rng.uniform(0, 1, n_sample)
    region_pts = {}
    for x, y in zip(xs, ys):
        bits = []
        for i in range(4):
            ex = x - (cx + offsets_x[i])
            ey = y - (cy + offsets_y[i])
            theta = np.deg2rad(-rotations[i])
            rx = ex * np.cos(theta) - ey * np.sin(theta)
            ry = ex * np.sin(theta) + ey * np.cos(theta)
            inside = (rx / (w/2))**2 + (ry / (h/2))**2 <= 1.0
            bits.append("1" if inside else "0")
        bit = "".join(bits)
        region_pts.setdefault(bit, []).append((x, y))

    # Place the count label at the median x,y of points in that region
    placed = 0
    for bit, pts in region_pts.items():
        cnt = counts.get(bit, 0)
        if cnt == 0:
            continue
        if len(pts) < 15:
            continue
        xpos = float(np.median([p[0] for p in pts]))
        ypos = float(np.median([p[1] for p in pts]))
        nbits = bit.count("1")
        color = "black" if nbits <= 2 else "#2a0b0b"
        weight = "bold" if nbits >= 3 else "normal"
        bkg = None
        if nbits == 4:
            bkg = dict(boxstyle="round,pad=0.35", fc="#fff3a3",
                       ec="#b8860b", lw=2)
            color = "#5c4300"
        elif nbits == 3:
            bkg = dict(boxstyle="round,pad=0.3", fc="#fde9ef",
                       ec="#943158", lw=1.4)
            color = "#5e1a36"
        elif nbits == 2:
            bkg = dict(boxstyle="round,pad=0.25", fc="#e8f4f8",
                       ec="#2e6a80", lw=1)

        ax.text(xpos, ypos, f"n = {cnt:,}",
                ha="center", va="center", fontsize=10 if nbits < 4 else 11.5,
                fontweight=weight, color=color,
                bbox=bkg)
        placed += 1
    print(f"  Annotated counts on {placed} regions (non-zero n)")

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.legend(loc="upper center", fontsize=11,
              bbox_to_anchor=(0.5, 1.02),
              ncol=2, frameon=True, framealpha=0.95,
              edgecolor="#666666")

fig2, ax2 = plt.subplots(figsize=(16, 13))
draw_4ellipse_venn(gene_sets, order, organ_colors, ax2)
fig2.suptitle("4-Organ Venn — Shared Significant Fibrosis DEGs\n"
              "Kidney  •  Liver  •  Lung  •  Skin"
              "         |         Filter: adj_p < 0.05, |logFC| > 0.585",
              fontsize=18, fontweight="bold", y=0.985, color="#1a1a4a")

# Add a side info box with key intersections
stats_text = []
stats_text.append("KEY OVERLAPS:")
stats_text.append(f"  All 4 organs shared:    {counts.get('1111',0):>5,}")
stats_text.append(f"  Only K & Li & Lu (no S):  {counts.get('1110',0):>5,}")
stats_text.append(f"  Only K & Li & S (no Lu):  {counts.get('1101',0):>5,}")
stats_text.append(f"  Only K & Lu & S (no Li):  {counts.get('1011',0):>5,}")
stats_text.append(f"  Only Li & Lu & S (no K):  {counts.get('0111',0):>5,}")
stats_text.append(f"  Kidney-specific only:     {counts.get('1000',0):>5,}")
stats_text.append(f"  Liver-specific only:      {counts.get('0100',0):>5,}")
stats_text.append(f"  Lung-specific only:       {counts.get('0010',0):>5,}")
stats_text.append(f"  Skin-specific only:       {counts.get('0001',0):>5,}")

ax2.text(1.01, 0.5, "\n".join(stats_text),
         transform=ax2.transAxes, ha="left", va="center", fontsize=11.5,
         fontweight="bold", color="#1a1a4a",
         bbox=dict(boxstyle="round,pad=0.6", fc="#f5f7ff",
                   ec="#7a7abb", lw=1.8))

out2 = os.path.join(BASE, "venn_4organ_manual_ellipses.png")
fig2.savefig(out2, dpi=220, bbox_inches="tight", facecolor="white")
plt.close(fig2)
print(f"[SAVED] Manual 4-ellipse Venn: {out2}")

# ---------- SUMMARY PRINT ----------
print("\n" + "="*78)
print("4-ORGAN VENN SUMMARY")
print("="*78)
print(f"\n{'Set':<8} {'# genes':>10}")
print("-"*22)
for o in order:
    print(f"{o:<8} {len(gene_sets[o]):>10,}")

print("\nALL 4 SHARED (pan-fibrotic core):  {:,} genes".format(counts.get("1111", 0)))
three_way_total = sum(c for b, c in counts.items() if b.count("1") == 3)
two_way_total   = sum(c for b, c in counts.items() if b.count("1") == 2)
one_way_total   = sum(c for b, c in counts.items() if b.count("1") == 1)
print(f"Genes in exactly 3 of 4 organs:  {three_way_total:,}")
print(f"Genes in exactly 2 of 4 organs:  {two_way_total:,}")
print(f"Genes organ-specific (only 1):   {one_way_total:,}")
print("="*78)

# Also save region counts as CSV for reference
rows = []
for bits, cnt in sorted(counts.items()):
    present = [order[i] for i, b in enumerate(bits) if b == "1"]
    absent  = [order[i] for i, b in enumerate(bits) if b == "0"]
    rows.append({
        "bit_mask": bits,
        "count": cnt,
        "n_organs_shared": len(present),
        "present_in": ",".join(present) if present else "(none)",
        "absent_in": ",".join(absent) if absent else "(none)",
    })
regions_csv = os.path.join(BASE, "venn_4organ_region_counts.csv")
pd.DataFrame(rows).sort_values("n_organs_shared", ascending=False).to_csv(regions_csv, index=False)
print(f"[SAVED] All 16 region counts to CSV: {regions_csv}")
print("DONE.")
