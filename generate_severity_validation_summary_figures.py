#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Cross-Organ Severity Validation Summary & Honest Reporting Audit Panels
================================================================================
Figure 1: cross_organ_severity_validation_summary.png
  - Top: 2x2 Organ Grid (Liver, Kidney, Lung, Skin) showing metadata icon,
    genes tested (5), significant count, and top gene Spearman rho as a large number.
  - Bottom: Horizontal bar chart comparing % significant genes across organs.
  
Figure 2: severity_honest_reporting_audit_table.png
  - Honest Reporting / Limitations audit infographic table with traffic-light badges
    and 1-line explanatory notes per organ.
"""

import os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

warnings.filterwarnings("ignore")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR  = os.path.join(BASE_DIR, "validation_2")
os.makedirs(OUT_DIR, exist_ok=True)

# Dataset metrics
ORGANS = ["Liver", "Kidney", "Lung", "Skin"]

DATA = {
    "Liver": {
        "accession": "GSE162694",
        "metric_name": "METAVIR Stage (F0-F4)",
        "tested": 5,
        "full_sig": 5,
        "disease_sig": 0,
        "disease_trend": 2, # COL3A1 (rho=0.58, p=0.076), AEBP1 (rho=0.51, p=0.133)
        "top_gene": "COL3A1",
        "top_rho_full": 0.883,
        "top_rho_do": 0.585,
        "color": "#E07B54",
        "status_badge": "Tested — Positive Trend",
        "status_code": "yellow",
        "note": "Real METAVIR stages (F0-F4; n=143). Full-sample p<0.001 (rho=+0.88). Disease-only n=10 power-limited; COL3A1 shows strong trend (rho=+0.58, p=0.076).",
    },
    "Kidney": {
        "accession": "GSE66494",
        "metric_name": "%TIF Fibrosis (5%-75%)",
        "tested": 5,
        "full_sig": 5,
        "disease_sig": 0,
        "disease_trend": 0,
        "top_gene": "COL1A2",
        "top_rho_full": 0.841,
        "top_rho_do": 0.285,
        "color": "#5B8DB8",
        "status_badge": "Tested — Largely Null",
        "status_code": "yellow",
        "note": "Real %TIF (5%-75%; n=61). Full-sample p<0.001 (rho=+0.84). Disease-only n=10 power-limited with compressed 5-stage ordinal scale.",
    },
    "Lung": {
        "accession": "GSE47460",
        "metric_name": "FVC% Deficit (1%-87%)",
        "tested": 5,
        "full_sig": 5,
        "disease_sig": 1, # COL3A1 (rho=0.71, p=0.022)
        "disease_trend": 0,
        "top_gene": "COL3A1",
        "top_rho_full": 0.880,
        "top_rho_do": 0.709,
        "color": "#6BAE75",
        "status_badge": "Tested — Significant",
        "status_code": "green",
        "note": "Widest continuous severity range (1%-87% FVC% deficit; n=441). COL3A1 maintains disease-only significance (rho=+0.71, p=0.022).",
    },
    "Skin": {
        "accession": "GSE130955",
        "metric_name": "mRSS Skin Score (6-43)",
        "tested": 5,
        "full_sig": 5,
        "disease_sig": 0,
        "disease_trend": 1, # VWF (rho=0.55, p=0.098)
        "top_gene": "VWF",
        "top_rho_full": 0.870,
        "top_rho_do": 0.552,
        "color": "#A97DC9",
        "status_badge": "Tested — Largely Null",
        "status_code": "yellow",
        "note": "Real mRSS score (6-43; n=55). Full-sample p<0.001 (rho=+0.87). Disease-only n=10 power-limited; VWF shows positive trend (rho=+0.55, p=0.098).",
    },
}


# =============================================================================
# FIGURE 1: 2x2 Grid + Horizontal Bar Chart
# =============================================================================
def generate_figure_1():
    print("Generating Figure 1: cross_organ_severity_validation_summary.png ...")
    fig = plt.figure(figsize=(13, 10), facecolor="#F8F9FA")
    
    # Outer title
    fig.suptitle("Cross-Organ Severity/Dose-Response Validation Summary", 
                 fontsize=18, fontweight="bold", y=0.97, color="#111111")
    fig.text(0.5, 0.935, "Evaluation of 5 Pan-Fibrotic Core Genes Across Clinical Severity Metrics in 4 Organs",
             ha="center", fontsize=11, color="#555555", style="italic")

    # Grid Spec: 2 rows for 2x2 grid, 1 row for lower bar chart
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.95], hspace=0.38, wspace=0.25,
                         left=0.08, right=0.92, top=0.89, bottom=0.08)

    # ── Top 2x2 Grid ──────────────────────────────────────────────────────────
    coords = [(0, 0), (0, 1), (1, 0), (1, 1)] # Liver, Kidney, Lung, Skin
    for idx, organ in enumerate(ORGANS):
        r, c = coords[idx]
        ax = fig.add_subplot(gs[r, c])
        ax.set_facecolor("white")
        
        # Border box & subtle shadow background
        rect = patches.FancyBboxPatch((0.02, 0.04), 0.96, 0.92, boxstyle="round,pad=0.03,rounding_size=0.04",
                                    facecolor="white", edgecolor="#D0D7DE", linewidth=1.5)
        ax.add_patch(rect)
        
        # Header banner inside box
        banner = patches.FancyBboxPatch((0.02, 0.72), 0.96, 0.24, boxstyle="round,pad=0.01,rounding_size=0.03",
                                       facecolor=DATA[organ]["color"], alpha=0.15, edgecolor="none")
        ax.add_patch(banner)

        # Organ Name & Accession
        ax.text(0.06, 0.83, organ.upper(), fontsize=14, fontweight="bold", color=DATA[organ]["color"], va="center")
        ax.text(0.06, 0.74, f"{DATA[organ]['accession']} | {DATA[organ]['metric_name']}", fontsize=8.5, color="#444444", va="center")
        
        # Metadata Availability Icon (Green Checkmark Badge)
        check_box = patches.FancyBboxPatch((0.74, 0.74), 0.22, 0.18, boxstyle="round,pad=0.02,rounding_size=0.03",
                                          facecolor="#E8F5E9", edgecolor="#2E7D32", linewidth=1.2)
        ax.add_patch(check_box)
        ax.text(0.85, 0.83, "[✓] Metadata", fontsize=8, fontweight="bold", color="#1B5E20", ha="center", va="center")

        # Main Large Rho Number Display
        ax.text(0.06, 0.45, f"+{DATA[organ]['top_rho_full']:.3f}", fontsize=28, fontweight="bold", color="#111111", va="center")
        ax.text(0.06, 0.25, f"Top Gene: {DATA[organ]['top_gene']} (Full-Sample Spearman ρ)", fontsize=8.5, fontweight="bold", color="#333333", va="center")

        # Stats Breakdown Pills
        ax.text(0.58, 0.52, f"Genes Tested: {DATA[organ]['tested']}", fontsize=9, color="#555555", va="center")
        ax.text(0.58, 0.40, f"Full-Sample Sig: {DATA[organ]['full_sig']}/{DATA[organ]['tested']} (100%)", fontsize=9, fontweight="bold", color="#2E7D32", va="center")
        ax.text(0.58, 0.28, f"Disease-Only Sig: {DATA[organ]['disease_sig']}/{DATA[organ]['tested']}", fontsize=9, color="#D32F2F" if DATA[organ]['disease_sig']>0 else "#666666", va="center")

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    # ── Lower Horizontal Bar Chart ────────────────────────────────────────────
    ax_bar = fig.add_subplot(gs[2, :])
    ax_bar.set_facecolor("white")
    
    # We compare Full-Sample vs Disease-Only % Significant Genes
    y_positions = np.arange(len(ORGANS))
    bar_height = 0.32

    # Data for bars (% significant)
    full_sample_pct = [100.0, 100.0, 100.0, 100.0]
    disease_only_pct = [0.0, 0.0, 20.0, 0.0] # 20% for Lung (COL3A1)
    
    # Bar Plot 1: Full Sample (All 100% Red/Coral)
    bars_full = ax_bar.barh(y_positions + bar_height/2, full_sample_pct, height=bar_height, 
                            color="#E53935", alpha=0.9, edgecolor="#B71C1C", linewidth=1.2, label="Full-Sample Correlation (Controls Included, p < 0.05)")

    # Bar Plot 2: Disease Only (Muted Gray / Red for Sig)
    disease_colors = ["#9E9E9E", "#9E9E9E", "#E53935", "#9E9E9E"] # Lung is red (sig), others gray
    bars_do = ax_bar.barh(y_positions - bar_height/2, [max(p, 2.0) for p in disease_only_pct], height=bar_height, 
                          color=disease_colors, alpha=0.85, edgecolor="#424242", linewidth=1.2, label="Disease-Only Correlation (Controls Excluded, p < 0.05)")

    # Value Labels on Bars
    for i, p in enumerate(full_sample_pct):
        ax_bar.text(p + 1.5, y_positions[i] + bar_height/2, f"100% (5/5 genes, p < 0.001)", 
                    va="center", fontsize=8.5, fontweight="bold", color="#B71C1C")

    for i, p in enumerate(disease_only_pct):
        if p > 0:
            ax_bar.text(p + 1.5, y_positions[i] - bar_height/2, f"20% (1/5 genes: COL3A1 p=0.022)", 
                        va="center", fontsize=8.5, fontweight="bold", color="#B71C1C")
        else:
            note_str = "0% strict p<0.05 (COL3A1 ρ=+0.58 p=0.076 trend)" if ORGANS[i]=="Liver" else ("0% strict p<0.05 (VWF ρ=+0.55 p=0.098 trend)" if ORGANS[i]=="Skin" else "0% strict p<0.05 (n=10 power limited)")
            ax_bar.text(4.0, y_positions[i] - bar_height/2, note_str, 
                        va="center", fontsize=8, color="#555555", style="italic")

    ax_bar.set_yticks(y_positions)
    ax_bar.set_yticklabels([f"{o}" for o in ORGANS], fontsize=11, fontweight="bold")
    ax_bar.set_xlim(0, 125)
    ax_bar.set_xlabel("% of Tested Genes with Significant Severity Correlation (p < 0.05)", fontsize=10, fontweight="bold")
    ax_bar.set_title("Comparison of % Significant Severity Correlation (Full-Sample vs Disease-Only)", fontsize=11, fontweight="bold", pad=10)
    
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.grid(axis="x", linestyle=":", alpha=0.5)

    # Legend
    # Custom legend elements including hatched gray for "Metadata Not Available"
    import matplotlib.patches as mpatches
    leg_full = mpatches.Patch(facecolor="#E53935", edgecolor="#B71C1C", label="Significant Correlation (p < 0.05)")
    leg_null = mpatches.Patch(facecolor="#9E9E9E", edgecolor="#424242", label="Tested — Available but Not Significant")
    leg_na   = mpatches.Patch(facecolor="#E0E0E0", edgecolor="#666666", hatch="///", label="Metadata Not Available (Not Attempted)")
    
    ax_bar.legend(handles=[leg_full, leg_null, leg_na], loc="lower right", fontsize=8.5, frameon=True, facecolor="#F8F9FA", edgecolor="#D0D7DE")

    out_path = os.path.join(OUT_DIR, "cross_organ_severity_validation_summary.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved Figure 1 -> {out_path}")


# =============================================================================
# FIGURE 2: Honest Reporting Audit Table / Infographic Panel
# =============================================================================
def generate_figure_2():
    print("Generating Figure 2: severity_honest_reporting_audit_table.png ...")
    fig = plt.figure(figsize=(12, 7.5), facecolor="#F8F9FA")

    # Titles
    fig.suptitle("Severity/Dose-Response Testing: Honest Reporting Across Organs", 
                 fontsize=17, fontweight="bold", y=0.96, color="#111111")
    fig.text(0.5, 0.915, "Transparent Audit of Clinical Staging Metadata, Statistical Power Constraints, and Dose-Response Outcomes",
             ha="center", fontsize=10.5, color="#555555", style="italic")

    ax = fig.add_axes([0.05, 0.06, 0.90, 0.81])
    ax.set_facecolor("white")

    # Outer border for audit table
    outer_box = patches.FancyBboxPatch((0.01, 0.02), 0.98, 0.96, boxstyle="round,pad=0.02,rounding_size=0.03",
                                       facecolor="white", edgecolor="#D0D7DE", linewidth=1.5)
    ax.add_patch(outer_box)

    # Table Header Banner
    tbl_header = patches.FancyBboxPatch((0.01, 0.85), 0.98, 0.13, boxstyle="round,pad=0.01,rounding_size=0.02",
                                        facecolor="#ECEFF1", edgecolor="#CFD8DC", linewidth=1.0)
    ax.add_patch(tbl_header)
    
    ax.text(0.04, 0.915, "ORGAN", fontsize=10, fontweight="bold", color="#37474F", va="center")
    ax.text(0.24, 0.915, "TESTING STATUS & BADGE", fontsize=10, fontweight="bold", color="#37474F", va="center")
    ax.text(0.58, 0.915, "TRANSPARENT AUDIT NOTE & BIOLOGICAL RATIONALE", fontsize=10, fontweight="bold", color="#37474F", va="center")

    # 4 Rows for 4 Organs
    y_starts = [0.66, 0.46, 0.26, 0.06]
    row_height = 0.17

    badge_styles = {
        "green":  {"bg": "#E8F5E9", "fg": "#1B5E20", "border": "#2E7D32"},
        "yellow": {"bg": "#FFF8E1", "fg": "#F57F17", "border": "#FBC02D"},
        "gray":   {"bg": "#F5F5F5", "fg": "#616161", "border": "#9E9E9E"},
    }

    for idx, organ in enumerate(ORGANS):
        y = y_starts[idx]
        
        # Row divider line
        if idx < 3:
            ax.plot([0.03, 0.97], [y - 0.02, y - 0.02], color="#E0E0E0", linewidth=1.0, linestyle="-")

        # Organ Name & Metadata Metric
        ax.text(0.04, y + 0.09, organ, fontsize=13, fontweight="bold", color=DATA[organ]["color"], va="center")
        ax.text(0.04, y + 0.03, f"{DATA[organ]['accession']} ({DATA[organ]['metric_name']})", fontsize=8.5, color="#666666", va="center")

        # Status Badge Pill
        b_code = DATA[organ]["status_code"]
        b_style = badge_styles[b_code]
        
        badge_box = patches.FancyBboxPatch((0.24, y + 0.03), 0.28, 0.08, boxstyle="round,pad=0.02,rounding_size=0.03",
                                          facecolor=b_style["bg"], edgecolor=b_style["border"], linewidth=1.2)
        ax.add_patch(badge_box)
        ax.text(0.38, y + 0.07, DATA[organ]["status_badge"], fontsize=9, fontweight="bold", 
                color=b_style["fg"], ha="center", va="center")

        # Audit Note
        ax.text(0.56, y + 0.06, DATA[organ]["note"], fontsize=8.5, color="#222222", va="center", wrap=True)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    out_path = os.path.join(OUT_DIR, "severity_honest_reporting_audit_table.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved Figure 2 -> {out_path}")


def main():
    generate_figure_1()
    generate_figure_2()
    print("All validation summary figures generated successfully!")


if __name__ == "__main__":
    main()
