#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Disease-Only Spearman Correlation — Stricter Internal Check
============================================================
Re-runs Spearman correlation using FIBROTIC samples only (controls dropped)
for all 4 organs: kidney, liver, lung, skin.

If rho stays meaningfully positive and significant even without controls
anchoring the low end, this is the strongest possible claim.

Outputs:
  validation_2/<organ>/correlation/<organ>_disease_only_spearman.csv
  validation_2/<organ>/correlation/<organ>_disease_only_spearman.png
  validation_2/disease_only_spearman_summary.csv  (cross-organ table)
"""

import os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace") if hasattr(sys.stdout, "reconfigure") else None

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR   = os.path.join(BASE_DIR, "organ_validation2_data")
OUT_DIR     = os.path.join(BASE_DIR, "validation_2")
VAL1_DIR    = os.path.join(BASE_DIR, "validation_1_all_organs")

GENES = ["AEBP1", "COL1A1", "COL1A2", "COL3A1", "VWF"]   # 5-gene pan-fibrotic signature

PALETTE = {
    "kidney": "#5B8DB8",
    "liver":  "#E07B54",
    "lung":   "#6BAE75",
    "skin":   "#A97DC9",
}

# =============================================================================
def load_val1_genes(organ):
    p = os.path.join(VAL1_DIR, organ, f"{organ}_validated_signature.csv")
    return pd.read_csv(p)["gene"].astype(str).tolist() if os.path.exists(p) else GENES


def run_disease_only(organ):
    print(f"\n{'#'*60}\n# ORGAN: {organ.upper()} — DISEASE-ONLY SPEARMAN\n{'#'*60}")
    corr_dir = os.path.join(OUT_DIR, organ, "correlation")
    os.makedirs(corr_dir, exist_ok=True)

    sg_path   = os.path.join(INPUT_DIR, organ, "sample_groups.csv")
    expr_path = os.path.join(INPUT_DIR, organ, "expr_matrix.csv")
    sg        = pd.read_csv(sg_path)
    expr_df   = pd.read_csv(expr_path)

    # ── Drop controls — keep fibrotic only ───────────────────────────────────
    fib_sg = sg[sg["group"] == "fibrotic"].copy().reset_index(drop=True)
    if fib_sg["severity_numeric"].isna().all():
        print("  [SKIP] No severity_numeric in fibrotic samples.")
        return None

    fib_samples    = fib_sg["sample_id"].tolist()
    fib_severities = fib_sg["severity_numeric"].values.astype(float)
    n_fib          = len(fib_samples)

    print(f"  Fibrotic samples : {n_fib}")
    print(f"  Severity range   : {fib_severities.min():.2f} – {fib_severities.max():.2f}")

    # ── Gene expression subset ────────────────────────────────────────────────
    avail     = [g for g in GENES if g in expr_df["gene"].values]
    miss      = [g for g in GENES if g not in expr_df["gene"].values]
    if miss:
        print(f"  Genes not in matrix (skipped): {miss}")

    expr_sub  = expr_df[expr_df["gene"].isin(avail)].set_index("gene")[fib_samples]

    # ── Spearman per gene ─────────────────────────────────────────────────────
    rows, sig = [], []
    for g in expr_sub.index:
        gv            = expr_sub.loc[g].values.astype(float)
        rho, pval     = stats.spearmanr(gv, fib_severities)
        rows.append({"gene": g, "spearman_rho": round(rho, 6), "p_value": pval,
                     "significant": pval < 0.05})
        if pval < 0.05:
            sig.append(g)

    corr_df = (pd.DataFrame(rows)
               .sort_values("p_value")
               .reset_index(drop=True))

    csv_path = os.path.join(corr_dir, f"{organ}_disease_only_spearman.csv")
    corr_df.to_csv(csv_path, index=False)

    # ── Console table ─────────────────────────────────────────────────────────
    print(f"\n  {'Gene':<12} {'Rho':>8}  {'p-value':>14}  Sig?")
    print(f"  {'-'*48}")
    for _, r in corr_df.iterrows():
        flag = " YES" if r.p_value < 0.05 else "  no"
        print(f"  {r.gene:<12} {r.spearman_rho:>8.4f}  {r.p_value:>14.4e}  {flag}")
    print(f"\n  Significant (p<0.05): {sig if sig else 'None'}")
    print(f"  Saved: {csv_path}")

    # ── Scatter subplots ──────────────────────────────────────────────────────
    nc    = 3
    nr    = int(np.ceil(len(corr_df) / nc))
    color = PALETTE.get(organ, "#888")

    fig, axes = plt.subplots(nr, nc, figsize=(5*nc, 4*nr), constrained_layout=True)
    axf       = np.array(axes).flatten()

    for i, (_, row) in enumerate(corr_df.iterrows()):
        ax = axf[i]; g = row.gene
        gv = expr_sub.loc[g].values.astype(float)
        ax.scatter(fib_severities, gv, color=color, edgecolors="white",
                   linewidth=0.6, s=70, zorder=3, alpha=0.85)
        m, b = np.polyfit(fib_severities, gv, 1)
        xl   = np.linspace(fib_severities.min(), fib_severities.max(), 100)
        ax.plot(xl, m*xl+b, color="#222", linewidth=1.8, linestyle="--", zorder=2)

        bg = "#fff3f0" if row.p_value < 0.05 else "white"
        ax.set_facecolor(bg)
        ax.set_title(
            f"{g}\nrho={row.spearman_rho:+.3f},  p={row.p_value:.3e}"
            + (" [SIG]" if row.p_value < 0.05 else ""),
            fontsize=9, fontweight="bold"
        )
        ax.set_xlabel("Severity Numeric (disease only)", fontsize=8)
        ax.set_ylabel("Expression (log2)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.25, linestyle=":")

    for i in range(len(corr_df), len(axf)):
        axf[i].set_visible(False)

    fig.suptitle(
        f"{organ.capitalize()} — Disease-Only Spearman (controls excluded)\n"
        f"n={n_fib} fibrotic samples  |  {len(sig)}/{len(corr_df)} genes p<0.05",
        fontsize=12, fontweight="bold", color=color
    )

    png = os.path.join(corr_dir, f"{organ}_disease_only_spearman.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png}")

    corr_df["organ"] = organ
    return corr_df


# =============================================================================
def make_summary_heatmap(all_df, out_path):
    """Cross-organ rho heatmap (disease-only)."""
    pivot  = all_df.pivot(index="gene", columns="organ", values="spearman_rho")
    piv_p  = all_df.pivot(index="gene", columns="organ", values="p_value")

    # annotation: rho (p)
    annot  = pivot.copy().astype(str)
    for g in pivot.index:
        for o in pivot.columns:
            r = pivot.loc[g, o]
            p = piv_p.loc[g, o]
            star = "*" if p < 0.05 else ""
            annot.loc[g, o] = f"{r:.3f}{star}"

    fig, ax = plt.subplots(figsize=(7, 4))
    sns.heatmap(pivot, annot=annot, fmt="", cmap="RdYlGn", center=0,
                vmin=-1, vmax=1, linewidths=0.5,
                cbar_kws={"label": "Spearman rho", "shrink": 0.8},
                annot_kws={"fontsize": 11, "fontweight": "bold"}, ax=ax)
    ax.set_title("Disease-Only Spearman rho — 5 Genes × 4 Organs\n"
                 "(* = p<0.05; controls excluded)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Organ", fontsize=10); ax.set_ylabel("Gene", fontsize=10)
    ax.tick_params(axis="x", rotation=0); ax.tick_params(axis="y", rotation=0)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Cross-organ heatmap saved: {out_path}")


# =============================================================================
def make_dot_plot(all_df, out_path):
    """Dot plot: size = -log10(p), color = rho."""
    all_df = all_df.copy()
    all_df["-log10p"]  = -np.log10(all_df["p_value"].clip(lower=1e-15))
    all_df["sig_label"] = all_df["p_value"].apply(lambda x: "p<0.05" if x < 0.05 else "ns")

    organs = ["kidney", "liver", "lung", "skin"]
    fig, axes = plt.subplots(1, 4, figsize=(14, 5), sharey=True, constrained_layout=True)

    for ax, org in zip(axes, organs):
        sub = all_df[all_df["organ"] == org].sort_values("gene")
        sc  = ax.scatter(sub["spearman_rho"], sub["gene"],
                         s=sub["-log10p"]*60 + 40,
                         c=sub["spearman_rho"], cmap="RdYlGn",
                         vmin=-1, vmax=1, edgecolors="#555", linewidth=0.6,
                         zorder=3)
        ax.axvline(0, color="#bbb", linewidth=1, linestyle="--")
        ax.set_title(org.capitalize(), fontsize=11, fontweight="bold",
                     color=PALETTE.get(org, "#333"))
        ax.set_xlabel("Spearman rho", fontsize=9)
        ax.set_xlim(-0.1, 1.05)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", alpha=0.2)
        # label significance
        for _, r in sub.iterrows():
            if r.p_value < 0.05:
                ax.annotate("*", (r.spearman_rho + 0.02, r.gene),
                            fontsize=14, color="#c00", va="center")

    axes[0].set_ylabel("Gene", fontsize=10)
    fig.suptitle("Disease-Only Spearman: 5-Gene × 4-Organ — Controls Excluded\n"
                 "(dot size = -log10 p-value;  * = significant at p<0.05)",
                 fontsize=12, fontweight="bold")

    sm = plt.cm.ScalarMappable(cmap="RdYlGn", norm=plt.Normalize(-1, 1))
    sm.set_array([])
    fig.colorbar(sm, ax=axes, label="Spearman rho", shrink=0.6, pad=0.01)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Dot plot saved: {out_path}")


# =============================================================================
def main():
    print("\n" + "="*60)
    print("  DISEASE-ONLY SPEARMAN — 4 ORGANS")
    print("="*60)

    all_results = []

    for organ in ["kidney", "liver", "lung", "skin"]:
        df = run_disease_only(organ)
        if df is not None:
            all_results.append(df)

    if not all_results:
        print("\n[ERROR] No results collected."); return

    combined = pd.concat(all_results, ignore_index=True)

    # ── Cross-organ summary CSV ───────────────────────────────────────────────
    summary_csv = os.path.join(OUT_DIR, "disease_only_spearman_summary.csv")
    combined.to_csv(summary_csv, index=False)

    # ── Cross-organ pivot table ───────────────────────────────────────────────
    print(f"\n{'='*60}\nCROSS-ORGAN SUMMARY (disease-only Spearman)\n{'='*60}")
    pivot_rho = combined.pivot(index="gene", columns="organ", values="spearman_rho")
    pivot_p   = combined.pivot(index="gene", columns="organ", values="p_value")
    print("\n  Spearman rho (* = p<0.05):")
    print(f"  {'Gene':<10}", end="")
    for o in ["kidney","liver","lung","skin"]: print(f"  {o:>12}", end="")
    print()
    print("  " + "-"*58)
    for g in GENES:
        print(f"  {g:<10}", end="")
        for o in ["kidney","liver","lung","skin"]:
            rho  = pivot_rho.loc[g, o] if (g in pivot_rho.index and o in pivot_rho.columns) else float("nan")
            pval = pivot_p.loc[g, o]   if (g in pivot_p.index   and o in pivot_p.columns)   else 1.0
            star = "*" if pval < 0.05 else " "
            print(f"  {rho:>+10.4f}{star}", end="")
        print()

    # Count fully significant rows
    sig_across = []
    for g in GENES:
        all_sig = all(
            pivot_p.loc[g, o] < 0.05
            for o in ["kidney","liver","lung","skin"]
            if o in pivot_p.columns and g in pivot_p.index
        )
        if all_sig: sig_across.append(g)
    print(f"\n  Genes significant in ALL 4 organs (disease-only): {sig_across}")

    # ── Heatmap & dot plot ────────────────────────────────────────────────────
    heatmap_path = os.path.join(OUT_DIR, "disease_only_spearman_heatmap.png")
    dotplot_path = os.path.join(OUT_DIR, "disease_only_spearman_dotplot.png")
    make_summary_heatmap(combined, heatmap_path)
    make_dot_plot(combined, dotplot_path)

    # ── Final checklist ───────────────────────────────────────────────────────
    print(f"\n{'='*60}\nFINAL FILE CHECKLIST\n{'='*60}")
    files = [summary_csv, heatmap_path, dotplot_path]
    for organ in ["kidney","liver","lung","skin"]:
        cd = os.path.join(OUT_DIR, organ, "correlation")
        files += [
            os.path.join(cd, f"{organ}_disease_only_spearman.csv"),
            os.path.join(cd, f"{organ}_disease_only_spearman.png"),
        ]
    for p in files:
        ok = os.path.exists(p)
        status = "OK     " if ok else "MISSING"
        print(f"  [{status}]  {os.path.relpath(p, BASE_DIR)}")

    print("\n  DONE.\n")


if __name__ == "__main__":
    main()
