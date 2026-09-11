#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ECM vs Non-ECM — CORRECTED Visualization
=========================================
Uses the *corrected* independent discovery sets (validation GSEs excluded)
to regenerate the 6 comparison charts and final comparison statistics.

This replaces the inflated liver numbers from the duplicated-data version
with the honest independent validation results.
"""

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from itertools import combinations

# ---------------------------------------------------------------
# Theme (kept consistent with earlier viz for easy comparison)
# ---------------------------------------------------------------
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({
    "figure.facecolor": "#0B1020",
    "axes.facecolor":   "#111833",
    "axes.edgecolor":   "#2A3466",
    "axes.labelcolor":  "#E8ECFF",
    "axes.titlecolor":  "#FFFFFF",
    "text.color":       "#E8ECFF",
    "xtick.color":      "#C8CEEA",
    "ytick.color":      "#C8CEEA",
    "grid.color":       "#1E2750",
    "grid.alpha":       0.55,
    "legend.facecolor": "#111833",
    "legend.edgecolor": "#2A3466",
    "legend.labelcolor":"#E8ECFF",
    "font.family":      "DejaVu Sans",
    "axes.titleweight": "black",
    "axes.titlesize":   16,
    "axes.labelsize":   13,
    "savefig.facecolor":"#0B1020",
    "savefig.dpi":      220,
    "savefig.bbox":     "tight",
})
ECM_PALETTE = ["#22D3EE", "#F43F5E"]
ECM_DIV = sns.diverging_palette(220, 20, as_cmap=True, center="dark")

ORGANS = ["kidney", "liver", "lung", "skin"]
ADJ_P = 0.05

BASE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(BASE, "report", "figures")
os.makedirs(FIG, exist_ok=True)


def load():
    full = pd.read_csv(os.path.join(
        BASE, "results", "pan_fibrotic_core_genes_validated_corrected.csv"))
    full["gene"] = full["gene"].astype(str).str.upper().str.strip()
    # Ensure ECM flag is present — merge from corrected core raw if missing
    if "is_ECM_gene" not in full.columns or "matrisome_category" not in full.columns:
        raw = pd.read_csv(os.path.join(BASE, "pan_fibrotic_core_genes_corrected.csv"))
        raw["gene"] = raw["gene"].astype(str).str.upper().str.strip()
        full = full.merge(
            raw[["gene", "is_ECM_gene", "matrisome_category"]],
            on="gene", how="left", suffixes=("", "_raw"),
        )
        for col in ("is_ECM_gene", "matrisome_category"):
            if col not in full.columns or full[col].isna().all():
                full[col] = full.get(f"{col}_raw", pd.Series(dtype=object))
    ecm = full[full["is_ECM_gene"] == True].copy()
    non = full[full["is_ECM_gene"] == False].copy()
    print(f"Total corrected core genes: {len(full)}")
    print(f"  ECM genes:     {len(ecm)}")
    print(f"  Non-ECM genes: {len(non)}")
    return full, ecm, non


def fig1_per_organ_bar(ecm, non):
    rows = []
    for s_df, lab in [(ecm, "ECM"), (non, "Non-ECM")]:
        tot = len(s_df)
        for o in ORGANS:
            f = s_df[f"{o}_val_logFC"].notna().sum()
            c = int(s_df[f"{o}_validated"].fillna(False).sum())
            rows.append({
                "subset": lab, "organ": o,
                "pct_found_concordant": c / f * 100 if f else 0,
                "pct_total_concordant": c / tot * 100 if tot else 0,
                "n_concordant": c, "n_found": int(f), "n_total": tot,
            })
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(13, 7))
    sns.barplot(data=df, y="organ", x="pct_found_concordant",
                hue="subset", order=ORGANS, palette=ECM_PALETTE,
                edgecolor="#0B1020", linewidth=1.5, width=0.65, ax=ax)
    for p in ax.patches:
        w = p.get_width()
        if w > 0:
            ax.text(w + 1.2, p.get_y() + p.get_height() / 2,
                    f"{w:.1f}%", va="center", ha="left",
                    fontsize=12, fontweight="bold", color="#FFFFFF")
    ax.set_xlim(0, 110)
    ax.set_xlabel("% of Found Genes with Direction Concordance")
    ax.set_ylabel("")
    ax.set_title("Per-Organ Direction Concordance: ECM vs Non-ECM\n"
                 "(Corrected: Validation GEOs Excluded from Discovery)")
    ax.legend(title="Gene Subset", loc="lower right", frameon=True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    out = os.path.join(FIG, "A_per_organ_concordance_ecm_vs_non_ecm_CORRECTED.png")
    fig.savefig(out); plt.close(fig); print(f"  {out}")
    return df


def fig2_funnel(ecm, non):
    rows = []
    for s_df, lab in [(ecm, "ECM"), (non, "Non-ECM")]:
        tot = len(s_df)
        for crit, col in [("Dir. only", "n_organs_validated"),
                          ("Dir. + adj_p<0.05", "n_organs_validated_sig")]:
            for k in [0, 1, 2, 3, 4]:
                n = tot if k == 0 else int((s_df[col] >= k).sum())
                rows.append({"subset": lab, "criteria": crit,
                             "tier": k, "n": n,
                             "pct": round(n / tot * 100, 2) if tot else 0})
    fdf = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    for ax, crit in zip(axes, ["Dir. only", "Dir. + adj_p<0.05"]):
        sub = fdf[fdf["criteria"] == crit]
        tiers = [0, 1, 2, 3, 4]; x = np.arange(len(tiers)); w = 0.36
        ev = [sub[(sub["subset"] == "ECM") & (sub["tier"] == t)]["n"].values[0]
              if len(sub[(sub["subset"] == "ECM") & (sub["tier"] == t)]) else 0
              for t in tiers]
        nv = [sub[(sub["subset"] == "Non-ECM") & (sub["tier"] == t)]["n"].values[0]
              if len(sub[(sub["subset"] == "Non-ECM") & (sub["tier"] == t)]) else 0
              for t in tiers]
        ax.bar(x - w/2, ev, w, label=f"ECM ({len(ecm)})",
               color=ECM_PALETTE[0], edgecolor="#0B1020", linewidth=1.3)
        ax.bar(x + w/2, nv, w, label=f"Non-ECM ({len(non)})",
               color=ECM_PALETTE[1], edgecolor="#0B1020", linewidth=1.3)
        for xi, e, n_v in zip(x, ev, nv):
            if e > 0:
                ax.text(xi - w/2, e + max(ev + nv) * 0.015,
                        f"{e}", ha="center", va="bottom",
                        fontweight="bold", color="#E8ECFF", fontsize=11)
            if n_v > 0:
                ax.text(xi + w/2, n_v + max(ev + nv) * 0.015,
                        f"{n_v}", ha="center", va="bottom",
                        fontweight="bold", color="#E8ECFF", fontsize=11)
        ax.set_xticks(x); ax.set_xticklabels(["Total", ">=1", ">=2", ">=3", ">=4"])
        ax.set_xlabel("Organs Validated (Tier)"); ax.set_title(crit, fontsize=14)
        ax.legend(loc="upper right", frameon=True)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("Number of Genes")
    fig.suptitle("Validation Funnel: ECM vs Non-ECM (Corrected, Independent Cohorts)",
                 fontsize=17, fontweight="black", y=1.02)
    fig.tight_layout()
    out = os.path.join(FIG, "B_funnel_ecm_vs_non_ecm_CORRECTED.png")
    fig.savefig(out); plt.close(fig); print(f"  {out}")
    return fdf


def fig3_heatmap(full):
    f4 = full[full["n_organs_validated"] == 4].copy()
    if len(f4) == 0:
        print("  [SKIP] 0 genes validated in all 4.")
        return
    ecm_f4 = f4[f4["is_ECM_gene"] == True].sort_values("gene")
    non_f4 = f4[f4["is_ECM_gene"] == False].sort_values("gene")
    ord_df = pd.concat([ecm_f4, non_f4], ignore_index=True)
    cols, labels = [], []
    for o in ORGANS:
        cols += [f"{o}_logFC", f"{o}_val_logFC"]
        labels += [f"Discovery\n{o.capitalize()}", f"Validation\n{o.capitalize()}"]
    mat = ord_df[cols].values.astype(float); mat = np.where(np.isnan(mat), 0, mat)
    fig_h = max(6, 0.38 * len(ord_df) + 2.5)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    vmax = np.nanmax(np.abs(mat)) or 1
    sns.heatmap(mat, cmap=ECM_DIV, center=0, vmin=-vmax, vmax=vmax,
                xticklabels=labels, yticklabels=ord_df["gene"].values,
                ax=ax, cbar=True, linewidths=0.5, linecolor="#0B1020",
                cbar_kws={"label": "log2 Fold-Change", "shrink": 0.7, "pad": 0.02})
    ne = len(ecm_f4)
    if 0 < ne < len(ord_df):
        ax.hlines(ne, *ax.get_xlim(), colors="#FFFFFF",
                  linestyles="--", lw=2)
    ax.set_title(f"logFC Heatmap — {len(f4)} Genes Validated in All 4 Organs\n"
                 f"(Corrected: ECM={ne}, Non-ECM={len(f4) - ne})",
                 fontsize=15, pad=14)
    ax.set_ylabel(""); ax.tick_params(axis="y", labelsize=9)
    out = os.path.join(FIG, "C_logFC_heatmap_all4organs_CORRECTED.png")
    fig.savefig(out); plt.close(fig); print(f"  {out}")


def fig4_scatter(ecm, non):
    fig, axes = plt.subplots(2, 4, figsize=(20, 10), sharex=False, sharey=False)
    subsets = [(ecm, "ECM Genes", ECM_PALETTE[0], axes[0]),
               (non, "Non-ECM Genes", ECM_PALETTE[1], axes[1])]
    for s_df, s_lab, col, ax_row in subsets:
        for j, o in enumerate(ORGANS):
            ax = ax_row[j]
            xc, yc = f"{o}_logFC", f"{o}_val_logFC"
            tmp = s_df[[xc, yc]].dropna().copy()
            if len(tmp) == 0:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes); continue
            x, y = tmp[xc].values, tmp[yc].values
            c_mask = (x > 0) == (y > 0)
            ax.scatter(x[c_mask], y[c_mask], s=55, c=col,
                       edgecolors="#FFFFFF", linewidths=0.9, alpha=0.9,
                       label=f"Concordant ({c_mask.sum()})")
            d = ~c_mask
            if d.sum() > 0:
                ax.scatter(x[d], y[d], s=55, c="#6B7280",
                           edgecolors="#FFFFFF", linewidths=0.9, alpha=0.75,
                           marker="x", label=f"Discordant ({d.sum()})")
            lims = [min(x.min(), y.min()) * 1.15, max(x.max(), y.max()) * 1.15]
            ax.plot(lims, lims, "w--", lw=1.4, alpha=0.85)
            ax.axhline(0, color="#6B7280", lw=0.7, alpha=0.6)
            ax.axvline(0, color="#6B7280", lw=0.7, alpha=0.6)
            ax.set_xlim(*lims); ax.set_ylim(*lims)
            if len(tmp) >= 3:
                pr, _ = stats.pearsonr(x, y)
                ax.text(0.03, 0.97, f"r={pr:.2f}", transform=ax.transAxes,
                        va="top", ha="left", fontweight="bold", color="#FFFFFF",
                        bbox=dict(boxstyle="round,pad=0.25", fc="#0B1020",
                                  ec="#2A3466", alpha=0.85))
            ax.set_title(f"{o.capitalize()} (n={len(tmp)})", fontsize=13)
            if j == 0:
                ax.set_ylabel(f"{s_lab}\nValidation logFC", fontsize=11)
            if ax is ax_row[-1]:
                ax.legend(loc="lower right", fontsize=8, frameon=True)
    for ax in axes[-1, :]:
        ax.set_xlabel("Discovery logFC", fontsize=12)
    fig.suptitle("Discovery vs Validation logFC — ECM vs Non-ECM "
                 "(Corrected: Truly Independent Cohorts)",
                 fontsize=16, fontweight="black", y=1.015)
    fig.tight_layout()
    out = os.path.join(FIG, "D_scatter_ecm_vs_non_ecm_per_organ_CORRECTED.png")
    fig.savefig(out); plt.close(fig); print(f"  {out}")


def fig5_matrisome(full, ecm):
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey=True)
    all_ecm = ecm["matrisome_category"].fillna("Unknown")
    all_ecm = all_ecm[all_ecm != ""]
    ac = all_ecm.value_counts()

    f4 = full[full["n_organs_validated"] == 4]
    fe = f4[f4["is_ECM_gene"] == True]
    fs = fe["matrisome_category"].fillna("Unknown")
    fs = fs[fs != ""]
    fc = fs.value_counts()

    cats = list(ac.index)
    for c in fc.index:
        if c not in cats:
            cats.append(c)
    ac = ac.reindex(cats, fill_value=0); fc = fc.reindex(cats, fill_value=0)
    x = np.arange(len(cats))
    axes[0].bar(x, ac.values, color=ECM_PALETTE[0],
                edgecolor="#0B1020", lw=1.2)
    for xi, v in zip(x, ac.values):
        if v > 0:
            axes[0].text(xi, v + 0.25, str(v), ha="center", va="bottom",
                         fontweight="bold")
    axes[0].set_xticks(x); axes[0].set_xticklabels(cats, rotation=30, ha="right")
    axes[0].set_title(f"All {len(ecm)} ECM Discovery Genes\n(by Matrisome Category)")
    axes[0].set_ylabel("Gene Count")

    axes[1].bar(x, fc.values, color="#FACC15",
                edgecolor="#0B1020", lw=1.2)
    for xi, v in zip(x, fc.values):
        if v > 0:
            axes[1].text(xi, v + 0.2, str(v), ha="center", va="bottom",
                         fontweight="bold")
    axes[1].set_xticks(x); axes[1].set_xticklabels(cats, rotation=30, ha="right")
    n_fin_ecm = len(fe)
    axes[1].set_title(f"All-4-Organ Validated\nECM Genes (n={n_fin_ecm})")
    for ax in axes:
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    fig.suptitle("ECM Gene Breakdown — All Discovery vs Final Validated (Corrected)",
                 fontsize=15, fontweight="black", y=1.02)
    fig.tight_layout()
    out = os.path.join(FIG, "E_matrisome_category_breakdown_CORRECTED.png")
    fig.savefig(out); plt.close(fig); print(f"  {out}")


def fig6_combo(full):
    combos = []
    for k in range(1, 5):
        for combo in combinations(ORGANS, k):
            mask = pd.Series(True, index=full.index)
            for o in ORGANS:
                c = f"{o}_validated"
                if o in combo:
                    mask &= full[c].fillna(False)
                else:
                    mask &= ~(full[c].fillna(False))
            na = int(mask.sum())
            ne = int(((full["is_ECM_gene"] == True) & mask).sum())
            nn = na - ne
            if na > 0:
                combos.append({
                    "combo": "+".join(o[:3].upper() for o in combo),
                    "label": " + ".join(o.capitalize() for o in combo),
                    "n": na, "n_ecm": ne, "n_non": nn, "k": k,
                })
    if not combos:
        print("  [SKIP] combos empty."); return
    cdf = pd.DataFrame(combos).sort_values(["k", "n"], ascending=[True, False])
    fig, ax = plt.subplots(figsize=(14, 8))
    x = np.arange(len(cdf))
    ax.bar(x, cdf["n_ecm"].values, label="ECM", color=ECM_PALETTE[0],
           edgecolor="#0B1020", lw=1.1)
    ax.bar(x, cdf["n_non"].values, bottom=cdf["n_ecm"].values, label="Non-ECM",
           color=ECM_PALETTE[1], edgecolor="#0B1020", lw=1.1)
    for xi, e, nv in zip(x, cdf["n_ecm"], cdf["n_non"]):
        t = e + nv
        if t > 0:
            ax.text(xi, t + max(cdf["n"]) * 0.015,
                    str(t), ha="center", va="bottom",
                    fontweight="bold", color="#FFFFFF")
    ax.set_xticks(x)
    ax.set_xticklabels(cdf["label"].values, rotation=35, ha="right")
    ax.set_ylabel("Number of Genes")
    ax.set_xlabel("Validation Organ Combination (exact membership)")
    ax.set_title("Gene Validation Overlap by Organ Combination\n"
                 "(Corrected: Stacked ECM vs Non-ECM)")
    ax.legend(frameon=True, loc="upper right")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    out = os.path.join(FIG, "F_organ_combo_stacked_ecm_CORRECTED.png")
    fig.savefig(out); plt.close(fig); print(f"  {out}")


def summary_stats(full, ecm, non):
    print("\n================== CORRECTED SUMMARY STATISTICS ==================")
    print(f"Total corrected core genes: {len(full)}  (ECM={len(ecm)}, Non-ECM={len(non)})")
    for s_name, s_df in [("ECM   ", ecm), ("Non-ECM", non)]:
        t = len(s_df)
        a4 = int((s_df["n_organs_validated"] == 4).sum())
        g3 = int((s_df["n_organs_validated"] >= 3).sum())
        print(f"  {s_name} (n={t:>3}): >=3 organs = {g3:>3} ({g3/t*100:5.1f}%),  "
              f"ALL 4 = {a4:>3} ({a4/t*100:5.1f}%)")
    a = int((ecm["n_organs_validated"] == 4).sum())
    b = len(ecm) - a
    c = int((non["n_organs_validated"] == 4).sum())
    d = len(non) - c
    try:
        from scipy.stats import fisher_exact
        _, fp = fisher_exact([[a, b], [c, d]])
        print(f"Fisher exact test (ECM vs Non-ECM, all-4-organ rate): "
              f"ECM={a}/{len(ecm)}, Non-ECM={c}/{len(non)}, p={fp:.4f}")
    except Exception as e:
        print(f"Fisher test skipped: {e}")

    # Correlation summary (corrected)
    val_dir = os.path.join(BASE, "validation")
    vd = {}
    for o in ORGANS:
        p = os.path.join(val_dir, f"{o}_validation_DEGs.csv")
        dd = pd.read_csv(p); dd["gene"] = dd["gene"].astype(str).str.upper().str.strip()
        vd[o] = dd.set_index("gene")
    print("\n-- Pearson r (Discovery vs Validation logFC, corrected) --")
    for s_name, s_df in [("ECM", ecm), ("Non-ECM", non)]:
        print(f"  {s_name}:")
        for o in ORGANS:
            xc, yc = f"{o}_logFC", f"{o}_val_logFC"
            tmp = s_df[[xc, yc]].dropna()
            if len(tmp) < 3:
                r = float("nan")
            else:
                r, _ = stats.pearsonr(tmp[xc].values, tmp[yc].values)
            print(f"    {o:6s}: r={r:.3f}  (n={len(tmp)})")


def main():
    print("=" * 72)
    print("CORRECTED ECM vs Non-ECM VISUALIZATION")
    print("=" * 72)
    full, ecm, non = load()
    print("\n[PLOT] Generating charts ...")
    fig1_per_organ_bar(ecm, non)
    fig2_funnel(ecm, non)
    fig3_heatmap(full)
    fig4_scatter(ecm, non)
    fig5_matrisome(full, ecm)
    fig6_combo(full)
    summary_stats(full, ecm, non)
    print(f"\nDONE. All CORRECTED charts in: {FIG}")


if __name__ == "__main__":
    main()
