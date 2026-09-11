#!/usr/bin/env python3
"""
ECM vs Non-ECM Subset Validation — Comparative Visualization
=============================================================

Splits the 175 discovery pan-fibrotic genes into 41 ECM genes and
134 non-ECM genes, then runs the direction-concordance validation
pipeline *separately* on each subset so you can compare validation
rates between the two groups.

Generates publication-ready comparison charts (PNG) in ./report/figures/.
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

# ---------------------------------------------------------------------------
# Styling: professional "command center" dark theme per user preferences
# ---------------------------------------------------------------------------
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
    "font.weight":      "bold",
    "axes.titleweight": "black",
    "axes.titlesize":   17,
    "axes.labelsize":   13,
    "savefig.facecolor":"#0B1020",
    "savefig.dpi":      200,
    "savefig.bbox":     "tight",
})

ECM_PALETTE = ["#22D3EE", "#F43F5E"]  # Cyan = ECM, Rose = Non-ECM
ECM_CMAP    = sns.diverging_palette(220, 20, as_cmap=True, center="dark")

ORGANS = ["kidney", "liver", "lung", "skin"]
ADJ_P_THRESHOLD = 0.05


def setup_dirs(base_dir):
    fig_dir = os.path.join(base_dir, "report", "figures")
    os.makedirs(fig_dir, exist_ok=True)
    return {
        "base": base_dir,
        "figures": fig_dir,
        "core": os.path.join(base_dir, "results", "pan_fibrotic_core_genes_validated.csv"),
        "core_raw": os.path.join(base_dir, "pan_fibrotic_core_genes.csv"),
    }


def load_and_split(paths):
    """Load validated CSV and split into ECM / Non-ECM subsets."""
    df = pd.read_csv(paths["core"])
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()

    if "is_ECM_gene" not in df.columns:
        raw = pd.read_csv(paths["core_raw"])
        raw["gene"] = raw["gene"].astype(str).str.upper().str.strip()
        df = df.merge(raw[["gene", "is_ECM_gene", "matrisome_category"]],
                      on="gene", how="left", suffixes=("", "_raw"))

    ecm_df     = df[df["is_ECM_gene"] == True].copy()
    non_ecm_df = df[df["is_ECM_gene"] == False].copy()

    print(f"  Total core genes:    {len(df):>3}")
    print(f"  ECM genes subset:    {len(ecm_df):>3}  (is_ECM_gene=True)")
    print(f"  Non-ECM genes subset:{len(non_ecm_df):>3}  (is_ECM_gene=False)")
    return df, ecm_df, non_ecm_df


# ======================================================================
# Subset-level validation statistics
# ======================================================================

def per_subset_per_organ_rates(subset_df, label):
    """Return a per-organ validation rate table for a given subset."""
    rows = []
    total = len(subset_df)
    for organ in ORGANS:
        found    = subset_df[f"{organ}_val_logFC"].notna().sum()
        concord  = subset_df[f"{organ}_validated"].fillna(False).sum()
        sig      = subset_df[f"{organ}_validated_sig"].fillna(False).sum()
        rows.append({
            "subset": label,
            "organ": organ,
            "n_total": total,
            "n_found": int(found),
            "n_concordant": int(concord),
            "n_concordant_sig": int(sig),
            "pct_found_concordant":
                round(concord / found * 100, 2) if found else 0.0,
            "pct_total_concordant":
                round(concord / total * 100, 2) if total else 0.0,
            "pct_total_concordant_sig":
                round(sig / total * 100, 2) if total else 0.0,
        })
    return pd.DataFrame(rows)


def subset_funnel(subset_df, label):
    """Return funnel counts for a subset (direction-only & dir+sig tiers)."""
    total = len(subset_df)
    rows = []
    for criteria, col in [
        ("Direction only", "n_organs_validated"),
        ("Direction + adj_p<0.05", "n_organs_validated_sig"),
    ]:
        for k in [0, 1, 2, 3, 4]:
            if k == 0:
                n = total
            else:
                n = int((subset_df[col] >= k).sum())
            rows.append({
                "subset": label,
                "criteria": criteria,
                "tier": k,
                "gene_count": n,
                "percentage": round(n / total * 100, 2) if total else 0.0,
            })
    return pd.DataFrame(rows)


def subset_correlation(subset_df, val_degs, label):
    """Pearson/Spearman correlation table per organ for a subset."""
    rows = []
    for organ in ORGANS:
        col = f"{organ}_logFC"
        if col not in subset_df.columns:
            continue
        x = subset_df[col].values.astype(float)
        y_col = f"{organ}_val_logFC"
        y = subset_df[y_col].values.astype(float)
        mask = ~(np.isnan(x) | np.isnan(y))
        xc, yc = x[mask], y[mask]
        n = len(xc)
        if n < 3:
            pr = pp = sr = sp = np.nan
        else:
            pr, pp = stats.pearsonr(xc, yc)
            sr, sp = stats.spearmanr(xc, yc)
        rows.append({
            "subset": label, "organ": organ, "n_overlap": n,
            "pearson_r": round(pr, 4) if not pd.isna(pr) else np.nan,
            "pearson_p": pp,
            "spearman_r": round(sr, 4) if not pd.isna(sr) else np.nan,
            "spearman_p": sp,
        })
    return pd.DataFrame(rows)


def load_val_degs(base_dir):
    val_dir = os.path.join(base_dir, "validation")
    out = {}
    for o in ORGANS:
        p = os.path.join(val_dir, f"{o}_validation_DEGs.csv")
        d = pd.read_csv(p)
        d["gene"] = d["gene"].astype(str).str.upper().str.strip()
        out[o] = d.set_index("gene")
    return out


# ======================================================================
# Charts
# ======================================================================

def fig1_per_organ_concordance_bar(ecm_rates, non_rates, fig_dir):
    """
    Chart A: Grouped horizontal bar chart — % of genes direction-concordant
    per organ, ECM vs Non-ECM side-by-side.
    """
    frames = []
    for s_df, lab in [(ecm_rates, "ECM"), (non_rates, "Non-ECM")]:
        tmp = s_df.copy()
        tmp["subset"] = lab
        frames.append(tmp)
    plot_df = pd.concat(frames, ignore_index=True)

    fig, ax = plt.subplots(figsize=(13, 7))
    order = ORGANS
    sns.barplot(
        data=plot_df, y="organ", x="pct_found_concordant",
        hue="subset", order=order, palette=ECM_PALETTE, ax=ax,
        edgecolor="#0B1020", linewidth=1.5,
        width=0.65,
    )
    # Add value labels
    for p in ax.patches:
        w = p.get_width()
        if w > 0:
            ax.text(w + 1.2, p.get_y() + p.get_height() / 2,
                    f"{w:.1f}%", va="center", ha="left",
                    fontsize=12, fontweight="bold", color="#FFFFFF")

    ax.set_xlim(0, 110)
    ax.set_xlabel("% of Found Genes with Direction Concordance")
    ax.set_ylabel("")
    ax.set_title("Per-Organ Direction Concordance:\nECM vs Non-ECM Subsets")
    ax.legend(title="Gene Subset", loc="lower right", frameon=True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    out = os.path.join(fig_dir, "A_per_organ_concordance_ecm_vs_non_ecm.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"  [SAVED] {out}")
    return plot_df


def fig2_funnel_comparison(ecm_funnel, non_funnel, fig_dir):
    """
    Chart B: Two-panel funnel bar chart — validation tiers (0-4 organs)
    for ECM (left) and Non-ECM (right), both criteria rows shown.
    """
    ecm_funnel = ecm_funnel.copy(); ecm_funnel["subset"] = "ECM"
    non_funnel = non_funnel.copy(); non_funnel["subset"] = "Non-ECM"
    plot_df = pd.concat([ecm_funnel, non_funnel], ignore_index=True)

    criteria_list = ["Direction only", "Direction + adj_p<0.05"]
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)

    for ax, criteria in zip(axes, criteria_list):
        sub = plot_df[plot_df["criteria"] == criteria]
        tiers = [0, 1, 2, 3, 4]
        x = np.arange(len(tiers))
        w = 0.36
        ecm_vals = [
            sub[(sub["subset"] == "ECM") & (sub["tier"] == t)]["gene_count"].values[0]
            if len(sub[(sub["subset"] == "ECM") & (sub["tier"] == t)]) else 0
            for t in tiers
        ]
        non_vals = [
            sub[(sub["subset"] == "Non-ECM") & (sub["tier"] == t)]["gene_count"].values[0]
            if len(sub[(sub["subset"] == "Non-ECM") & (sub["tier"] == t)]) else 0
            for t in tiers
        ]
        ax.bar(x - w/2, ecm_vals, w, label="ECM (41)",
               color=ECM_PALETTE[0], edgecolor="#0B1020", linewidth=1.3)
        ax.bar(x + w/2, non_vals, w, label="Non-ECM (134)",
               color=ECM_PALETTE[1], edgecolor="#0B1020", linewidth=1.3)

        for xi, ev, nv in zip(x, ecm_vals, non_vals):
            if ev > 0:
                ax.text(xi - w/2, ev + max(ecm_vals)*0.015,
                        f"{ev}", ha="center", va="bottom",
                        fontweight="bold", color="#E8ECFF", fontsize=11)
            if nv > 0:
                ax.text(xi + w/2, nv + max(non_vals)*0.015,
                        f"{nv}", ha="center", va="bottom",
                        fontweight="bold", color="#E8ECFF", fontsize=11)

        ax.set_xticks(x)
        ax.set_xticklabels(["Total", "≥1", "≥2", "≥3", "≥4"])
        ax.set_xlabel("Organs Validated (Tier)")
        ax.set_title(criteria, fontsize=14)
        ax.legend(loc="upper right", frameon=True)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)

    axes[0].set_ylabel("Number of Genes")
    fig.suptitle("Validation Funnel — ECM vs Non-ECM Subsets",
                 fontsize=18, fontweight="black", y=1.02)
    fig.tight_layout()
    out = os.path.join(fig_dir, "B_funnel_ecm_vs_non_ecm.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"  [SAVED] {out}")
    return plot_df


def fig3_all4organs_heatmap(final_df, fig_dir):
    """
    Chart C: Diverging heatmap of discovery vs validation logFC for the
    51 genes validated in all 4 organs, grouped by ECM status.
    """
    if len(final_df) == 0:
        print("  [SKIP] No genes validated in all 4 organs; skipping heatmap.")
        return

    # Reorder: ECM genes first, then non-ECM
    ecm_final = final_df[final_df["is_ECM_gene"] == True].sort_values("gene")
    non_final = final_df[final_df["is_ECM_gene"] == False].sort_values("gene")
    ordered = pd.concat([ecm_final, non_final], ignore_index=True)

    cols = []
    col_labels = []
    for o in ORGANS:
        cols.append(f"{o}_logFC")
        col_labels.append(f"Discovery\n{o.capitalize()}")
        cols.append(f"{o}_val_logFC")
        col_labels.append(f"Validation\n{o.capitalize()}")

    mat = ordered[cols].values.astype(float)
    mat = np.where(np.isnan(mat), 0, mat)

    fig_h = max(6, 0.32 * len(ordered) + 2.5)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    vmax = np.nanmax(np.abs(mat))
    if vmax == 0:
        vmax = 1

    sns.heatmap(
        mat,
        cmap=ECM_CMAP, center=0, vmin=-vmax, vmax=vmax,
        xticklabels=col_labels,
        yticklabels=ordered["gene"].values,
        ax=ax,
        annot=False, cbar=True,
        linewidths=0.5, linecolor="#0B1020",
        cbar_kws={"label": "log2 Fold-Change", "shrink": 0.7,
                  "pad": 0.02},
    )
    # ECM separator
    n_ecm = len(ecm_final)
    if n_ecm > 0 and n_ecm < len(ordered):
        ax.hlines(n_ecm, *ax.get_xlim(), colors="#FFFFFF",
                  linestyles="--", linewidths=2, label="ECM boundary")

    ax.set_title(f"logFC Heatmap — {len(final_df)} Genes Validated in All 4 Organs\n"
                 f"(ECM: {n_ecm}   Non-ECM: {len(final_df)-n_ecm})",
                 fontsize=15, pad=14)
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=9)
    out = os.path.join(fig_dir, "C_logFC_heatmap_all4organs.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"  [SAVED] {out}")


def fig4_scatter_disc_vs_val_per_organ(ecm_df, non_ecm_df, fig_dir):
    """
    Chart D: 2x4 grid of scatter plots, discovery logFC vs validation
    logFC. Rows: ECM (top), Non-ECM (bottom). Cols: 4 organs.
    Diagonal reference line + concordance coloring.
    """
    fig, axes = plt.subplots(2, 4, figsize=(20, 10), sharex=False, sharey=False)
    subsets = [
        (ecm_df, "ECM Genes (41)", ECM_PALETTE[0], axes[0]),
        (non_ecm_df, "Non-ECM Genes (134)", ECM_PALETTE[1], axes[1]),
    ]

    for s_df, s_label, color, ax_row in subsets:
        for j, organ in enumerate(ORGANS):
            ax = ax_row[j]
            x_col = f"{organ}_logFC"
            y_col = f"{organ}_val_logFC"
            tmp = s_df[[x_col, y_col]].dropna().copy()
            if len(tmp) == 0:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            x = tmp[x_col].values
            y = tmp[y_col].values
            concord = (x > 0) == (y > 0)

            ax.scatter(x[concord], y[concord], s=48,
                       c=color, edgecolors="#FFFFFF", linewidths=0.8,
                       alpha=0.88, label=f"Concordant ({concord.sum()})")
            discord = ~concord
            if discord.sum() > 0:
                ax.scatter(x[discord], y[discord], s=48,
                           c="#6B7280", edgecolors="#FFFFFF", linewidths=0.8,
                           alpha=0.75, marker="x",
                           label=f"Discordant ({discord.sum()})")

            lims = [min(x.min(), y.min()) * 1.15, max(x.max(), y.max()) * 1.15]
            ax.plot(lims, lims, "w--", linewidth=1.4, alpha=0.85,
                    label="y=x (ref)")
            ax.axhline(0, color="#6B7280", linewidth=0.7, alpha=0.6)
            ax.axvline(0, color="#6B7280", linewidth=0.7, alpha=0.6)
            ax.set_xlim(*lims); ax.set_ylim(*lims)

            if len(tmp) >= 3:
                pr, _ = stats.pearsonr(x, y)
                ax.text(0.03, 0.97, f"r = {pr:.2f}",
                        transform=ax.transAxes, va="top", ha="left",
                        fontweight="bold", color="#FFFFFF",
                        bbox=dict(boxstyle="round,pad=0.25", fc="#0B1020",
                                  ec="#2A3466", alpha=0.85))

            ax.set_title(f"{organ.capitalize()}  (n={len(tmp)})",
                         fontsize=13)
            if j == 0:
                ax.set_ylabel(f"{s_label}\nValidation logFC", fontsize=12)
            if ax is ax_row[-1]:
                ax.legend(loc="lower right", fontsize=8, frameon=True)

    for ax in axes[-1, :]:
        ax.set_xlabel("Discovery logFC", fontsize=12)

    fig.suptitle("Discovery vs Validation logFC — ECM vs Non-ECM "
                 "(per organ scatter with concordance coloring)",
                 fontsize=17, fontweight="black", y=1.015)
    fig.tight_layout()
    out = os.path.join(fig_dir, "D_scatter_ecm_vs_non_ecm_per_organ.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"  [SAVED] {out}")


def fig5_matrisome_breakdown(final_df, ecm_df, fig_dir):
    """
    Chart E: Bar chart of matrisome category breakdown for the ECM genes
    in the final all-4-organ validated list vs all ECM discovery genes.
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey=True)

    all_ecm = ecm_df["matrisome_category"].fillna("Unknown")
    all_ecm = all_ecm[all_ecm != ""]
    all_counts = all_ecm.value_counts()

    if len(final_df) > 0:
        fin_ecm = final_df[final_df["is_ECM_gene"] == True]
        fin_series = fin_ecm["matrisome_category"].fillna("Unknown")
        fin_series = fin_series[fin_series != ""]
        fin_counts = fin_series.value_counts()
    else:
        fin_counts = pd.Series(dtype=int)

    # Align categories
    all_cats = list(all_counts.index)
    for c in fin_counts.index:
        if c not in all_cats:
            all_cats.append(c)
    all_counts = all_counts.reindex(all_cats, fill_value=0)
    fin_counts = fin_counts.reindex(all_cats, fill_value=0)

    x = np.arange(len(all_cats))
    axes[0].bar(x, all_counts.values, color=ECM_PALETTE[0],
                edgecolor="#0B1020", linewidth=1.2)
    for xi, v in zip(x, all_counts.values):
        if v > 0:
            axes[0].text(xi, v + 0.25, str(v), ha="center", va="bottom",
                         fontweight="bold")
    axes[0].set_xticks(x); axes[0].set_xticklabels(all_cats, rotation=30, ha="right")
    axes[0].set_title(f"All {len(ecm_df)} ECM Discovery Genes\n(by Matrisome Category)")
    axes[0].set_ylabel("Gene Count")

    axes[1].bar(x, fin_counts.values, color="#FACC15",
                edgecolor="#0B1020", linewidth=1.2)
    for xi, v in zip(x, fin_counts.values):
        if v > 0:
            axes[1].text(xi, v + 0.2, str(v), ha="center", va="bottom",
                         fontweight="bold")
    axes[1].set_xticks(x); axes[1].set_xticklabels(all_cats, rotation=30, ha="right")
    n_fin = len(final_df[final_df["is_ECM_gene"] == True]) if len(final_df) > 0 else 0
    axes[1].set_title(f"Final All-4-Organ Validated\nECM Genes (n={n_fin})")

    for ax in axes:
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    fig.suptitle("ECM Gene Breakdown — All Discovery vs Final Validated",
                 fontsize=16, fontweight="black", y=1.02)
    fig.tight_layout()
    out = os.path.join(fig_dir, "E_matrisome_category_breakdown.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"  [SAVED] {out}")


def fig6_upset_ecm_vs_non_ecm(df, fig_dir):
    """
    Chart F: Stacked validation overlap (organ combination) breakdown
    colored by ECM vs non-ECM status, using a bar chart of combination
    sizes (since upsetplot mixing is complex; this is the standard
    readable alt for 2 subsets).
    """
    from itertools import combinations
    val_cols = [f"{o}_validated" for o in ORGANS]
    combos = []
    for k in range(1, 5):
        for combo in combinations(ORGANS, k):
            mask = pd.Series(True, index=df.index)
            for o in ORGANS:
                col = f"{o}_validated"
                if o in combo:
                    mask &= df[col].fillna(False)
                else:
                    mask &= ~(df[col].fillna(False))
            n_all = int(mask.sum())
            n_ecm = int(((df["is_ECM_gene"] == True) & mask).sum())
            n_non = n_all - n_ecm
            if n_all > 0:
                combos.append({
                    "combo": "+".join([o[:3].upper() for o in combo]),
                    "full_label": " + ".join([o.capitalize() for o in combo]),
                    "n": n_all, "n_ecm": n_ecm, "n_non": n_non,
                    "k": k,
                })
    if not combos:
        print("  [SKIP] No combos found.")
        return
    combos_df = pd.DataFrame(combos).sort_values(["k", "n"], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(14, 8))
    x = np.arange(len(combos_df))
    ax.bar(x, combos_df["n_ecm"].values, label="ECM",
           color=ECM_PALETTE[0], edgecolor="#0B1020", linewidth=1.1)
    ax.bar(x, combos_df["n_non"].values, bottom=combos_df["n_ecm"].values,
           label="Non-ECM", color=ECM_PALETTE[1],
           edgecolor="#0B1020", linewidth=1.1)

    for xi, e, n in zip(x, combos_df["n_ecm"], combos_df["n_non"]):
        tot = e + n
        if tot > 0:
            ax.text(xi, tot + max(combos_df["n"]) * 0.015,
                    str(tot), ha="center", va="bottom",
                    fontweight="bold", color="#FFFFFF")

    ax.set_xticks(x)
    ax.set_xticklabels(combos_df["full_label"].values, rotation=35, ha="right")
    ax.set_ylabel("Number of Genes")
    ax.set_xlabel("Validation Organ Combination (exact membership)")
    ax.set_title("Gene Validation Overlap by Organ Combination\n"
                 "(Stacked: ECM vs Non-ECM)")
    ax.legend(frameon=True, loc="upper right")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    out = os.path.join(fig_dir, "F_organ_combo_stacked_ecm.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"  [SAVED] {out}")


# ======================================================================
# Main
# ======================================================================

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    paths = setup_dirs(base_dir)

    print("=" * 72)
    print("ECM vs NON-ECM SUBSET VALIDATION — COMPARATIVE ANALYSIS")
    print("=" * 72)

    # ---------------- Load & split ----------------
    print("\n[LOAD] Splitting 175 core genes into ECM / Non-ECM subsets ...")
    full_df, ecm_df, non_ecm_df = load_and_split(paths)
    val_degs = load_val_degs(base_dir)

    # ---------------- Subset statistics ----------------
    print("\n[COMPUTE] Per-organ validation rates per subset ...")
    ecm_rates = per_subset_per_organ_rates(ecm_df, "ECM")
    non_rates = per_subset_per_organ_rates(non_ecm_df, "Non-ECM")
    rates_all = pd.concat([ecm_rates, non_rates], ignore_index=True)
    print("\n  -- Per-Organ Direction Concordance Rates --")
    print(rates_all[["subset", "organ", "n_total", "n_found",
                     "n_concordant", "pct_found_concordant",
                     "pct_total_concordant_sig"]].to_string(index=False))

    print("\n[COMPUTE] Validation funnels per subset ...")
    ecm_funnel = subset_funnel(ecm_df, "ECM")
    non_funnel = subset_funnel(non_ecm_df, "Non-ECM")

    print("\n[COMPUTE] Correlations per subset ...")
    ecm_corr = subset_correlation(ecm_df, val_degs, "ECM")
    non_corr = subset_correlation(non_ecm_df, val_degs, "Non-ECM")
    corr_all = pd.concat([ecm_corr, non_corr], ignore_index=True)
    print("\n  -- Pearson Correlation (Discovery vs Validation logFC) --")
    print(corr_all[["subset", "organ", "n_overlap", "pearson_r",
                    "pearson_p"]].to_string(index=False, na_rep="N/A"))

    # Save comparison tables
    rates_csv = os.path.join(paths["base"], "results",
                             "ecm_vs_non_ecm_per_organ_rates.csv")
    funnel_csv = os.path.join(paths["base"], "results",
                              "ecm_vs_non_ecm_funnel.csv")
    corr_csv = os.path.join(paths["base"], "results",
                            "ecm_vs_non_ecm_correlations.csv")
    rates_all.to_csv(rates_csv, index=False)
    pd.concat([ecm_funnel, non_funnel], ignore_index=True).to_csv(
        funnel_csv, index=False)
    corr_all.to_csv(corr_csv, index=False)
    print(f"\n  [SAVED] {rates_csv}")
    print(f"  [SAVED] {funnel_csv}")
    print(f"  [SAVED] {corr_csv}")

    # ---------------- Charts ----------------
    print("\n[PLOT] Generating comparison charts ...")
    fig1_per_organ_concordance_bar(ecm_rates, non_rates, paths["figures"])
    fig2_funnel_comparison(ecm_funnel, non_funnel, paths["figures"])

    final_all4 = full_df[full_df["n_organs_validated"] == 4].copy()
    fig3_all4organs_heatmap(final_all4, paths["figures"])
    fig4_scatter_disc_vs_val_per_organ(ecm_df, non_ecm_df, paths["figures"])
    fig5_matrisome_breakdown(final_all4, ecm_df, paths["figures"])
    fig6_upset_ecm_vs_non_ecm(full_df, paths["figures"])

    # ---------------- Printed summary of key findings ----------------
    print("\n" + "=" * 72)
    print("SUMMARY OF KEY FINDINGS (ECM vs NON-ECM)")
    print("=" * 72)
    for subset_name, s_df in [("ECM   (n=41)", ecm_df),
                              ("Non-ECM (n=134)", non_ecm_df)]:
        all4 = int((s_df["n_organs_validated"] == 4).sum())
        ge3  = int((s_df["n_organs_validated"] >= 3).sum())
        ge1  = int((s_df["n_organs_validated"] >= 1).sum())
        tot  = len(s_df)
        print(f"\n  {subset_name}:")
        print(f"    Validated in >= 1 organ : {ge1:>3} / {tot:>3}  "
              f"({ge1/tot*100:5.1f}%)")
        print(f"    Validated in >= 3 organs: {ge3:>3} / {tot:>3}  "
              f"({ge3/tot*100:5.1f}%)")
        print(f"    Validated in ALL 4     : {all4:>3} / {tot:>3}  "
              f"({all4/tot*100:5.1f}%)")

    # Statistical test: Fisher exact / chi-square for ECM vs non-ECM in
    # "all 4 organs validated"
    a = int((ecm_df["n_organs_validated"] == 4).sum())
    b = len(ecm_df) - a
    c = int((non_ecm_df["n_organs_validated"] == 4).sum())
    d = len(non_ecm_df) - c
    table = np.array([[a, b], [c, d]])
    try:
        from scipy.stats import fisher_exact
        _, fp = fisher_exact(table)
        print(f"\n  Fisher exact test (ECM vs Non-ECM, all-4-organ rate):")
        print(f"    ECM     : {a}/{len(ecm_df)} = {a/len(ecm_df)*100:.1f}%")
        print(f"    Non-ECM : {c}/{len(non_ecm_df)} = {c/len(non_ecm_df)*100:.1f}%")
        print(f"    p-value : {fp:.4f}")
    except Exception as e:
        print(f"  Fisher test skipped: {e}")

    print("\nDONE. All charts saved to ./report/figures/")


if __name__ == "__main__":
    main()
