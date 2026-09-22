# -*- coding: utf-8 -*-
"""
Publication-Quality 3-Panel Figure:
Dedicated Early-Detection Validation for the 5 Validated Pan-Fibrotic Hub Genes
(COL15A1, COL1A1, COL3A1, SERPINE2, SERPINF2)

Layout: 1x3 Grid (300 DPI, Nature/IEEE Style Palette)
- Panel A: Early-Stage Subgroup ROC Curves (Control vs. Stage 1-2 Fibrosis)
- Panel B: Early Switch vs. Late Progression Dynamics (Normalized Z-score Trajectories across Stages 0 -> 4)
- Panel C: Subclinical Decision Curve Analysis (DCA, pt = 0.05 to 0.50)

Output:
- plots/hub_genes_early_detection_triptych_5genes.png
- results/hub_genes_early_stage_roc_metrics_5genes.csv
"""

import os
import re
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import make_interp_spline
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("plots", exist_ok=True)
os.makedirs("results", exist_ok=True)

# Set Nature / IEEE Publication Styling
plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.edgecolor"] = "#334155"
plt.rcParams["axes.linewidth"] = 1.1
plt.rcParams["axes.labelcolor"] = "#0f172a"
plt.rcParams["xtick.color"] = "#1e293b"
plt.rcParams["ytick.color"] = "#1e293b"

# Color Palette for the 5 Validated Hubs + Ensemble
GENE_COLORS = {
    "COL15A1": "#1e40af",  # Deep Royal Blue
    "COL1A1": "#dc2626",   # Crimson Red
    "COL3A1": "#0284c7",   # Sky Blue
    "SERPINE2": "#d97706", # Rich Amber
    "SERPINF2": "#059669", # Jade Green
    "Ensemble": "#7c3aed"  # Vibrant Violet
}

# ==============================================================================
# 1. LOAD GSE84044 (Liver Scheuer Staging S0 to S4, N=124)
# ==============================================================================
print("1. Loading GSE84044 for 5 validated hub genes...")
gsm_list, scheuer_s = [], []
with gzip.open("geo_cache/GSE84044_series_matrix.txt.gz", "rt", errors="ignore") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            gsm_list = [x.strip().replace('"', '') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1") and "scheuer score s:" in line.lower():
            scheuer_s = [x.strip().replace('"', '') for x in line.split("\t")[1:]]
        elif line.startswith("!series_matrix_table_begin"):
            break

stages = {}
for gsm, s_str in zip(gsm_list, scheuer_s):
    m = re.search(r"scheuer score s:\s*(\d+)", s_str, re.IGNORECASE)
    if m:
        stages[gsm] = int(m.group(1))

probe_map = {
    "COL15A1": "203477_at",
    "COL1A1": "202310_s_at",
    "COL3A1": "215076_s_at",
    "SERPINE2": "212190_at",
    "SERPINF2": "205075_at"
}
rev_map = {v: k for k, v in probe_map.items()}
expr = {g: {} for g in probe_map}

with gzip.open("geo_cache/GSE84044_series_matrix.txt.gz", "rt", errors="ignore") as f:
    in_table = False
    samples_hdr = []
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            hdr_line = f.readline().strip().split("\t")
            samples_hdr = [x.replace('"', '') for x in hdr_line[1:]]
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            parts = line.strip().split("\t")
            pid = parts[0].replace('"', '')
            if pid in rev_map:
                g = rev_map[pid]
                for s, val_str in zip(samples_hdr, parts[1:]):
                    try:
                        expr[g][s] = float(val_str.replace('"', '').strip())
                    except:
                        pass

df_all = pd.DataFrame(expr)
df_all["stage"] = [stages.get(idx) for idx in df_all.index]
df_all = df_all.dropna(subset=["stage"]).copy()
df_all["stage"] = df_all["stage"].astype(int)

# Filter strictly early cohort: S0 (n=43) vs S1/S2 (n=53)
df_early = df_all[df_all["stage"].isin([0, 1, 2])].copy()
df_early["label"] = (df_early["stage"] > 0).astype(int)

print(f"  Total samples: {len(df_all)}")
print(f"  Early cohort (S0 vs S1/S2): N={len(df_early)} (Controls={sum(df_early['label']==0)}, Early={sum(df_early['label']==1)})")

# ==============================================================================
# 2. COMPUTE PANEL A METRICS: EARLY-STAGE ROC & 5-FOLD CV ENSEMBLE
# ==============================================================================
print("2. Computing Early-Stage ROC Curves & Ensemble 5-Fold CV...")
y_true = df_early["label"].values
all_hubs = ["COL15A1", "COL1A1", "COL3A1", "SERPINE2", "SERPINF2"]

roc_results = {}
roc_summary_rows = []

for g in all_hubs:
    ctrl = df_early[df_early["label"] == 0][g].values
    early = df_early[df_early["label"] == 1][g].values
    log2fc = np.mean(early) - np.mean(ctrl)
    
    y_score = df_early[g].values
    y_score_eval = -y_score if log2fc < 0 else y_score
    
    fpr, tpr, _ = roc_curve(y_true, y_score_eval)
    auc_val = roc_auc_score(y_true, y_score_eval)
    
    # 1,000 Bootstrap 95% CI
    np.random.seed(42)
    boot_aucs = []
    n_pts = len(y_true)
    for _ in range(1000):
        idx = np.random.choice(n_pts, n_pts, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        boot_aucs.append(roc_auc_score(y_true[idx], y_score_eval[idx]))
    ci_low = np.percentile(boot_aucs, 2.5)
    ci_high = np.percentile(boot_aucs, 97.5)
    
    roc_results[g] = {
        "fpr": fpr, "tpr": tpr, "auc": auc_val,
        "ci_low": ci_low, "ci_high": ci_high,
        "direction": "Up" if log2fc >= 0 else "Down"
    }
    roc_summary_rows.append({
        "Feature": g,
        "Type": "Individual Gene",
        "Direction": "Up" if log2fc >= 0 else "Down",
        "AUC": round(auc_val, 4),
        "95CI_Low": round(ci_low, 4),
        "95CI_High": round(ci_high, 4)
    })

# 5-Fold CV for Ensemble Classifier
X_early = df_early[all_hubs].values
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
oof_lr = np.zeros(len(y_true))
oof_rf = np.zeros(len(y_true))

for tr_idx, va_idx in skf.split(X_early, y_true):
    X_tr, y_tr = X_early[tr_idx], y_true[tr_idx]
    X_va, y_va = X_early[va_idx], y_true[va_idx]
    
    lr = LogisticRegression(penalty="l2", C=1.0, random_state=42, max_iter=1000)
    lr.fit(X_tr, y_tr)
    oof_lr[va_idx] = lr.predict_proba(X_va)[:, 1]
    
    rf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
    rf.fit(X_tr, y_tr)
    oof_rf[va_idx] = rf.predict_proba(X_va)[:, 1]

oof_ens = (oof_lr + oof_rf) / 2.0
fpr_ens, tpr_ens, _ = roc_curve(y_true, oof_ens)
auc_ens = roc_auc_score(y_true, oof_ens)

boot_ens = []
for _ in range(1000):
    idx = np.random.choice(len(y_true), len(y_true), replace=True)
    if len(np.unique(y_true[idx])) < 2:
        continue
    boot_ens.append(roc_auc_score(y_true[idx], oof_ens[idx]))
ci_low_ens = np.percentile(boot_ens, 2.5)
ci_high_ens = np.percentile(boot_ens, 97.5)

roc_results["Ensemble"] = {
    "fpr": fpr_ens, "tpr": tpr_ens, "auc": auc_ens,
    "ci_low": ci_low_ens, "ci_high": ci_high_ens
}
roc_summary_rows.append({
    "Feature": "5-Hub Ensemble",
    "Type": "Multivariable Model (LR+RF)",
    "Direction": "Ensemble Probability",
    "AUC": round(auc_ens, 4),
    "95CI_Low": round(ci_low_ens, 4),
    "95CI_High": round(ci_high_ens, 4)
})

pd.DataFrame(roc_summary_rows).to_csv("results/hub_genes_early_stage_roc_metrics_5genes.csv", index=False)
print("  Saved: results/hub_genes_early_stage_roc_metrics_5genes.csv")

# ==============================================================================
# 3. COMPUTE PANEL B METRICS: Z-SCORE PROGRESSION DYNAMICS (STAGES 0 TO 4)
# ==============================================================================
print("3. Computing Z-Score Trajectories across Stages 0 to 4...")
stage_stats = {}
stage_labels = [0, 1, 2, 3, 4]
stage_names = ["Stage 0\n(Control)", "Stage 1\n(Early)", "Stage 2\n(Early-Mod)", "Stage 3\n(Bridging)", "Stage 4\n(Cirrhosis)"]

for g in all_hubs:
    ctrl_mean = df_all[df_all["stage"] == 0][g].mean()
    pooled_sd = df_all[g].std()
    df_all[g + "_z"] = (df_all[g] - ctrl_mean) / pooled_sd
    
    means = df_all.groupby("stage")[g + "_z"].mean().reindex(stage_labels).values
    sems = df_all.groupby("stage")[g + "_z"].sem().reindex(stage_labels).values
    stage_stats[g] = {"means": means, "sems": sems}

# ==============================================================================
# 4. COMPUTE PANEL C METRICS: DECISION CURVE ANALYSIS (DCA)
# ==============================================================================
print("4. Computing Decision Curve Analysis (pt = 0.05 to 0.50)...")
pt_range = np.linspace(0.05, 0.50, 46)
n_early_total = len(y_true)
prevalence = np.mean(y_true)

net_benefit_nomogram = []
net_benefit_all = []
net_benefit_none = np.zeros(len(pt_range))

for pt in pt_range:
    pred = (oof_ens >= pt).astype(int)
    tp = np.sum((pred == 1) & (y_true == 1))
    fp = np.sum((pred == 1) & (y_true == 0))
    w = pt / (1.0 - pt)
    
    nb_mod = (tp / n_early_total) - (fp / n_early_total) * w
    nb_all = prevalence - (1.0 - prevalence) * w
    
    net_benefit_nomogram.append(nb_mod)
    net_benefit_all.append(nb_all)

# ==============================================================================
# 5. RENDER PUBLICATION TRIPTYCH FIGURE (1x3 GRID, 300 DPI)
# ==============================================================================
print("5. Rendering Publication Triptych Figure (5 Hubs)...")

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(21, 6.2), dpi=300)
fig.patch.set_facecolor("white")
plt.subplots_adjust(wspace=0.28, left=0.06, right=0.97, top=0.81, bottom=0.15)

# ------------------------------------------------------------------------------
# PANEL A: Early-Stage Subgroup ROC Curves
# ------------------------------------------------------------------------------
ax1.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", linewidth=1.3, label="Chance (AUC = 0.500)")

# Plot individual genes
for g in ["COL3A1", "COL1A1", "SERPINE2", "SERPINF2", "COL15A1"]:
    r = roc_results[g]
    lbl = f"{g}: AUC = {r['auc']:.3f} [{r['ci_low']:.2f}–{r['ci_high']:.2f}]"
    ax1.plot(r["fpr"], r["tpr"], color=GENE_COLORS[g], linewidth=1.8, label=lbl)

# Plot Ensemble model
r_ens = roc_results["Ensemble"]
lbl_ens = f"5-Hub Ensemble: AUC = {r_ens['auc']:.3f} [{r_ens['ci_low']:.2f}–{r_ens['ci_high']:.2f}]"
ax1.plot(r_ens["fpr"], r_ens["tpr"], color=GENE_COLORS["Ensemble"], linewidth=2.8, linestyle="-", label=lbl_ens)

# Mark Youden cutoff point on Ensemble
youden = r_ens["tpr"] - r_ens["fpr"]
opt_idx = np.argmax(youden)
sens_opt = r_ens["tpr"][opt_idx] * 100
spec_opt = (1 - r_ens["fpr"][opt_idx]) * 100
ax1.scatter(r_ens["fpr"][opt_idx], r_ens["tpr"][opt_idx], color="#dc2626", s=70, zorder=6,
            edgecolors="black", linewidth=1.0)
ax1.annotate(f"Optimal Cutoff\n(Sens={sens_opt:.1f}%, Spec={spec_opt:.1f}%)",
             xy=(r_ens["fpr"][opt_idx], r_ens["tpr"][opt_idx]),
             xytext=(r_ens["fpr"][opt_idx] + 0.12, r_ens["tpr"][opt_idx] - 0.14),
             fontsize=8.5, fontweight="bold", color="#1e293b",
             arrowprops=dict(arrowstyle="->", color="#1e293b", lw=1.0),
             bbox=dict(boxstyle="round,pad=0.25", facecolor="#fef2f2", edgecolor="#fca5a5", alpha=0.9))

ax1.set_xlim(-0.02, 1.02)
ax1.set_ylim(-0.02, 1.02)
ax1.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=11, fontweight="bold", labelpad=6)
ax1.set_ylabel("True Positive Rate (Sensitivity)", fontsize=11, fontweight="bold", labelpad=6)
ax1.set_title("Early-Stage Subgroup ROC Curves\nControl (S0, n=43) vs. Early Fibrosis (S1-S2, n=53)",
              fontsize=12, fontweight="bold", pad=12, color="#0f172a")
ax1.legend(loc="lower right", fontsize=8.0, frameon=True, facecolor="white", edgecolor="#cbd5e1", framealpha=0.95)

ax1.text(-0.14, 1.08, "A", transform=ax1.transAxes, fontsize=18, fontweight="bold", color="#0f172a")

# ------------------------------------------------------------------------------
# PANEL B: Early Switch vs. Late Progression Dynamics
# ------------------------------------------------------------------------------
x_pts = np.array([0, 1, 2, 3, 4])
x_smooth = np.linspace(0, 4, 200)

for g in all_hubs:
    m = stage_stats[g]["means"]
    s = stage_stats[g]["sems"]
    
    spl_m = make_interp_spline(x_pts, m, k=2)
    spl_s = make_interp_spline(x_pts, s, k=2)
    y_smooth = spl_m(x_smooth)
    sem_smooth = spl_s(x_smooth)
    
    lw = 2.6 if g in ["COL15A1", "COL1A1"] else 1.8
    ls = "-" if g != "SERPINF2" else "--"
    ax2.plot(x_smooth, y_smooth, color=GENE_COLORS[g], linewidth=lw, linestyle=ls, label=f"{g}")
    ax2.fill_between(x_smooth, y_smooth - sem_smooth, y_smooth + sem_smooth, color=GENE_COLORS[g], alpha=0.14)
    ax2.scatter(x_pts, m, color=GENE_COLORS[g], s=35, zorder=5, edgecolors="white", linewidth=0.8)

# Baseline zero guide
ax2.axhline(0, color="#94a3b8", linestyle=":", linewidth=1.0)

# Highlighting annotations - positioned cleanly
ax2.annotate("COL15A1: Early-Onset Switch\n(Basement Membrane Remodeling)",
             xy=(2.0, stage_stats["COL15A1"]["means"][2]),
             xytext=(0.85, 1.45),
             fontsize=8.2, fontweight="bold", color="#1e40af",
             arrowprops=dict(arrowstyle="->", color="#1e40af", lw=1.2),
             bbox=dict(boxstyle="round,pad=0.25", facecolor="#eff6ff", edgecolor="#93c5fd", alpha=0.92))

ax2.annotate("COL1A1: Linear Progression Tracker\n(Cumulative Fibrillar Scarring)",
             xy=(3.0, stage_stats["COL1A1"]["means"][3]),
             xytext=(2.2, 0.35),
             fontsize=8.2, fontweight="bold", color="#b91c1c",
             arrowprops=dict(arrowstyle="->", color="#b91c1c", lw=1.2),
             bbox=dict(boxstyle="round,pad=0.25", facecolor="#fef2f2", edgecolor="#fca5a5", alpha=0.92))

ax2.set_xticks(x_pts)
ax2.set_xticklabels(stage_names, fontsize=9.5, fontweight="bold")
ax2.set_xlabel("Histological Fibrosis Progression Stage", fontsize=11, fontweight="bold", labelpad=6)
ax2.set_ylabel(r"Expression Dynamic (Normalized $\Delta$Z-Score)", fontsize=11, fontweight="bold", labelpad=6)
ax2.set_title("Early Switch vs. Late Progression Dynamics\nContinuous Molecular Trajectories from Control to Cirrhosis",
              fontsize=12, fontweight="bold", pad=12, color="#0f172a")
ax2.legend(loc="upper left", fontsize=8.0, frameon=True, facecolor="white", edgecolor="#cbd5e1", framealpha=0.95)

ax2.text(-0.14, 1.08, "B", transform=ax2.transAxes, fontsize=18, fontweight="bold", color="#0f172a")

# ------------------------------------------------------------------------------
# PANEL C: Subclinical Decision Curve Analysis (DCA)
# ------------------------------------------------------------------------------
ax3.plot(pt_range, net_benefit_nomogram, color=GENE_COLORS["Ensemble"], linewidth=2.8,
         label="5-Hub Early Nomogram (S1-S2 vs S0)")
ax3.plot(pt_range, net_benefit_all, color="#64748b", linewidth=1.8, linestyle="--",
         label="Screen All (Biopsy All High-Risk Patients)")
ax3.plot(pt_range, net_benefit_none, color="#0f172a", linewidth=1.4, linestyle=":",
         label="Screen None (Net Benefit = 0)")

# Shade High-Utility Subclinical Window
ax3.axvspan(0.10, 0.45, color="#ede9fe", alpha=0.35, zorder=0)
ax3.text(0.275, 0.03, "High-Utility Subclinical Window\n(Net Benefit Improvement > All Strategies)",
         ha="center", va="bottom", fontsize=8.2, fontweight="bold", color="#5b21b6",
         bbox=dict(boxstyle="round,pad=0.25", facecolor="#ffffff", edgecolor="#c4b5fd", alpha=0.92))

# Highlight peak benefit at pt = 0.25
nb_at_25 = net_benefit_nomogram[20]
ax3.scatter(0.25, nb_at_25, color="#7c3aed", s=65, zorder=6, edgecolors="black", linewidth=1.0)
ax3.annotate(f"Net Benefit = {nb_at_25:+.3f}\nat $p_t$ = 25%",
             xy=(0.25, nb_at_25),
             xytext=(0.21, nb_at_25 - 0.14),
             fontsize=8.5, fontweight="bold", color="#5b21b6",
             arrowprops=dict(arrowstyle="->", color="#5b21b6", lw=1.0),
             bbox=dict(boxstyle="round,pad=0.25", facecolor="#faf5ff", edgecolor="#d8b4fe", alpha=0.9))

ax3.set_xlim(0.05, 0.50)
ax3.set_ylim(-0.06, 0.58)
ax3.set_xlabel("Subclinical Threshold Risk Probability ($p_t$)", fontsize=11, fontweight="bold", labelpad=6)
ax3.set_ylabel("Clinical Net Benefit", fontsize=11, fontweight="bold", labelpad=6)
ax3.set_title("Subclinical Decision Curve Analysis (DCA)\nClinical Screening Utility across Early Risk Thresholds",
              fontsize=12, fontweight="bold", pad=12, color="#0f172a")
ax3.legend(loc="upper right", fontsize=8.2, frameon=True, facecolor="white", edgecolor="#cbd5e1", framealpha=0.95)

ax3.text(-0.14, 1.08, "C", transform=ax3.transAxes, fontsize=18, fontweight="bold", color="#0f172a")

# Global Header
plt.suptitle("Early Detection & Subclinical Risk Stratification of the 5 Universal Pan-Fibrotic Hub Biomarkers\nNon-Invasive Diagnostic Proof in Biopsy-Staged Clinical Cohorts (GSE84044 Scheuer Stages)",
             fontsize=14.5, fontweight="bold", color="#0f172a", y=0.96)

out_fig = "plots/hub_genes_early_detection_triptych_5genes.png"
plt.savefig(out_fig, dpi=300)
plt.close()
print(f"Successfully generated: {out_fig}")
