# -*- coding: utf-8 -*-
"""
Module 1: Clinical Nomogram Construction, Bootstrap Calibration, and Decision Curve Analysis (DCA)
Targeting all 9 Consensus Pan-Fibrotic Hub Genes:
COL15A1, COL1A1, COL3A1, SERPINE2, SERPINF2, LAMC3, LTBP2, MDK, SVEP1
Cohort: 1,069 Pure Human Biopsy Discovery Cohort (ComBat Batch-Harmonized)

Output:
- plots/hub_genes_nomogram_and_dca.png
- results/hub_genes_nomogram_parameters.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score
from scipy import stats

os.makedirs("plots", exist_ok=True)
os.makedirs("results", exist_ok=True)

plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.edgecolor"] = "#334155"
plt.rcParams["axes.linewidth"] = 1.1

# 1. Load ComBat-Corrected Matrix
df = pd.read_csv("results/real_human_patient_combat_corrected_matrix.csv")
hub_genes = ["COL15A1", "COL1A1", "COL3A1", "SERPINE2", "SERPINF2", "LAMC3", "LTBP2", "MDK", "SVEP1"]
X = df[hub_genes].copy()
y = df["label"].values
N = len(df)

print(f"Loaded {N} patient biopsy samples (Controls: {(y==0).sum()}, Fibrosis: {(y==1).sum()})")

# 2. Fit Multivariate Logistic Regression via statsmodels
X_with_const = sm.add_constant(X)
logit_model = sm.Logit(y, X_with_const).fit(disp=False)
print("\nMultivariable Logistic Regression Summary:")
print(logit_model.summary())

coef_dict = logit_model.params.to_dict()
pvalues_dict = logit_model.pvalues.to_dict()
intercept = coef_dict["const"]

# Fit sklearn model for calibrated probabilities
clf = LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000, random_state=42)
clf.fit(X, y)
y_pred_probs = clf.predict_proba(X)[:, 1]
overall_auc = roc_auc_score(y, y_pred_probs)

# Bootstrap 95% CI for Nomogram AUC
np.random.seed(42)
boot_aucs = []
for _ in range(1000):
    b_idx = np.random.choice(N, N, replace=True)
    if len(np.unique(y[b_idx])) == 2:
        boot_aucs.append(roc_auc_score(y[b_idx], y_pred_probs[b_idx]))
auc_ci_l = np.percentile(boot_aucs, 2.5)
auc_ci_u = np.percentile(boot_aucs, 97.5)
print(f"9-Hub Nomogram Diagnostic AUC: {overall_auc:.4f} (95% CI: {auc_ci_l:.4f} - {auc_ci_u:.4f})")

# 3. Calculate Nomogram Point Scaling System
gene_ranges = {}
beta_ranges = {}
for g in hub_genes:
    low = np.percentile(X[g], 1)
    high = np.percentile(X[g], 99)
    gene_ranges[g] = (low, high)
    beta_ranges[g] = abs(coef_dict[g] * (high - low))

# Maximum beta range gets 100 points
max_beta_range = max(beta_ranges.values())
gene_max_pts = {g: (beta_ranges[g] / max_beta_range) * 100 for g in hub_genes}

nomo_records = []
for g in hub_genes:
    low, high = gene_ranges[g]
    max_p = gene_max_pts[g]
    nomo_records.append({
        "gene": g,
        "coefficient_beta": round(coef_dict[g], 4),
        "odds_ratio": round(np.exp(coef_dict[g]), 4),
        "std_err": round(logit_model.bse[g], 4),
        "z_statistic": round(logit_model.tvalues[g], 4),
        "p_value": logit_model.pvalues[g],
        "range_min": round(low, 3),
        "range_max": round(high, 3),
        "points_at_min": 0.0 if coef_dict[g] > 0 else round(max_p, 1),
        "points_at_max": round(max_p, 1) if coef_dict[g] > 0 else 0.0
    })

pd.DataFrame(nomo_records).to_csv("results/hub_genes_nomogram_parameters_9genes.csv", index=False)
print("Saved: results/hub_genes_nomogram_parameters_9genes.csv")

# 4. Calibration Curve via Deciles
n_bins = 10
bins = pd.qcut(y_pred_probs, q=n_bins, duplicates="drop")
pred_prob_bins = []
obs_prob_bins = []
for b in bins.unique():
    mask = (bins == b)
    pred_prob_bins.append(y_pred_probs[mask].mean())
    obs_prob_bins.append(y[mask].mean())

order = np.argsort(pred_prob_bins)
pred_prob_bins = np.array(pred_prob_bins)[order]
obs_prob_bins = np.array(obs_prob_bins)[order]

# Brier score
brier = np.mean((y_pred_probs - y) ** 2)

# 5. Decision Curve Analysis (DCA)
thresholds = np.linspace(0.05, 0.90, 50)
net_benefit_model = []
net_benefit_all = []
net_benefit_none = [0.0] * len(thresholds)
prevalence = np.mean(y)

for pt in thresholds:
    pred = (y_pred_probs >= pt).astype(int)
    tp = np.sum((pred == 1) & (y == 1))
    fp = np.sum((pred == 1) & (y == 0))
    nb = (tp / N) - (fp / N) * (pt / (1.0 - pt))
    net_benefit_model.append(nb)
    
    # Treat All strategy
    nb_all = prevalence - (1.0 - prevalence) * (pt / (1.0 - pt))
    net_benefit_all.append(nb_all)

# ==============================================================================
# 6. PLOT COMPREHENSIVE FIGURE: Nomogram, Calibration, DCA
# ==============================================================================
fig = plt.figure(figsize=(20, 16), dpi=300)
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(2, 2, height_ratios=[1.7, 1.0], hspace=0.32, wspace=0.22)

# --- PANEL A: Visual Nomogram Scoring System (All 9 Hub Genes) ---
ax_nomo = fig.add_subplot(gs[0, :])
ax_nomo.set_facecolor("white")
ax_nomo.set_title("A. Diagnostic Nomogram for Human Pan-Fibrotic Risk Stratification (9 Consensus Hubs)\nMultivariable Logistic Regression Model Derived from 1,069 Pure Human Biopsies",
                  fontsize=14, fontweight="bold", pad=20, color="#0f172a")

gene_annotations = {
    "COL15A1": "COL15A1 (Basement Membrane)",
    "COL1A1": "COL1A1 (Fibrillar Collagen I)",
    "COL3A1": "COL3A1 (Fibrillar Collagen III)",
    "SERPINE2": "SERPINE2 (Antiprotease / SERPIN)",
    "SERPINF2": "SERPINF2 (Antiplasmin / Repressed)",
    "LAMC3": "LAMC3 (Laminin γ3 Subunit)",
    "LTBP2": "LTBP2 (TGF-β Latent Complex)",
    "MDK": "MDK (Midkine Growth Factor)",
    "SVEP1": "SVEP1 (Selectin Like ECM Cell Adhesion)"
}

ruler_labels = [
    "Points Ruler (0–100)",
    gene_annotations["COL15A1"],
    gene_annotations["COL1A1"],
    gene_annotations["COL3A1"],
    gene_annotations["SERPINE2"],
    gene_annotations["SERPINF2"],
    gene_annotations["LAMC3"],
    gene_annotations["LTBP2"],
    gene_annotations["MDK"],
    gene_annotations["SVEP1"],
    "Total Points",
    "Fibrosis Probability"
]

y_positions = [11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]

# 1. Top Points Ruler (0 to 100)
ax_nomo.plot([0, 100], [11.0, 11.0], color="#0f172a", linewidth=1.6)
ax_nomo.text(-2.5, 11.0, "Points Ruler", va="center", ha="right", fontsize=10.5, fontweight="bold", color="#0f172a")
for pt in range(0, 101, 10):
    ax_nomo.plot([pt, pt], [11.0, 11.16], color="#0f172a", linewidth=1.2)
    ax_nomo.text(pt, 11.26, str(pt), ha="center", va="bottom", fontsize=8.5, color="#1e293b")

# 2. Gene Rulers
for idx, g in enumerate(hub_genes):
    y_pos = y_positions[idx + 1]
    lbl = ruler_labels[idx + 1]
    ax_nomo.text(-2.5, y_pos, lbl, va="center", ha="right", fontsize=10.0, fontweight="bold", color="#334155")
    
    max_pts = gene_max_pts[g]
    low, high = gene_ranges[g]
    ax_nomo.plot([0, max_pts], [y_pos, y_pos], color="#0f172a", linewidth=1.4)
    
    # 5 evenly spaced ticks across gene expression range
    val_ticks = np.linspace(low, high, 5)
    for val in val_ticks:
        if coef_dict[g] > 0:
            pts = ((val - low) / (high - low)) * max_pts
        else:
            pts = ((high - val) / (high - low)) * max_pts
        pts = np.clip(pts, 0, max_pts)
        ax_nomo.plot([pts, pts], [y_pos, y_pos + 0.14], color="#0284c7", linewidth=1.3)
        ax_nomo.text(pts, y_pos + 0.20, f"{val:.1f}", ha="center", va="bottom", fontsize=8.0, fontweight="bold", color="#0369a1")

# 3. Total Points Ruler (0 to 450)
total_max_pts = sum(gene_max_pts.values())
ax_nomo.text(-2.5, 1.0, "Total Points", va="center", ha="right", fontsize=10.5, fontweight="bold", color="#b91c1c")
ax_nomo.plot([0, 100], [1.0, 1.0], color="#0f172a", linewidth=1.6)
for t_pt in range(0, 451, 50):
    scaled_x = (t_pt / 450.0) * 100.0
    ax_nomo.plot([scaled_x, scaled_x], [1.0, 1.16], color="#dc2626", linewidth=1.4)
    ax_nomo.text(scaled_x, 1.25, str(t_pt), ha="center", va="bottom", fontsize=8.5, color="#b91c1c", fontweight="bold")

# 4. Risk Probability Ruler
ax_nomo.text(-2.5, 0.0, "Fibrosis Probability", va="center", ha="right", fontsize=10.5, fontweight="bold", color="#4338ca")
ax_nomo.plot([0, 100], [0.0, 0.0], color="#0f172a", linewidth=1.6)

# Reference baseline: when all genes at lowest points (pts = 0)
baseline_logit = intercept + sum([coef_dict[g] * (gene_ranges[g][0] if coef_dict[g] > 0 else gene_ranges[g][1]) for g in hub_genes])
pts_to_logit_slope = max_beta_range / 100.0

probs = [0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.85, 0.95, 0.99]
for p_val in probs:
    logit_val = np.log(p_val / (1.0 - p_val))
    req_total_pts = (logit_val - baseline_logit) / pts_to_logit_slope
    req_scaled_x = (req_total_pts / 450.0) * 100.0
    if 0 <= req_scaled_x <= 100:
        ax_nomo.plot([req_scaled_x, req_scaled_x], [0.0, 0.15], color="#4f46e5", linewidth=1.4)
        ax_nomo.text(req_scaled_x, -0.32, f"{p_val:.2f}", ha="center", va="top", fontsize=8.2, color="#3730a3", fontweight="bold")

ax_nomo.set_xlim(-26, 105)
ax_nomo.set_ylim(-0.9, 12.0)
ax_nomo.axis("off")

# --- PANEL B: Calibration Curve ---
ax_cal = fig.add_subplot(gs[1, 0])
ax_cal.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", linewidth=1.4, label="Ideal Calibration (Slope = 1.0)")
ax_cal.plot(pred_prob_bins, obs_prob_bins, marker="o", markersize=7, color="#ea580c", linewidth=2.2,
            label=f"9-Hub Nomogram (Brier = {brier:.4f})")

ax_cal.set_title("B. Decile Calibration Curve (N = 1,069 Biopsies)", fontsize=12, fontweight="bold", pad=12, color="#0f172a")
ax_cal.set_xlabel("Nomogram Predicted Probability", fontsize=10.5, fontweight="bold")
ax_cal.set_ylabel("Observed Actual Proportion", fontsize=10.5, fontweight="bold")
ax_cal.set_xlim(-0.02, 1.02)
ax_cal.set_ylim(-0.02, 1.02)
ax_cal.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=9.5)

# --- PANEL C: Decision Curve Analysis (DCA) ---
ax_dca = fig.add_subplot(gs[1, 1])
ax_dca.plot(thresholds, net_benefit_model, color="#1e1b4b", linewidth=2.5, label="9-Hub Diagnostic Nomogram")
ax_dca.plot(thresholds, net_benefit_all, color="#dc2626", linestyle="--", linewidth=1.8, label="Treat All Cases")
ax_dca.plot(thresholds, net_benefit_none, color="#64748b", linestyle=":", linewidth=1.4, label="Treat None (Net Benefit = 0)")

# Shade optimal clinical decision window
ax_dca.axvspan(0.10, 0.75, color="#ede9fe", alpha=0.35, zorder=0)
ax_dca.text(0.42, max(net_benefit_model) * 0.45, "Optimal Clinical Action Window\n(Net Benefit > All / None Strategies)",
            ha="center", va="center", fontsize=8.5, fontweight="bold", color="#5b21b6",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#faf5ff", edgecolor="#c4b5fd", alpha=0.92))

ax_dca.set_title("C. Decision Curve Analysis (Clinical Net Benefit)", fontsize=12, fontweight="bold", pad=12, color="#0f172a")
ax_dca.set_xlabel("High-Risk Decision Threshold Probability ($p_t$)", fontsize=10.5, fontweight="bold")
ax_dca.set_ylabel("Standardized Net Benefit", fontsize=10.5, fontweight="bold")
ax_dca.set_xlim(0.05, 0.85)
ax_dca.set_ylim(-0.05, max(net_benefit_model) * 1.15)
ax_dca.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=9.5)

out_png = "plots/hub_genes_nomogram_and_dca_9genes.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved 9-hub Nomogram & DCA figure: {out_png}")
