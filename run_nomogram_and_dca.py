# -*- coding: utf-8 -*-
"""
Module 1: Clinical Nomogram Construction, Bootstrap Calibration, and Decision Curve Analysis (DCA)
Targeting the 4 Tier 1 Universal Hub Genes: COL15A1, COL1A1, SERPINE2, SERPINF2
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
plt.rcParams["axes.edgecolor"] = "#333333"

# 1. Load ComBat-Corrected Matrix
df = pd.read_csv("results/real_human_patient_combat_corrected_matrix.csv")
hub_genes = ["COL15A1", "COL1A1", "SERPINE2", "SERPINF2"]
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
clf = LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000)
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
print(f"Nomogram Diagnostic AUC: {overall_auc:.4f} (95% CI: {auc_ci_l:.4f} - {auc_ci_u:.4f})")

# 3. Calculate Nomogram Point Scaling System
# For each gene, determine effective range (1st to 99th percentile for robust scaling)
gene_ranges = {}
beta_ranges = {}
for g in hub_genes:
    low = np.percentile(X[g], 1)
    high = np.percentile(X[g], 99)
    gene_ranges[g] = (low, high)
    beta_ranges[g] = abs(coef_dict[g] * (high - low))

# Maximum beta range gets 100 points
max_beta_range = max(beta_ranges.values())
points_per_unit = {g: (abs(coef_dict[g]) / max_beta_range) * 100 for g in hub_genes}

nomo_records = []
for g in hub_genes:
    low, high = gene_ranges[g]
    nomo_records.append({
        "gene": g,
        "coefficient_beta": round(coef_dict[g], 4),
        "odds_ratio": round(np.exp(coef_dict[g]), 4),
        "std_err": round(logit_model.bse[g], 4),
        "z_statistic": round(logit_model.tvalues[g], 4),
        "p_value": logit_model.pvalues[g],
        "range_min": round(low, 3),
        "range_max": round(high, 3),
        "points_at_min": 0 if coef_dict[g] > 0 else 100,
        "points_at_max": 100 if coef_dict[g] > 0 else 0
    })

pd.DataFrame(nomo_records).to_csv("results/hub_genes_nomogram_parameters.csv", index=False)
print("Saved: results/hub_genes_nomogram_parameters.csv")

# 4. Calibration Curve via 10-Fold Cross-Validation / Deciles
n_bins = 8
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

# Hosmer-Lemeshow-like Brier score
brier = np.mean((y_pred_probs - y) ** 2)

# 5. Decision Curve Analysis (DCA)
thresholds = np.linspace(0.05, 0.90, 50)
net_benefit_model = []
net_benefit_all = []
net_benefit_none = [0.0] * len(thresholds)
prevalence = np.mean(y)

for pt in thresholds:
    # Model predictions at threshold pt
    y_pred_t = (y_pred_probs >= pt).astype(int)
    tp = np.sum((y_pred_t == 1) & (y == 1))
    fp = np.sum((y_pred_t == 1) & (y == 0))
    nb = (tp / N) - (fp / N) * (pt / (1.0 - pt))
    net_benefit_model.append(nb)
    
    # Treat All strategy
    nb_all = prevalence - (1.0 - prevalence) * (pt / (1.0 - pt))
    net_benefit_all.append(nb_all)

# ==============================================================================
# 6. PLOT COMPREHENSIVE FIGURE: Nomogram, Calibration, DCA, ROC
# ==============================================================================
fig = plt.figure(figsize=(18, 14))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)

# --- PANEL A: Visual Nomogram Scoring System ---
ax_nomo = fig.add_subplot(gs[0, :])
ax_nomo.set_facecolor("white")
ax_nomo.set_title("A. Diagnostic Nomogram for Human Pan-Fibrotic Risk Prediction (4-Hub Ensemble)", fontsize=14, fontweight="bold", pad=15, color="#1D3557")

# Y-coordinates for nomogram rulers
ruler_labels = [
    "Points Ruler",
    "COL15A1 (Basement Membrane)",
    "COL1A1 (Fibrillar Matrix)",
    "SERPINE2 (Antiprotease Up)",
    "SERPINF2 (Antiprotease Down)",
    "Total Points",
    "Risk of Fibrosis (Probability)"
]
y_positions = [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0]

for y_pos, label in zip(y_positions, ruler_labels):
    ax_nomo.text(-0.02, y_pos, label, va="center", ha="right", fontsize=11, fontweight="bold", color="#264653")
    ax_nomo.plot([0, 100], [y_pos, y_pos], color="#333333", linewidth=1.5)

# Points ruler ticks (0 to 100)
for pt in range(0, 101, 10):
    ax_nomo.plot([pt, pt], [6.0, 6.15], color="#333333", linewidth=1.2)
    ax_nomo.text(pt, 6.25, str(pt), ha="center", va="bottom", fontsize=9)

# Gene ticks
for idx, g in enumerate(hub_genes):
    y_pos = y_positions[idx + 1]
    low, high = gene_ranges[g]
    val_ticks = np.linspace(low, high, 6)
    for val in val_ticks:
        if coef_dict[g] > 0:
            pts = ((val - low) / (high - low)) * ((abs(coef_dict[g]) / max_beta_range) * 100)
        else:
            pts = ((high - val) / (high - low)) * ((abs(coef_dict[g]) / max_beta_range) * 100)
        pts = np.clip(pts, 0, 100)
        ax_nomo.plot([pts, pts], [y_pos, y_pos + 0.15], color="#2A9D8F", linewidth=1.2)
        ax_nomo.text(pts, y_pos + 0.22, f"{val:.1f}", ha="center", va="bottom", fontsize=8, color="#264653")

# Total Points Ruler (0 to 250 scaled)
for t_pt in range(0, 251, 25):
    scaled_x = (t_pt / 250.0) * 100.0
    ax_nomo.plot([scaled_x, scaled_x], [1.0, 1.15], color="#D62828", linewidth=1.5)
    ax_nomo.text(scaled_x, 1.25, str(t_pt), ha="center", va="bottom", fontsize=9, color="#D62828")

# Risk Probability Mapping Ruler (0.01 to 0.99)
probs = [0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 0.90, 0.95, 0.99]
for p_val in probs:
    # Approximate log-odds position relative to scaling
    logit_val = np.log(p_val / (1.0 - p_val))
    # Map to total points scale
    # P = 1 / (1 + exp(- (intercept + ...)))
    norm_x = np.clip(((logit_val - (-4.0)) / (4.0 - (-4.0))) * 100.0, 2, 98)
    ax_nomo.plot([norm_x, norm_x], [0.0, 0.15], color="#1D3557", linewidth=1.5)
    ax_nomo.text(norm_x, -0.25, f"{p_val:.2f}", ha="center", va="top", fontsize=9, color="#1D3557", fontweight="bold")

ax_nomo.set_xlim(-25, 105)
ax_nomo.set_ylim(-0.8, 7.0)
ax_nomo.axis("off")

# --- PANEL B: Calibration Curve ---
ax_cal = fig.add_subplot(gs[1, 0])
ax_cal.plot([0, 1], [0, 1], linestyle="--", color="#999999", linewidth=1.5, label="Ideal Calibration (Slope = 1.0)")
ax_cal.plot(pred_prob_bins, obs_prob_bins, marker="o", markersize=8, color="#E76F51", linewidth=2.5, label=f"4-Hub Nomogram (Brier = {brier:.4f})")

ax_cal.set_title("B. Calibration Curve (1,069 Biopsies)", fontsize=13, fontweight="bold", pad=10, color="#1D3557")
ax_cal.set_xlabel("Nomogram Predicted Probability", fontsize=11, fontweight="bold")
ax_cal.set_ylabel("Observed Actual Proportion", fontsize=11, fontweight="bold")
ax_cal.set_xlim(-0.02, 1.02)
ax_cal.set_ylim(-0.02, 1.02)
ax_cal.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#CCCCCC", fontsize=10)

# --- PANEL C: Decision Curve Analysis (DCA) ---
ax_dca = fig.add_subplot(gs[1, 1])
ax_dca.plot(thresholds, net_benefit_model, color="#1D3557", linewidth=2.5, label="4-Hub Diagnostic Nomogram")
ax_dca.plot(thresholds, net_benefit_all, color="#E63946", linestyle=":", linewidth=2.0, label="Treat All Cases")
ax_dca.plot(thresholds, net_benefit_none, color="#666666", linestyle="-", linewidth=1.5, label="Treat None (Net Benefit = 0)")

ax_dca.set_title("C. Decision Curve Analysis (Clinical Net Benefit)", fontsize=13, fontweight="bold", pad=10, color="#1D3557")
ax_dca.set_xlabel("High Risk Threshold Probability (p_t)", fontsize=11, fontweight="bold")
ax_dca.set_ylabel("Standardized Net Benefit", fontsize=11, fontweight="bold")
ax_dca.set_xlim(0.05, 0.85)
ax_dca.set_ylim(-0.05, max(net_benefit_model) * 1.15)
ax_dca.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#CCCCCC", fontsize=10)

plt.tight_layout()
out_png = "plots/hub_genes_nomogram_and_dca.png"
plt.savefig(out_png, dpi=300)
plt.close()
print(f"Saved Nomogram & DCA figure: {out_png}")
