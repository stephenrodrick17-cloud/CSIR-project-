# -*- coding: utf-8 -*-
import os
import gzip
import glob
import json
import pandas as pd
import numpy as np
from pycombat import Combat
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.utils.class_weight import compute_sample_weight
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import discovery_config

# 1. Enforce Cohort Isolation Guard
discovery_config.verify_cohort_isolation()

base_dir = r"D:\CSIR"
results_dir = os.path.join(base_dir, "results")
plots_dir = os.path.join(base_dir, "plots")
geo_cache = os.path.join(base_dir, "geo_cache")
os.makedirs(results_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

# 2. Target 24 Clean Core ECM Genes
ecm_df = pd.read_csv(os.path.join(base_dir, "ecm_clean_genes.csv"))
ecm_genes = sorted(ecm_df["gene"].tolist())
print(f"Target Feature Set ({len(ecm_genes)} Clean Core ECM Genes): {ecm_genes}\n")

# 3. Load Real Patient Expression Matrix (799 Samples across 8 Studies)
raw_df = pd.read_csv(os.path.join(results_dir, "real_human_patient_ml_training_matrix.csv"))
print(f"Loaded Real Patient Matrix: {len(raw_df)} samples across {raw_df['study'].nunique()} studies.")

# 4. Missing Value & Pre-Scaling Variance Audit
print("\n" + "=" * 80)
print("MISSING VALUE & PRE-SCALING VARIANCE AUDIT (FGF14, MDK vs COL1A1, AEBP1)")
print("=" * 80)
for g in ["FGF14", "MDK", "COL1A1", "AEBP1"]:
    n_miss = raw_df[g].isna().sum()
    print(f"Gene: {g} | Missing: {n_miss}/{len(raw_df)} | Overall Mean: {raw_df[g].mean():.3f} | Var: {raw_df[g].var():.3f}")
    grp = raw_df.groupby("study")[g].agg(["count", "mean", "var"]).round(3)
    print(grp.to_string())
    print()

# 5. Impute probe gaps using study-level median
X_raw = raw_df[ecm_genes].copy()
for col in X_raw.columns:
    if X_raw[col].isna().any():
        X_raw[col] = X_raw[col].fillna(raw_df.groupby("study")[col].transform("median")).fillna(X_raw[col].median())

# 6. Apply ComBat Batch Correction across 8 GEO Studies
print("\nApplying ComBat Batch Correction across 8 GEO studies...")
combat = Combat()
X_combat_vals = combat.fit_transform(X_raw.values, raw_df["study"].values)
X_combat = pd.DataFrame(X_combat_vals, columns=ecm_genes, index=raw_df.index)

# Save ComBat-corrected matrix
combat_df = X_combat.copy()
combat_df["GSM"] = raw_df["GSM"]
combat_df["study"] = raw_df["study"]
combat_df["organ"] = raw_df["organ"]
combat_df["label"] = raw_df["label"]
combat_df.to_csv(os.path.join(results_dir, "real_human_patient_combat_corrected_matrix.csv"), index=False)
print(f"Saved ComBat-corrected matrix to {os.path.join(results_dir, 'real_human_patient_combat_corrected_matrix.csv')}")

# 7. Multi-Organ Stratified Balancing
# Equalize organ sample contributions so no organ dominates the fit
np.random.seed(42)
balanced_indices = []
for org in ["Kidney", "Liver", "Lungs", "Skin"]:
    org_df = raw_df[raw_df["organ"] == org]
    ctrl_idx = org_df[org_df["label"] == 0].index
    dis_idx = org_df[org_df["label"] == 1].index
    
    n_ctrl = min(len(ctrl_idx), 18)
    n_dis = min(len(dis_idx), 42)
    
    c_chosen = np.random.choice(ctrl_idx, n_ctrl, replace=False)
    d_chosen = np.random.choice(dis_idx, n_dis, replace=False)
    balanced_indices.extend(list(c_chosen) + list(d_chosen))

df_bal = combat_df.loc[balanced_indices].copy()
X_bal = X_combat.loc[balanced_indices].copy()
y_bal = df_bal["label"].values

print(f"\nBalanced 4-Organ Cohort constructed: {len(df_bal)} samples")
print(df_bal["organ"].value_counts().to_string())
n_ctrl = int((df_bal["label"] == 0).sum())
n_dis = int((df_bal["label"] == 1).sum())
print(f"Disease breakdown: {n_ctrl} Controls vs {n_dis} Disease cases")

# 8. Standardize Features
scaler = StandardScaler()
X_bal_scaled = pd.DataFrame(scaler.fit_transform(X_bal), columns=ecm_genes)

# --- Algorithm 1: LASSO (L1 Regularization, 10-fold CV) ---
print("\n[1/4] Running LASSO with 10-fold CV on ComBat-corrected balanced data...")
lasso = LogisticRegressionCV(
    Cs=100, cv=10, penalty="l1", solver="saga", scoring="roc_auc",
    class_weight="balanced", random_state=42, max_iter=5000, n_jobs=-1
)
lasso.fit(X_bal_scaled, y_bal)
lasso_coefs = pd.Series(lasso.coef_[0], index=ecm_genes)
lasso_selected = set(lasso_coefs[lasso_coefs.abs() > 1e-4].index)
print(f"  -> LASSO selected {len(lasso_selected)} / {len(ecm_genes)} features (optimal C = {lasso.C_[0]:.4f})")

# --- Algorithm 2: SVM-RFE (Linear Kernel, class-weighted) ---
print("[2/4] Running SVM-RFE on ComBat-corrected balanced data...")
svc_linear = SVC(kernel="linear", C=1.0, class_weight="balanced", random_state=42)
rfe = RFE(estimator=svc_linear, n_features_to_select=10, step=1)
rfe.fit(X_bal_scaled, y_bal)
svm_ranking = pd.Series(rfe.ranking_, index=ecm_genes)
svm_selected = set(svm_ranking[svm_ranking == 1].index)
print(f"  -> SVM-RFE selected {len(svm_selected)} top features")

# --- Algorithm 3: Random Forest (500 trees, class-weighted) ---
print("[3/4] Running Random Forest on ComBat-corrected balanced data...")
rf = RandomForestClassifier(n_estimators=500, max_depth=5, class_weight="balanced", random_state=42, n_jobs=-1)
rf.fit(X_bal_scaled, y_bal)
rf_importance = pd.Series(rf.feature_importances_, index=ecm_genes)
rf_threshold = 1.0 / len(ecm_genes)
rf_selected = set(rf_importance[rf_importance >= rf_threshold].index)
print(f"  -> Random Forest selected {len(rf_selected)} top features (importance >= {rf_threshold:.4f})")

# --- Algorithm 4: XGBoost (500 trees, class-weighted) ---
print("[4/4] Running XGBoost on ComBat-corrected balanced data...")
sample_weights = compute_sample_weight("balanced", y_bal)
xgb = XGBClassifier(n_estimators=500, max_depth=3, learning_rate=0.05, random_state=42, eval_metric="logloss")
xgb.fit(X_bal_scaled, y_bal, sample_weight=sample_weights)
xgb_importance = pd.Series(xgb.feature_importances_, index=ecm_genes)
xgb_threshold = 1.0 / len(ecm_genes)
xgb_selected = set(xgb_importance[xgb_importance >= xgb_threshold].index)
print(f"  -> XGBoost selected {len(xgb_selected)} dominant features")

# 9. Consensus Feature Table
results = []
for g in ecm_genes:
    v_lasso = int(g in lasso_selected)
    v_svm = int(g in svm_selected)
    v_rf = int(g in rf_selected)
    v_xgb = int(g in xgb_selected)
    total_votes = v_lasso + v_svm + v_rf + v_xgb
    
    auc_ind = roc_auc_score(y_bal, X_bal[g])
    if auc_ind < 0.5:
        auc_ind = 1.0 - auc_ind
        
    results.append({
        "Gene": g,
        "LASSO_Selected": v_lasso,
        "SVM_RFE_Selected": v_svm,
        "RandomForest_Selected": v_rf,
        "XGBoost_Selected": v_xgb,
        "Total_ML_Votes": total_votes,
        "Diagnostic_ROC_AUC": round(auc_ind, 4),
        "LASSO_Coef": round(lasso_coefs.get(g, 0.0), 4),
        "RF_Importance": round(rf_importance.get(g, 0.0), 4),
        "XGB_Importance": round(xgb_importance.get(g, 0.0), 4),
        "SVM_Rank": int(svm_ranking.get(g, 99))
    })

res_df = pd.DataFrame(results).sort_values(by=["Total_ML_Votes", "Diagnostic_ROC_AUC"], ascending=False)
res_df.to_csv(os.path.join(results_dir, "ml_4model_24ecm_hub_biomarkers.csv"), index=False)

print("\n" + "=" * 80)
print("FINAL COMBAT-CORRECTED 4-MODEL CONSENSUS FEATURE SELECTION")
print("=" * 80)
print(res_df.to_string(index=False))

# 10. Generate High-Resolution Consensus Plot
plt.figure(figsize=(14, 8))
colors = ["#1b4965" if v == 4 else "#2a9d8f" if v == 3 else "#e76f51" if v == 2 else "#adb5bd" for v in res_df["Total_ML_Votes"]]
bars = plt.barh(res_df["Gene"][::-1], res_df["Total_ML_Votes"][::-1], color=colors[::-1], edgecolor="black", linewidth=0.8)
plt.xlabel("ML Consensus Votes (Out of 4 Algorithms: LASSO, SVM-RFE, RF, XGBoost)", fontsize=12, fontweight="bold")
plt.title("Pan-Fibrotic Consensus Hub Biomarkers (ComBat-Corrected & Organ-Balanced Real Data)", fontsize=14, fontweight="bold", pad=15)
plt.xlim(0, 4.5)
plt.xticks([0, 1, 2, 3, 4], ["0", "1", "2", "3 (Consensus Hub)", "4 (Unanimous Core)"], fontsize=11)
plt.axvline(x=3, color="#e63946", linestyle="--", linewidth=1.5, label="Consensus Threshold (>= 3 Votes)")

for bar in bars:
    w = bar.get_width()
    plt.text(w + 0.08, bar.get_y() + bar.get_height()/2.0, f"{int(w)}/4", va="center", ha="left", fontsize=10, fontweight="bold")

plt.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.9)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "ml_4model_consensus_hub_biomarkers.png"), dpi=300)
plt.close()
print(f"\nSaved consensus plot to {os.path.join(plots_dir, 'ml_4model_consensus_hub_biomarkers.png')}")
