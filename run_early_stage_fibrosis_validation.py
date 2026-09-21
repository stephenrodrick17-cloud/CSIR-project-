# -*- coding: utf-8 -*-
"""
Dedicated Early-Detection Validation for the 4 Pan-Fibrotic Tier 1 Hub Genes
(COL15A1, COL1A1, SERPINE2, SERPINF2)

Clinical Severity Cohorts:
1. GSE84044 (Chronic Hepatitis B Microarray, N=124 total -> N=96 early cohort: 43 S0 vs 53 S1/S2)
2. GSE135251 (NASH/NAFLD RNA-seq, N=216 total -> N=148 early cohort: 46 F0 vs 102 F1/F2)

Steps:
1. Stage Subset Filtering: Strict S0/F0 (Control) vs S1/S2 or F1/F2 (Early-Stage Fibrosis). Exclude late stage (S3/S4 or F3/F4).
2. Individual Early ROC AUC: Mann-Whitney U test, log2FC, and 1,000-bootstrap 95% CI AUCs.
3. Multi-Gene Early Classifiers: 5-fold Stratified CV (Logistic Regression, Random Forest, Ensemble) -> Out-of-fold AUC, Sensitivity, Specificity.
4. Decision Curve Analysis (DCA): Clinical Net Benefit across threshold probabilities pt in [0.05, 0.50].
5. Outputs:
   - plots/hub_genes_early_stage_validation.png
   - results/hub_genes_early_stage_validation_metrics.csv
   - results/hub_genes_early_stage_dca_curve.csv
"""

import os
import re
import gzip
import tarfile
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs("results", exist_ok=True)
os.makedirs("plots", exist_ok=True)

plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.9

HUB_GENES = ["COL15A1", "COL1A1", "SERPINE2", "SERPINF2"]
GENE_COLORS = {
    "COL15A1": "#1e3a8a",  # Deep Blue
    "COL1A1": "#b91c1c",   # Crimson
    "SERPINE2": "#d97706", # Amber/Orange
    "SERPINF2": "#047857"  # Emerald Green
}

# ==============================================================================
# 1. PARSE GSE84044 (Liver Microarray, Scheuer Staging S0 vs S1/S2)
# ==============================================================================
print("1. Parsing GSE84044 (Scheuer Staging)...")
gsm_84, scheuer_s = [], []
with gzip.open("geo_cache/GSE84044_series_matrix.txt.gz", "rt", errors="ignore") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            gsm_84 = [x.strip().replace('"', '') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1") and "scheuer score s:" in line.lower():
            scheuer_s = [x.strip().replace('"', '') for x in line.split("\t")[1:]]
        elif line.startswith("!series_matrix_table_begin"):
            break

stages_84 = {}
for gsm, s_str in zip(gsm_84, scheuer_s):
    m = re.search(r"scheuer score s:\s*(\d+)", s_str, re.IGNORECASE)
    stages_84[gsm] = int(m.group(1)) if m else np.nan

probe_map_84 = {
    "COL15A1": "203477_at",
    "COL1A1": "202310_s_at",
    "SERPINE2": "212190_at",
    "SERPINF2": "205075_at"
}
rev_map_84 = {v: k for k, v in probe_map_84.items()}
expr_84 = {g: {} for g in probe_map_84}

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
            if pid in rev_map_84:
                g = rev_map_84[pid]
                for s, val_str in zip(samples_hdr, parts[1:]):
                    try:
                        expr_84[g][s] = float(val_str.replace('"', '').strip())
                    except:
                        pass

df_84 = pd.DataFrame(expr_84)
df_84["stage"] = [stages_84.get(idx) for idx in df_84.index]

# Filter strictly S0 vs (S1 or S2), exclude S3 and S4
df_early_84 = df_84[df_84["stage"].isin([0, 1, 2])].copy()
df_early_84["label"] = (df_early_84["stage"] > 0).astype(int)

n_ctrl_84 = (df_early_84["label"] == 0).sum()
n_early_84 = (df_early_84["label"] == 1).sum()
print(f"  GSE84044 Early Cohort: N = {len(df_early_84)} (Controls S0: {n_ctrl_84}, Early S1/S2: {n_early_84})")

# ==============================================================================
# 2. PARSE GSE135251 (Liver RNA-seq, Kleiner Staging F0 vs F1/F2)
# ==============================================================================
print("2. Parsing GSE135251 (Kleiner Staging)...")
cached_135_path = "results/gse135251_early_stage_cache.csv"

if os.path.exists(cached_135_path):
    print("  Loading GSE135251 from local cache...")
    df_early_135 = pd.read_csv(cached_135_path, index_col=0)
else:
    gsm_135, fib_list_135 = [], []
    with gzip.open("geo_cache/GSE135251_series_matrix.txt.gz", "rt", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                gsm_135 = [x.strip().replace('"', '') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and any("fibrosis stage:" in v.lower() for v in line.split("\t")[1:]):
                fib_list_135 = [x.strip().replace('"', '') for x in line.split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                break

    stages_135 = {}
    for gsm, s_str in zip(gsm_135, fib_list_135):
        m = re.search(r"fibrosis stage:\s*(\d+)", s_str, re.IGNORECASE)
        stages_135[gsm] = int(m.group(1)) if m else np.nan

    ens_map_135 = {
        "COL15A1": "ENSG00000204262",
        "COL1A1": "ENSG00000108821",
        "SERPINE2": "ENSG00000135914",
        "SERPINF2": "ENSG00000167711"
    }
    rev_map_135 = {v: k for k, v in ens_map_135.items()}
    sample_expr_135 = {g: {} for g in ens_map_135}

    with tarfile.open("geo_cache/GSE135251_RAW.tar", "r") as tar:
        for member in tar.getmembers():
            m_gsm = re.match(r"(GSM\d+)", member.name)
            if not m_gsm:
                continue
            gsm = m_gsm.group(1)
            f = tar.extractfile(member)
            with gzip.open(f, "rt") as gz:
                sample_c = {}
                tot = 0.0
                for line in gz:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        ens = parts[0]
                        try:
                            c = float(parts[1])
                            if ens in rev_map_135:
                                sample_c[rev_map_135[ens]] = c
                            if not ens.startswith("__"):
                                tot += c
                        except:
                            pass
                for g in ens_map_135:
                    cpm = (sample_c.get(g, 0.0) / tot) * 1e6 if tot > 0 else 0.0
                    sample_expr_135[g][gsm] = np.log2(cpm + 1.0)

    df_135 = pd.DataFrame(sample_expr_135)
    df_135["stage"] = [stages_135.get(idx) for idx in df_135.index]
    df_early_135 = df_135[df_135["stage"].isin([0, 1, 2])].copy()
    df_early_135["label"] = (df_early_135["stage"] > 0).astype(int)
    df_early_135.to_csv(cached_135_path)

n_ctrl_135 = (df_early_135["label"] == 0).sum()
n_early_135 = (df_early_135["label"] == 1).sum()
print(f"  GSE135251 Early Cohort: N = {len(df_early_135)} (Controls F0: {n_ctrl_135}, Early F1/F2: {n_early_135})")

# ==============================================================================
# 3. INDIVIDUAL EARLY ROC AUC & STATISTICAL METRICS (1,000 BOOTSTRAP CIs)
# ==============================================================================
print("3. Computing individual early-detection statistics & bootstrap 95% CIs...")

def compute_individual_metrics(df_sub, cohort_name):
    metrics_list = []
    p_raw_list = []
    roc_curves_dict = {}
    
    y_true = df_sub["label"].values
    
    for g in HUB_GENES:
        ctrl_vals = df_sub[df_sub["label"] == 0][g].values
        early_vals = df_sub[df_sub["label"] == 1][g].values
        
        # Log2 Fold Change (mean difference on log2 scale)
        log2fc = np.mean(early_vals) - np.mean(ctrl_vals)
        
        # Two-sided Mann-Whitney U test
        u_stat, p_raw = stats.mannwhitneyu(early_vals, ctrl_vals, alternative="two-sided")
        p_raw_list.append(p_raw)
        
        # ROC Score (inverted for downregulated markers like SERPINF2)
        y_score = df_sub[g].values
        if log2fc < 0:
            y_score_eval = -y_score
        else:
            y_score_eval = y_score
            
        fpr, tpr, _ = roc_curve(y_true, y_score_eval)
        auc_val = roc_auc_score(y_true, y_score_eval)
        
        # 1,000 Bootstrap 95% CI
        np.random.seed(42)
        boot_aucs = []
        n_samples = len(y_true)
        for _ in range(1000):
            b_idx = np.random.choice(n_samples, n_samples, replace=True)
            if len(np.unique(y_true[b_idx])) < 2:
                continue
            boot_aucs.append(roc_auc_score(y_true[b_idx], y_score_eval[b_idx]))
            
        ci_lower = np.percentile(boot_aucs, 2.5)
        ci_upper = np.percentile(boot_aucs, 97.5)
        
        roc_curves_dict[g] = {
            "fpr": fpr, "tpr": tpr, "auc": auc_val,
            "ci_lower": ci_lower, "ci_upper": ci_upper,
            "log2fc": log2fc, "u_stat": u_stat, "p_raw": p_raw
        }
        
        metrics_list.append({
            "Cohort": cohort_name,
            "Gene": g,
            "Log2FC": log2fc,
            "MWU_U": u_stat,
            "MWU_p_raw": p_raw,
            "ROC_AUC": auc_val,
            "AUC_95CI_Lower": ci_lower,
            "AUC_95CI_Upper": ci_upper
        })
        
    p_adj_list = stats.false_discovery_control(p_raw_list)
    for i, p_adj in enumerate(p_adj_list):
        metrics_list[i]["MWU_p_adj"] = p_adj
        roc_curves_dict[HUB_GENES[i]]["p_adj"] = p_adj
        
    return pd.DataFrame(metrics_list), roc_curves_dict

df_metrics_84, rocs_84 = compute_individual_metrics(df_early_84, "GSE84044 (Scheuer S0 vs S1/S2)")
df_metrics_135, rocs_135 = compute_individual_metrics(df_early_135, "GSE135251 (Kleiner F0 vs F1/F2)")

# ==============================================================================
# 4. MULTI-GENE 5-FOLD STRATIFIED CV CLASSIFIERS (LR, RF, ENSEMBLE)
# ==============================================================================
print("4. Training 5-fold cross-validated early classifiers...")

def train_cv_classifiers(df_sub):
    X = df_sub[HUB_GENES].values
    y = df_sub["label"].values
    
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof_lr = np.zeros(len(y))
    oof_rf = np.zeros(len(y))
    
    for train_idx, val_idx in skf.split(X, y):
        X_tr, y_tr = X[train_idx], y[train_idx]
        X_va, y_va = X[val_idx], y[val_idx]
        
        lr = LogisticRegression(penalty="l2", C=1.0, random_state=42)
        lr.fit(X_tr, y_tr)
        oof_lr[val_idx] = lr.predict_proba(X_va)[:, 1]
        
        rf = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
        rf.fit(X_tr, y_tr)
        oof_rf[val_idx] = rf.predict_proba(X_va)[:, 1]
        
    oof_ens = (oof_lr + oof_rf) / 2.0
    
    # Evaluate Out-of-Fold ROC
    fpr_lr, tpr_lr, _ = roc_curve(y, oof_lr)
    auc_lr = roc_auc_score(y, oof_lr)
    
    fpr_rf, tpr_rf, _ = roc_curve(y, oof_rf)
    auc_rf = roc_auc_score(y, oof_rf)
    
    fpr_ens, tpr_ens, thresh_ens = roc_curve(y, oof_ens)
    auc_ens = roc_auc_score(y, oof_ens)
    
    # 1,000 Bootstrap CI for Ensemble
    np.random.seed(42)
    boot_ens = []
    for _ in range(1000):
        b_idx = np.random.choice(len(y), len(y), replace=True)
        if len(np.unique(y[b_idx])) < 2:
            continue
        boot_ens.append(roc_auc_score(y[b_idx], oof_ens[b_idx]))
    ci_lower_ens = np.percentile(boot_ens, 2.5)
    ci_upper_ens = np.percentile(boot_ens, 97.5)
    
    # Optimal threshold via Youden's Index
    youden = tpr_ens - fpr_ens
    best_idx = np.argmax(youden)
    opt_th = thresh_ens[best_idx]
    pred_ens = (oof_ens >= opt_th).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred_ens).ravel()
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    
    results = {
        "oof_lr": oof_lr, "auc_lr": auc_lr, "fpr_lr": fpr_lr, "tpr_lr": tpr_lr,
        "oof_rf": oof_rf, "auc_rf": auc_rf, "fpr_rf": fpr_rf, "tpr_rf": tpr_rf,
        "oof_ens": oof_ens, "auc_ens": auc_ens, "ci_lower_ens": ci_lower_ens, "ci_upper_ens": ci_upper_ens,
        "fpr_ens": fpr_ens, "tpr_ens": tpr_ens,
        "opt_th": opt_th, "sens": sens, "spec": spec,
        "y": y
    }
    return results

cv_res_84 = train_cv_classifiers(df_early_84)
cv_res_135 = train_cv_classifiers(df_early_135)

print(f"  GSE84044 5-Fold CV Ensemble AUC = {cv_res_84['auc_ens']:.4f} (95% CI: {cv_res_84['ci_lower_ens']:.3f}-{cv_res_84['ci_upper_ens']:.3f})")
print(f"    Sensitivity = {cv_res_84['sens']:.3f}, Specificity = {cv_res_84['spec']:.3f} (at optimal threshold {cv_res_84['opt_th']:.3f})")
print(f"  GSE135251 5-Fold CV Ensemble AUC = {cv_res_135['auc_ens']:.4f} (95% CI: {cv_res_135['ci_lower_ens']:.3f}-{cv_res_135['ci_upper_ens']:.3f})")
print(f"    Sensitivity = {cv_res_135['sens']:.3f}, Specificity = {cv_res_135['spec']:.3f} (at optimal threshold {cv_res_135['opt_th']:.3f})")

# ==============================================================================
# 5. DECISION CURVE ANALYSIS (DCA) FOR SUBCLINICAL SCREENING (pt = 0.05 to 0.50)
# ==============================================================================
print("5. Computing Decision Curve Analysis (DCA)...")

thresholds = np.linspace(0.05, 0.50, 46)

def compute_dca(y_true, y_prob, pt_range):
    n = len(y_true)
    prev = np.mean(y_true)
    net_benefit_model = []
    net_benefit_all = []
    net_benefit_none = np.zeros(len(pt_range))
    
    for pt in pt_range:
        # Model predictions at threshold pt
        y_pred = (y_prob >= pt).astype(int)
        tp = np.sum((y_pred == 1) & (y_true == 1))
        fp = np.sum((y_pred == 1) & (y_true == 0))
        
        # Net Benefit formula: TP/N - FP/N * (pt / (1 - pt))
        weight = pt / (1.0 - pt)
        nb_mod = (tp / n) - (fp / n) * weight
        net_benefit_model.append(nb_mod)
        
        # Treat All formula: Prev - (1 - Prev) * (pt / (1 - pt))
        nb_all = prev - (1.0 - prev) * weight
        net_benefit_all.append(nb_all)
        
    return pd.DataFrame({
        "threshold": pt_range,
        "net_benefit_model": net_benefit_model,
        "net_benefit_all": net_benefit_all,
        "net_benefit_none": net_benefit_none
    })

df_dca_84 = compute_dca(cv_res_84["y"], cv_res_84["oof_ens"], thresholds)
df_dca_84["cohort"] = "GSE84044 (Scheuer S0 vs S1/S2)"
df_dca_84.to_csv("results/hub_genes_early_stage_dca_curve.csv", index=False)

# ==============================================================================
# 6. OUTPUT PUBLICATION-READY SUMMARY TABLE
# ==============================================================================
print("6. Formatting publication summary tables...")

summary_rows = []
for idx, row in df_metrics_84.iterrows():
    g = row["Gene"]
    summary_rows.append({
        "Gene Symbol": g,
        "Cohort": "GSE84044 (Scheuer S0 vs S1/S2)",
        "Sample Size": f"N={len(df_early_84)} (Ctrl={n_ctrl_84}, Early={n_early_84})",
        "Early vs Control Log2FC": f"{row['Log2FC']:+.3f}",
        "Mann-Whitney U": f"{row['MWU_U']:.1f}",
        "Mann-Whitney p-raw": f"{row['MWU_p_raw']:.2e}",
        "Mann-Whitney p-adj": f"{row['MWU_p_adj']:.2e}",
        "Early ROC AUC (95% CI)": f"{row['ROC_AUC']:.3f} ({row['AUC_95CI_Lower']:.3f}–{row['AUC_95CI_Upper']:.3f})",
        "Ensemble Early AUC (5-Fold CV)": f"{cv_res_84['auc_ens']:.3f} ({cv_res_84['ci_lower_ens']:.3f}–{cv_res_84['ci_upper_ens']:.3f})",
        "Ensemble Sensitivity": f"{cv_res_84['sens']:.1%}",
        "Ensemble Specificity": f"{cv_res_84['spec']:.1%}"
    })

# Add comparative cohort GSE135251
for idx, row in df_metrics_135.iterrows():
    g = row["Gene"]
    summary_rows.append({
        "Gene Symbol": g,
        "Cohort": "GSE135251 (Kleiner F0 vs F1/F2)",
        "Sample Size": f"N={len(df_early_135)} (Ctrl={n_ctrl_135}, Early={n_early_135})",
        "Early vs Control Log2FC": f"{row['Log2FC']:+.3f}",
        "Mann-Whitney U": f"{row['MWU_U']:.1f}",
        "Mann-Whitney p-raw": f"{row['MWU_p_raw']:.2e}",
        "Mann-Whitney p-adj": f"{row['MWU_p_adj']:.2e}",
        "Early ROC AUC (95% CI)": f"{row['ROC_AUC']:.3f} ({row['AUC_95CI_Lower']:.3f}–{row['AUC_95CI_Upper']:.3f})",
        "Ensemble Early AUC (5-Fold CV)": f"{cv_res_135['auc_ens']:.3f} ({cv_res_135['ci_lower_ens']:.3f}–{cv_res_135['ci_upper_ens']:.3f})",
        "Ensemble Sensitivity": f"{cv_res_135['sens']:.1%}",
        "Ensemble Specificity": f"{cv_res_135['spec']:.1%}"
    })

df_summary_pub = pd.DataFrame(summary_rows)
out_summary_csv = "results/hub_genes_early_stage_validation_metrics.csv"
df_summary_pub.to_csv(out_summary_csv, index=False)
print(f"Saved: {out_summary_csv}")

# ==============================================================================
# 7. GENERATE MULTI-PANEL PUBLICATION FIGURE (300 DPI)
# ==============================================================================
print("7. Generating publication figure...")

fig = plt.figure(figsize=(19, 13), dpi=300)
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1.2, 1, 1], hspace=0.32, wspace=0.28)

# ------------------------------------------------------------------------------
# PANEL A: Expression Boxplots (S0 vs S1/S2 in GSE84044)
# ------------------------------------------------------------------------------
ax_box = fig.add_subplot(gs[0, 0])
melt_data = []
for g in HUB_GENES:
    for _, r in df_early_84.iterrows():
        melt_data.append({
            "Gene": g,
            "Expression": r[g],
            "Group": "Control (S0, n=43)" if r["label"] == 0 else "Early Fibrosis (S1/S2, n=53)"
        })
df_melt = pd.DataFrame(melt_data)

palette_box = {"Control (S0, n=43)": "#3b82f6", "Early Fibrosis (S1/S2, n=53)": "#ef4444"}
sns.boxplot(data=df_melt, x="Gene", y="Expression", hue="Group", palette=palette_box,
            ax=ax_box, width=0.55, fliersize=0, boxprops=dict(alpha=0.85, linewidth=1.2))
sns.stripplot(data=df_melt, x="Gene", y="Expression", hue="Group", palette=palette_box,
              ax=ax_box, dodge=True, alpha=0.55, size=4.5, jitter=0.2, edgecolor="black", linewidth=0.4)

handles, labels = ax_box.get_legend_handles_labels()
ax_box.legend(handles[:2], labels[:2], loc="upper right", frameon=True, fontsize=9.5, facecolor="white")

ax_box.set_title("A. Early-Stage Expression Profile (GSE84044)\nStrictly S0 Controls vs S1/S2 Early Disease",
                 fontsize=12, fontweight="bold", pad=8)
ax_box.set_ylabel("Microarray log2 Signal Intensity", fontsize=10.5, fontweight="bold")
ax_box.set_xlabel("4 Pan-Fibrotic Tier-1 Hub Genes", fontsize=10.5, fontweight="bold")

y_tops = {"COL15A1": 9.2, "COL1A1": 13.6, "SERPINE2": 11.2, "SERPINF2": 12.0}
for i, g in enumerate(HUB_GENES):
    p_adj_val = rocs_84[g]["p_adj"]
    fc_val = rocs_84[g]["log2fc"]
    sig_label = f"p={p_adj_val:.2e}\nLog2FC={fc_val:+.2f}"
    ax_box.text(i, y_tops[g], sig_label, ha="center", va="bottom", fontsize=8.2,
                fontweight="bold", color="#1e293b", bbox=dict(boxstyle="round,pad=0.2", facecolor="#f8fafc", edgecolor="#cbd5e1", alpha=0.85))

ax_box.set_ylim(4.0, 15.0)

# ------------------------------------------------------------------------------
# PANEL B: Individual Early ROC Curves (GSE84044)
# ------------------------------------------------------------------------------
ax_roc_ind = fig.add_subplot(gs[0, 1])
ax_roc_ind.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", linewidth=1.2, label="Chance (AUC = 0.500)")

for g in HUB_GENES:
    info = rocs_84[g]
    lbl = f"{g} (AUC = {info['auc']:.3f}, 95% CI {info['ci_lower']:.3f}–{info['ci_upper']:.3f})"
    ax_roc_ind.plot(info["fpr"], info["tpr"], color=GENE_COLORS[g], linewidth=2.0, label=lbl)

ax_roc_ind.set_title("B. Individual Early ROC Discriminatory Power\nGSE84044 (S0 Controls vs S1/S2 Early Fibrosis)",
                     fontsize=12, fontweight="bold", pad=8)
ax_roc_ind.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=10.5, fontweight="bold")
ax_roc_ind.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10.5, fontweight="bold")
ax_roc_ind.legend(loc="lower right", fontsize=8.2, frameon=True, facecolor="white")
ax_roc_ind.set_xlim(-0.02, 1.02)
ax_roc_ind.set_ylim(-0.02, 1.02)

# ------------------------------------------------------------------------------
# PANEL C: 5-Fold Cross-Validated Multi-Gene Classifier ROC
# ------------------------------------------------------------------------------
ax_roc_cv = fig.add_subplot(gs[0, 2])
ax_roc_cv.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", linewidth=1.2, label="Chance (AUC = 0.500)")

ax_roc_cv.plot(cv_res_84["fpr_lr"], cv_res_84["tpr_lr"], color="#3b82f6", linewidth=1.8, linestyle="-.",
               label=f"Logistic Regression (AUC = {cv_res_84['auc_lr']:.3f})")
ax_roc_cv.plot(cv_res_84["fpr_rf"], cv_res_84["tpr_rf"], color="#10b981", linewidth=1.8, linestyle=":",
               label=f"Random Forest (AUC = {cv_res_84['auc_rf']:.3f})")
ax_roc_cv.plot(cv_res_84["fpr_ens"], cv_res_84["tpr_ens"], color="#7c3aed", linewidth=2.4,
               label=f"Ensemble 4-Hub (AUC = {cv_res_84['auc_ens']:.3f})\n95% CI: [{cv_res_84['ci_lower_ens']:.3f}–{cv_res_84['ci_upper_ens']:.3f}]")

opt_idx = np.argmax(cv_res_84["tpr_ens"] - cv_res_84["fpr_ens"])
ax_roc_cv.scatter(cv_res_84["fpr_ens"][opt_idx], cv_res_84["tpr_ens"][opt_idx],
                  color="#dc2626", s=75, zorder=5, edgecolors="black", linewidth=1.2,
                  label=f"Optimal Threshold (Sens={cv_res_84['sens']:.1%}, Spec={cv_res_84['spec']:.1%})")

ax_roc_cv.set_title("C. Multi-Gene 5-Fold Cross-Validation\nEarly Subclinical Detection (S0 vs S1/S2)",
                    fontsize=12, fontweight="bold", pad=8)
ax_roc_cv.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=10.5, fontweight="bold")
ax_roc_cv.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10.5, fontweight="bold")
ax_roc_cv.legend(loc="lower right", fontsize=8.0, frameon=True, facecolor="white")
ax_roc_cv.set_xlim(-0.02, 1.02)
ax_roc_cv.set_ylim(-0.02, 1.02)

# ------------------------------------------------------------------------------
# PANEL D: Decision Curve Analysis (DCA) Net Benefit
# ------------------------------------------------------------------------------
ax_dca = fig.add_subplot(gs[1, 0])
ax_dca.plot(df_dca_84["threshold"], df_dca_84["net_benefit_model"], color="#7c3aed", linewidth=2.4,
            label=f"4-Hub Ensemble Early Model (GSE84044)")
ax_dca.plot(df_dca_84["threshold"], df_dca_84["net_benefit_all"], color="#94a3b8", linewidth=1.5, linestyle="--",
            label="Treat All (Screen All High-Risk Patients)")
ax_dca.plot(df_dca_84["threshold"], df_dca_84["net_benefit_none"], color="#334155", linewidth=1.5, linestyle=":",
            label="Treat None (Net Benefit = 0)")

ax_dca.set_title("D. Decision Curve Analysis (DCA)\nSubclinical Screening Thresholds ($p_t$ = 0.05 to 0.50)",
                 fontsize=12, fontweight="bold", pad=8)
ax_dca.set_xlabel("Threshold Risk Probability ($p_t$)", fontsize=10.5, fontweight="bold")
ax_dca.set_ylabel("Clinical Net Benefit", fontsize=10.5, fontweight="bold")
ax_dca.legend(loc="upper right", fontsize=8.5, frameon=True, facecolor="white")
ax_dca.set_xlim(0.05, 0.50)
ax_dca.set_ylim(-0.05, 0.60)
ax_dca.axvspan(0.10, 0.45, color="#ede9fe", alpha=0.35, label="High-Utility Window")

# ------------------------------------------------------------------------------
# PANEL E: Independent External Replication (GSE135251 RNA-seq, F0 vs F1/F2)
# ------------------------------------------------------------------------------
ax_rep = fig.add_subplot(gs[1, 1])
ax_rep.plot([0, 1], [0, 1], linestyle="--", color="#94a3b8", linewidth=1.2, label="Chance (AUC = 0.500)")

for g in HUB_GENES:
    info = rocs_135[g]
    lbl = f"{g} (AUC = {info['auc']:.3f}, 95% CI {info['ci_lower']:.3f}–{info['ci_upper']:.3f})"
    ax_rep.plot(info["fpr"], info["tpr"], color=GENE_COLORS[g], linewidth=1.8, label=lbl)

ax_rep.plot(cv_res_135["fpr_ens"], cv_res_135["tpr_ens"], color="#7c3aed", linewidth=2.4, linestyle="-",
            label=f"4-Hub Ensemble (AUC = {cv_res_135['auc_ens']:.3f})")

ax_rep.set_title("E. External Early Cohort Replication\nGSE135251 RNA-seq (F0 Controls vs F1/F2 Early Fibrosis, N=148)",
                 fontsize=12, fontweight="bold", pad=8)
ax_rep.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=10.5, fontweight="bold")
ax_rep.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10.5, fontweight="bold")
ax_rep.legend(loc="lower right", fontsize=8.0, frameon=True, facecolor="white")
ax_rep.set_xlim(-0.02, 1.02)
ax_rep.set_ylim(-0.02, 1.02)

# ------------------------------------------------------------------------------
# PANEL F: Summary Metrics Heatmap / Comparison Table
# ------------------------------------------------------------------------------
ax_tbl = fig.add_subplot(gs[1, 2])
ax_tbl.axis("off")

table_data = [
    ["Gene", "GSE84044\nAUC (95% CI)", "GSE84044\np-adj", "GSE135251\nAUC (95% CI)", "GSE135251\np-adj"],
    ["COL15A1", f"{rocs_84['COL15A1']['auc']:.3f} ({rocs_84['COL15A1']['ci_lower']:.2f}–{rocs_84['COL15A1']['ci_upper']:.2f})", f"{rocs_84['COL15A1']['p_adj']:.2e}", f"{rocs_135['COL15A1']['auc']:.3f} ({rocs_135['COL15A1']['ci_lower']:.2f}–{rocs_135['COL15A1']['ci_upper']:.2f})", f"{rocs_135['COL15A1']['p_adj']:.2e}"],
    ["COL1A1", f"{rocs_84['COL1A1']['auc']:.3f} ({rocs_84['COL1A1']['ci_lower']:.2f}–{rocs_84['COL1A1']['ci_upper']:.2f})", f"{rocs_84['COL1A1']['p_adj']:.2e}", f"{rocs_135['COL1A1']['auc']:.3f} ({rocs_135['COL1A1']['ci_lower']:.2f}–{rocs_135['COL1A1']['ci_upper']:.2f})", f"{rocs_135['COL1A1']['p_adj']:.2e}"],
    ["SERPINE2", f"{rocs_84['SERPINE2']['auc']:.3f} ({rocs_84['SERPINE2']['ci_lower']:.2f}–{rocs_84['SERPINE2']['ci_upper']:.2f})", f"{rocs_84['SERPINE2']['p_adj']:.2e}", f"{rocs_135['SERPINE2']['auc']:.3f} ({rocs_135['SERPINE2']['ci_lower']:.2f}–{rocs_135['SERPINE2']['ci_upper']:.2f})", f"{rocs_135['SERPINE2']['p_adj']:.2e}"],
    ["SERPINF2", f"{rocs_84['SERPINF2']['auc']:.3f} ({rocs_84['SERPINF2']['ci_lower']:.2f}–{rocs_84['SERPINF2']['ci_upper']:.2f})", f"{rocs_84['SERPINF2']['p_adj']:.2e}", f"{rocs_135['SERPINF2']['auc']:.3f} ({rocs_135['SERPINF2']['ci_lower']:.2f}–{rocs_135['SERPINF2']['ci_upper']:.2f})", f"{rocs_135['SERPINF2']['p_adj']:.2e}"],
    ["4-Hub Ens", f"{cv_res_84['auc_ens']:.3f} ({cv_res_84['ci_lower_ens']:.2f}–{cv_res_84['ci_upper_ens']:.2f})", "5-Fold CV", f"{cv_res_135['auc_ens']:.3f} ({cv_res_135['ci_lower_ens']:.2f}–{cv_res_135['ci_upper_ens']:.2f})", "5-Fold CV"]
]

tab = ax_tbl.table(cellText=table_data, loc="center", cellLoc="center")
tab.auto_set_font_size(False)
tab.set_fontsize(8.8)
tab.scale(1.05, 1.85)

for (r_idx, c_idx), cell in tab.get_celld().items():
    cell.set_edgecolor("#94a3b8")
    cell.set_linewidth(0.8)
    if r_idx == 0:
        cell.set_facecolor("#1e3a8a")
        cell.set_text_props(weight="bold", color="white")
    elif r_idx == 5:
        cell.set_facecolor("#ede9fe")
        cell.set_text_props(weight="bold", color="#5b21b6")
    elif r_idx % 2 == 1:
        cell.set_facecolor("#f8fafc")

ax_tbl.set_title("F. Cross-Cohort Early Detection Summary\nConsensus Performance Across Microarray & RNA-seq",
                 fontsize=12, fontweight="bold", pad=12)

plt.suptitle("Dedicated Subclinical & Early-Stage Fibrosis Validation of 4 Tier-1 Universal Hub Biomarkers\nStrict Elimination of Late-Stage Fibrosis / Cirrhosis (Healthy Controls S0/F0 vs Early-Stage S1/S2 or F1/F2 Only)",
             fontsize=15, fontweight="bold", y=0.98)

out_fig_path = "plots/hub_genes_early_stage_validation.png"
plt.savefig(out_fig_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Generated Figure: {out_fig_path}")
print("EARLY-STAGE VALIDATION PIPELINE COMPLETED SUCCESSFULLY.")
