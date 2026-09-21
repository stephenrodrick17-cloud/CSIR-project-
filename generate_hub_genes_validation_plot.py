# -*- coding: utf-8 -*-
"""
Generate Canonical Validation Figure for the 4 Tier 1 Universal Hub Genes:
COL15A1, COL1A1, SERPINE2, SERPINF2

Layout: 4 Rows x 3 Columns (or 4 Rows x 2 Columns with Mann-Whitney U and Spearman Severity)
- Column 1: Held-Out Validation 2 Two-Group Differential Expression (Mann-Whitney U Test, Normal vs Cirrhosis, GSE14323)
- Column 2: Clinical Severity Dose-Response Correlation (Spearman rho, Histological Fibrosis Stage)
- Column 3: Diagnostic ROC Analysis (True Positive Rate vs False Positive Rate, AUC with 95% CI)

Output: plots/hub_genes_validation_mwu_spearman.png
"""
import gzip
import re
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_curve, roc_auc_score

os.makedirs("plots", exist_ok=True)
os.makedirs("results", exist_ok=True)

# Set high-grade styling
plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 1.0

# 4 Tier 1 Hub Genes
hub_genes = ["COL15A1", "COL1A1", "SERPINE2", "SERPINF2"]
gene_subtitles = {
    "COL15A1": "#1 Universal Matrix Hub (Basement Membrane Organizer)",
    "COL1A1": "Canonical Fibrillar Collagen (Mechanical Stiffness Driver)",
    "SERPINE2": "Antiprotease Hub (Upregulated Fibrinolysis Inhibitor)",
    "SERPINF2": "Antiprotease Hub (Conserved Downregulated Serpin)"
}

# Microarray probe mapping (GSE14323, GSE84044)
probe_map = {
    "COL15A1": "203477_at",
    "COL1A1": "202310_s_at",
    "SERPINE2": "212190_at",
    "SERPINF2": "205075_at"
}
rev_probe_map = {v: k for k, v in probe_map.items()}

# Ensembl mapping (GSE162694 RNA-seq)
ens_map = {
    "COL15A1": "ENSG00000204291",
    "COL1A1": "ENSG00000108821",
    "SERPINE2": "ENSG00000135914",
    "SERPINF2": "ENSG00000167711"
}
rev_ens_map = {v: k for k, v in ens_map.items()}

# ==============================================================================
# 1. PARSE GSE14323 (Held-Out Liver Val 2, n=60: 19 Normal, 41 Cirrhosis)
# ==============================================================================
print("1. Parsing GSE14323 for Mann-Whitney U and ROC analysis...")
sample_expr_val2 = {g: {"normal": [], "cirrhosis": []} for g in hub_genes}

with gzip.open("geo_cache/GSE14323_family.soft.gz", "rt", encoding="utf-8", errors="ignore") as f:
    curr_type = None
    in_table = False
    for line in f:
        if line.startswith("^SAMPLE = "):
            in_table = False
            curr_type = None
        elif line.startswith("!Sample_source_name_ch1 = "):
            src = line.strip().split(" = ")[1].lower()
            if "cirrhosis" in src and "hcc" not in src:
                curr_type = "cirrhosis"
            elif "normal" in src:
                curr_type = "normal"
        elif line.startswith("!sample_table_begin"):
            if curr_type in ("cirrhosis", "normal"):
                in_table = True
        elif line.startswith("!sample_table_end"):
            in_table = False
        elif in_table:
            parts = line.strip().split("\t")
            if len(parts) >= 2 and parts[0] in rev_probe_map:
                g = rev_probe_map[parts[0]]
                sample_expr_val2[g][curr_type].append(float(parts[1]))

mwu_results = {}
roc_results = {}
for g in hub_genes:
    ctrl = sample_expr_val2[g]["normal"]
    dis = sample_expr_val2[g]["cirrhosis"]
    u_stat, p_raw = stats.mannwhitneyu(dis, ctrl, alternative="two-sided")
    
    # ROC
    y_true = [0] * len(ctrl) + [1] * len(dis)
    y_score = ctrl + dis
    # if downregulated in fibrosis (like SERPINF2), invert for ROC
    if np.median(dis) < np.median(ctrl):
        y_score_roc = [-x for x in y_score]
    else:
        y_score_roc = y_score
    fpr, tpr, _ = roc_curve(y_true, y_score_roc)
    auc = roc_auc_score(y_true, y_score_roc)
    
    # 1000-bootstrap 95% CI for AUC
    np.random.seed(42)
    boot_aucs = []
    y_true_arr = np.array(y_true)
    y_score_arr = np.array(y_score_roc)
    for _ in range(1000):
        idx = np.random.choice(len(y_true), len(y_true), replace=True)
        if len(np.unique(y_true_arr[idx])) < 2:
            continue
        boot_aucs.append(roc_auc_score(y_true_arr[idx], y_score_arr[idx]))
    ci_lower = np.percentile(boot_aucs, 2.5)
    ci_upper = np.percentile(boot_aucs, 97.5)
    
    mwu_results[g] = {
        "u": u_stat, "p_raw": p_raw,
        "ctrl_median": np.median(ctrl), "dis_median": np.median(dis),
        "ctrl": ctrl, "dis": dis
    }
    roc_results[g] = {
        "fpr": fpr, "tpr": tpr, "auc": auc,
        "ci_lower": ci_lower, "ci_upper": ci_upper
    }

# Benjamini-Hochberg FDR for MWU
p_raws_mwu = [mwu_results[g]["p_raw"] for g in hub_genes]
p_adjs_mwu = stats.false_discovery_control(p_raws_mwu)
for g, p_adj in zip(hub_genes, p_adjs_mwu):
    mwu_results[g]["p_adj"] = p_adj

# ==============================================================================
# 2. PARSE GSE84044 (Liver Scheuer Staging S0-S4, n=124 total, n=87 disease S1-S4)
# ==============================================================================
print("2. Parsing GSE84044 for Spearman clinical severity tracking...")
with gzip.open("geo_cache/GSE84044_series_matrix.txt.gz", "rt", encoding="utf-8", errors="ignore") as f:
    sample_ids_84044 = []
    chars_84044 = []
    expr_dict_84044 = {}
    in_table = False
    for line in f:
        line = line.rstrip("\r\n")
        if line.startswith("!Sample_geo_accession"):
            sample_ids_84044 = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            chars_84044.append([x.strip('"') for x in line.split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            in_table = True
        elif line.startswith("!series_matrix_table_end"):
            in_table = False
        elif in_table:
            parts = line.split("\t")
            probe_id = parts[0].strip('"')
            if probe_id in rev_probe_map:
                gene = rev_probe_map[probe_id]
                expr_dict_84044[gene] = [float(x) if x not in ("null", "NA", "") else np.nan for x in parts[1:]]

stages_84044 = []
for idx in range(len(sample_ids_84044)):
    sample_chars = [chars_84044[c_idx][idx] for c_idx in range(len(chars_84044))]
    stage_val = np.nan
    for c in sample_chars:
        m = re.search(r"scheuer.*?([0-4])", c, re.IGNORECASE)
        if m:
            stage_val = float(m.group(1))
            break
        elif "stage:" in c.lower() or "fibrosis:" in c.lower():
            m2 = re.search(r"([0-4])", c)
            if m2:
                stage_val = float(m2.group(1))
                break
    stages_84044.append(stage_val)

df_84044 = pd.DataFrame(expr_dict_84044, index=sample_ids_84044)
df_84044["stage"] = stages_84044
df_dis_84044 = df_84044[df_84044["stage"] > 0].dropna(subset=["stage"])

spearman_results = {}
p_raws_spearman = []
for g in hub_genes:
    sub = df_dis_84044.dropna(subset=[g])
    rho, p_val = stats.spearmanr(sub[g], sub["stage"])
    spearman_results[g] = {
        "rho": rho, "p_raw": p_val, "n": len(sub),
        "expr": sub[g].values, "stage": sub["stage"].values
    }
    p_raws_spearman.append(p_val)

p_adjs_spearman = stats.false_discovery_control(p_raws_spearman)
for g, p_adj in zip(hub_genes, p_adjs_spearman):
    spearman_results[g]["p_adj"] = p_adj

# ==============================================================================
# 3. GENERATE COMPREHENSIVE 4-ROW X 3-COLUMN VALIDATION FIGURE
# ==============================================================================
print("3. Generating high-resolution 4x3 validation figure...")
fig, axes = plt.subplots(4, 3, figsize=(18, 20))
fig.patch.set_facecolor("white")

palette_groups = {"Control (n=19)": "#2A9D8F", "Cirrhosis (n=41)": "#E76F51"}

for row_idx, g in enumerate(hub_genes):
    # -------------------------------------------------------------------------
    # COLUMN 1: MANN-WHITNEY U TEST BOXPLOT (GSE14323)
    # -------------------------------------------------------------------------
    ax_box = axes[row_idx, 0]
    ctrl_vals = mwu_results[g]["ctrl"]
    dis_vals = mwu_results[g]["dis"]
    
    df_plot_box = pd.DataFrame({
        "Expression": ctrl_vals + dis_vals,
        "Group": ["Control (n=19)"] * len(ctrl_vals) + ["Cirrhosis (n=41)"] * len(dis_vals)
    })
    
    sns.boxplot(
        data=df_plot_box, x="Group", y="Expression", hue="Group", ax=ax_box,
        palette=palette_groups, width=0.42, boxprops=dict(alpha=0.8, edgecolor="#264653", linewidth=1.5),
        whiskerprops=dict(linewidth=1.2, color="#264653"), capprops=dict(linewidth=1.2, color="#264653"),
        medianprops=dict(color="#D62828", linewidth=2.5), legend=False,
        showmeans=True, meanprops={"marker":"D", "markerfacecolor":"#FFD166", "markeredgecolor":"#264653", "markersize":7}
    )
    sns.stripplot(
        data=df_plot_box, x="Group", y="Expression", ax=ax_box,
        color="#1D3557", alpha=0.6, jitter=0.2, size=6, edgecolor="white", linewidth=0.5
    )
    
    u_val = mwu_results[g]["u"]
    p_adj_mwu = mwu_results[g]["p_adj"]
    log2fc = np.median(dis_vals) - np.median(ctrl_vals)
    dir_str = "Upregulated" if log2fc > 0 else "Downregulated"
    
    ax_box.set_title(
        f"{g}: Held-Out Val 2 ({dir_str})\nMann-Whitney U = {u_val:.1f} | p_adj = {p_adj_mwu:.2e}",
        fontsize=13, fontweight="bold", pad=10, color="#1D3557"
    )
    ax_box.set_xlabel("")
    ax_box.set_ylabel(f"{g} Log2 Intensity", fontsize=11, fontweight="bold")
    ax_box.tick_params(labelsize=10)
    
    # -------------------------------------------------------------------------
    # COLUMN 2: SPEARMAN CLINICAL SEVERITY REGRESSION (GSE84044 Scheuer Stages)
    # -------------------------------------------------------------------------
    ax_scat = axes[row_idx, 1]
    sub_expr = spearman_results[g]["expr"]
    sub_stage = spearman_results[g]["stage"]
    
    df_scat = pd.DataFrame({"Stage": sub_stage, "Expression": sub_expr})
    
    # Add subtle jitter to stages for visibility
    np.random.seed(42 + row_idx)
    jitter = np.random.normal(0, 0.08, size=len(sub_stage))
    df_scat["Jittered_Stage"] = df_scat["Stage"] + jitter
    
    sns.regplot(
        data=df_scat, x="Jittered_Stage", y="Expression", ax=ax_scat,
        color="#457B9D",
        scatter_kws={"alpha": 0.65, "s": 45, "color": "#1D3557", "edgecolor": "white", "linewidths": 0.5},
        line_kws={"color": "#E63946", "linewidth": 2.5}
    )
    
    rho_val = spearman_results[g]["rho"]
    p_adj_spear = spearman_results[g]["p_adj"]
    n_pts = spearman_results[g]["n"]
    
    ax_scat.set_title(
        f"{g}: Fibrosis Severity Dose-Response\nSpearman rho = {rho_val:+.3f} | p_adj = {p_adj_spear:.2e} (n={n_pts})",
        fontsize=13, fontweight="bold", pad=10, color="#1D3557"
    )
    ax_scat.set_xlabel("Histological Fibrosis Stage (Scheuer S1 - S4)", fontsize=11, fontweight="bold")
    ax_scat.set_ylabel(f"{g} Log2 Intensity", fontsize=11, fontweight="bold")
    ax_scat.set_xticks([1, 2, 3, 4])
    ax_scat.set_xticklabels(["S1 (Portal)", "S2 (Periportal)", "S3 (Bridging)", "S4 (Cirrhosis)"], fontsize=10)
    ax_scat.tick_params(labelsize=10)
    
    # -------------------------------------------------------------------------
    # COLUMN 3: DIAGNOSTIC ROC CURVE (GSE14323)
    # -------------------------------------------------------------------------
    ax_roc = axes[row_idx, 2]
    fpr = roc_results[g]["fpr"]
    tpr = roc_results[g]["tpr"]
    auc_val = roc_results[g]["auc"]
    ci_l = roc_results[g]["ci_lower"]
    ci_u = roc_results[g]["ci_upper"]
    
    ax_roc.plot(fpr, tpr, color="#E76F51", linewidth=2.5, label=f"ROC Curve (AUC = {auc_val:.3f})")
    ax_roc.plot([0, 1], [0, 1], color="#999999", linestyle="--", linewidth=1.5, label="Random Chance (AUC = 0.50)")
    
    ax_roc.set_title(
        f"{g}: Held-Out Diagnostic ROC\nAUC = {auc_val:.3f} (95% CI: {ci_l:.3f} - {ci_u:.3f})",
        fontsize=13, fontweight="bold", pad=10, color="#1D3557"
    )
    ax_roc.set_xlabel("1 - Specificity (False Positive Rate)", fontsize=11, fontweight="bold")
    ax_roc.set_ylabel("Sensitivity (True Positive Rate)", fontsize=11, fontweight="bold")
    ax_roc.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#CCCCCC", fontsize=10)
    ax_roc.set_xlim(-0.02, 1.02)
    ax_roc.set_ylim(-0.02, 1.02)
    ax_roc.tick_params(labelsize=10)

plt.suptitle(
    "Universal Pan-Fibrotic Hub Biomarkers Validation (N=4 Tier 1 Hubs)\nHeld-Out Differential Expression (Mann-Whitney U), Severity Dose-Response (Spearman rho), and Diagnostic ROC Analysis",
    fontsize=16, fontweight="bold", y=0.995, color="#1D3557"
)

plt.tight_layout(rect=[0, 0.02, 1, 0.98])
out_plot_path = os.path.join("plots", "hub_genes_validation_mwu_spearman.png")
plt.savefig(out_plot_path, dpi=300)
plt.close()
print(f"Saved canonical validation figure: {out_plot_path}")

# ==============================================================================
# 4. SAVE CANONICAL NUMERICAL RESULTS CSV
# ==============================================================================
out_csv_path = os.path.join("results", "hub_genes_validation_mwu_spearman_metrics.csv")
records = []
for g in hub_genes:
    records.append({
        "gene": g,
        "biological_role": gene_subtitles[g],
        "val2_cohort": "GSE14323 (Liver Microarray)",
        "val2_n_ctrl": len(mwu_results[g]["ctrl"]),
        "val2_n_dis": len(mwu_results[g]["dis"]),
        "val2_ctrl_median": round(mwu_results[g]["ctrl_median"], 3),
        "val2_dis_median": round(mwu_results[g]["dis_median"], 3),
        "val2_mwu_u": mwu_results[g]["u"],
        "val2_mwu_p_raw": mwu_results[g]["p_raw"],
        "val2_mwu_p_adj": mwu_results[g]["p_adj"],
        "val2_diagnostic_roc_auc": round(roc_results[g]["auc"], 4),
        "val2_roc_ci_lower": round(roc_results[g]["ci_lower"], 4),
        "val2_roc_ci_upper": round(roc_results[g]["ci_upper"], 4),
        "severity_cohort": "GSE84044 (Scheuer Staging S1-S4)",
        "severity_n_dis": spearman_results[g]["n"],
        "severity_spearman_rho": round(spearman_results[g]["rho"], 4),
        "severity_spearman_p_raw": spearman_results[g]["p_raw"],
        "severity_spearman_p_adj": spearman_results[g]["p_adj"]
    })

pd.DataFrame(records).to_csv(out_csv_path, index=False)
print(f"Saved numerical validation table: {out_csv_path}")
