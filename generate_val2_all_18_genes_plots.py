# -*- coding: utf-8 -*-
"""
Generate Complete 18-Gene Validation 2 Figures:
1. results/val2_all_18_genes_severity_spearman.csv:
   Full Spearman correlation statistics across all 18 surviving genes for GSE162694 (METAVIR) and GSE84044 (Scheuer).
2. plots/val2_all_18_genes_mannwhitney_boxplots.png:
   6x3 grid of Mann-Whitney U Disease vs Control boxplots for all 18 genes.
3. plots/val2_severity_spearman_correlations.png:
   6x3 grid of Spearman clinical severity dose-response scatterplots for all 18 genes.
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

os.makedirs("results", exist_ok=True)
os.makedirs("plots", exist_ok=True)

# 18 Surviving Genes in Validation 2
genes_18 = [
    'COL15A1', 'AEBP1', 'COL1A2', 'COL3A1', 'SERPINE2', 'VWF',
    'CCL19', 'CCL2', 'SERPINH1', 'COL1A1', 'SERPINF2', 'CCL5',
    'LTBP2', 'SVEP1', 'CLEC2D', 'LAMC3', 'MFAP4', 'PDGFD'
]

# Probe mappings for Affymetrix microarray cohorts
probe_map = {
    'COL15A1': '203477_at',
    'AEBP1': '201792_at',
    'COL1A2': '202403_s_at',
    'COL3A1': '215076_s_at',
    'SERPINE2': '212190_at',
    'VWF': '202112_at',
    'CCL19': '210072_at',
    'CCL2': '216598_s_at',
    'SERPINH1': '207714_s_at',
    'COL1A1': '202310_s_at',
    'SERPINF2': '205075_at',
    'CCL5': '1405_i_at',
    'LTBP2': '204682_at',
    'SVEP1': '213247_at',
    'CLEC2D': '220132_s_at',
    'LAMC3': '219407_s_at',
    'MFAP4': '212713_at',
    'PDGFD': '219304_s_at'
}
rev_probe_map = {v: k for k, v in probe_map.items()}

# Ensembl mappings for RNA-seq GSE162694
ensembl_map = {
    'COL15A1': 'ENSG00000204291',
    'AEBP1': 'ENSG00000106624',
    'COL1A2': 'ENSG00000164692',
    'COL3A1': 'ENSG00000168542',
    'SERPINE2': 'ENSG00000135914',
    'VWF': 'ENSG00000110799',
    'CCL19': 'ENSG00000172216',
    'CCL2': 'ENSG00000108691',
    'SERPINH1': 'ENSG00000149257',
    'COL1A1': 'ENSG00000108821',
    'SERPINF2': 'ENSG00000167711',
    'CCL5': 'ENSG00000161570',
    'LTBP2': 'ENSG00000119314',
    'SVEP1': 'ENSG00000165124',
    'CLEC2D': 'ENSG00000139187',
    'LAMC3': 'ENSG00000136286',
    'MFAP4': 'ENSG00000166482',
    'PDGFD': 'ENSG00000156299'
}
rev_ens_map = {v: k for k, v in ensembl_map.items()}

# ==============================================================================
# 1. PARSE GSE14323 (Liver Val 2) FOR MANN-WHITNEY U
# ==============================================================================
sample_data_val2 = {g: {'cirrhosis': [], 'normal': []} for g in probe_map}
with gzip.open('geo_cache/GSE14323_family.soft.gz', 'rt', encoding='utf-8', errors='ignore') as f:
    curr_type = None
    in_table = False
    for line in f:
        if line.startswith('^SAMPLE = '):
            in_table = False
            curr_type = None
        elif line.startswith('!Sample_source_name_ch1 = '):
            src = line.strip().split(' = ')[1].lower()
            if 'cirrhosis' in src and 'hcc' not in src:
                curr_type = 'cirrhosis'
            elif 'normal' in src:
                curr_type = 'normal'
        elif line.startswith('!sample_table_begin'):
            if curr_type in ('cirrhosis', 'normal'):
                in_table = True
        elif line.startswith('!sample_table_end'):
            in_table = False
        elif in_table:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and parts[0] in rev_probe_map:
                g = rev_probe_map[parts[0]]
                sample_data_val2[g][curr_type].append(float(parts[1]))

mwu_stats = {}
p_raws = []
for g in genes_18:
    c = sample_data_val2[g]['normal']
    d = sample_data_val2[g]['cirrhosis']
    u, p = stats.mannwhitneyu(d, c, alternative='two-sided')
    mwu_stats[g] = {
        'u': u, 'p_raw': p,
        'med_c': np.median(c), 'med_d': np.median(d)
    }
    p_raws.append(p)

p_adjs = stats.false_discovery_control(p_raws)
for g, p_adj in zip(genes_18, p_adjs):
    mwu_stats[g]['p_adj'] = p_adj

# ==============================================================================
# 2. PARSE GSE162694 (Liver METAVIR Staging F1-F4, n=77)
# ==============================================================================
titles, accs, fib_chars = [], [], []
with gzip.open('geo_cache/GSE162694_series_matrix.txt.gz', 'rt', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            titles = [x.strip('"') for x in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_geo_accession'):
            accs = [x.strip('"') for x in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            if 'fibrosis' in line.lower():
                fib_chars = [x.strip('"') for x in line.strip().split('\t')[1:]]
        elif line.startswith('!series_matrix_table_begin'):
            break

col_to_stage_162694 = {}
for t in titles:
    col = t.split()[-1]
    m = re.search(r'_(F[0-4]|N)', t)
    if m:
        st_str = m.group(1)
        col_to_stage_162694[col] = 0.0 if st_str in ('N', 'F0') else float(st_str[1])

if fib_chars:
    for t, fc in zip(titles, fib_chars):
        col = t.split()[-1]
        m = re.search(r'([0-4])', fc)
        if m:
            col_to_stage_162694[col] = float(m.group(1))

with gzip.open('geo_cache/GSE162694_raw_counts.csv.gz', 'rt') as f:
    cols = [c.strip('"') for c in f.readline().strip().split(',')][1:]
    gene_counts_162694 = {}
    total_counts_162694 = np.zeros(len(cols))
    for line in f:
        parts = line.strip().split(',')
        ens = parts[0].strip('"')
        vals = np.array([float(x) for x in parts[1:]])
        total_counts_162694 += vals
        if ens in rev_ens_map:
            gene = rev_ens_map[ens]
            gene_counts_162694[gene] = vals

cpm_162694 = {}
for g, counts in gene_counts_162694.items():
    cpm_162694[g] = np.log2((counts / total_counts_162694) * 1e6 + 1)

dis_idx_162694 = [i for i, c in enumerate(cols) if col_to_stage_162694.get(c, 0) > 0]
dis_stages_162694 = np.array([col_to_stage_162694[cols[i]] for i in dis_idx_162694])

spearman_162694 = {}
p_raws_162694 = []
for g in genes_18:
    expr = [cpm_162694[g][i] for i in dis_idx_162694]
    rho, p = stats.spearmanr(expr, dis_stages_162694)
    spearman_162694[g] = {'rho': rho, 'p_raw': p, 'n': len(expr), 'expr': expr}
    p_raws_162694.append(p)

p_adjs_162694 = stats.false_discovery_control(p_raws_162694)
for g, p_adj in zip(genes_18, p_adjs_162694):
    spearman_162694[g]['p_adj'] = p_adj

# ==============================================================================
# 3. PARSE GSE84044 (Liver Scheuer Staging S1-S4, n=87)
# ==============================================================================
with gzip.open('geo_cache/GSE84044_series_matrix.txt.gz', 'rt', encoding='utf-8', errors='ignore') as f:
    sample_ids_84044 = []
    chars_84044 = []
    expr_dict_84044 = {}
    in_table = False
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith('!Sample_geo_accession'):
            sample_ids_84044 = [x.strip('"') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            chars_84044.append([x.strip('"') for x in line.split('\t')[1:]])
        elif line.startswith('!series_matrix_table_begin'):
            in_table = True
        elif line.startswith('!series_matrix_table_end'):
            in_table = False
        elif in_table:
            parts = line.split('\t')
            probe_id = parts[0].strip('"')
            if probe_id in rev_probe_map:
                gene = rev_probe_map[probe_id]
                expr_dict_84044[gene] = [float(x) if x not in ('null', 'NA', '') else np.nan for x in parts[1:]]

stages_84044 = []
for idx in range(len(sample_ids_84044)):
    sample_chars = [chars_84044[c_idx][idx] for c_idx in range(len(chars_84044))]
    stage_val = np.nan
    for c in sample_chars:
        m = re.search(r'scheuer.*?([0-4])', c, re.IGNORECASE)
        if m:
            stage_val = float(m.group(1))
            break
        elif 'stage:' in c.lower() or 'fibrosis:' in c.lower():
            m2 = re.search(r'([0-4])', c)
            if m2:
                stage_val = float(m2.group(1))
                break
    stages_84044.append(stage_val)

df_84044 = pd.DataFrame(expr_dict_84044, index=sample_ids_84044)
df_84044['stage'] = stages_84044
df_dis_84044 = df_84044[df_84044['stage'] > 0].dropna(subset=['stage'])

spearman_84044 = {}
p_raws_84044 = []
for g in genes_18:
    sub = df_dis_84044.dropna(subset=[g])
    rho, p = stats.spearmanr(sub[g], sub['stage'])
    spearman_84044[g] = {'rho': rho, 'p_raw': p, 'n': len(sub), 'expr': sub[g].values, 'stages': sub['stage'].values}
    p_raws_84044.append(p)

p_adjs_84044 = stats.false_discovery_control(p_raws_84044)
for g, p_adj in zip(genes_18, p_adjs_84044):
    spearman_84044[g]['p_adj'] = p_adj

# ==============================================================================
# 4. SAVE SUMMARY CSV FOR ALL 18 GENES
# ==============================================================================
summary_rows = []
for g in genes_18:
    summary_rows.append({
        'gene': g,
        'val2_mwu_u': mwu_stats[g]['u'],
        'val2_mwu_p_raw': mwu_stats[g]['p_raw'],
        'val2_mwu_p_adj': mwu_stats[g]['p_adj'],
        'gse162694_metavir_n_dis': spearman_162694[g]['n'],
        'gse162694_metavir_rho_dis': spearman_162694[g]['rho'],
        'gse162694_metavir_p_raw': spearman_162694[g]['p_raw'],
        'gse162694_metavir_p_adj': spearman_162694[g]['p_adj'],
        'gse84044_scheuer_n_dis': spearman_84044[g]['n'],
        'gse84044_scheuer_rho_dis': spearman_84044[g]['rho'],
        'gse84044_scheuer_p_raw': spearman_84044[g]['p_raw'],
        'gse84044_scheuer_p_adj': spearman_84044[g]['p_adj']
    })

df_summary = pd.DataFrame(summary_rows)
out_csv = "results/val2_all_18_genes_severity_spearman.csv"
df_summary.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")

# ==============================================================================
# 5. GENERATE PLOT 1: 6x3 MANN-WHITNEY U BOXPLOT GRID (ALL 18 GENES)
# ==============================================================================
plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

fig, axes = plt.subplots(6, 3, figsize=(18, 26), dpi=300)
fig.patch.set_facecolor("white")
palette_box = {"Control (n=19)": "#3b82f6", "Cirrhosis (n=41)": "#ef4444"}

for idx, g in enumerate(genes_18):
    row, col = idx // 3, idx % 3
    ax = axes[row, col]
    
    vals_ctrl = sample_data_val2[g]["normal"]
    vals_dis = sample_data_val2[g]["cirrhosis"]
    
    plot_df = pd.DataFrame({
        "Expression (log2)": vals_ctrl + vals_dis,
        "Group": ["Control (n=19)"] * len(vals_ctrl) + ["Cirrhosis (n=41)"] * len(vals_dis)
    })
    
    sns.boxplot(
        data=plot_df, x="Group", y="Expression (log2)", hue="Group", ax=ax,
        palette=palette_box, width=0.45, boxprops=dict(alpha=0.75), legend=False,
        showmeans=True, meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":5}
    )
    sns.stripplot(
        data=plot_df, x="Group", y="Expression (log2)", ax=ax,
        color="black", alpha=0.55, jitter=0.2, size=4.5
    )
    
    st = mwu_stats[g]
    ax.set_title(f"{g} (GSE14323 Held-Out Liver Val2)", fontsize=11.5, fontweight="bold", pad=8)
    ax.set_ylabel("log2 Expression", fontsize=10, fontweight="bold")
    ax.set_xlabel("")
    
    annotation = (f"Mann-Whitney U = {st['u']:.1f}\n"
                  f"Raw p = {st['p_raw']:.2e}\n"
                  f"FDR p = {st['p_adj']:.2e}\n"
                  f"Ctrl Med: {st['med_c']:.2f} | Dis Med: {st['med_d']:.2f}")
    
    ax.text(0.04, 0.96, annotation, transform=ax.transAxes, fontsize=8.5,
            verticalalignment="top", bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#cbd5e1", alpha=0.9))

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.suptitle("Validation 2: Non-Parametric Mann-Whitney U Testing Across All 18 Surviving Pan-Fibrotic Core ECM Genes\nHeld-Out Human Liver Fibrosis Cohort (GSE14323: 19 Normal Controls vs. 41 HCV Cirrhosis Biopsies)",
             fontsize=14.5, fontweight="bold", y=0.99)

out_mwu_grid = "plots/val2_all_18_genes_mannwhitney_boxplots.png"
plt.savefig(out_mwu_grid, dpi=300, bbox_inches="tight")
plt.close()
print(f"Generated: {out_mwu_grid}")

# ==============================================================================
# 6. GENERATE PLOT 2: 6x3 SPEARMAN CLINICAL SEVERITY REGRESSION GRID (ALL 18 GENES)
# ==============================================================================
fig, axes = plt.subplots(6, 3, figsize=(18, 26), dpi=300)
fig.patch.set_facecolor("white")
np.random.seed(42)

for idx, g in enumerate(genes_18):
    row, col = idx // 3, idx % 3
    ax = axes[row, col]
    
    # Use GSE162694 METAVIR F1-F4 real data
    y_expr = np.array(spearman_162694[g]['expr'])
    x_stages = dis_stages_162694
    r_val = spearman_162694[g]['rho']
    p_val = spearman_162694[g]['p_raw']
    p_adj = spearman_162694[g]['p_adj']
    
    # Add gentle jitter for discrete stage visibility
    jitter_x = x_stages + np.random.normal(0, 0.08, size=len(x_stages))
    
    # Color points based on significance
    pt_color = "#dc2626" if p_adj < 0.05 else "#64748b"
    ax.scatter(jitter_x, y_expr, color=pt_color, edgecolors="black", linewidths=0.5, alpha=0.75, s=45)
    
    # Trendline
    m, b = np.polyfit(x_stages, y_expr, 1)
    x_line = np.linspace(1, 4, 100)
    ax.plot(x_line, m*x_line + b, color="#1e293b", linewidth=1.8, linestyle="--")
    
    ax.set_title(f"{g} — METAVIR Severity Tracking\nGSE162694 Disease Biopsies (n=77)", fontsize=11.5, fontweight="bold", pad=8)
    ax.set_xlabel("Clinical Fibrosis Stage (METAVIR F1–F4)", fontsize=10, fontweight="bold")
    ax.set_ylabel("log2(CPM + 1) Expression", fontsize=10, fontweight="bold")
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(["F1\n(Portal)", "F2\n(Periportal)", "F3\n(Bridging)", "F4\n(Cirrhosis)"], fontsize=9)
    
    sig_str = "FDR p < 0.05 (Sig)" if p_adj < 0.05 else "Non-Sig"
    box_fc = "#fef2f2" if p_adj < 0.05 else "#f8fafc"
    box_ec = "#fca5a5" if p_adj < 0.05 else "#cbd5e1"
    
    annot_text = (f"Spearman rho = {r_val:+.3f}\n"
                  f"Raw p = {p_val:.2e}\n"
                  f"FDR p = {p_adj:.2e} ({sig_str})\n"
                  f"Scheuer rho (GSE84044): {spearman_84044[g]['rho']:+.3f}")
    
    ax.text(0.04, 0.96, annot_text, transform=ax.transAxes, fontsize=8.5,
            verticalalignment="top", bbox=dict(boxstyle="round,pad=0.35", facecolor=box_fc, edgecolor=box_ec, alpha=0.9))

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.suptitle("Clinical Disease Severity Tracking Across All 18 Surviving Pan-Fibrotic Core ECM Genes\nSpearman Rank Correlation Regressions Against Histological Fibrosis Stages (GSE162694, n=77 Biopsies)",
             fontsize=14.5, fontweight="bold", y=0.99)

out_sev_grid = "plots/val2_severity_spearman_correlations.png"
plt.savefig(out_sev_grid, dpi=300, bbox_inches="tight")
plt.close()
print(f"Generated: {out_sev_grid}")
print("ALL 18 GENES PROCESSED AND VISUALIZED SUCCESSFULLY.")
