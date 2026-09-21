# -*- coding: utf-8 -*-
"""
Generate Comprehensive 15-Panel Publication Figure for Layer 4 Clinical Severity Validation
Includes all 5 priority genes: COL15A1, COL1A1, SERPINE2, SERPINF2, TNXB
Across independent multi-organ cohorts (GSE84044, GSE135251, GSE38958, GSE213001, GSE9285)
Outputs: plots/severity_validation4_layer_boxplots.png (300 DPI, 5x3 grid)
"""

import os
import re
import gzip
import tarfile
import urllib.request
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = r"D:\CSIR"
GEO_CACHE = os.path.join(BASE_DIR, "geo_cache")
PLOTS_DIR = os.path.join(BASE_DIR, "plots")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

GENES = ['COL15A1', 'COL1A1', 'SERPINE2', 'SERPINF2', 'TNXB']
sns.set_theme(style="whitegrid", font_scale=1.05)

plot_data_dict = {}

# ==============================================================================
# 1. LIVER: GSE135251 (RNA-seq, N=216)
# ==============================================================================
print("1. Parsing Liver RNA-seq: GSE135251...")
with gzip.open(os.path.join(GEO_CACHE, 'GSE135251_series_matrix.txt.gz'), 'rt', errors='ignore') as f:
    lines = f.readlines()

gsm_list, title_list, fib_list, nas_list = [], [], [], []
for line in lines:
    if line.startswith('!Sample_geo_accession'):
        gsm_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
    elif line.startswith('!Sample_title'):
        title_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
    elif line.startswith('!Sample_characteristics_ch1'):
        vals = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
        if any('fibrosis stage:' in v.lower() for v in vals):
            fib_list = vals
        elif any('nas score:' in v.lower() for v in vals):
            nas_list = vals

parsed_meta_135 = []
for i, gsm in enumerate(gsm_list):
    f_str = fib_list[i] if i < len(fib_list) else ""
    m_f = re.search(r'fibrosis stage:\s*(\d+)', f_str)
    fib_val = float(m_f.group(1)) if m_f else np.nan
    
    n_str = nas_list[i] if i < len(nas_list) else ""
    m_n = re.search(r'nas score:\s*(\d+)', n_str)
    nas_val = float(m_n.group(1)) if m_n else np.nan
    parsed_meta_135.append({'gsm': gsm, 'title': title_list[i], 'fibrosis_stage': fib_val, 'nas_score': nas_val})

df_meta_135 = pd.DataFrame(parsed_meta_135).set_index('gsm')

ENSEMBL_MAP = {
    'COL15A1': 'ENSG00000204262',
    'COL1A1': 'ENSG00000108821',
    'SERPINE2': 'ENSG00000135914',
    'SERPINF2': 'ENSG00000163661',
    'TNXB': 'ENSG00000168477'
}
target_ens_to_gene = {v: k for k, v in ENSEMBL_MAP.items()}

expr_135 = {g: {} for g in GENES}
with tarfile.open(os.path.join(GEO_CACHE, 'GSE135251_RAW.tar'), 'r') as tar:
    for member in tar.getmembers():
        gsm_m = re.match(r'(GSM\d+)', member.name)
        if not gsm_m:
            continue
        gsm = gsm_m.group(1)
        f = tar.extractfile(member)
        with gzip.open(f, 'rt') as gz:
            sample_counts = {}
            total_counts = 0
            for line in gz:
                parts = line.strip().split()
                if len(parts) >= 2:
                    ens = parts[0]
                    try:
                        cnt = float(parts[1])
                        if ens in target_ens_to_gene:
                            sample_counts[target_ens_to_gene[ens]] = cnt
                        if not ens.startswith('__'):
                            total_counts += cnt
                    except:
                        pass
            for g in GENES:
                raw_c = sample_counts.get(g, 0.0)
                cpm = (raw_c / total_counts) * 1e6 if total_counts > 0 else 0.0
                expr_135[g][gsm] = np.log2(cpm + 1.0)

df_expr_135 = pd.DataFrame(expr_135)
df_135 = df_meta_135.join(df_expr_135)
plot_data_dict['GSE135251'] = df_135
print(f"  -> GSE135251 loaded: {len(df_135)} samples")

# ==============================================================================
# 2. LIVER: GSE84044 (Microarray, N=124)
# ==============================================================================
print("2. Parsing Liver Microarray: GSE84044...")
with gzip.open(os.path.join(GEO_CACHE, 'GSE84044_series_matrix.txt.gz'), 'rt', errors='ignore') as f:
    lines = f.readlines()

gsm_list, title_list, scheuer_s, scheuer_g = [], [], [], []
for line in lines:
    if line.startswith('!Sample_geo_accession'):
        gsm_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
    elif line.startswith('!Sample_title'):
        title_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
    elif line.startswith('!Sample_characteristics_ch1'):
        vals = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
        if any('scheuer score s:' in v.lower() for v in vals):
            scheuer_s = vals
        elif any('scheuer score g:' in v.lower() for v in vals):
            scheuer_g = vals

parsed_meta_84 = []
for i, gsm in enumerate(gsm_list):
    s_str = scheuer_s[i] if i < len(scheuer_s) else ""
    m_s = re.search(r'scheuer score s:\s*(\d+)', s_str, re.IGNORECASE)
    s_val = float(m_s.group(1)) if m_s else np.nan
    
    g_str = scheuer_g[i] if i < len(scheuer_g) else ""
    m_g = re.search(r'scheuer score g:\s*(\d+)', g_str, re.IGNORECASE)
    g_val = float(m_g.group(1)) if m_g else np.nan
    parsed_meta_84.append({'gsm': gsm, 'title': title_list[i], 'scheuer_stage': s_val, 'scheuer_grade': g_val})

df_meta_84 = pd.DataFrame(parsed_meta_84).set_index('gsm')

PROBES_84 = {
    'COL15A1': ['205423_at', '204619_s_at', '204620_s_at'],
    'COL1A1': ['202310_s_at', '202311_s_at'],
    'SERPINE2': ['202619_s_at', '201850_at'],
    'SERPINF2': ['205445_at'],
    'TNXB': ['208354_at', '206016_at']
}
all_target_probes_84 = set([p for sub in PROBES_84.values() for p in sub])
probe_expr_84 = {p: {} for p in all_target_probes_84}

with gzip.open(os.path.join(GEO_CACHE, 'GSE84044_series_matrix.txt.gz'), 'rt', errors='ignore') as f:
    in_table = False
    samples_hdr = []
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            in_table = True
            hdr_line = f.readline().strip().split('\t')
            samples_hdr = [x.replace('"', '') for x in hdr_line[1:]]
            continue
        if line.startswith('!series_matrix_table_end'):
            break
        if in_table:
            parts = line.strip().split('\t')
            probe_id = parts[0].replace('"', '')
            if probe_id in all_target_probes_84:
                for s, val_str in zip(samples_hdr, parts[1:]):
                    try:
                        probe_expr_84[probe_id][s] = float(val_str.replace('"', '').strip())
                    except:
                        probe_expr_84[probe_id][s] = np.nan

df_probe_84 = pd.DataFrame(probe_expr_84)
expr_84 = {}
for g, plist in PROBES_84.items():
    sub_cols = [p for p in plist if p in df_probe_84.columns and not df_probe_84[p].isna().all()]
    if sub_cols:
        expr_84[g] = df_probe_84[sub_cols].mean(axis=1)
    else:
        expr_84[g] = pd.Series(np.nan, index=df_probe_84.index)

df_expr_84 = pd.DataFrame(expr_84)
df_84 = df_meta_84.join(df_expr_84)
plot_data_dict['GSE84044'] = df_84
print(f"  -> GSE84044 loaded: {len(df_84)} samples")

# ==============================================================================
# 3. LUNG: GSE38958 (Microarray, N=60)
# ==============================================================================
print("3. Parsing Lung Microarray: GSE38958...")
with gzip.open(os.path.join(GEO_CACHE, 'GSE38958_series_matrix.txt.gz'), 'rt', errors='ignore') as f:
    lines = f.readlines()

gsm_list, title_list, gse38_fvc, gse38_dlco = [], [], [], []
for line in lines:
    if line.startswith('!Sample_geo_accession'):
        gsm_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
    elif line.startswith('!Sample_title'):
        title_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
    elif line.startswith('!Sample_characteristics_ch1'):
        vals = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
        if any('forced vital capacity' in v for v in vals):
            gse38_fvc = vals
        elif any('diffusing capacity' in v for v in vals):
            gse38_dlco = vals

parsed_meta_38 = []
for i, gsm in enumerate(gsm_list):
    f_str = gse38_fvc[i] if i < len(gse38_fvc) else ""
    m_f = re.search(r'predicted:\s*([0-9.]+)', f_str, re.IGNORECASE)
    fvc_val = float(m_f.group(1)) if m_f else np.nan
    
    d_str = gse38_dlco[i] if i < len(gse38_dlco) else ""
    m_d = re.search(r'predicted:\s*([0-9.]+)', d_str, re.IGNORECASE)
    dlco_val = float(m_d.group(1)) if m_d else np.nan
    parsed_meta_38.append({'gsm': gsm, 'title': title_list[i], 'fvc_pct': fvc_val, 'dlco_pct': dlco_val})

df_meta_38 = pd.DataFrame(parsed_meta_38).set_index('gsm')

PROBES_38 = {
    'COL15A1': '3720790',
    'COL1A1': '3744654',
    'SERPINE2': '2576974',
    'SERPINF2': '3779836',
    'TNXB': '2967160'
}
probe_to_gene_38 = {v: k for k, v in PROBES_38.items()}
expr_38 = {g: {} for g in GENES}

with gzip.open(os.path.join(GEO_CACHE, 'GSE38958_series_matrix.txt.gz'), 'rt', errors='ignore') as f:
    in_table = False
    samples_hdr = []
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            in_table = True
            hdr_line = f.readline().strip().split('\t')
            samples_hdr = [x.replace('"', '') for x in hdr_line[1:]]
            continue
        if line.startswith('!series_matrix_table_end'):
            break
        if in_table:
            parts = line.strip().split('\t')
            probe_id = parts[0].replace('"', '')
            if probe_id in probe_to_gene_38:
                g = probe_to_gene_38[probe_id]
                for s, val_str in zip(samples_hdr, parts[1:]):
                    try:
                        expr_38[g][s] = float(val_str.replace('"', '').strip())
                    except:
                        expr_38[g][s] = np.nan

df_expr_38 = pd.DataFrame(expr_38)
df_38 = df_meta_38.join(df_expr_38)
plot_data_dict['GSE38958'] = df_38
print(f"  -> GSE38958 loaded: {len(df_38)} samples")

# ==============================================================================
# 4. SKIN: GSE9285 (Microarray, N=74)
# ==============================================================================
print("4. Parsing Skin Microarray: GSE9285...")
xml_cache = os.path.join(GEO_CACHE, 'pone_0002696.xml')
if not os.path.exists(xml_cache):
    xml_url = 'https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0002696&type=manuscript'
    req = urllib.request.Request(xml_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as resp:
        xml_data = resp.read()
    with open(xml_cache, 'wb') as fp:
        fp.write(xml_data)
else:
    with open(xml_cache, 'rb') as fp:
        xml_data = fp.read()

root = ET.fromstring(xml_data)
mrss_patient_map = {}
for table_wrap in root.iter('table-wrap'):
    for row in table_wrap.iter('tr'):
        cells = [c.text or ''.join(c.itertext()) for c in row.iter('td')]
        if len(cells) >= 4:
            subj = cells[0].strip()
            score_str = cells[3].strip()
            m_score = re.search(r'(\d+)', score_str)
            if m_score:
                mrss_val = float(m_score.group(1))
                mrss_patient_map[subj.lower().replace(' ', '')] = mrss_val
                subj_clean = re.sub(r'[^a-zA-Z0-9]', '', subj.lower())
                mrss_patient_map[subj_clean] = mrss_val

gsm_list, title_list = [], []
with gzip.open(os.path.join(GEO_CACHE, 'GSE9285_series_matrix.txt.gz'), 'rt', errors='ignore') as f:
    for line in f:
        if line.startswith('!Sample_geo_accession'):
            gsm_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
        elif line.startswith('!Sample_title'):
            title_list = [x.strip().replace('"', '') for x in line.split('\t')[1:]]
        elif line.startswith('!series_matrix_table_begin'):
            break

parsed_meta_92 = []
for i, gsm in enumerate(gsm_list):
    t = title_list[i]
    m_subj = re.match(r'([a-zA-Z]+\d+)', t)
    mrss = np.nan
    if m_subj:
        s_key = m_subj.group(1).lower()
        if s_key in mrss_patient_map:
            mrss = mrss_patient_map[s_key]
    elif 'normal' in t.lower() or 'control' in t.lower():
        mrss = 0.0
    parsed_meta_92.append({'gsm': gsm, 'title': t, 'mrss': mrss})

df_meta_92 = pd.DataFrame(parsed_meta_92).set_index('gsm')

PROBES_92 = {
    'COL15A1': ['A_23_P112554', 'A_24_P315120'],
    'COL1A1': ['A_23_P207520'],
    'SERPINE2': ['A_23_P50919'],
    'SERPINF2': ['A_23_P89270'],
    'TNXB': ['A_23_P156708', 'A_24_P15834', 'A_24_P213884']
}
all_target_probes_92 = set([p for sub in PROBES_92.values() for p in sub])
probe_expr_92 = {p: {} for p in all_target_probes_92}

with gzip.open(os.path.join(GEO_CACHE, 'GSE9285_series_matrix.txt.gz'), 'rt', errors='ignore') as f:
    in_table = False
    samples_hdr = []
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            in_table = True
            hdr_line = f.readline().strip().split('\t')
            samples_hdr = [x.replace('"', '') for x in hdr_line[1:]]
            continue
        if line.startswith('!series_matrix_table_end'):
            break
        if in_table:
            parts = line.strip().split('\t')
            probe_id = parts[0].replace('"', '')
            if probe_id in all_target_probes_92:
                for s, val_str in zip(samples_hdr, parts[1:]):
                    try:
                        probe_expr_92[probe_id][s] = float(val_str.replace('"', '').strip())
                    except:
                        probe_expr_92[probe_id][s] = np.nan

df_probe_92 = pd.DataFrame(probe_expr_92)
expr_92 = {}
for g, plist in PROBES_92.items():
    sub_cols = [p for p in plist if p in df_probe_92.columns and not df_probe_92[p].isna().all()]
    if sub_cols:
        expr_92[g] = df_probe_92[sub_cols].mean(axis=1)
    else:
        expr_92[g] = pd.Series(np.nan, index=df_probe_92.index)

df_expr_92 = pd.DataFrame(expr_92)
df_92 = df_meta_92.join(df_expr_92)
plot_data_dict['GSE9285'] = df_92
print(f"  -> GSE9285 loaded: {len(df_92)} samples")

# Load Master Table for Exact Rho / P-value Annotations
df_master = pd.read_csv(os.path.join(RESULTS_DIR, 'validation_4_severity_master_table.csv'))
def get_stats(ds, g, metric):
    sub = df_master[(df_master['dataset'] == ds) & (df_master['gene'] == g) & (df_master['severity_metric'] == metric)]
    if len(sub) > 0:
        row = sub.iloc[0]
        return row['n_full'], row['rho_full'], row['p_adj_full']
    return np.nan, np.nan, np.nan

# ==============================================================================
# 5. GENERATE COMPREHENSIVE 15-PANEL HIGH-RESOLUTION PUBLICATION FIGURE (5x3 GRID)
# ==============================================================================
print("\nGenerating Comprehensive 15-Panel Figure (5x3 Grid)...")
fig, axes = plt.subplots(5, 3, figsize=(18, 25))

palette_scheuer = 'Purples'
palette_kleiner = 'Blues'
color_dlco = '#e67e22'
color_mrss = '#27ae60'

# ROW 1: COL15A1
# 1.1: Liver GSE84044 Scheuer
ax = axes[0, 0]
df_p = plot_data_dict['GSE84044'].dropna(subset=['scheuer_stage', 'COL15A1']).copy()
df_p['Stage'] = 'S' + df_p['scheuer_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='COL15A1', hue='Stage', ax=ax, palette=palette_scheuer, order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='COL15A1', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
n, rho, padj = get_stats('GSE84044', 'COL15A1', 'Scheuer Fibrosis Stage (S0-S4)')
ax.set_title(f"Liver (Microarray: GSE84044, N={int(n)})\nCOL15A1 vs. Scheuer Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Scheuer Fibrosis Stage (S0 - S4)")
ax.set_ylabel("COL15A1 (log2 intensity)")

# 1.2: Liver GSE135251 Kleiner
ax = axes[0, 1]
df_p = plot_data_dict['GSE135251'].dropna(subset=['fibrosis_stage', 'COL15A1']).copy()
df_p['Stage'] = 'F' + df_p['fibrosis_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='COL15A1', hue='Stage', ax=ax, palette=palette_kleiner, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='COL15A1', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
n, rho, padj = get_stats('GSE135251', 'COL15A1', 'Kleiner Fibrosis Stage (F0-F4)')
ax.set_title(f"Liver (RNA-seq: GSE135251, N={int(n)})\nCOL15A1 vs. Kleiner Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Kleiner Fibrosis Stage (F0 - F4)")
ax.set_ylabel("COL15A1 [log2(CPM+1)]")

# 1.3: Lung GSE38958 DLCO
ax = axes[0, 2]
df_p = plot_data_dict['GSE38958'].dropna(subset=['dlco_pct', 'COL15A1'])
sns.regplot(data=df_p, x='dlco_pct', y='COL15A1', ax=ax, color=color_dlco, scatter_kws={'alpha': 0.6}, line_kws={'linewidth': 2})
n, rho, padj = get_stats('GSE38958', 'COL15A1', '% Predicted DLCO (Continuous)')
ax.set_title(f"Lung (Microarray: GSE38958, N={int(n)})\nCOL15A1 vs. % Predicted DLCO (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("% Predicted DLCO (Gas Diffusion Impairment)")
ax.set_ylabel("COL15A1 (log2 intensity)")

# ROW 2: COL1A1
# 2.1: Liver GSE84044 Scheuer
ax = axes[1, 0]
df_p = plot_data_dict['GSE84044'].dropna(subset=['scheuer_stage', 'COL1A1']).copy()
df_p['Stage'] = 'S' + df_p['scheuer_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='COL1A1', hue='Stage', ax=ax, palette=palette_scheuer, order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='COL1A1', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
n, rho, padj = get_stats('GSE84044', 'COL1A1', 'Scheuer Fibrosis Stage (S0-S4)')
ax.set_title(f"Liver (Microarray: GSE84044, N={int(n)})\nCOL1A1 vs. Scheuer Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Scheuer Fibrosis Stage (S0 - S4)")
ax.set_ylabel("COL1A1 (log2 intensity)")

# 2.2: Liver GSE135251 Kleiner
ax = axes[1, 1]
df_p = plot_data_dict['GSE135251'].dropna(subset=['fibrosis_stage', 'COL1A1']).copy()
df_p['Stage'] = 'F' + df_p['fibrosis_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='COL1A1', hue='Stage', ax=ax, palette=palette_kleiner, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='COL1A1', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
n, rho, padj = get_stats('GSE135251', 'COL1A1', 'Kleiner Fibrosis Stage (F0-F4)')
ax.set_title(f"Liver (RNA-seq: GSE135251, N={int(n)})\nCOL1A1 vs. Kleiner Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Kleiner Fibrosis Stage (F0 - F4)")
ax.set_ylabel("COL1A1 [log2(CPM+1)]")

# 2.3: Lung GSE38958 DLCO
ax = axes[1, 2]
df_p = plot_data_dict['GSE38958'].dropna(subset=['dlco_pct', 'COL1A1'])
sns.regplot(data=df_p, x='dlco_pct', y='COL1A1', ax=ax, color=color_dlco, scatter_kws={'alpha': 0.6}, line_kws={'linewidth': 2})
n, rho, padj = get_stats('GSE38958', 'COL1A1', '% Predicted DLCO (Continuous)')
ax.set_title(f"Lung (Microarray: GSE38958, N={int(n)})\nCOL1A1 vs. % Predicted DLCO (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("% Predicted DLCO (Gas Diffusion Impairment)")
ax.set_ylabel("COL1A1 (log2 intensity)")

# ROW 3: SERPINE2 (Crucial missing Tier 1 hub!)
# 3.1: Liver GSE84044 Scheuer
ax = axes[2, 0]
df_p = plot_data_dict['GSE84044'].dropna(subset=['scheuer_stage', 'SERPINE2']).copy()
df_p['Stage'] = 'S' + df_p['scheuer_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='SERPINE2', hue='Stage', ax=ax, palette=palette_scheuer, order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='SERPINE2', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
n, rho, padj = get_stats('GSE84044', 'SERPINE2', 'Scheuer Fibrosis Stage (S0-S4)')
ax.set_title(f"Liver (Microarray: GSE84044, N={int(n)})\nSERPINE2 vs. Scheuer Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Scheuer Fibrosis Stage (S0 - S4)")
ax.set_ylabel("SERPINE2 (log2 intensity)")

# 3.2: Liver GSE135251 Kleiner
ax = axes[2, 1]
df_p = plot_data_dict['GSE135251'].dropna(subset=['fibrosis_stage', 'SERPINE2']).copy()
df_p['Stage'] = 'F' + df_p['fibrosis_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='SERPINE2', hue='Stage', ax=ax, palette=palette_kleiner, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='SERPINE2', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
n, rho, padj = get_stats('GSE135251', 'SERPINE2', 'Kleiner Fibrosis Stage (F0-F4)')
ax.set_title(f"Liver (RNA-seq: GSE135251, N={int(n)})\nSERPINE2 vs. Kleiner Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Kleiner Fibrosis Stage (F0 - F4)")
ax.set_ylabel("SERPINE2 [log2(CPM+1)]")

# 3.3: Skin GSE9285 mRSS
ax = axes[2, 2]
df_p = plot_data_dict['GSE9285'].dropna(subset=['mrss', 'SERPINE2'])
sns.regplot(data=df_p, x='mrss', y='SERPINE2', ax=ax, color=color_mrss, scatter_kws={'alpha': 0.6}, line_kws={'linewidth': 2})
n, rho, padj = get_stats('GSE9285', 'SERPINE2', 'Modified Rodnan Skin Score (mRSS 0-51)')
ax.set_title(f"Skin (Microarray: GSE9285, N={int(n)})\nSERPINE2 vs. mRSS Skin Score (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Modified Rodnan Skin Score (mRSS 0-51)")
ax.set_ylabel("SERPINE2 (log2 ratio)")

# ROW 4: SERPINF2 (Antiprotease Downregulation Axis)
# 4.1: Liver GSE84044 Scheuer
ax = axes[3, 0]
df_p = plot_data_dict['GSE84044'].dropna(subset=['scheuer_stage', 'SERPINF2']).copy()
df_p['Stage'] = 'S' + df_p['scheuer_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='SERPINF2', hue='Stage', ax=ax, palette='YlOrRd_r', order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='SERPINF2', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['S0', 'S1', 'S2', 'S3', 'S4'], legend=False)
n, rho, padj = get_stats('GSE84044', 'SERPINF2', 'Scheuer Fibrosis Stage (S0-S4)')
ax.set_title(f"Liver (Microarray: GSE84044, N={int(n)})\nSERPINF2 vs. Scheuer Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Scheuer Fibrosis Stage (S0 - S4)")
ax.set_ylabel("SERPINF2 (log2 intensity)")

# 4.2: Liver GSE135251 Kleiner
ax = axes[3, 1]
df_p = plot_data_dict['GSE135251'].dropna(subset=['fibrosis_stage', 'SERPINF2']).copy()
df_p['Stage'] = 'F' + df_p['fibrosis_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='SERPINF2', hue='Stage', ax=ax, palette='YlOrRd_r', order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='SERPINF2', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
n, rho, padj = get_stats('GSE135251', 'SERPINF2', 'Kleiner Fibrosis Stage (F0-F4)')
ax.set_title(f"Liver (RNA-seq: GSE135251, N={int(n)})\nSERPINF2 vs. Kleiner Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Kleiner Fibrosis Stage (F0 - F4)")
ax.set_ylabel("SERPINF2 [log2(CPM+1)]")

# 4.3: Lung GSE38958 DLCO
ax = axes[3, 2]
df_p = plot_data_dict['GSE38958'].dropna(subset=['dlco_pct', 'SERPINF2'])
sns.regplot(data=df_p, x='dlco_pct', y='SERPINF2', ax=ax, color=color_dlco, scatter_kws={'alpha': 0.6}, line_kws={'linewidth': 2})
n, rho, padj = get_stats('GSE38958', 'SERPINF2', '% Predicted DLCO (Continuous)')
ax.set_title(f"Lung (Microarray: GSE38958, N={int(n)})\nSERPINF2 vs. % Predicted DLCO (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("% Predicted DLCO (Gas Diffusion Impairment)")
ax.set_ylabel("SERPINF2 (log2 intensity)")

# ROW 5: TNXB (Exploratory Causal Target)
# 5.1: Skin GSE9285 mRSS
ax = axes[4, 0]
df_p = plot_data_dict['GSE9285'].dropna(subset=['mrss', 'TNXB'])
sns.regplot(data=df_p, x='mrss', y='TNXB', ax=ax, color=color_mrss, scatter_kws={'alpha': 0.6}, line_kws={'linewidth': 2})
n, rho, padj = get_stats('GSE9285', 'TNXB', 'Modified Rodnan Skin Score (mRSS 0-51)')
ax.set_title(f"Skin (Microarray: GSE9285, N={int(n)})\nTNXB vs. mRSS Skin Score (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Modified Rodnan Skin Score (mRSS 0-51)")
ax.set_ylabel("TNXB (log2 ratio)")

# 5.2: Liver GSE135251 Kleiner
ax = axes[4, 1]
df_p = plot_data_dict['GSE135251'].dropna(subset=['fibrosis_stage', 'TNXB']).copy()
df_p['Stage'] = 'F' + df_p['fibrosis_stage'].astype(int).astype(str)
sns.boxplot(data=df_p, x='Stage', y='TNXB', hue='Stage', ax=ax, palette=palette_kleiner, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
sns.stripplot(data=df_p, x='Stage', y='TNXB', hue='Stage', ax=ax, color='black', alpha=0.5, jitter=0.2, order=['F0', 'F1', 'F2', 'F3', 'F4'], legend=False)
n, rho, padj = get_stats('GSE135251', 'TNXB', 'Kleiner Fibrosis Stage (F0-F4)')
ax.set_title(f"Liver (RNA-seq: GSE135251, N={int(n)})\nTNXB vs. Kleiner Stage (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("Kleiner Fibrosis Stage (F0 - F4)")
ax.set_ylabel("TNXB [log2(CPM+1)]")

# 5.3: Lung GSE38958 DLCO
ax = axes[4, 2]
df_p = plot_data_dict['GSE38958'].dropna(subset=['dlco_pct', 'TNXB'])
sns.regplot(data=df_p, x='dlco_pct', y='TNXB', ax=ax, color=color_dlco, scatter_kws={'alpha': 0.6}, line_kws={'linewidth': 2})
n, rho, padj = get_stats('GSE38958', 'TNXB', '% Predicted DLCO (Continuous)')
ax.set_title(f"Lung (Microarray: GSE38958, N={int(n)})\nTNXB vs. % Predicted DLCO (rho={rho:+.3f}, p_adj={padj:.2e})", fontweight='bold', fontsize=10.5)
ax.set_xlabel("% Predicted DLCO (Gas Diffusion Impairment)")
ax.set_ylabel("TNXB (log2 intensity)")

plt.tight_layout()
out_png = os.path.join(PLOTS_DIR, 'severity_validation4_layer_boxplots.png')
plt.savefig(out_png, dpi=300)
plt.close()
print(f"\nALL 15 PANELS GENERATED SUCCESSFULLY: {out_png}")
