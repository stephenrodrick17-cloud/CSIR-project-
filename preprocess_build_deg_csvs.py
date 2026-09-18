# Preprocessing Script: Build Per-Organ DEG CSVs from Locked Discovery Cohorts
import os
import glob
import json
import pandas as pd
import numpy as np
import discovery_config


def load_biodbnet_mappings(organ_dir):
    pattern = os.path.join(organ_dir, '**', 'bioDBnet*.txt')
    files = glob.glob(pattern, recursive=True)
    mapping = {}
    for fpath in files:
        try:
            with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
            for line in lines:
                parts = line.strip().split('	')
                if len(parts) >= 2:
                    gene_id, symbol = parts[0].strip(), parts[1].strip().upper()
                    if gene_id and symbol and symbol not in ('-', 'NONE', 'NAN', 'NA'):
                        mapping[gene_id] = symbol
        except Exception as e:
            print(f'    Warning reading {fpath}: {e}')
    return mapping


def load_all_mappings(base_dir):
    global_map = {}
    for organ in ['Kidney', 'Liver', 'Lungs', 'Skin']:
        organ_path = os.path.join(base_dir, organ)
        if os.path.isdir(organ_path):
            sub_map = load_biodbnet_mappings(organ_path)
            global_map.update(sub_map)

    # Add GenBank accession map if available
    acc_map_file = os.path.join(base_dir, 'geo_cache', 'gb_acc_to_symbol_map.json')
    if os.path.exists(acc_map_file):
        try:
            with open(acc_map_file, 'r', encoding='utf-8') as f:
                acc_map = json.load(f)
                for k, v in acc_map.items():
                    if v and str(v).strip().upper() not in ('-', 'NONE', 'NAN', 'NA'):
                        global_map[str(k).strip()] = str(v).strip().upper()
                        global_map[str(k).strip().upper()] = str(v).strip().upper()
            print(f'  Loaded {len(acc_map)} GenBank accession mappings')
        except Exception as e:
            print(f'  Warning loading acc map: {e}')

    return global_map


def process_dataset_vectorized(tsv_path, bioDBnet_map):
    try:
        df = pd.read_csv(tsv_path, sep='	')
    except Exception as e:
        print(f'    Failed to read {tsv_path}: {e}')
        return pd.DataFrame()

    df.columns = [c.strip() for c in df.columns]

    # Find columns case-insensitively
    cols_lower = {c.lower(): c for c in df.columns}
    
    adj_p_col = next((cols_lower[c] for c in ['adj.p.val', 'padj', 'adj_p_value', 'adj_pval', 'adj.p.value', 'fdr'] if c in cols_lower), None)
    p_col = next((cols_lower[c] for c in ['p.value', 'pvalue', 'p_value', 'pvalue', 'pval'] if c in cols_lower), None)
    lfc_col = next((cols_lower[c] for c in ['logfc', 'log2foldchange', 'lfc'] if c in cols_lower), None)
    id_col = next((cols_lower[c] for c in ['id', 'geneid', 'probeid', 'probe_id'] if c in cols_lower), None)
    sym_col = next((cols_lower[c] for c in ['symbol', 'gene.symbol', 'gene_symbol', 'genesymbol'] if c in cols_lower), None)
    gi_col = cols_lower.get('gi', None)
    gb_col = next((cols_lower[c] for c in ['gb_acc', 'gb_list', 'genbank', 'acc'] if c in cols_lower), None)

    if not all([adj_p_col, p_col, lfc_col, id_col]):
        print(f'    Warning: skipping {os.path.basename(tsv_path)}, missing required cols. Found: {list(df.columns)}')
        return pd.DataFrame()

    work = pd.DataFrame({
        'raw_id': df[id_col].astype(str),
        'adj_p_value': pd.to_numeric(df[adj_p_col], errors='coerce'),
        'p_value': pd.to_numeric(df[p_col], errors='coerce'),
        'logFC': pd.to_numeric(df[lfc_col], errors='coerce'),
    })

    if sym_col is not None:
        work['sym_direct'] = df[sym_col].astype(str).str.strip().str.upper()
        work.loc[work['sym_direct'].isin(['', 'NAN', 'NA', '-', 'NONE', 'NULL']), 'sym_direct'] = np.nan
    else:
        work['sym_direct'] = np.nan

    def safe_map(series, mapping):
        mapped = series.map(mapping)
        return mapped.astype('object').where(mapped.isna(), mapped.astype(str).str.upper())

    # Map GI
    if gi_col is not None:
        work['gi_str'] = df[gi_col].astype(str).str.strip().replace(['', 'NAN', 'NA', '-', 'NONE'], np.nan)
        work['sym_gi'] = safe_map(work['gi_str'], bioDBnet_map)
    else:
        work['sym_gi'] = np.nan

    # Map GB_ACC
    if gb_col is not None:
        work['gb_str'] = df[gb_col].astype(str).str.strip().replace(['', 'NAN', 'NA', '-', 'NONE'], np.nan)
        work['sym_gb'] = safe_map(work['gb_str'], bioDBnet_map)
    else:
        work['sym_gb'] = np.nan

    # Map probe ID
    work['id_clean'] = work['raw_id'].astype(str).str.strip()
    work['id_no_at'] = work['id_clean'].where(~work['id_clean'].str.endswith('_at'), work['id_clean'].str[:-3])
    work['sym_at'] = safe_map(work['id_no_at'], bioDBnet_map)
    work['sym_id_str'] = safe_map(work['id_clean'], bioDBnet_map)

    # Hierarchical symbol selection
    work['gene'] = work['sym_direct']
    work['gene'] = work['gene'].fillna(work['sym_gb'])
    work['gene'] = work['gene'].fillna(work['sym_gi'])
    work['gene'] = work['gene'].fillna(work['sym_at'])
    work['gene'] = work['gene'].fillna(work['sym_id_str'])

    work = work.dropna(subset=['gene', 'adj_p_value', 'logFC'])
    if len(work) == 0:
        return pd.DataFrame()

    work['gene'] = work['gene'].astype(str).str.strip()
    work = work[~work['gene'].isin(['', 'NAN', 'NA', '-', 'NONE', 'NULL'])]

    # Deduplicate within dataset keeping the best adj_p_value
    work = work.sort_values('adj_p_value', ascending=True)
    work = work.drop_duplicates(subset='gene', keep='first')

    return work[['gene', 'logFC', 'p_value', 'adj_p_value']].reset_index(drop=True)


def aggregate_organ_datasets(organ_name, organ_base_dir, bioDBnet_map):
    print(f'\nProcessing {organ_name} Discovery Datasets...')
    all_dfs = []

    pattern = os.path.join(organ_base_dir, '**', '*.top.table.tsv')
    all_tsv_files = sorted(glob.glob(pattern, recursive=True))

    tsv_files = []
    excludes = discovery_config.EXCLUDE_FROM_DISCOVERY.get(organ_name, set())
    for f in all_tsv_files:
        if any(exc in f for exc in excludes):
            print(f'  [HELD OUT FOR VALIDATION] {os.path.basename(f)}')
        else:
            tsv_files.append(f)

    print(f'  Found {len(tsv_files)} Discovery dataset TSV files:')
    for tsv_path in tsv_files:
        ds_name = os.path.basename(os.path.dirname(tsv_path))
        print(f'    - {ds_name}/{os.path.basename(tsv_path)}')
        ds_df = process_dataset_vectorized(tsv_path, bioDBnet_map)
        if len(ds_df) > 0:
            print(f'      -> {len(ds_df)} valid genes mapped')
            all_dfs.append(ds_df)
        else:
            print(f'      -> 0 valid genes mapped')

    if not all_dfs:
        print(f'  [ERROR] No valid discovery data for {organ_name}')
        return pd.DataFrame()

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f'  Combined across {len(all_dfs)} discovery cohorts: {len(combined)} records')

    # De-duplicate across cohorts within organ (best adj_p_value per gene)
    combined = combined.sort_values('adj_p_value', ascending=True)
    combined = combined.drop_duplicates(subset='gene', keep='first').reset_index(drop=True)

    print(f'  Final unique discovery genes for {organ_name}: {len(combined)}')
    return combined


def main():
    base_dir = r'D:\CSIR'
    results_dir = os.path.join(base_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    organs = [
        ('Kidney', 'kidney_DEGs.csv'),
        ('Liver',  'liver_DEGs.csv'),
        ('Lungs',  'lung_DEGs.csv'),
        ('Skin',   'skin_DEGs.csv'),
    ]

    print('=' * 75)
    print('LOCKED DISCOVERY PIPELINE: GENERATING PER-ORGAN DEG CSVS')
    print('=' * 75)

    global_bio_map = load_all_mappings(base_dir)
    print(f'Total global mapping dictionary entries: {len(global_bio_map)}')

    organ_results = {}
    for organ_dirname, output_csv in organs:
        organ_path = os.path.join(base_dir, organ_dirname)
        if not os.path.isdir(organ_path):
            continue

        organ_df = aggregate_organ_datasets(organ_dirname, organ_path, global_bio_map)
        if len(organ_df) > 0:
            out_path = os.path.join(results_dir, output_csv)
            organ_df.to_csv(out_path, index=False)
            print(f'  [SAVED] {out_path} ({len(organ_df)} unique genes)')
            organ_results[organ_dirname] = organ_df

    print('\n' + '=' * 75)
    print('DISCOVERY SUMMARY (SIGNIFICANT AT adj_p < 0.05 & |logFC| >= 0.585)')
    print('=' * 75)
    for name, df in organ_results.items():
        sig = df[(df['adj_p_value'] < 0.05) & (df['logFC'].abs() >= 0.585)]
        print(f'  {name:8s}: {len(df):6d} total genes, {len(sig):6d} significant DEGs')


if __name__ == '__main__':
    main()
