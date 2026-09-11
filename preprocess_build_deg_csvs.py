#!/usr/bin/env python3
"""
Preprocessing script (VECTORIZED) to build per-organ DEG CSV files from raw GEO datasets.

Handles:
- Affymetrix arrays (IDs with _at suffix -> bioDBnet Gene ID -> Gene Symbol)
- Illumina arrays (IDs with ILMN_ prefix; GI column + Gene.symbol column)
- DESeq2-style output with Symbol column
- Aggregation: for duplicate genes, keep the most significant (smallest adj_p_value)
"""

import os
import glob
import pandas as pd
import numpy as np

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def load_biodbnet_mappings(organ_dir):
    """
    Load all bioDBnet mapping files in a given organ directory tree.
    Returns a dict: {Gene ID (str) -> Gene Symbol (str)}
    Also includes integer-keyed versions for flexible matching.
    """
    mapping = {}
    pattern = os.path.join(organ_dir, "**", "bioDBnet*.txt")
    for fpath in glob.glob(pattern, recursive=True):
        try:
            df = pd.read_csv(fpath, sep="\t")
            df.columns = [c.strip() for c in df.columns]
            if "Gene ID" in df.columns and "Gene Symbol" in df.columns:
                df = df.dropna(subset=["Gene Symbol"])
                df["Gene Symbol"] = df["Gene Symbol"].astype(str).str.strip()
                df["Gene ID_str"] = df["Gene ID"].astype(str).str.strip()
                valid = df[
                    (df["Gene Symbol"] != "") &
                    (df["Gene Symbol"] != "-") &
                    (df["Gene Symbol"].str.lower() != "nan")
                ]
                for gid_str, sym in zip(valid["Gene ID_str"], valid["Gene Symbol"]):
                    mapping[gid_str] = sym
                try:
                    df["Gene ID_int"] = df["Gene ID"].astype(float).astype(int).astype(str)
                    for gid_int, sym in zip(df["Gene ID_int"], df["Gene Symbol"]):
                        mapping.setdefault(gid_int, sym)
                except (ValueError, TypeError):
                    pass
        except Exception as e:
            print(f"    Warning: could not parse {os.path.basename(fpath)}: {e}")
    return mapping


def process_dataset_vectorized(tsv_path, bioDBnet_map):
    """
    Process a single top.table.tsv using VECTORIZED operations.
    Returns a DataFrame: gene, logFC, p_value, adj_p_value  (one row per unique gene)
    """
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        print(f"    Failed to read {tsv_path}: {e}")
        return pd.DataFrame()

    df.columns = [c.strip() for c in df.columns]

    # --- Identify required columns ---
    adj_p_col = next((c for c in ["adj.P.Val", "padj", "adj_p_value", "adj_pval", "adj.P.Value", "FDR"] if c in df.columns), None)
    p_col = next((c for c in ["P.Value", "pvalue", "p_value", "PValue", "pval"] if c in df.columns), None)
    lfc_col = next((c for c in ["logFC", "log2FoldChange", "LFC", "lfc"] if c in df.columns), None)
    id_col = next((c for c in ["ID", "GeneID", "ProbeID", "probe_id", "id"] if c in df.columns), None)
    sym_col = next((c for c in ["Symbol", "Gene.symbol", "Gene_symbol", "gene_symbol"] if c in df.columns), None)
    gi_col = "GI" if "GI" in df.columns else None

    if not all([adj_p_col, p_col, lfc_col, id_col]):
        print(f"    Warning: skipping {os.path.basename(tsv_path)}, missing cols. Found: {list(df.columns)}")
        return pd.DataFrame()

    work = pd.DataFrame({
        "raw_id": df[id_col].astype(str),
        "adj_p_value": pd.to_numeric(df[adj_p_col], errors="coerce"),
        "p_value": pd.to_numeric(df[p_col], errors="coerce"),
        "logFC": pd.to_numeric(df[lfc_col], errors="coerce"),
    })

    if sym_col is not None:
        work["sym_direct"] = df[sym_col].astype(str).str.strip().str.upper()
        work.loc[work["sym_direct"].isin(["", "NAN", "NA", "-", "NONE"]), "sym_direct"] = np.nan
    else:
        work["sym_direct"] = np.nan

    def safe_map_upper(series, mapping):
        mapped = series.map(mapping)
        if mapped.dtype == object:
            return mapped.str.upper()
        return mapped.astype("object").where(mapped.isna(), mapped.astype(str).str.upper())

    if gi_col is not None:
        work["gi_str"] = df[gi_col].astype(str).str.strip()
        work["gi_str"] = work["gi_str"].replace(["", "NAN", "NA", "-", "NONE"], np.nan)
        work["sym_gi"] = safe_map_upper(work["gi_str"], bioDBnet_map)
        try:
            work["gi_int"] = pd.to_numeric(df[gi_col], errors="coerce").astype("Int64").astype(str)
            work.loc[work["gi_int"] == "<NA>", "gi_int"] = np.nan
            mask = work["sym_gi"].isna() & work["gi_int"].notna()
            work.loc[mask, "sym_gi"] = safe_map_upper(work.loc[mask, "gi_int"], bioDBnet_map)
        except Exception:
            pass
    else:
        work["sym_gi"] = np.nan

    work["id_clean"] = work["raw_id"].astype(str).str.strip()

    work["id_no_at"] = work["id_clean"].where(~work["id_clean"].str.endswith("_at"), work["id_clean"].str[:-3])
    work["sym_at"] = safe_map_upper(work["id_no_at"], bioDBnet_map)
    try:
        work["id_int"] = pd.to_numeric(work["id_no_at"], errors="coerce").astype("Int64").astype(str)
        work.loc[work["id_int"] == "<NA>", "id_int"] = np.nan
        mask = work["sym_at"].isna() & work["id_int"].notna()
        work.loc[mask, "sym_at"] = safe_map_upper(work.loc[mask, "id_int"], bioDBnet_map)
    except Exception:
        pass

    work["sym_id_str"] = safe_map_upper(work["id_clean"], bioDBnet_map)
    try:
        work["id_generic_int"] = pd.to_numeric(work["id_clean"], errors="coerce").astype("Int64").astype(str)
        work.loc[work["id_generic_int"] == "<NA>", "id_generic_int"] = np.nan
        mask = work["sym_id_str"].isna() & work["id_generic_int"].notna()
        work.loc[mask, "sym_id_str"] = safe_map_upper(work.loc[mask, "id_generic_int"], bioDBnet_map)
    except Exception:
        pass

    work["gene"] = work["sym_direct"]
    work["gene"] = work["gene"].fillna(work["sym_gi"])
    work["gene"] = work["gene"].fillna(work["sym_at"])
    work["gene"] = work["gene"].fillna(work["sym_id_str"])

    work = work.dropna(subset=["gene", "adj_p_value", "logFC"])

    if len(work) == 0:
        return pd.DataFrame()

    work["gene"] = work["gene"].astype(str).str.strip()
    work = work[work["gene"] != ""]

    work = work.sort_values("adj_p_value", ascending=True)
    work = work.drop_duplicates(subset="gene", keep="first")

    return work[["gene", "logFC", "p_value", "adj_p_value"]].reset_index(drop=True)


def aggregate_organ_datasets(organ_name, organ_base_dir, bioDBnet_map):
    """
    Process all datasets for a single organ, aggregate (best adj_p per gene).
    """
    print(f"\nProcessing {organ_name}...")
    all_dfs = []

    pattern = os.path.join(organ_base_dir, "**", "*.top.table.tsv")
    tsv_files = sorted(glob.glob(pattern, recursive=True))
    print(f"  Found {len(tsv_files)} dataset TSV files")

    for tsv_path in tsv_files:
        ds_name = os.path.basename(os.path.dirname(tsv_path))
        print(f"    Dataset {ds_name}/{os.path.basename(tsv_path)}")
        ds_df = process_dataset_vectorized(tsv_path, bioDBnet_map)
        if len(ds_df) > 0:
            print(f"      -> {len(ds_df)} genes with valid symbols")
            all_dfs.append(ds_df)
        else:
            print(f"      -> no valid genes (check gene symbol mapping)")

    if not all_dfs:
        print(f"  [ERROR] No valid data for {organ_name}")
        return pd.DataFrame()

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"  Combined: {len(combined)} gene records across all datasets")

    combined = combined.sort_values("adj_p_value", ascending=True)
    combined = combined.drop_duplicates(subset="gene", keep="first")
    combined = combined.sort_values("adj_p_value", ascending=True).reset_index(drop=True)

    print(f"  After de-duplication (best p-value per gene): {len(combined)} unique genes")
    return combined


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(base_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    organs = [
        ("Kidney", "kidney_DEGs.csv"),
        ("Liver", "liver_DEGs.csv"),
        ("Lungs", "lung_DEGs.csv"),
        ("Skin", "skin_DEGs.csv"),
    ]

    print("=" * 70)
    print("BUILDING PER-ORGAN DEG CSV FILES FROM RAW GEO DATASETS")
    print("=" * 70)

    print("\nLoading GLOBAL bioDBnet mappings (across all organs)...")
    global_bio_map = {}
    for organ_dirname, _ in organs:
        organ_path = os.path.join(base_dir, organ_dirname)
        if os.path.isdir(organ_path):
            sub_map = load_biodbnet_mappings(organ_path)
            global_bio_map.update(sub_map)
    print(f"  Total loaded: {len(global_bio_map)} Gene ID -> Symbol mappings")

    organ_results = {}
    for organ_dirname, output_csv in organs:
        organ_path = os.path.join(base_dir, organ_dirname)
        if not os.path.isdir(organ_path):
            print(f"\n[SKIP] Organ directory not found: {organ_path}")
            continue

        print(f"\nUsing global bioDBnet map for {organ_dirname}")
        organ_df = aggregate_organ_datasets(organ_dirname, organ_path, global_bio_map)

        if len(organ_df) > 0:
            out_path = os.path.join(results_dir, output_csv)
            organ_df.to_csv(out_path, index=False)
            print(f"  [SAVED] {out_path}  ({len(organ_df)} genes)")
            organ_results[organ_dirname] = organ_df
        else:
            print(f"  [SKIP] No genes to save for {organ_dirname}")

    print("\n" + "=" * 70)
    print("PREPROCESSING COMPLETE")
    print("=" * 70)
    for name, df in organ_results.items():
        sig = df[(df["adj_p_value"] < 0.05) & (df["logFC"].abs() > 0.585)]
        print(f"  {name}: {len(df)} total genes, {len(sig)} significant "
              f"(adj_p<0.05 & |logFC|>0.585)")


if __name__ == "__main__":
    main()
