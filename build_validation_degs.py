#!/usr/bin/env python3
"""
Generate Validation DEG CSV Files
==================================
For each organ, build {organ}_validation_DEGs.csv from the specified
independent GEO validation cohort (NOT used in discovery):

  kidney = GSE200818   (Affymetrix + bioDBnet mapping)
  liver  = GSE162694   (DESeq2 output, has Symbol directly)
  lung   = GSE24206    (Affymetrix + Gene.symbol column directly)
  skin   = GSE58095    (Illumina ILMN IDs + GI column + Gene.symbol)

Applies the SAME probe-to-gene aggregation logic as the discovery pipeline:
  * map probes to gene symbols via direct columns or bioDBnet
  * for a gene with multiple probes, keep the row with the smallest adj_p_value
  * output columns: gene, logFC, p_value, adj_p_value
"""

import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "validation")
os.makedirs(OUTDIR, exist_ok=True)

# Load the GLOBAL bioDBnet map (same approach as discovery pipeline)
def load_biodbnet_all():
    mapping = {}
    for root, _, files in os.walk(BASE):
        for fn in files:
            if fn.startswith("bioDBnet") and fn.endswith(".txt"):
                fpath = os.path.join(root, fn)
                try:
                    df = pd.read_csv(fpath, sep="\t")
                    df.columns = [c.strip() for c in df.columns]
                    if "Gene ID" in df.columns and "Gene Symbol" in df.columns:
                        df = df.dropna(subset=["Gene Symbol"])
                        df["Gene Symbol"] = df["Gene Symbol"].astype(str).str.strip()
                        filt = df[
                            (df["Gene Symbol"] != "") &
                            (df["Gene Symbol"] != "-") &
                            (df["Gene Symbol"].str.lower() != "nan")
                        ]
                        for gid, sym in zip(filt["Gene ID"].astype(str).str.strip(),
                                            filt["Gene Symbol"]):
                            mapping[gid] = sym
                        try:
                            gi = filt["Gene ID"].astype(float).astype(int).astype(str)
                            for g, s in zip(gi, filt["Gene Symbol"]):
                                mapping.setdefault(g, s)
                        except Exception:
                            pass
                except Exception as e:
                    pass
    return mapping

bio_map = load_biodbnet_all()
print(f"Loaded global bioDBnet map: {len(bio_map)} gene ID -> symbol entries")


def safe_map_upper(series, mapping):
    mapped = series.map(mapping)
    if mapped.dtype == object:
        return mapped.str.upper()
    return mapped.astype("object").where(mapped.isna(), mapped.astype(str).str.upper())


def process_dataset(tsv_path):
    """
    Generic DEG table processor: reads a .top.table.tsv (or .tsv) file,
    maps probe/gene IDs to uppercase gene symbols, and returns a DataFrame
    with columns: gene, logFC, p_value, adj_p_value (one row per unique gene,
    best adj_p kept).
    """
    print(f"\n  Reading {os.path.basename(tsv_path)} ...")
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        print(f"    FAILED to read: {e}")
        return pd.DataFrame()

    df.columns = [c.strip() for c in df.columns]

    # Detect columns
    adj_p_col = next((c for c in ["adj.P.Val", "padj", "adj_p_value", "adj_pval",
                                   "adj.P.Value", "FDR"] if c in df.columns), None)
    p_col     = next((c for c in ["P.Value", "pvalue", "p_value", "PValue",
                                   "pval"] if c in df.columns), None)
    lfc_col   = next((c for c in ["logFC", "log2FoldChange", "LFC", "lfc"]
                      if c in df.columns), None)
    id_col    = next((c for c in ["ID", "GeneID", "Gene ID", "ProbeID", "probe_id",
                                   "id"] if c in df.columns), None)
    sym_col   = next((c for c in ["Symbol", "Gene.symbol", "Gene_symbol",
                                   "gene_symbol", "Gene Symbol"]
                      if c in df.columns), None)
    gi_col    = "GI" if "GI" in df.columns else None

    print(f"    Columns detected: adj_p={adj_p_col}, p={p_col}, logFC={lfc_col}, "
          f"ID={id_col}, Symbol={sym_col}, GI={gi_col}")

    missing = [n for n, c in [("adj_p", adj_p_col), ("p", p_col),
                                ("logFC", lfc_col), ("ID", id_col)] if c is None]
    if missing:
        print(f"    SKIP: missing required columns: {missing}")
        return pd.DataFrame()

    work = pd.DataFrame({
        "raw_id":      df[id_col].astype(str),
        "adj_p_value": pd.to_numeric(df[adj_p_col], errors="coerce"),
        "p_value":     pd.to_numeric(df[p_col],     errors="coerce"),
        "logFC":       pd.to_numeric(df[lfc_col],   errors="coerce"),
    })

    # --- 1) Direct symbol column ---
    if sym_col is not None:
        work["sym_direct"] = (df[sym_col].astype(str).str.strip()
                                .str.upper()
                                .replace(["", "NAN", "NA", "-", "NONE"], np.nan))
    else:
        work["sym_direct"] = np.nan

    # --- 2) GI column (Illumina) ---
    if gi_col is not None:
        work["gi_str"] = df[gi_col].astype(str).str.strip().replace(
            ["", "NAN", "NA", "-", "NONE"], np.nan)
        work["sym_gi"] = safe_map_upper(work["gi_str"], bio_map)
        try:
            work["gi_int"] = (pd.to_numeric(df[gi_col], errors="coerce")
                                .astype("Int64").astype(str)
                                .replace("<NA>", np.nan))
            mask = work["sym_gi"].isna() & work["gi_int"].notna()
            work.loc[mask, "sym_gi"] = safe_map_upper(work.loc[mask, "gi_int"], bio_map)
        except Exception:
            pass
    else:
        work["sym_gi"] = np.nan

    # --- 3) ID with _at stripped (Affymetrix) ---
    work["id_clean"] = work["raw_id"].astype(str).str.strip()
    work["id_no_at"] = work["id_clean"].where(
        ~work["id_clean"].str.endswith("_at"), work["id_clean"].str[:-3])
    work["sym_at"] = safe_map_upper(work["id_no_at"], bio_map)
    try:
        work["id_int"] = (pd.to_numeric(work["id_no_at"], errors="coerce")
                            .astype("Int64").astype(str)
                            .replace("<NA>", np.nan))
        mask = work["sym_at"].isna() & work["id_int"].notna()
        work.loc[mask, "sym_at"] = safe_map_upper(work.loc[mask, "id_int"], bio_map)
    except Exception:
        pass

    # --- 4) General ID mapping (string + int) ---
    work["sym_id_str"] = safe_map_upper(work["id_clean"], bio_map)
    try:
        work["id_generic_int"] = (pd.to_numeric(work["id_clean"], errors="coerce")
                                    .astype("Int64").astype(str)
                                    .replace("<NA>", np.nan))
        mask = work["sym_id_str"].isna() & work["id_generic_int"].notna()
        work.loc[mask, "sym_id_str"] = safe_map_upper(
            work.loc[mask, "id_generic_int"], bio_map)
    except Exception:
        pass

    # Hierarchical fallback
    work["gene"] = work["sym_direct"]
    work["gene"] = work["gene"].fillna(work["sym_gi"])
    work["gene"] = work["gene"].fillna(work["sym_at"])
    work["gene"] = work["gene"].fillna(work["sym_id_str"])

    work = work.dropna(subset=["gene", "adj_p_value", "logFC"])
    if len(work) == 0:
        print("    No valid genes after symbol mapping!")
        return pd.DataFrame()

    work["gene"] = work["gene"].astype(str).str.strip()
    work = work[work["gene"] != ""]

    # Keep most significant row per gene
    work = work.sort_values("adj_p_value", ascending=True)
    work = work.drop_duplicates(subset="gene", keep="first")
    work = work[["gene", "logFC", "p_value", "adj_p_value"]].reset_index(drop=True)
    print(f"    -> {len(work):,} unique genes retained")
    return work


# ---- Validation dataset spec ----
validation_specs = [
    ("kidney", os.path.join(BASE, "Kidney", "18", "GSE200818.top.table.tsv")),
    ("liver",  os.path.join(BASE, "Liver",  "694", "GSE162694.top.table.tsv")),
    ("lung",   os.path.join(BASE, "Lungs",  "206",
                            "GSE24206.top.table (2).tsv")),
    ("skin",   os.path.join(BASE, "Skin",   "095", "GSE58095.top.table.tsv")),
]

print("\n" + "="*72)
print("GENERATING VALIDATION DEG CSVs")
print("="*72)

for organ_name, tsv_path in validation_specs:
    print(f"\n=== {organ_name.upper()} (validation) ===")
    if not os.path.isfile(tsv_path):
        print(f"  [ERROR] File not found: {tsv_path}")
        continue
    deg_df = process_dataset(tsv_path)
    if len(deg_df) == 0:
        continue
    out_csv = os.path.join(OUTDIR, f"{organ_name}_validation_DEGs.csv")
    deg_df.to_csv(out_csv, index=False)
    print(f"  [SAVED] {out_csv}  ({len(deg_df):,} rows)")
    sig = deg_df[(deg_df["adj_p_value"] < 0.05) & (deg_df["logFC"].abs() > 0.585)]
    print(f"  Significant (adj_p<0.05 & |logFC|>0.585): {len(sig):,} genes")

print("\nDONE.")
