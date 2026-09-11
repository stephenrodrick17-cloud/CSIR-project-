#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fix liver discovery DEGs by EXCLUDING the validation dataset GSE162694.
Then re-run pan-fibrotic core gene discovery, validation pipeline,
and ECM vs Non-ECM comparison charts with the CORRECTED independent data.
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import glob

# Import helpers directly from the preprocessing script (so we use identical logic)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from preprocess_build_deg_csvs import (
    load_biodbnet_mappings,
    process_dataset_vectorized,
)


BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "results")
os.makedirs(RESULTS, exist_ok=True)

# Validation accessions to EXCLUDE from discovery
VAL_EXCLUDES = {
    "kidney": {"GSE200818"},
    "liver":  {"GSE162694"},
    "lung":   {"GSE24206"},
    "skin":   {"GSE58095"},
}

ORGAN_DIRS = {
    "kidney": "Kidney",
    "liver":  "Liver",
    "lung":   "Lungs",
    "skin":   "Skin",
}


def build_organ_discovery_clean(organ_name, exclude_accessions, global_bio_map):
    """Build discovery DEG table for one organ, excluding specified GSE IDs."""
    organ_dirname = ORGAN_DIRS[organ_name]
    organ_path = os.path.join(BASE, organ_dirname)

    print(f"\n=== {organ_name.upper()} (clean discovery) ===")
    print(f"  Excluding accessions: {exclude_accessions}")

    pattern = os.path.join(organ_path, "**", "*.top.table.tsv")
    all_tsvs = sorted(glob.glob(pattern, recursive=True))

    kept_tsvs = []
    for t in all_tsvs:
        excluded = False
        for acc in exclude_accessions:
            if acc in os.path.basename(t):
                excluded = True
                print(f"  [EXCLUDE] {t}")
                break
        if not excluded:
            kept_tsvs.append(t)

    print(f"  Kept {len(kept_tsvs)} discovery TSVs (out of {len(all_tsvs)} total)")
    all_dfs = []
    for tsv_path in kept_tsvs:
        ds_name = os.path.basename(os.path.dirname(tsv_path))
        print(f"    Dataset {ds_name}/{os.path.basename(tsv_path)}")
        ds_df = process_dataset_vectorized(tsv_path, global_bio_map)
        if len(ds_df) > 0:
            print(f"      -> {len(ds_df)} genes")
            all_dfs.append(ds_df)
        else:
            print(f"      -> no valid genes")

    if not all_dfs:
        print(f"  [ERROR] No valid data for {organ_name} after exclusions!")
        return pd.DataFrame()

    combined = pd.concat(all_dfs, ignore_index=True)
    combined = combined.sort_values("adj_p_value", ascending=True)
    combined = combined.drop_duplicates(subset="gene", keep="first")
    combined = combined.sort_values("adj_p_value", ascending=True).reset_index(drop=True)

    print(f"  Combined {len(combined)} unique genes (best adj_p per gene)")
    return combined


def main():
    # ---------------------------------------------------------------
    # Step 1. Rebuild per-organ discovery CSVs, excluding validation GSEs
    # ---------------------------------------------------------------
    print("=" * 72)
    print("REBUILDING DISCOVERY DEGs (excluding validation GEO accessions)")
    print("=" * 72)

    print("\nLoading global bioDBnet map ...")
    global_bio_map = {}
    for organ_dirname in ORGAN_DIRS.values():
        organ_path = os.path.join(BASE, organ_dirname)
        if os.path.isdir(organ_path):
            sub_map = load_biodbnet_mappings(organ_path)
            global_bio_map.update(sub_map)
    print(f"  {len(global_bio_map)} Gene ID -> Symbol mappings loaded")

    new_degs = {}
    for organ_name in ["kidney", "liver", "lung", "skin"]:
        degs = build_organ_discovery_clean(
            organ_name, VAL_EXCLUDES[organ_name], global_bio_map
        )
        if len(degs) == 0:
            print(f"[FATAL] {organ_name} returned 0 genes. Quitting.")
            sys.exit(1)
        out = os.path.join(RESULTS, f"{organ_name}_DEGs.csv")
        degs.to_csv(out, index=False)
        new_degs[organ_name] = degs
        sig = degs[(degs["adj_p_value"] < 0.05) & (degs["logFC"].abs() > 0.585)]
        print(f"  [SAVED] {out}: {len(degs):,} genes, {len(sig):,} significant")

    # Quick sanity: confirm liver discovery vs validation no longer identical
    print("\nQuick post-fix check — liver disc vs val correlation:")
    new_liver = pd.read_csv(os.path.join(RESULTS, "liver_DEGs.csv"))
    val_liver = pd.read_csv(os.path.join(BASE, "validation", "liver_validation_DEGs.csv"))
    merged = new_liver.merge(val_liver, on="gene", suffixes=("_disc", "_val"))
    from scipy import stats as sstats
    r, p = sstats.pearsonr(merged["logFC_disc"].astype(float),
                           merged["logFC_val"].astype(float))
    same_lfc = np.isclose(merged["logFC_disc"].astype(float),
                          merged["logFC_val"].astype(float),
                          atol=1e-6).sum()
    print(f"  Overlap genes: {len(merged)}")
    print(f"  Identical logFC: {same_lfc} / {len(merged)}  ({same_lfc/len(merged)*100:.2f}%)")
    print(f"  Pearson r: {r:.4f}  (was 0.94 before fix — expect ~0.3-0.7 now)")
    if r < 0.9:
        print("  -> Liver independence restored ✓")
    else:
        print("  -> WARNING: liver r still suspiciously high. Check exclusions manually.")


if __name__ == "__main__":
    main()
