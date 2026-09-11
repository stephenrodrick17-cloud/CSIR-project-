#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liver dataset independence check.

Given the near-perfect discovery-vs-validation correlation seen in the liver
results (r = 0.97 for ECM genes, points sitting almost exactly on y = x),
this script checks for two possible causes:

  1. METADATA OVERLAP  — the discovery and validation GEO series actually
     share the same patients/samples (common when a lab deposits the same
     cohort under two different accessions, e.g. one for the array data and
     one for a later reanalysis).

  2. DATA-LEVEL DUPLICATION — a coding/loading bug caused the same
     expression matrix (or the same subset of samples) to be used for both
     "discovery" and "validation" without you intending it.

Fill in the two accession numbers you actually used for liver discovery
and liver validation below, then run this script.
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np

try:
    import GEOparse
    HAS_GEO = True
except ImportError:
    HAS_GEO = False
    print("[WARN] GEOparse not installed. Will run offline data-only check "
          "against your local DEG CSVs instead.")

# ------------------------------------------------------------------
# EDIT THESE TWO LINES to match what you actually used for liver
# ------------------------------------------------------------------
DISCOVERY_ACCESSION = "GSE77627"      # <-- change if different
VALIDATION_ACCESSION = "GSE162694"    # <-- change if different
DEST_DIR = "./data/liver"
# ------------------------------------------------------------------

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEST_DIR_ABS = os.path.join(BASE_DIR, "data", "liver")
os.makedirs(DEST_DIR_ABS, exist_ok=True)


def load_series(accession, destdir):
    if not HAS_GEO:
        return None
    print(f"\nLoading {accession} via GEOparse (may take a minute on first download) ...")
    try:
        gse = GEOparse.get_GEO(geo=accession, destdir=destdir)
        print(f"  -> OK: {len(gse.gsms)} samples loaded")
        return gse
    except Exception as e:
        print(f"  -> GEOparse FAILED for {accession}: {e}")
        print("     (proceeding with offline CSV check instead)")
        return None


def get_sample_metadata(gse):
    """Pull every metadata field per sample into one dataframe, so we can
    compare titles, characteristics, and any patient/donor identifiers."""
    rows = []
    if gse is None:
        return pd.DataFrame()
    for gsm_name, gsm in gse.gsms.items():
        meta_flat = {}
        for key, val in gsm.metadata.items():
            meta_flat[key] = " | ".join(val) if isinstance(val, list) else str(val)
        meta_flat["sample_id"] = gsm_name
        rows.append(meta_flat)
    return pd.DataFrame(rows)


def check_metadata_overlap(meta_disc, meta_val):
    print("\n" + "=" * 70)
    print("STEP 1: METADATA-LEVEL OVERLAP CHECK")
    print("=" * 70)

    if len(meta_disc) == 0 or len(meta_val) == 0:
        print("  [SKIP] Metadata not available (GEOparse download failed). "
              "See offline CSV check in Step 2 below.")
        return

    disc_ids = set(meta_disc["sample_id"])
    val_ids = set(meta_val["sample_id"])
    id_overlap = disc_ids & val_ids
    print(f"Directly shared GSM sample IDs: {len(id_overlap)}")
    if id_overlap:
        print("  -> WARNING: identical sample accessions in both sets:", id_overlap)

    candidate_fields = [c for c in meta_disc.columns
                        if any(k in c.lower() for k in
                                ["subject", "patient", "donor", "individual", "characteristics", "title"])]
    print(f"\nChecking possible patient/donor identifier fields: {candidate_fields}")

    for field in candidate_fields:
        if field not in meta_val.columns:
            continue
        disc_vals = set(meta_disc[field].dropna().astype(str))
        val_vals = set(meta_val[field].dropna().astype(str))
        overlap = disc_vals & val_vals
        if overlap and len(overlap) > 0:
            frac = len(overlap) / min(len(disc_vals), len(val_vals))
            flag = "  <-- POSSIBLE OVERLAP" if frac > 0.05 else ""
            print(f"  Field '{field}': {len(overlap)} shared values "
                  f"({frac*100:.1f}% of smaller set){flag}")

    print("\nFirst 5 discovery sample titles:")
    if "title" in meta_disc.columns:
        print(meta_disc[["sample_id", "title"]].head(5).to_string(index=False))
    print("\nFirst 5 validation sample titles:")
    if "title" in meta_val.columns:
        print(meta_val[["sample_id", "title"]].head(5).to_string(index=False))
    print("\n--> Manually confirm these look like genuinely different patients/cohorts,")
    print("    not the same subjects renamed under a new accession.")


def get_expression_matrix(gse):
    if gse is None:
        return None
    frames = []
    for gsm_name, gsm in gse.gsms.items():
        if gsm.table is not None and "VALUE" in gsm.table.columns:
            frames.append(gsm.table.set_index(gsm.table.columns[0])["VALUE"].rename(gsm_name))
    if not frames:
        return None
    return pd.concat(frames, axis=1)


def check_data_level_duplication(expr_disc, expr_val):
    print("\n" + "=" * 70)
    print("STEP 2: DATA-LEVEL DUPLICATION CHECK")
    print("=" * 70)

    # ---- ALWAYS: offline CSV check (the smoking-gun we already found) ----
    print("\n[2A] Offline CSV comparison: liver_DEGs.csv (discovery aggregate)")
    print("     vs liver_validation_DEGs.csv (GSE162694 standalone) ...")
    disc_csv = os.path.join(BASE_DIR, "results", "liver_DEGs.csv")
    val_csv  = os.path.join(BASE_DIR, "validation", "liver_validation_DEGs.csv")
    if os.path.isfile(disc_csv) and os.path.isfile(val_csv):
        d = pd.read_csv(disc_csv)
        v = pd.read_csv(val_csv)
        d["gene"] = d["gene"].astype(str).str.upper().str.strip()
        v["gene"] = v["gene"].astype(str).str.upper().str.strip()
        m = d.merge(v, on="gene", suffixes=("_disc", "_val"))
        print(f"  Genes in discovery liver DEGs: {len(d):,}")
        print(f"  Genes in validation liver DEGs: {len(v):,}")
        print(f"  Overlapping genes:             {len(m):,}")

        # Exact match check
        lfc_same = (np.isclose(m["logFC_disc"].values.astype(float),
                               m["logFC_val"].values.astype(float),
                               atol=1e-6, rtol=0)).sum()
        adjp_same = (np.isclose(m["adj_p_value_disc"].values.astype(float),
                                m["adj_p_value_val"].values.astype(float),
                                atol=1e-10, rtol=0)).sum()
        print(f"  Genes with IDENTICAL logFC:    {lfc_same:,} / {len(m):,}  "
              f"({lfc_same/len(m)*100:.2f}%)")
        print(f"  Genes with IDENTICAL adj_p:    {adjp_same:,} / {len(m):,}  "
              f"({adjp_same/len(m)*100:.2f}%)")

        if lfc_same / len(m) > 0.95:
            print("\n  *** CONFIRMED DATA-LEVEL DUPLICATION ***")
            print("  -> The liver discovery aggregate contains the EXACT same")
            print("     GSE162694 rows that are used for validation. This is a BUG")
            print("     in the discovery pipeline: the glob pattern picks up ALL")
            print("     Liver/**/*.top.table.tsv files including the GSE162694 that")
            print("     you reserved for validation. Fix: re-build discovery DEGs")
            print("     excluding Liver/694/GSE162694.top.table.tsv from the aggregate.")
        else:
            print("  No mass exact-match duplication found via CSV check.")

        # Correlation check of logFC
        if len(m) >= 10:
            from scipy import stats as sstats
            r, p = sstats.pearsonr(m["logFC_disc"].values.astype(float),
                                   m["logFC_val"].values.astype(float))
            print(f"  Pearson r (logFC) between disc & val: {r:.4f}  (p={p:.2e})")
            if r > 0.95:
                print("  -> Near-1.0 correlation is NOT expected between genuinely")
                print("     independent liver fibrosis cohorts. This reinforces the")
                print("     exact-match finding above.")
    else:
        print("  [SKIP] Required CSVs not found on disk.")

    # ---- Optional: GEOparse VALUE-matrix based sample cross-correlation ----
    print("\n[2B] GEOparse raw-expression cross-sample check ...")
    if expr_disc is None or expr_val is None:
        print("  [SKIP] Raw VALUE matrices not available from GEOparse. "
              "Relied on Step 2A instead.")
        return

    common_genes = expr_disc.index.intersection(expr_val.index)
    print(f"  Common probes/genes between the two series: {len(common_genes)}")
    if len(common_genes) < 10:
        print("  Too few common identifiers to compare reliably.")
        return

    disc_sub = expr_disc.loc[common_genes].astype(float)
    val_sub = expr_val.loc[common_genes].astype(float)

    corr_matrix = pd.DataFrame(
        np.corrcoef(disc_sub.T.values, val_sub.T.values)[:disc_sub.shape[1], disc_sub.shape[1]:],
        index=disc_sub.columns,
        columns=val_sub.columns,
    )

    max_corr = corr_matrix.values.max()
    max_loc = np.unravel_index(np.argmax(corr_matrix.values), corr_matrix.shape)
    disc_sample = corr_matrix.index[max_loc[0]]
    val_sample = corr_matrix.columns[max_loc[1]]

    print(f"\n  Highest single-sample cross-correlation: {max_corr:.4f}")
    print(f"    Between '{disc_sample}' (disc) and '{val_sample}' (val)")

    if max_corr > 0.98:
        print("    -> STRONG WARNING: duplicated/reused sample.")
    elif max_corr > 0.9:
        print("    -> Worth closer manual look (technical replicates?).")
    else:
        print("    -> No individual sample pair looks duplicated. Reassuring.")

    print(f"  Average cross-correlation: {corr_matrix.values.mean():.4f}")
    print(f"  Median  cross-correlation: {np.median(corr_matrix.values):.4f}")


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    dest_dir = os.path.join(base, "data", "liver")
    os.makedirs(dest_dir, exist_ok=True)

    print("=" * 70)
    print(f"LIVER DATASET INDEPENDENCE CHECK")
    print(f"  Discovery: {DISCOVERY_ACCESSION}  (plus GSE164760, GSE89377 in aggregate)")
    print(f"  Validation: {VALIDATION_ACCESSION}")
    print("=" * 70)

    # ---- GEOparse metadata check ----
    gse_disc = load_series(DISCOVERY_ACCESSION, dest_dir)
    gse_val  = load_series(VALIDATION_ACCESSION, dest_dir)
    meta_disc = get_sample_metadata(gse_disc)
    meta_val  = get_sample_metadata(gse_val)
    check_metadata_overlap(meta_disc, meta_val)

    # ---- Data-level duplication check (ALSO includes offline CSV check) ----
    expr_disc = get_expression_matrix(gse_disc)
    expr_val  = get_expression_matrix(gse_val)
    check_data_level_duplication(expr_disc, expr_val)

    print("\n" + "=" * 70)
    print("SUMMARY / WHAT TO DO NEXT")
    print("=" * 70)
    print("""
1. If Step 1 flagged shared patient/donor IDs, or Step 2 found a
   near-1.0 correlation pair: your two liver datasets are NOT independent.
   Do not use this pair as discovery+validation. Replace the validation
   set with a genuinely separate GEO accession before trusting the liver
   validation numbers.

2. If both steps come back clean but the r=0.97 pattern persists: check
   your own DEG script for a bug — e.g. accidentally loading the discovery
   CSV twice, or a copy-paste error where the validation DEG loop still
   points at the discovery dataframe/variable name.

3. Either way, do not proceed to Mendelian randomization or the final
   biomarker shortlist using the current liver validation numbers until
   this is resolved, since it inflates confidence in your overall 51-gene
   pan-fibrotic result.
""")


if __name__ == "__main__":
    main()
