#!/usr/bin/env python3
"""
Cross-Organ Fibrosis Biomarker Validation Pipeline (Steps 1-2)
==============================================================

Performs direction concordance validation of the 175-gene pan-fibrotic
discovery signature against independent GEO validation cohorts for
kidney, liver, lung, and skin.

Outputs:
  - results/pan_fibrotic_core_genes_validated.csv  (full annotated table)
  - results/final_validated_pan_fibrotic_genes.csv (all-4-organ validated)
  - Console summary of validation funnel statistics
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats


ORGANS = ["kidney", "liver", "lung", "skin"]
ADJ_P_THRESHOLD = 0.05


# =============================================================================
# I/O Helpers
# =============================================================================

def setup_paths(base_dir):
    """Return dictionary of standard project paths."""
    paths = {
        "base": base_dir,
        "results": os.path.join(base_dir, "results"),
        "validation": os.path.join(base_dir, "validation"),
        "core_genes": os.path.join(base_dir, "pan_fibrotic_core_genes.csv"),
        "out_validated": os.path.join(
            base_dir, "results", "pan_fibrotic_core_genes_validated.csv"
        ),
        "out_final": os.path.join(
            base_dir, "results", "final_validated_pan_fibrotic_genes.csv"
        ),
    }
    for key in ("results", "validation"):
        os.makedirs(paths[key], exist_ok=True)
    return paths


def load_core_genes(core_csv_path):
    """Load the 175 discovery pan-fibrotic core genes table."""
    df = pd.read_csv(core_csv_path)
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    print(f"  [LOAD] {len(df)} discovery core genes from "
          f"{os.path.basename(core_csv_path)}")
    return df


def load_validation_degs(validation_dir):
    """
    Load all four validation DEG tables into a dict of DataFrames.
    Keys: organ name -> DataFrame indexed by uppercase gene symbol.
    """
    degs = {}
    for organ in ORGANS:
        csv_path = os.path.join(validation_dir, f"{organ}_validation_DEGs.csv")
        if not os.path.isfile(csv_path):
            raise FileNotFoundError(f"Missing validation CSV: {csv_path}")
        df = pd.read_csv(csv_path)
        df["gene"] = df["gene"].astype(str).str.upper().str.strip()
        df = df.drop_duplicates(subset="gene", keep="first")
        df = df.set_index("gene")
        degs[organ] = df
        print(f"  [LOAD] {organ:6s}: {len(df):>6,} genes")
    return degs


# =============================================================================
# Step 1 — Direction Concordance Check
# =============================================================================

def direction_concordance(disc_logfc, val_logfc):
    """
    Return True if both values have the same sign (both positive or both
    negative), False if opposite sign, and NaN if either is missing.
    """
    if pd.isna(disc_logfc) or pd.isna(val_logfc):
        return np.nan
    if disc_logfc == 0 or val_logfc == 0:
        return False
    return (disc_logfc > 0) == (val_logfc > 0)


def annotate_validation_status(core_df, val_degs):
    """
    For each organ, add columns:
      - {organ}_val_logFC        : validation logFC (NaN if gene absent)
      - {organ}_val_adj_p_value  : validation adj_p_value (NaN if absent)
      - {organ}_validated        : direction concordant? (True/False/NA)
      - {organ}_validated_sig    : direction concordant AND adj_p < 0.05
    """
    df = core_df.copy()

    for organ in ORGANS:
        val_df = val_degs[organ]
        disc_lfc_col = f"{organ}_logFC"

        # Lookup validation values for every discovery gene
        genes = df["gene"].values
        present = [g in val_df.index for g in genes]
        val_logfc = np.where(
            present,
            val_df.reindex(genes)["logFC"].values,
            np.nan,
        )
        val_adjp = np.where(
            present,
            val_df.reindex(genes)["adj_p_value"].values,
            np.nan,
        )

        df[f"{organ}_val_logFC"] = val_logfc
        df[f"{organ}_val_adj_p_value"] = val_adjp

        # Direction concordance
        validated = []
        validated_sig = []
        for i, gene in enumerate(genes):
            disc_lfc = df.iloc[i][disc_lfc_col]
            v_lfc = val_logfc[i]
            v_adjp = val_adjp[i]
            concord = direction_concordance(disc_lfc, v_lfc)
            validated.append(concord)
            if pd.isna(concord) or not concord:
                validated_sig.append(False if not pd.isna(concord) else np.nan)
            else:
                validated_sig.append(
                    (not pd.isna(v_adjp)) and (v_adjp < ADJ_P_THRESHOLD)
                )

        df[f"{organ}_validated"] = validated
        df[f"{organ}_validated_sig"] = validated_sig

    # Derived counts
    val_cols = [f"{o}_validated" for o in ORGANS]
    sig_cols = [f"{o}_validated_sig" for o in ORGANS]

    df["n_organs_validated"] = df[val_cols].fillna(False).sum(axis=1).astype(int)
    df["n_organs_validated_sig"] = (
        df[sig_cols].fillna(False).sum(axis=1).astype(int)
    )

    n_found = df[[f"{o}_val_logFC" for o in ORGANS]].notna().sum(axis=1)
    df["n_organs_found_in_val"] = n_found.astype(int)

    return df


# =============================================================================
# Step 2 — Strict Validation Funnel
# =============================================================================

def build_funnel_table(df):
    """
    Return a DataFrame describing, step by step, how many genes pass each
    validation tier (both for direction-only and direction + significance).
    """
    total = len(df)
    rows = []

    for label, n_col in [
        ("Direction concordance only", "n_organs_validated"),
        ("Direction + adj_p < 0.05", "n_organs_validated_sig"),
    ]:
        for k in [4, 3, 2, 1, 0]:
            if k == 0:
                n_k = (df[n_col] >= 0).sum()  # always true; placeholder total
            else:
                n_k = (df[n_col] >= k).sum()
            pct = n_k / total * 100 if total else 0.0
            rows.append({
                "criteria": label,
                "tier": f">= {k} organ(s)",
                "gene_count": int(n_k) if k != 0 else int(total),
                "percentage": round(pct, 2) if k != 0 else 100.00,
            })
        # Replace the k=0 row with the correct total entry
        # (find last row for this label and overwrite)
    # Build cleaner rows
    rows = []
    total = len(df)
    for label, n_col in [
        ("Direction concordance only", "n_organs_validated"),
        ("Direction + adj_p < 0.05", "n_organs_validated_sig"),
    ]:
        rows.append({
            "criteria": label,
            "tier": "Total discovery genes",
            "gene_count": total,
            "percentage": 100.00,
        })
        for k in [1, 2, 3, 4]:
            n_k = int((df[n_col] >= k).sum())
            pct = round(n_k / total * 100, 2) if total else 0.0
            rows.append({
                "criteria": label,
                "tier": f"Validated in >= {k} organ(s)",
                "gene_count": n_k,
                "percentage": pct,
            })
        rows.append({
            "criteria": label,
            "tier": f"Validated in ALL 4 organs",
            "gene_count": int((df[n_col] == 4).sum()),
            "percentage": round((df[n_col] == 4).sum() / total * 100, 2)
            if total else 0.0,
        })

    funnel_df = pd.DataFrame(rows)
    return funnel_df


def per_organ_validation_rates(df):
    """Return per-organ validation statistics (concordance, significance)."""
    rows = []
    total = len(df)
    for organ in ORGANS:
        found = df[f"{organ}_val_logFC"].notna().sum()
        concordant = df[f"{organ}_validated"].fillna(False).sum()
        concordant_sig = df[f"{organ}_validated_sig"].fillna(False).sum()

        # Of genes found in validation dataset, what % concordant?
        pct_of_found = concordant / found * 100 if found else 0.0
        # Of all discovery genes (175)?
        pct_of_total = concordant / total * 100 if total else 0.0

        rows.append({
            "organ": organ,
            "total_discovery_genes": total,
            "genes_found_in_val": int(found),
            "direction_concordant": int(concordant),
            "concordant_and_sig_adjp": int(concordant_sig),
            "pct_of_found_concordant": round(pct_of_found, 2),
            "pct_of_all_discovery_concordant": round(pct_of_total, 2),
            "pct_of_all_discovery_concordant_sig": round(
                concordant_sig / total * 100, 2
            ) if total else 0.0,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Step 4 — Statistical Summary (correlation)
# =============================================================================

def compute_correlation_summary(core_df, val_degs):
    """
    For each organ compute Pearson & Spearman correlation between discovery
    logFC and validation logFC (across the 175 core genes that are present).
    Returns a DataFrame.
    """
    rows = []
    for organ in ORGANS:
        disc_lfc_col = f"{organ}_logFC"
        disc = core_df[["gene", disc_lfc_col]].copy().set_index("gene")
        val = val_degs[organ][["logFC"]].rename(
            columns={"logFC": "val_logFC"}
        )
        merged = disc.join(val, how="inner").dropna()

        n = len(merged)
        if n < 3:
            pearson_r = pearson_p = spearman_r = spearman_p = np.nan
        else:
            x = merged[disc_lfc_col].values.astype(float)
            y = merged["val_logFC"].values.astype(float)
            pearson_r, pearson_p = stats.pearsonr(x, y)
            spearman_r, spearman_p = stats.spearmanr(x, y)

        rows.append({
            "organ": organ,
            "n_genes_overlap": n,
            "pearson_r": round(pearson_r, 4) if not pd.isna(pearson_r) else np.nan,
            "pearson_p": pearson_p,
            "spearman_r": round(spearman_r, 4) if not pd.isna(spearman_r) else np.nan,
            "spearman_p": spearman_p,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Main
# =============================================================================

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    paths = setup_paths(base_dir)

    print("=" * 72)
    print("CROSS-ORGAN FIBROSIS BIOMARKER VALIDATION PIPELINE")
    print("                  (Steps 1-2 + Correlation)")
    print("=" * 72)

    # ------------------------------------------------------------------
    print("\n[STEP 0] Loading input data ...")
    # ------------------------------------------------------------------
    core_df = load_core_genes(paths["core_genes"])
    val_degs = load_validation_degs(paths["validation"])

    # ------------------------------------------------------------------
    print("\n[STEP 1] Direction concordance check ...")
    # ------------------------------------------------------------------
    annotated_df = annotate_validation_status(core_df, val_degs)
    print("  [OK] Added per-organ validation status columns")

    for organ in ORGANS:
        vc = annotated_df[f"{organ}_validated"].value_counts(dropna=False)
        found = int(annotated_df[f"{organ}_val_logFC"].notna().sum())
        concord = int(annotated_df[f"{organ}_validated"].fillna(False).sum())
        sig = int(annotated_df[f"{organ}_validated_sig"].fillna(False).sum())
        print(f"    {organ:6s}: found={found:>4},  concordant={concord:>4},  "
              f"concordant+sig(adj_p<0.05)={sig:>4}")

    # ------------------------------------------------------------------
    print("\n[STEP 2] Strict validation funnel ...")
    # ------------------------------------------------------------------
    funnel_df = build_funnel_table(annotated_df)
    per_organ_df = per_organ_validation_rates(annotated_df)

    print("\n  -- Validation Funnel Summary --")
    for criteria_name in funnel_df["criteria"].unique():
        print(f"\n  {criteria_name}:")
        sub = funnel_df[funnel_df["criteria"] == criteria_name]
        for _, row in sub.iterrows():
            print(f"    {row['tier']:<35s}  "
                  f"{row['gene_count']:>5,} genes  ({row['percentage']:>5.2f}%)")

    print("\n  -- Per-Organ Validation Rates --")
    print(per_organ_df.to_string(index=False))

    # ------------------------------------------------------------------
    print("\n[STEP 4] Discovery vs Validation correlation ...")
    # ------------------------------------------------------------------
    corr_df = compute_correlation_summary(core_df, val_degs)
    print("\n  -- Pearson & Spearman Correlation (discovery logFC vs validation logFC) --")
    print(corr_df.to_string(index=False, na_rep="N/A"))

    # ------------------------------------------------------------------
    print("\n[STEP 5] Saving output files ...")
    # ------------------------------------------------------------------
    annotated_df.to_csv(paths["out_validated"], index=False)
    print(f"  [SAVED] {paths['out_validated']}  ({len(annotated_df):,} rows)")

    all4 = annotated_df[annotated_df["n_organs_validated"] == 4].copy()
    all4_sorted = all4.sort_values("n_organs_validated_sig", ascending=False)
    all4_sorted.to_csv(paths["out_final"], index=False)
    print(f"  [SAVED] {paths['out_final']}  ({len(all4_sorted):,} genes)")

    # Also save funnel + per-organ tables as CSVs for later reference
    funnel_csv = os.path.join(paths["results"], "validation_funnel_summary.csv")
    per_organ_csv = os.path.join(paths["results"], "validation_per_organ_rates.csv")
    corr_csv = os.path.join(paths["results"], "validation_correlation_summary.csv")
    funnel_df.to_csv(funnel_csv, index=False)
    per_organ_df.to_csv(per_organ_csv, index=False)
    corr_df.to_csv(corr_csv, index=False)
    print(f"  [SAVED] {funnel_csv}")
    print(f"  [SAVED] {per_organ_csv}")
    print(f"  [SAVED] {corr_csv}")

    print("\n" + "=" * 72)
    print("VALIDATION PIPELINE (STEPS 1-2) COMPLETE")
    print("=" * 72)
    print(f"\n  Final validated genes (all 4 organs, direction-concordant): "
          f"{len(all4_sorted):,}")
    if len(all4_sorted) > 0:
        print("  Gene list:")
        for i, g in enumerate(all4_sorted["gene"].tolist(), 1):
            sig_n = all4_sorted.iloc[i-1]["n_organs_validated_sig"]
            print(f"    {i:>3}. {g:<15s}  [sig in {sig_n}/4 organs]")


if __name__ == "__main__":
    main()
