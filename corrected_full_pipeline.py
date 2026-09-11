#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Corrected full pipeline:
  1. Rebuild per-organ discovery DEGs EXCLUDING validation GEOs
     (with patched gene-symbol column detection for Gene.Symbol / GI variants)
  2. Re-run pan_fibrotic_analysis to get corrected core-genes intersection
  3. Re-run validation pipeline
  4. Re-generate ECM vs Non-ECM charts
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import glob

# ---------------------------------------------------------------------------
# Monkey-patch the column detection BEFORE importing the helper module,
# so that Gene.Symbol (capital S) and a broader set are recognised.
# Also patch safe_map_upper for robustness.
# ---------------------------------------------------------------------------
import preprocess_build_deg_csvs as _pp

# ---------------------------------------------------------------------------
# Patched process_dataset_vectorized: more flexible symbol + GI detection
# ---------------------------------------------------------------------------
_ORIG_PROCESS = _pp.process_dataset_vectorized


def safe_map_upper(series, mapping):
    mapped = series.map(mapping)
    if mapped.dtype == object:
        return mapped.str.upper()
    return mapped.astype("object").where(mapped.isna(),
                                          mapped.astype(str).str.upper())


def process_dataset_patched(tsv_path, bioDBnet_map):
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        print(f"    Failed to read {tsv_path}: {e}")
        return pd.DataFrame()

    df.columns = [c.strip() for c in df.columns]

    # Broader column name search (case-insensitive match + aliases)
    def find_col(candidates, df_cols):
        lower_map = {c.lower(): c for c in df_cols}
        for cand in candidates:
            if cand in df_cols:
                return cand
            if cand.lower() in lower_map:
                return lower_map[cand.lower()]
        return None

    adj_p_col = find_col(
        ["adj.P.Val", "padj", "adj_p_value", "adj_pval",
         "adj.P.Value", "FDR", "adj p-value", "adj.Pval", "p.adj"],
        df.columns)
    p_col = find_col(
        ["P.Value", "pvalue", "p_value", "PValue", "pval",
         "p-value", "P.Value", "raw p-value"],
        df.columns)
    lfc_col = find_col(
        ["logFC", "log2FoldChange", "LFC", "lfc", "log2fc",
         "Log Fold Change"],
        df.columns)
    id_col = find_col(
        ["ID", "GeneID", "ProbeID", "probe_id", "id",
         "Probe_Set_ID", "ProbeSetID"],
        df.columns)
    sym_col = find_col(
        ["Symbol", "Gene.symbol", "Gene_symbol", "gene_symbol",
         "Gene Symbol", "Gene.Symbol", "GeneSymbol", "Symbols",
         "geneName", "Gene Name", "Associated Gene Name"],
        df.columns)
    gi_col = find_col(["GI", "gi", "GeneIndex", "GENEINDEX"], df.columns)

    if None in (adj_p_col, p_col, lfc_col, id_col):
        print(f"    Warning: skipping {os.path.basename(tsv_path)}, "
              f"missing cols (adj_p={adj_p_col}, p={p_col}, "
              f"logFC={lfc_col}, ID={id_col})")
        print(f"    Available: {list(df.columns)}")
        return pd.DataFrame()

    work = pd.DataFrame({
        "raw_id": df[id_col].astype(str),
        "adj_p_value": pd.to_numeric(df[adj_p_col], errors="coerce"),
        "p_value": pd.to_numeric(df[p_col], errors="coerce"),
        "logFC": pd.to_numeric(df[lfc_col], errors="coerce"),
    })

    if sym_col is not None:
        work["sym_direct"] = (
            df[sym_col].astype(str).str.strip().str.upper()
              .replace(["", "NAN", "NA", "-", "NONE", "NULL"], np.nan)
        )
    else:
        work["sym_direct"] = np.nan

    # GI / ILMN probe mapping (fallbacks)
    if gi_col is not None:
        work["gi_str"] = (
            df[gi_col].astype(str).str.strip()
              .replace(["", "NAN", "NA", "-", "NONE", "NULL"], np.nan)
        )
        work["sym_gi"] = safe_map_upper(work["gi_str"], bioDBnet_map)
        try:
            work["gi_int"] = (
                pd.to_numeric(df[gi_col], errors="coerce")
                  .astype("Int64").astype(str).replace("<NA>", np.nan)
            )
            mask = work["sym_gi"].isna() & work["gi_int"].notna()
            work.loc[mask, "sym_gi"] = safe_map_upper(
                work.loc[mask, "gi_int"], bioDBnet_map)
        except Exception:
            pass
    else:
        work["sym_gi"] = np.nan

    # ID -> _at stripping (Affymetrix)
    work["id_clean"] = work["raw_id"].astype(str).str.strip()
    work["id_no_at"] = work["id_clean"].where(
        ~work["id_clean"].str.endswith("_at"), work["id_clean"].str[:-3])
    work["sym_at"] = safe_map_upper(work["id_no_at"], bioDBnet_map)
    try:
        work["id_int"] = (
            pd.to_numeric(work["id_no_at"], errors="coerce")
              .astype("Int64").astype(str).replace("<NA>", np.nan)
        )
        mask = work["sym_at"].isna() & work["id_int"].notna()
        work.loc[mask, "sym_at"] = safe_map_upper(
            work.loc[mask, "id_int"], bioDBnet_map)
    except Exception:
        pass

    # General ID mapping (string + int) — also covers ILMN_ prefix IDs via map
    work["sym_id_str"] = safe_map_upper(work["id_clean"], bioDBnet_map)
    try:
        work["id_generic_int"] = (
            pd.to_numeric(work["id_clean"], errors="coerce")
              .astype("Int64").astype(str).replace("<NA>", np.nan)
        )
        mask = work["sym_id_str"].isna() & work["id_generic_int"].notna()
        work.loc[mask, "sym_id_str"] = safe_map_upper(
            work.loc[mask, "id_generic_int"], bioDBnet_map)
    except Exception:
        pass

    # Extra: try GB_ACC (GenBank accession) mapping via bioDBnet (some files)
    gb_col = find_col(["GB_ACC", "GB_LIST", "GB_ACC_LIST", "gb_acc"], df.columns)
    work["sym_gb"] = np.nan
    if gb_col is not None:
        work["gb_clean"] = (
            df[gb_col].astype(str).str.strip()
              .replace(["", "NAN", "NA", "-", "NONE"], np.nan)
        )
        work["sym_gb"] = safe_map_upper(work["gb_clean"], bioDBnet_map)

    # Hierarchical fallback order (most direct -> most indirect)
    work["gene"] = work["sym_direct"]
    work["gene"] = work["gene"].fillna(work["sym_gb"])  # GenBank if available
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


# Swap in the patched version for the duration of this script
_pp.process_dataset_vectorized = process_dataset_patched


BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "results")
os.makedirs(RESULTS, exist_ok=True)

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


def build_clean_discovery():
    print("=" * 72)
    print("STEP 1: REBUILD DISCOVERY DEGs (excluding validation GEOs)")
    print("         (with patched Gene.Symbol + GI column detection)")
    print("=" * 72)

    print("\nLoading global bioDBnet map ...")
    global_bio_map = {}
    for organ_dirname in ORGAN_DIRS.values():
        organ_path = os.path.join(BASE, organ_dirname)
        if os.path.isdir(organ_path):
            sub_map = _pp.load_biodbnet_mappings(organ_path)
            global_bio_map.update(sub_map)
    print(f"  {len(global_bio_map)} entries loaded")

    new_degs = {}
    for organ_name in ["kidney", "liver", "lung", "skin"]:
        organ_dirname = ORGAN_DIRS[organ_name]
        organ_path = os.path.join(BASE, organ_dirname)
        excludes = VAL_EXCLUDES[organ_name]

        print(f"\n=== {organ_name.upper()} ===")
        print(f"  Excluding GSE accessions: {excludes}")

        pattern = os.path.join(organ_path, "**", "*.top.table.tsv")
        all_tsvs = sorted(glob.glob(pattern, recursive=True))

        kept = []
        for t in all_tsvs:
            if any(acc in os.path.basename(t) for acc in excludes):
                print(f"  [EXCLUDE] {os.path.relpath(t, BASE)}")
                continue
            kept.append(t)

        print(f"  Kept {len(kept)} / {len(all_tsvs)} TSV files")
        all_dfs = []
        for tsv_path in kept:
            ds = f"{os.path.basename(os.path.dirname(tsv_path))}/{os.path.basename(tsv_path)}"
            print(f"    {ds}")
            ds_df = process_dataset_patched(tsv_path, global_bio_map)
            if len(ds_df) > 0:
                print(f"      -> {len(ds_df):,} genes with valid symbols")
                all_dfs.append(ds_df)
            else:
                print(f"      -> !! NO VALID GENES (check mapping) !!")
                # Try the original processor as a fallback too
                fb = _ORIG_PROCESS(tsv_path, global_bio_map)
                if len(fb) > 0:
                    print(f"         (fallback orig processor gave {len(fb)} — using it)")
                    all_dfs.append(fb)

        if not all_dfs:
            print(f"  [FATAL] No genes for {organ_name}! Abort.")
            sys.exit(1)

        combined = pd.concat(all_dfs, ignore_index=True)
        combined = combined.sort_values("adj_p_value", ascending=True)
        combined = combined.drop_duplicates(subset="gene", keep="first")
        combined = combined.sort_values("adj_p_value", ascending=True).reset_index(drop=True)

        out = os.path.join(RESULTS, f"{organ_name}_DEGs.csv")
        combined.to_csv(out, index=False)
        new_degs[organ_name] = combined
        sig = combined[(combined["adj_p_value"] < 0.05) &
                       (combined["logFC"].abs() > 0.585)]
        print(f"  -> {len(combined):,} total genes; "
              f"{len(sig):,} significant (adj_p<0.05 & |logFC|>0.585)")
        print(f"  [SAVED] {out}")

    # Quick sanity on liver
    print("\nLiver sanity check (disc vs GSE162694 validation):")
    liv_disc = pd.read_csv(os.path.join(RESULTS, "liver_DEGs.csv"))
    liv_val  = pd.read_csv(os.path.join(BASE, "validation",
                                        "liver_validation_DEGs.csv"))
    m = liv_disc.merge(liv_val, on="gene", suffixes=("_d", "_v"))
    from scipy import stats as sstats
    r, _ = sstats.pearsonr(m["logFC_d"].astype(float),
                           m["logFC_v"].astype(float))
    same = np.isclose(m["logFC_d"].astype(float),
                      m["logFC_v"].astype(float), atol=1e-6).sum()
    print(f"  Overlap: {len(m)}  Identical logFC: {same}  Pearson r: {r:.4f}")
    if r < 0.90 and same / len(m) < 0.10:
        print("  -> Liver independence CONFIRMED. OK")
    return new_degs


def run_pan_fibrotic_discovery():
    """Replicate pan_fibrotic_analysis.py using fresh discovery DEGs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_contents

    print("\n" + "=" * 72)
    print("STEP 2: Re-compute pan-fibrotic core genes (4-organ intersection)")
    print("=" * 72)

    deg_data_raw = {}
    csv_files = {
        "kidney": os.path.join(RESULTS, "kidney_DEGs.csv"),
        "liver":  os.path.join(RESULTS, "liver_DEGs.csv"),
        "lung":   os.path.join(RESULTS, "lung_DEGs.csv"),
        "skin":   os.path.join(RESULTS, "skin_DEGs.csv"),
    }
    for organ, fp in csv_files.items():
        deg_data_raw[organ] = pd.read_csv(fp)
        print(f"  {organ:6s}: {len(deg_data_raw[organ]):,} genes")

    ADJ_P = 0.05
    LFC = 0.585
    significant_genes = {}
    gene_sets = {}
    for organ, df in deg_data_raw.items():
        filt = df[(df["adj_p_value"] < ADJ_P) & (df["logFC"].abs() > LFC)].copy()
        filt["gene"] = filt["gene"].astype(str).str.upper()
        significant_genes[organ] = filt
        gene_sets[organ] = set(filt["gene"].values)
        print(f"  {organ:6s}: {len(filt):,} significant "
              f"({len(filt)/len(df)*100:.1f}%)")

    # UpSet plot
    print("\n  Re-generating UpSet plot ...")
    upset_data = from_contents(gene_sets)
    fig = plt.figure(figsize=(12, 8))
    upset = UpSet(upset_data, subset_size="count", show_counts=True,
                  sort_by="cardinality", show_percentages=True,
                  element_size=40)
    upset.plot(fig=fig)
    fig.suptitle(
        "Overlap of Significant DEGs Across 4 Organs (Corrected: Validation Excluded)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    upset_out = os.path.join(BASE, "upset_plot_4organs_corrected.png")
    plt.savefig(upset_out, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"  [SAVED] {upset_out}")

    organs_list = list(gene_sets.keys())
    core_genes_set = gene_sets[organs_list[0]].copy()
    for o in organs_list[1:]:
        core_genes_set &= gene_sets[o]
    core_genes_list = sorted(core_genes_set)
    print(f"\n  -> {len(core_genes_list)} pan-fibrotic core genes "
          f"(shared across all 4 organs AFTER correction)")

    core_df = pd.DataFrame({"gene": core_genes_list})
    for o in organs_list:
        organ_df = significant_genes[o][["gene", "logFC", "adj_p_value"]].copy()
        organ_df = organ_df.rename(columns={
            "logFC": f"{o}_logFC", "adj_p_value": f"{o}_adj_p_value",
        })
        core_df = core_df.merge(organ_df, on="gene", how="left")

    # ECM annotation
    ecm_xlsx = os.path.join(BASE, "reference", "ECM_genes_all.xlsx")
    sheet = "Hs_ECM_Masterlist"
    if os.path.exists(ecm_xlsx):
        ecm_df = pd.read_excel(ecm_xlsx, sheet_name=sheet, header=1,
                               engine="openpyxl")
        ecm_gene_col = "Gene Symbol"
        ecm_df[ecm_gene_col] = ecm_df[ecm_gene_col].astype(str).str.upper()
        mat_candidates = ["Matrisome Category", "matrisome_category",
                          "Category", "Matrisome_category"]
        mat_col = next((c for c in mat_candidates if c in ecm_df.columns), None)
        mat_lookup = {}
        if mat_col:
            mat_lookup = dict(zip(ecm_df[ecm_gene_col],
                                  ecm_df[mat_col].fillna("Unknown")))
        ecm_set = set(ecm_df[ecm_gene_col].values)
        core_df["is_ECM_gene"] = core_df["gene"].isin(ecm_set)
        core_df["matrisome_category"] = core_df["gene"].map(mat_lookup).fillna("")
        n_ecm = core_df["is_ECM_gene"].sum()
        print(f"  ECM genes in new core: {n_ecm} / {len(core_df)} "
              f"({n_ecm/len(core_df)*100:.1f}%)")
    else:
        core_df["is_ECM_gene"] = False
        core_df["matrisome_category"] = ""

    core_csv = os.path.join(BASE, "pan_fibrotic_core_genes_corrected.csv")
    core_df.to_csv(core_csv, index=False)
    print(f"  [SAVED] {core_csv}  ({len(core_df)} genes)")
    return core_df


def run_validation_pipeline(core_df):
    """Run validation pipeline against corrected core genes."""
    print("\n" + "=" * 72)
    print("STEP 3: Direction-concordance validation (corrected core genes)")
    print("=" * 72)

    val_dir = os.path.join(BASE, "validation")
    val_degs = {}
    for o in ORGAN_DIRS:
        p = os.path.join(val_dir, f"{o}_validation_DEGs.csv")
        d = pd.read_csv(p)
        d["gene"] = d["gene"].astype(str).str.upper().str.strip()
        val_degs[o] = d.set_index("gene")
        print(f"  {o:6s} validation: {len(val_degs[o]):,} genes")

    core = core_df.copy()
    core["gene"] = core["gene"].astype(str).str.upper().str.strip()
    ORGANS = list(ORGAN_DIRS.keys())

    def dir_concord(d, v):
        if pd.isna(d) or pd.isna(v) or d == 0 or v == 0:
            if pd.isna(d) or pd.isna(v):
                return np.nan
            return False
        return (d > 0) == (v > 0)

    for o in ORGANS:
        genes = core["gene"].values
        vdf = val_degs[o]
        found = [g in vdf.index for g in genes]
        v_lfc = np.where(found, vdf.reindex(genes)["logFC"].values, np.nan)
        v_adj = np.where(found, vdf.reindex(genes)["adj_p_value"].values, np.nan)

        core[f"{o}_val_logFC"] = v_lfc
        core[f"{o}_val_adj_p_value"] = v_adj

        v_list, vs_list = [], []
        d_lfc_col = f"{o}_logFC"
        for i, g in enumerate(genes):
            d_lfc = core.iloc[i][d_lfc_col]
            vl, va = v_lfc[i], v_adj[i]
            dc = dir_concord(d_lfc, vl)
            v_list.append(dc)
            if pd.isna(dc) or not dc:
                vs_list.append(np.nan if pd.isna(dc) else False)
            else:
                vs_list.append((not pd.isna(va)) and (va < 0.05))
        core[f"{o}_validated"] = v_list
        core[f"{o}_validated_sig"] = vs_list

    val_cols = [f"{o}_validated" for o in ORGANS]
    sig_cols = [f"{o}_validated_sig" for o in ORGANS]
    core["n_organs_validated"] = core[val_cols].fillna(False).sum(axis=1).astype(int)
    core["n_organs_validated_sig"] = core[sig_cols].fillna(False).sum(axis=1).astype(int)

    out_full = os.path.join(RESULTS,
                            "pan_fibrotic_core_genes_validated_corrected.csv")
    out_final = os.path.join(RESULTS,
                             "final_validated_pan_fibrotic_genes_corrected.csv")
    core.to_csv(out_full, index=False)
    all4 = core[core["n_organs_validated"] == 4].copy()
    all4_sorted = all4.sort_values("n_organs_validated_sig", ascending=False)
    all4_sorted.to_csv(out_final, index=False)
    print(f"  [SAVED] {out_full}  ({len(core)} rows)")
    print(f"  [SAVED] {out_final}  ({len(all4_sorted)} all-4-organ validated)")

    # Funnel table
    print("\n  Validation funnel (direction only):")
    total = len(core)
    for k in [1, 2, 3, 4]:
        n = (core["n_organs_validated"] >= k).sum()
        print(f"    >= {k} organ(s): {n:>4} / {total:>4}  ({n/total*100:5.1f}%)")
    print(f"    ALL 4 organs : {(core['n_organs_validated']==4).sum():>4} / {total:>4}")
    return core


def main():
    build_clean_discovery()
    core_df = run_pan_fibrotic_discovery()
    validated = run_validation_pipeline(core_df)
    print("\n" + "=" * 72)
    print("CORRECTED PIPELINE COMPLETE — now re-run ecm_vs_non_ecm_validation_viz.py")
    print("with paths updated to *_corrected.csv files for the final charts.")
    print("=" * 72)


if __name__ == "__main__":
    main()
