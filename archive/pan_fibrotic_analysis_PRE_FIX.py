#!/usr/bin/env python3
"""
Pan-Fibrotic Core Gene Analysis Script

This script performs cross-organ differential gene expression analysis
across kidney, liver, lung, and skin fibrosis datasets, identifying
shared (pan-fibrotic) core genes and annotating them with ECM reference data.
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents

# =============================================================================
# 1. Load all four CSV files
# =============================================================================
print("=" * 70)
print("STEP 1: Loading DEG CSV files from /results directory")
print("=" * 70)

# Define file paths for each organ's DEG results
base_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(base_dir, "results")

csv_files = {
    "kidney": os.path.join(results_dir, "kidney_DEGs.csv"),
    "liver": os.path.join(results_dir, "liver_DEGs.csv"),
    "lung": os.path.join(results_dir, "lung_DEGs.csv"),
    "skin": os.path.join(results_dir, "skin_DEGs.csv"),
}

# Dictionary to store raw DEG DataFrames per organ
deg_data_raw = {}
for organ, filepath in csv_files.items():
    if os.path.exists(filepath):
        deg_data_raw[organ] = pd.read_csv(filepath)
        print(f"  [OK] Loaded {organ}: {len(deg_data_raw[organ])} genes total")
    else:
        print(f"  [MISSING] File not found: {filepath}")
        raise FileNotFoundError(f"Required file not found: {filepath}")

# =============================================================================
# 2. Filter each to significant genes only
#    (adj_p_value < 0.05 AND abs(logFC) > 0.585)
# =============================================================================
print("\n" + "=" * 70)
print("STEP 2: Filtering for significant genes")
print("  Criteria: adj_p_value < 0.05 AND |logFC| > 0.585")
print("=" * 70)

# Thresholds for significance
ADJ_P_THRESHOLD = 0.05
LOGFC_THRESHOLD = 0.585

# Dictionary to store filtered significant gene sets per organ
significant_genes = {}
for organ, df in deg_data_raw.items():
    # Apply significance filters
    filtered = df[
        (df["adj_p_value"] < ADJ_P_THRESHOLD) &
        (df["logFC"].abs() > LOGFC_THRESHOLD)
    ].copy()

    # Ensure gene names are uppercase for consistent matching
    filtered["gene"] = filtered["gene"].astype(str).str.upper()
    significant_genes[organ] = filtered

    n_total = len(df)
    n_sig = len(filtered)
    pct = (n_sig / n_total * 100) if n_total > 0 else 0
    print(f"  {organ}: {n_sig} / {n_total} significant genes ({pct:.1f}%)")

# Convert to sets for set operations (gene symbols only)
gene_sets = {
    organ: set(df["gene"].values) for organ, df in significant_genes.items()
}

# =============================================================================
# 3. Create an UpSet plot for all 4 sets
#    (matplotlib_venn doesn't support 4 circles, so UpSet is the standard approach)
# =============================================================================
print("\n" + "=" * 70)
print("STEP 3: Generating UpSet plot for 4-set overlap visualization")
print("=" * 70)

# Create the UpSet input structure from the gene sets
# from_contents expects a dict of {set_name: list/set_of_members}
upset_data = from_contents(gene_sets)

# Create figure and UpSet plot
fig = plt.figure(figsize=(12, 8))
upset = UpSet(
    upset_data,
    subset_size="count",
    show_counts=True,
    sort_by="cardinality",
    show_percentages=True,
    element_size=40,
)

# Add the intersection matrix and bar plots to the figure
upset.plot(fig=fig)

# Add title
fig.suptitle(
    "Overlap of Significant DEGs Across Kidney, Liver, Lung, and Skin",
    fontsize=14,
    fontweight="bold",
    y=1.02,
)

# Save the plot to a PNG file
upset_plot_path = os.path.join(base_dir, "upset_plot_4organs.png")
plt.savefig(upset_plot_path, dpi=300, bbox_inches="tight")
print(f"  [SAVED] UpSet plot saved to: {upset_plot_path}")
plt.close(fig)

# =============================================================================
# 4. Find genes shared across ALL 4 organs (intersection)
#    and save to pan_fibrotic_core_genes.csv
# =============================================================================
print("\n" + "=" * 70)
print("STEP 4: Identifying pan-fibrotic core genes (shared across all 4 organs)")
print("=" * 70)

# Compute the intersection of all 4 gene sets
organs_list = list(gene_sets.keys())
core_genes_set = gene_sets[organs_list[0]].copy()
for organ in organs_list[1:]:
    core_genes_set &= gene_sets[organ]

core_genes_list = sorted(core_genes_set)
print(f"  Found {len(core_genes_list)} pan-fibrotic core genes shared across all 4 organs")

# Build a DataFrame with additional info (logFC and p-values per organ)
core_genes_df = pd.DataFrame({"gene": core_genes_list})

# Merge in logFC and adj_p_value from each organ for reference
for organ in organs_list:
    organ_df = significant_genes[organ][["gene", "logFC", "adj_p_value"]].copy()
    organ_df = organ_df.rename(columns={
        "logFC": f"{organ}_logFC",
        "adj_p_value": f"{organ}_adj_p_value"
    })
    core_genes_df = core_genes_df.merge(organ_df, on="gene", how="left")

# Print the list to console
print("\n  Pan-fibrotic core gene list:")
for i, gene in enumerate(core_genes_list, 1):
    print(f"    {i:3d}. {gene}")

# Save to CSV
core_genes_csv_path = os.path.join(base_dir, "pan_fibrotic_core_genes.csv")
core_genes_df.to_csv(core_genes_csv_path, index=False)
print(f"\n  [SAVED] Core genes saved to: {core_genes_csv_path}")

# =============================================================================
# 5. Load ECM reference list and annotate core genes
#    Check which core genes are in the ECM masterlist
# =============================================================================
print("\n" + "=" * 70)
print("STEP 5: Annotating core genes with ECM reference data")
print("=" * 70)

# Path to the ECM reference Excel file
reference_dir = os.path.join(base_dir, "reference")
ecm_xlsx_path = os.path.join(reference_dir, "ECM_genes_all.xlsx")
ecm_sheet_name = "Hs_ECM_Masterlist"

if os.path.exists(ecm_xlsx_path):
    # Load the ECM masterlist (header on row 2, so header=1 in 0-indexed)
    ecm_df = pd.read_excel(
        ecm_xlsx_path,
        sheet_name=ecm_sheet_name,
        header=1,
        engine="openpyxl",
    )
    print(f"  [OK] Loaded ECM reference: {len(ecm_df)} entries")

    # Standardize the gene symbol column to uppercase for matching
    ecm_gene_col = "Gene Symbol"
    if ecm_gene_col not in ecm_df.columns:
        raise ValueError(
            f"Column '{ecm_gene_col}' not found in ECM sheet. "
            f"Available columns: {list(ecm_df.columns)}"
        )
    ecm_df[ecm_gene_col] = ecm_df[ecm_gene_col].astype(str).str.upper()

    # Identify the matrisome category column (different versions may name it differently)
    matrisome_candidates = [
        "Matrisome Category",
        "matrisome_category",
        "Category",
        "Matrisome_category",
    ]
    matrisome_col = None
    for candidate in matrisome_candidates:
        if candidate in ecm_df.columns:
            matrisome_col = candidate
            break

    if matrisome_col is None:
        print(f"  [WARNING] No matrisome category column found. "
              f"Available columns: {list(ecm_df.columns)}")
        matrisome_lookup = {}
    else:
        print(f"  Using matrisome category column: '{matrisome_col}'")
        # Create lookup dict: gene -> matrisome category
        matrisome_lookup = dict(zip(
            ecm_df[ecm_gene_col],
            ecm_df[matrisome_col].fillna("Unknown")
        ))

    # Build the set of ECM genes for fast membership testing
    ecm_gene_set = set(ecm_df[ecm_gene_col].values)

    # Add is_ECM_gene and matrisome_category columns to the core genes DataFrame
    core_genes_df["is_ECM_gene"] = core_genes_df["gene"].isin(ecm_gene_set)
    core_genes_df["matrisome_category"] = core_genes_df["gene"].map(matrisome_lookup)

    # Fill NaN matrisome_category with empty string for non-ECM genes
    core_genes_df["matrisome_category"] = core_genes_df["matrisome_category"].fillna("")

    print("  [OK] ECM annotation complete")
else:
    print(f"  [MISSING] ECM reference file not found at: {ecm_xlsx_path}")
    print("  Adding empty annotation columns so the script completes.")
    core_genes_df["is_ECM_gene"] = False
    core_genes_df["matrisome_category"] = ""

# Re-save the annotated CSV (overwriting the previous version)
core_genes_df.to_csv(core_genes_csv_path, index=False)
print(f"  [SAVED] Annotated core genes saved to: {core_genes_csv_path}")

# =============================================================================
# 6. Print summary statistics
# =============================================================================
print("\n" + "=" * 70)
print("STEP 6: Summary Statistics")
print("=" * 70)

total_core = len(core_genes_df)
n_ecm = core_genes_df["is_ECM_gene"].sum()
pct_ecm = (n_ecm / total_core * 100) if total_core > 0 else 0.0

print(f"  Total pan-fibrotic core genes (shared across all 4 organs): {total_core}")
print(f"  Number of core genes that are ECM genes:                   {n_ecm}")
print(f"  Percentage of core genes that are ECM genes:                {pct_ecm:.1f}%")

# Breakdown by matrisome category if available
if len(core_genes_df[core_genes_df["is_ECM_gene"]]) > 0 and \
        "matrisome_category" in core_genes_df.columns:
    category_counts = (
        core_genes_df[core_genes_df["is_ECM_gene"]]["matrisome_category"]
        .value_counts()
    )
    if len(category_counts) > 0:
        print("\n  ECM gene breakdown by matrisome category:")
        for cat, count in category_counts.items():
            print(f"    {cat}: {count}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
