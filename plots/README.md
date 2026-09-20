# Canonical Publication Figures Directory (`plots/`)

This directory contains the final, publication-grade figures for the **CSIR Pan-Fibrotic Multi-Organ Biomarker Study**. Every figure here has been verified against the locked dataset and source files.

---

### Core Numbers Enforced Across All Figures
- **Discovery Cohorts**: $N = 14$ independent cohorts across 4 human organs (Kidney: 1, Liver: 3, Lung: 5, Skin: 5).
- **Discovery Sample Size**: $N = 1,069$ human patient tissue samples (170 non-fibrotic controls, 899 fibrotic cases).
- **Core Significance Filter**: eBayes moderated $t$-test $|\log_2\text{FC}| \ge 0.585$ (1.5-fold), Benjamini-Hochberg FDR $p < 0.05$.
- **Pan-Fibrotic Conserved Core**: **86 DEGs** shared across all 4 organs.
- **Clean Core ECM Program**: **24 genes** overlapping the Human Matrisome Masterlist ($N=1,027$).
- **Validation 1 (Internal)**: 24/24 direction-concordant across all 4 organs; 12/24 statistically significant in $\ge 2/4$ organs.
- **Validation 2 (Held-Out Multi-Center Cohorts)**: 18/24 replicated in $\ge 2/4$ organs; 6/24 replicated in $1/4$ organ.
- **Ensemble ML**: 4-model (LASSO, Random Forest, XGBoost, SVM-RFE) $\times$ 5-seed consensus on 1,069 pure discovery samples.
- **Evidence Tiers**:
  - **Tier 1 (Full Spectrum, $N=4$)**: `COL15A1`, `SERPINE2`, `COL1A1`, `SERPINF2` (Val2 $\ge 2/4$ + ML stability $\ge 4/5$).
  - **Tier 2 (Biological Validation-Only, $N=14$)**: `AEBP1`, `COL3A1`, `COL1A2`, `VWF`, `CCL2`, `SERPINH1`, `CCL19`, `LTBP2`, `CCL5`, `SVEP1`, `LAMC3`, `CLEC2D`, `PDGFD`, `MFAP4`.
  - **Not Supported ($N=6$)**: `MDK`, `FGF14`, `COLEC11`, `SPARCL1`, `BMP1`, `CCL21` (Val2 replicated in only 1/4 organ).

---

### Figure Index

| File | Description | Source Script / Input Data |
| :--- | :--- | :--- |
| `study_design_funnel_corrected.png` | Comprehensive 5-stage discovery and validation funnel diagram tracking cohort counts, sample sizes, and gene attrition at each layer. | `generate_final_figures.py` |
| `venn_4organ_manual_ellipses.png` | 4-Organ Venn diagram showing DEGs for Kidney (11,758), Liver (2,268), Lung (7,803), Skin (2,979), and the conserved 86-gene core intersection. | `venn_4organs.py` |
| `upset_plot_4organs_corrected.png` | UpSet plot showing all 15 subset intersections across the 4 organs; 100% mathematically and cardinality identical to the 4-organ Venn. | `generate_upset_plot.py` |
| `venn_organ_vs_ecm_compendium.png` | Multi-panel Venn compendium: 2x2 grid of organ DEGs vs. ECM Masterlist (1,027 genes), 3-way ECM intersection, and 4-way pan-fibrotic core (86 DEGs, 24 ECM). | `venn_organ_vs_ecm.py` |
| `clean_ecm_validation_survival_barchart.png` | Survival analysis of all 24 Clean Core ECM genes across Validation 1 (Internal) and Validation 2 (Held-out), grouped into Tier 1, Tier 2, and Not Supported. | `generate_clean_survival_barchart.py` |
| `ml_4model_consensus_hub_biomarkers.png` | Machine learning ranking of all 24 Clean Core ECM genes by mean diagnostic ROC AUC (5-seed average) across LASSO, Random Forest, XGBoost, and SVM-RFE. | `generate_final_figures.py` |
| `val2_platform_stratified_comparison.png` | Platform-stratified (Microarray: Kidney GSE30529, Liver GSE14323 vs. RNA-seq: Lung GSE83717, Skin GSE125362) logFC concordance heatmap across all 24 genes. | `results/val2_platform_stratified_table.csv` |
| `val2_mannwhitney_spearman_combined.png` | Validation 2 dual-significant hub biomarkers (`AEBP1`, `COL1A1`, `COL1A2`, `VWF`): Mann-Whitney U disease vs. control boxplots + Spearman severity correlation regressions. | `validation_2/` held-out datasets |
| `severity_validation4_layer_boxplots.png` | Clinical severity stage boxplots across independent cohorts for `COL15A1`, `COL1A1`, and `SERPINF2`. | Clinical staging data |
