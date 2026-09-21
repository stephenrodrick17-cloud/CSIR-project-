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
| `val2_mannwhitney_spearman_combined.png` | Validation 2 dual-significant hub biomarkers (`AEBP1`, `COL1A1`, `COL1A2`, `VWF`): Mann-Whitney U disease vs. control boxplots + Spearman severity correlation regressions. | Held-out Validation 2 clinical cohorts (`GSE14323`, `GSE162694`) |
| `severity_validation4_layer_boxplots.png` | **Supplementary External Severity Exploration**: Comprehensive $5 \times 3$ grid (15 panels) of clinical severity stage boxplots and regressions across external cohorts (`GSE84044`, `GSE135251`, `GSE38958`, `GSE9285`) for 5 priority genes (`COL15A1`, `COL1A1`, `SERPINE2`, `SERPINF2`, `TNXB`). Completely separate from the locked 14+4+4 pipeline. | `generate_layer4_comprehensive_plots.py` |
| `val2_all_18_genes_mannwhitney_boxplots.png` | Comprehensive $6 \times 3$ grid of Mann-Whitney U Disease vs. Control boxplots + strip points for all 18 surviving ECM genes in held-out Liver Val2 (`GSE14323`). | `generate_val2_all_18_genes_plots.py` |
| `val2_severity_spearman_correlations.png` | Comprehensive $6 \times 3$ grid of Spearman clinical disease severity tracking regressions across all 18 Validation 2 surviving genes staged in held-out Liver biopsy cohorts (METAVIR F1–F4 in GSE162694; Scheuer S1–S4 in GSE84044). | `generate_val2_all_18_genes_plots.py` |
| `hub_genes_validation_mwu_spearman.png` | **Tier 1 Universal Hub Biomarkers Validation**: Comprehensive $4 \times 3$ grid of Mann-Whitney U test (Held-Out Val 2 GSE14323), Spearman clinical severity regression (Scheuer staging GSE84044), and diagnostic ROC curves for the 4 Tier 1 Hub genes (`COL15A1`, `COL1A1`, `SERPINE2`, `SERPINF2`). | `generate_hub_genes_validation_plot.py` |
| `hub_genes_nomogram_and_dca.png` | **Diagnostic Nomogram & Clinical Decision Analysis**: Multivariable logistic regression nomogram, 1,000-bootstrap calibration curve, and Decision Curve Analysis (DCA) for the 4 Hub genes across 1,069 human biopsies (AUC = 0.901). | `run_nomogram_and_dca.py` |
| `hub_genes_ppi_network.png` | **Protein-Protein Interaction (PPI) Network**: STRING v12 high-confidence interaction architecture connecting the 4 Hub genes with master upstream regulators (`TGFB1`, `SMAD3`, `CTGF`), proteases (`MMP1/2`, `TIMP1`, `PLG`), and matrix scaffolding (`FN1`, `ITGB1`). | `run_ppi_network_analysis.py` |
| `hub_genes_gsea_hallmark_pathways.png` | **Gene Set Enrichment Analysis (GSEA)**: Multi-panel running enrichment plots across MSigDB Hallmark and KEGG pathways (EMT, ECM-Receptor, TGF-beta, Focal Adhesion, Coagulation/Serpin axis) stratified by 4-Hub biomarker score. | `run_gsea_enrichment_analysis.py` |
| `hub_genes_immune_infiltration.png` | **Immune Infiltration & Microenvironment Crosstalk**: Spearman correlation heatmap and subset comparisons across 12 immune/stromal populations in 1,069 biopsies, demonstrating coupling to M2 macrophages, Tregs, and vascular rarefaction. | `run_immune_infiltration_analysis.py` |
| `hub_genes_candidate_drugs.png` | **Candidate Drug Repurposing Network**: Bipartite pharmacological target network linking the 4 Hub genes to approved clinical standards (Pirfenidone, Nintedanib), serine protease inhibitors (Camostat, Nafamostat, Gabexate), and collagen synthesis modulators (Tranilast, Halofuginone). | `run_drug_repurposing_analysis.py` |
| `hub_genes_hpa_ihc_summary.png` | **In Silico IHC Protein Validation**: Human Protein Atlas (HPA v23) pathology staining profiles verifying protein-level upregulation of COL15A1, COL1A1, and SERPINE2 and loss of SERPINF2 in human kidney, liver, lung, and skin fibrosis. | `run_hpa_ihc_validation.py` |
| `hub_genes_early_stage_validation.png` | **Dedicated Subclinical Early-Stage Validation**: 6-panel comprehensive figure strictly comparing Healthy Controls (S0/F0) vs. Early-Stage Fibrosis (S1/S2 or F1/F2) across microarray (`GSE84044`, N=96) and RNA-seq (`GSE135251`, N=148), featuring individual ROC curves, 5-fold cross-validated multi-gene ensemble classifiers, and subclinical screening Decision Curve Analysis (DCA). | `run_early_stage_fibrosis_validation.py` |
| `hub_genes_early_detection_triptych.png` | **Early Detection & Subclinical Risk Stratification (1x3 Triptych)**: Nature/IEEE-style publication 3-panel figure featuring (A) Early-Stage ROC Curves (S0 vs S1-S2), (B) Early-Onset Switch vs Linear Progression Dynamics across histological stages, and (C) Subclinical Decision Curve Analysis across risk thresholds ($p_t = 0.05$ to $0.50$). | `generate_early_detection_triptych.py` |
| `hub_genes_docking_and_orthogonal_validation.png` | **Molecular Docking & Multi-Method Orthogonal Validation**: Publication 4-panel figure proving authenticity via (A) In Silico Molecular Docking Binding Energies ($\Delta G$ in kcal/mol) against PDB structures (`4D7N`, `2R9Y`, `1BKV`), (B) Structural Binding Pocket and residue contact maps, (C) Single-Cell RNA-Seq (scRNA-Seq) cellular localization, and (D) 6-Tier Orthogonal Evidence Pyramid. | `generate_docking_and_orthogonal_validation_figure.py` |





