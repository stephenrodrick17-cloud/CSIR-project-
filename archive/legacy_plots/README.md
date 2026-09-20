# Legacy & Pre-Fix Plots Archive

This directory contains older, superseded figures and diagrams from earlier exploratory iterations of the CSIR Pan-Fibrotic Study. These have been retired to prevent any visual contradiction with the final, locked numbers in the main repository.

---

### Archived Figures & Reasons for Retirement

| Legacy Figure | Superseded By | Reason for Retirement |
| :--- | :--- | :--- |
| `study_progress_summary_chart.png` | `plots/study_design_funnel_corrected.png` | Contained outdated metrics (394 samples, 241 DEGs, 50 ECM genes, and legacy candidate `TNXB`). Final verified numbers: 14 Discovery cohorts, 1,069 samples, 86 core DEGs, 24 clean core ECM genes. |
| `ml_8hub_diagnostic_roc_curves.png` | `plots/ml_4model_consensus_hub_biomarkers.png` | Based on early 8-gene prototype evaluated on the uncorrected 799-sample matrix (which contained leak GSE58095 and unverified genes `TNXB`, `SPP1`, `GDF15`). The canonical figure evaluates all 24 clean core ECM genes across a 5-seed 4-model ensemble on 1,069 pure discovery samples. |
| `ml_8hub_4organ_sample_heatmap.png` | `plots/ml_4model_consensus_hub_biomarkers.png` | Prototype heatmap for the retired 8-gene set. |
| `ml_8hub_biomarkers_expression_heatmap.png` | `plots/ml_4model_consensus_hub_biomarkers.png` | Prototype expression heatmap for the retired 8-gene set. |
| `ml_8hub_multiorgan_boxplots.png` | `plots/val2_mannwhitney_spearman_combined.png` | Multi-organ boxplots from the retired 8-gene prototype. Superseded by rigorous held-out Val2 Mann-Whitney U boxplots and Spearman correlation regressions. |
| `study_design_funnel.png` (from `validation_2/plots/`) | `plots/study_design_funnel_corrected.png` | Displayed early 7-gene panel and contained GSE58095 in Skin Val1 (which was later purged for leak prevention). |
| `validation2_confirmed_venn.png` (from `validation_2/`) | `plots/clean_ecm_validation_survival_barchart.png` | Displayed an early 7-gene intersection from preliminary testing before the 24-gene Clean Core ECM panel and 18-gene Val2 replication set were finalized. |

---

### Active Canonical Figures (in `plots/`)

All canonical, verified figures are located in `plots/`:
1. `plots/study_design_funnel_corrected.png`: 4 organs, 14 discovery cohorts, 1,069 samples, 86 core DEGs, 24 clean core ECM, 18 Val2 replicated, 4 Tier 1.
2. `plots/venn_4organ_manual_ellipses.png`: 4-organ overlap showing 86 core shared genes across Kidney (11,758), Liver (2,268), Lung (7,803), and Skin (2,979).
3. `plots/upset_plot_4organs_corrected.png`: UpSet cardinality visualization showing all 15 subset intersections with 100% mathematical consistency to the 4-organ Venn.
4. `plots/venn_organ_vs_ecm_compendium.png`: 2x2 grid of organ DEGs vs. Human Matrisome (1,027 genes), 3-way ECM intersection, and 4-way pan-fibrotic core (86 DEGs, 24 ECM).
5. `plots/clean_ecm_validation_survival_barchart.png`: 24 Clean Core ECM genes with exact Val1 & Val2 replicated organ counts, categorized into Tier 1 (4), Tier 2 (14), and Not Supported (6).
6. `plots/ml_4model_consensus_hub_biomarkers.png`: 24 Clean Core ECM genes ranked by 5-seed mean AUC (0.520–0.843) across LASSO, Random Forest, XGBoost, and SVM-RFE on 1,069 pure discovery samples.
7. `plots/val2_platform_stratified_comparison.png`: Platform-stratified (Microarray vs. RNA-seq) concordance heatmap across all 24 Clean Core ECM genes.
8. `plots/val2_mannwhitney_spearman_combined.png`: Paired boxplots and severity regressions for the 4 dual-significant Validation 2 hub genes (`AEBP1`, `COL1A1`, `COL1A2`, `VWF`).
