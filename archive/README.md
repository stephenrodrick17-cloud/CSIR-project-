# Archive: Superseded & Pre-Correction Artifacts

These files represent earlier exploratory analyses, pre-correction artifacts, intermediate outputs, and duplicate files preserved for complete historical transparency. None of these files are active or referenced in the final study results.

## Archived Superseded Scripts
- `fix_real_severity_pipeline.py`: Exploratory severity calculation containing preliminary interpolation; superseded by the Nakagawa et al. audit (`kidney_severity_status.csv`) and validated disease-only Spearman correlations in `results/within_vs_pooled_correlation_check.csv`.
- `disease_only_spearman.py`: Early correlation script configured for the deleted `organ_validation2_data/` directory; superseded by `results/within_vs_pooled_correlation_check.csv` and `run_master_consolidation.py`.
- `run_master_validation2_pipeline.py`: Early pipeline script targeting deleted 98-gene ECM files and placeholder paths; superseded by `run_master_validation_pipeline.py`.
- `corrected_full_pipeline.py`: Monolithic early pipeline script; superseded by modular scripts (`discovery_config.py`, `preprocess_build_deg_csvs.py`, `run_master_validation_pipeline.py`, and `run_master_consolidation.py`).
- `pan_fibrotic_analysis_PRE_FIX.py`: Pre-fix analysis script before formal 3-tier cohort isolation.

## Archived Duplicate & Superseded Data Tables
- `ecm_clean_genes.csv`: Root duplicate; canonical active file is `results/ecm_clean_genes.csv`.
- `pan_fibrotic_core_genes_corrected.csv`: Root duplicate; canonical active file is `results/pan_fibrotic_core_genes_corrected.csv`.
- `within_vs_pooled_correlation_check.csv`: Root duplicate; canonical active file is `results/within_vs_pooled_correlation_check.csv`.
- `ml_training_matrix_clean.csv`: Early unnormalized ML matrix draft; superseded by ComBat batch-corrected matrix `results/real_human_patient_combat_corrected_matrix.csv`.
- `pooled_vs_original_severity.csv`: Early pooled (case+control) correlation table; superseded by strict disease-only analysis in `results/within_vs_pooled_correlation_check.csv`.
- `venn_4organ_region_counts.csv`: Intermediate 16-region Venn count table; superseded by canonical Venn outputs.
- `kidney_dataset_search_audit.csv`: Negative-screening audit table documenting public kidney cohorts lacking per-sample continuous fibrosis metadata.
- `kidney_severity_status.csv`: Detailed audit of Nakagawa et al. 2015 histological staging showing underpowered sample counts ($n=5$).
- `pan_fibrotic_core_genes_PRE_FIX_175genes.csv`: Pre-fix 175-gene candidate list before rigorous data-leak elimination.
- `venn_4organ_region_counts_PRE_FIX.csv`: Pre-fix Venn count table.

## Archived Duplicate Visualizations
- `clean_ecm_validation_survival_barchart.png`: Root duplicate; canonical figure is `plots/clean_ecm_validation_survival_barchart.png`.
- `study_design_funnel_corrected.png`: Root duplicate; canonical figure is `plots/study_design_funnel_corrected.png`.
- `upset_plot_4organs_corrected.png`: Root duplicate; canonical figure is `plots/upset_plot_4organs_corrected.png`.
- `venn_4organ_manual_ellipses.png`: Root duplicate; canonical figure is `plots/venn_4organ_manual_ellipses.png`.
- `upset_plot_4organs_PRE_FIX.png`: Pre-fix UpSet diagram based on early 175-gene set.
