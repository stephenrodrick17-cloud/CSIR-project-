# Archive: Superseded & Pre-Correction Artifacts

These files represent earlier exploratory analyses, pre-correction artifacts, intermediate outputs, duplicate files, and superseded gene count tables preserved for complete historical transparency. None of these files are active or referenced in the current live study results (which strictly use the locked **86 Pan-Fibrotic Core DEGs** and **24 Clean Core ECM Genes**).

## Archived Superseded Directories
- `archive/validation_2/`: Entire directory containing early exploratory Validation 2 tables and plots generated from superseded 98-gene / 50-gene ECM panels (e.g., `validation2_real_master_summary.csv` [392 tests = 4 organs $\times$ 98 genes], `validation2_master_summary.csv`, `all_organs_validation2_summary.csv` [51 genes], per-organ 51-gene results, and 7-gene `final_confirmed_panel.csv`). Superseded by canonical 24-gene held-out validation results in `results/validation2_ecm_core_results.csv` and `results/final_master_evidence_table.csv`.
- `archive/validation_1_all_organs/`: Contains early 15-gene organ-validated signatures (`kidney_validated_signature.csv`, `liver_validated_signature.csv`, `lung_validated_signature.csv`, `skin_validated_signature.csv`). Superseded by canonical 24-gene Validation 1 results in `results/validation1_layer2_all_results.csv`.
- `archive/legacy_plots/`: Historical plot archive containing early drafts with superseded metrics.
- `archive/legacy_pre_fix/`: Intermediate pre-fix DEG and Venn scripts.

## Archived Superseded Scripts
- `fix_real_severity_pipeline.py`: Exploratory severity calculation containing preliminary interpolation; superseded by the Nakagawa et al. audit (`kidney_severity_status.csv`) and validated disease-only Spearman correlations in `results/within_vs_pooled_correlation_check.csv`.
- `disease_only_spearman.py`: Early correlation script configured for the deleted `organ_validation2_data/` directory; superseded by `results/within_vs_pooled_correlation_check.csv` and `run_master_consolidation.py`.
- `run_master_validation2_pipeline.py`: Early pipeline script targeting deleted 98-gene ECM files and placeholder paths; superseded by `run_master_validation_pipeline.py`.
- `corrected_full_pipeline.py`: Monolithic early pipeline script; superseded by modular scripts (`discovery_config.py`, `preprocess_build_deg_csvs.py`, `run_master_validation_pipeline.py`, and `run_master_consolidation.py`).
- `pan_fibrotic_analysis_PRE_FIX.py`: Pre-fix analysis script before formal 3-tier cohort isolation.

## Archived Superseded Data Tables
- `forensic_repository_audit_report_241_50_superseded.csv`: Pre-audit snapshot report referencing early 241-gene core, 50-gene clean ECM, and 49/50 replication metrics.
- `validation2_full_core_results_241genes_superseded.csv`: Validation 2 results evaluated across the early 241-gene core. Superseded by canonical 24-gene `results/validation2_ecm_core_results.csv`.
- `ml_4model_hub_biomarkers_50genes_superseded.csv`: Machine learning hub results evaluated across the early 50-gene ECM set (including preliminary single-seed runs with leaked `GSE58095`). Superseded by 5-seed ComBat-corrected `results/ml_4model_24ecm_hub_biomarkers.csv`.
- `ml_8hub_diagnostic_roc_auc_metrics_superseded.csv`: Single-seed ROC AUC metrics for 8 early consensus hubs. Superseded by canonical 5-seed 24-gene AUC evaluations in `results/final_master_evidence_table.csv`.
- `final_validated_pan_fibrotic_genes_15genes_superseded.csv`: Early 15-gene signature table from `corrected_full_pipeline.py`.
- `pan_fibrotic_core_genes_validated_72genes_superseded.csv`: Intermediate 72-gene table from earlier pipeline iteration.
- `master_validation_all_layers_results_72genes_superseded.csv`: Intermediate 72-gene cross-layer validation table.
- `ecm_clean_genes.csv`: Root duplicate; canonical active file is `results/ecm_clean_genes.csv`.
- `pan_fibrotic_core_genes_corrected.csv`: Root duplicate; canonical active file is `results/pan_fibrotic_core_genes_corrected.csv`.
- `within_vs_pooled_correlation_check.csv`: Root duplicate; canonical active file is `results/within_vs_pooled_correlation_check.csv`.
- `ml_training_matrix_clean.csv`: Early unnormalized ML matrix draft; superseded by ComBat batch-corrected matrix `results/real_human_patient_combat_corrected_matrix.csv`.
- `pooled_vs_original_severity.csv`: Early pooled (case+control) correlation table; superseded by strict disease-only analysis in `results/within_vs_pooled_correlation_check.csv`.
- `venn_4organ_region_counts.csv`: Intermediate 16-region Venn count table; superseded by canonical Venn outputs.
- `kidney_dataset_search_audit.csv`: Negative-screening audit table documenting public kidney cohorts lacking per-sample continuous fibrosis metadata.
- `kidney_severity_status.csv`: Detailed audit of Nakagawa et al. 2015 histological staging showing underpowered sample counts ($n=5$).
- `pan_fibrotic_core_genes_PRE_FIX_175genes.csv`: Pre-fix 175-gene candidate list before rigorous data-leak elimination.
- `pan_fibrotic_core_genes_49_c37eeee_superseded.csv`: Intermediate 49-gene DEG output from commit `c37eeee` containing non-ECM bystanders. Superseded by canonical 86-core DEGs (`results/pan_fibrotic_core_genes_corrected.csv`) and 24 Clean Core ECM genes (`results/ecm_clean_genes.csv`).
- `venn_4organ_region_counts_PRE_FIX.csv`: Pre-fix Venn count table.

## Archived Visualizations
- `clean_ecm_validation_survival_barchart.png`: Root duplicate; canonical figure is `plots/clean_ecm_validation_survival_barchart.png`.
- `study_design_funnel_corrected.png`: Root duplicate; canonical figure is `plots/study_design_funnel_corrected.png`.
- `upset_plot_4organs_corrected.png`: Root duplicate; canonical figure is `plots/upset_plot_4organs_corrected.png`.
- `venn_4organ_manual_ellipses.png`: Root duplicate; canonical figure is `plots/venn_4organ_manual_ellipses.png`.
- `venn_organ_vs_ecm_compendium.png`: Root duplicate; canonical figure is `plots/venn_organ_vs_ecm_compendium.png`.
- `upset_plot_4organs_PRE_FIX.png`: Pre-fix UpSet diagram based on early 175-gene set.
