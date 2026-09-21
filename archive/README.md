# Archive: Superseded & Pre-Correction Artifacts

These files represent earlier exploratory analyses, pre-correction artifacts, intermediate outputs, duplicate files, and superseded gene count tables preserved for complete historical transparency. None of these files are active or referenced in the current live study results (which strictly use the locked **86 Pan-Fibrotic Core DEGs** and **24 Clean Core ECM Genes**).

## Archived Superseded Directories
- `archive/validation_2/`: Entire directory containing early exploratory Validation 2 tables and plots generated from superseded 98-gene / 50-gene ECM panels (e.g., `validation2_real_master_summary.csv` [392 tests = 4 organs $\times$ 98 genes], `validation2_master_summary.csv`, `all_organs_validation2_summary.csv` [51 genes], per-organ 51-gene results, and 7-gene `final_confirmed_panel.csv`). Superseded by canonical 24-gene held-out validation results in `results/validation2_ecm_core_results.csv` and `results/final_master_evidence_table.csv`.
- `archive/validation_1_all_organs/`: Contains early 15-gene organ-validated signatures (`kidney_validated_signature.csv`, `liver_validated_signature.csv`, `lung_validated_signature.csv`, `skin_validated_signature.csv`). Superseded by canonical 24-gene Validation 1 results in `results/validation1_layer2_all_results.csv`.
- `archive/legacy_plots/`: Historical plot archive containing early drafts with superseded metrics.
- `archive/legacy_pre_fix/`: Intermediate pre-fix DEG and Venn scripts.

## Archived Historical Reference Scripts & Supporting Data
- `fix_real_severity_pipeline.py`: Exploratory severity calculation; superseded by strict disease-only analysis in `results/within_vs_pooled_correlation_check.csv`.
- `corrected_full_pipeline.py`: Monolithic early pipeline script; superseded by modular audited pipeline.
- `kidney_dataset_search_audit.csv`: Negative-screening audit table documenting public kidney cohorts lacking per-sample continuous fibrosis metadata.
- `kidney_severity_status.csv`: Detailed audit of Nakagawa et al. 2015 histological staging showing underpowered sample counts ($n=5$).
- `within_vs_pooled_correlation_check.csv`: Historical reference table confirming disease-only vs pooled correlation stability.
- `pooled_vs_original_severity.csv`: Reference comparison table.

## Notice on Purged Conflicting Files
All superseded, unverified, intermediate gene candidate lists (including `pan_fibrotic_core_genes_49_c37eeee.csv`, `pan_fibrotic_core_genes_PRE_FIX_175genes.csv`, `master_validation_all_layers_results_72genes_superseded.csv`, `final_validated_pan_fibrotic_genes_15genes_superseded.csv`, `validation2_full_core_results_241genes_superseded.csv`, `ml_4model_hub_biomarkers_50genes_superseded.csv`, and `ml_8hub_diagnostic_roc_auc_metrics_superseded.csv`) have been **permanently deleted** from the repository to eliminate any possible conflict with the canonical locked study results (**86 Pan-Fibrotic Core DEGs**, **24 Clean Core ECM Genes**, and **4 Tier 1 Hub Biomarkers**).
