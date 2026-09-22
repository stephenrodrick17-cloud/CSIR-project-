"""
Consensus Machine Learning & WGCNA Discovery Pipeline Runner
Executes the audited 5-method multi-seed consensus pipeline (LASSO, SVM-RFE, Random Forest, XGBoost, and TOM Co-expression Clustering / WGCNA)
across the locked N=1,069 discovery matrix (9 GEO cohorts with paired controls) using 5 random seeds.

Canonical implementation: run_master_consolidation.py
Outputs:
  - results/real_human_patient_ml_training_matrix.csv
  - results/real_human_patient_combat_corrected_matrix.csv
  - results/stage2_hub_gene_vote_table.csv
  - results/final_master_evidence_table.csv
"""
import subprocess
import sys

print("Executing Audited 5-Method Multi-Seed Consensus Machine Learning & WGCNA Pipeline...")
subprocess.run([sys.executable, "run_master_consolidation.py"], check=True)

print("\nMachine learning and WGCNA pipeline execution complete.")
