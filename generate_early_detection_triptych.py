"""
Master runner for Module 2: Early Detection & Progression Trajectory Analysis
Executes both 5-gene unanimous consensus and 9-gene multi-model consensus suites on GSE84044 Scheuer histological stages.
Outputs:
  - plots/hub_genes_early_detection_triptych_5genes.png
  - results/hub_genes_early_stage_roc_metrics_5genes.csv
  - plots/hub_genes_early_detection_triptych_9genes.png
  - results/hub_genes_early_stage_roc_metrics_9genes.csv
"""
import subprocess
import sys

print("Executing 5-Gene Early Detection Triptych Analysis...")
subprocess.run([sys.executable, "generate_early_detection_triptych_5genes.py"], check=True)

print("\nExecuting 9-Gene Early Detection Triptych Analysis...")
subprocess.run([sys.executable, "generate_early_detection_triptych_9genes.py"], check=True)

print("\nModule 2 early detection triptych execution complete.")
