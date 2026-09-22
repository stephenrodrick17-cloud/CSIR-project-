"""
Master runner for Module 4: Immune & Stromal Microenvironment Infiltration Analysis
Executes both 5-gene unanimous consensus and 9-gene multi-model consensus suites on GSE84044 (N=124 human liver biopsies).
Outputs:
  - plots/hub_genes_immune_infiltration_5genes.png
  - results/hub_genes_immune_correlations_5genes.csv
  - plots/hub_genes_immune_infiltration_9genes.png
  - results/hub_genes_immune_correlations_9genes.csv
"""
import subprocess
import sys

print("Executing 5-Gene Immune & Stromal Infiltration Analysis...")
subprocess.run([sys.executable, "run_immune_infiltration_analysis_5genes.py"], check=True)

print("\nExecuting 9-Gene Immune & Stromal Infiltration Analysis...")
subprocess.run([sys.executable, "run_immune_infiltration_analysis_9genes.py"], check=True)

print("\nModule 4 immune infiltration analysis execution complete.")
