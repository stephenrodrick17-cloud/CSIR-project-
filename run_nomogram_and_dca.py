"""
Master runner for Module 3: Clinical Diagnostic Nomogram & Decision Curve Analysis (DCA)
Executes both 5-gene unanimous consensus and 9-gene multi-model consensus suites on N=1,069 discovery biopsies.
Outputs:
  - plots/hub_genes_nomogram_and_dca_5genes.png
  - results/hub_genes_nomogram_parameters_5genes.csv
  - plots/hub_genes_nomogram_and_dca_9genes.png
  - results/hub_genes_nomogram_parameters_9genes.csv
"""
import subprocess
import sys

print("Executing 5-Gene Diagnostic Nomogram & DCA...")
subprocess.run([sys.executable, "run_nomogram_and_dca_5genes.py"], check=True)

print("\nExecuting 9-Gene Diagnostic Nomogram & DCA...")
subprocess.run([sys.executable, "run_nomogram_and_dca_9genes.py"], check=True)

print("\nModule 3 diagnostic nomogram and DCA execution complete.")
