"""
Master runner for Module 1: Protein-Protein Interaction (PPI) Network Analysis (STRING v12)
Executes both 5-gene unanimous consensus and 9-gene multi-model consensus suites querying official STRING v12 REST API.
Outputs:
  - plots/hub_genes_ppi_network_5genes.png
  - results/hub_genes_ppi_network_nodes_5genes.csv
  - results/hub_genes_ppi_network_edges_5genes.csv
  - plots/hub_genes_ppi_network_9genes.png
  - results/hub_genes_ppi_network_nodes_9genes.csv
  - results/hub_genes_ppi_network_edges_9genes.csv
"""
import subprocess
import sys

print("Executing 5-Gene PPI Network Analysis (STRING v12)...")
subprocess.run([sys.executable, "run_ppi_network_analysis_5genes.py"], check=True)

print("\nExecuting 9-Gene PPI Network Analysis (STRING v12)...")
subprocess.run([sys.executable, "run_ppi_network_analysis_9genes.py"], check=True)

print("\nModule 1 PPI network analysis execution complete.")
