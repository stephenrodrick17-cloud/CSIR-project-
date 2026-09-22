"""
Master runner for Module 5: Human Protein Atlas (HPA v23) Baseline IHC & Patient RNA Characterization
Runs both 5-gene and 9-gene suites using genuine empirical HPA normal baseline data and patient biopsy RNA fold changes.
"""
import subprocess
import sys

print("Executing 5-Gene HPA IHC & RNA Characterization...")
subprocess.run([sys.executable, "run_hpa_ihc_validation_5genes.py"], check=True)

print("Executing 9-Gene HPA IHC & RNA Characterization...")
subprocess.run([sys.executable, "run_hpa_ihc_validation_9genes.py"], check=True)

print("Module 5 execution complete.")
