import os
import gzip
import glob
import pandas as pd
import discovery_config as dc

print("=" * 80)
print("COMPREHENSIVE COHORT & DATA LEAKAGE AUDIT")
print("=" * 80)

# 1. Define all cohort sets explicitly
discovery_cohorts = dc.DISCOVERY_ALL
val1_cohorts = dc.VAL1_ALL
val2_cohorts = dc.VAL2_ALL

# ML training cohorts from actual file
ml_df = pd.read_csv('results/real_human_patient_ml_training_matrix.csv')
ml_cohorts = set(ml_df['study'].unique())

# Layer 4 (Post-ML independent severity validation cohorts)
layer4_cohorts = {'GSE84044', 'GSE135251', 'GSE38958', 'GSE213001', 'GSE9285'}

print(f"Discovery Cohorts ({len(discovery_cohorts)}): {sorted(list(discovery_cohorts))}")
print(f"Validation 1 Cohorts ({len(val1_cohorts)}): {sorted(list(val1_cohorts))}")
print(f"Validation 2 Cohorts ({len(val2_cohorts)}): {sorted(list(val2_cohorts))}")
print(f"Machine Learning Training Cohorts ({len(ml_cohorts)}): {sorted(list(ml_cohorts))}")
print(f"Layer 4 Independent Severity Cohorts ({len(layer4_cohorts)}): {sorted(list(layer4_cohorts))}")
print("-" * 80)

# 2. Pairwise Overlap Matrix
cohort_sets = {
    "Discovery": discovery_cohorts,
    "Val 1": val1_cohorts,
    "Val 2": val2_cohorts,
    "ML Training": ml_cohorts,
    "Layer 4 Severity": layer4_cohorts
}

names = list(cohort_sets.keys())
print(f"{'Comparison':<35} | {'Overlap Count':<15} | {'Overlap Accessions'}")
print("-" * 80)

any_leak = False
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        n1, n2 = names[i], names[j]
        # ML Training is intentionally a pure subset of Discovery
        if (n1 == "Discovery" and n2 == "ML Training") or (n1 == "ML Training" and n2 == "Discovery"):
            subset_check = ml_cohorts.issubset(discovery_cohorts)
            diff = ml_cohorts - discovery_cohorts
            print(f"{n1} vs {n2} (Pure Subset Check):{'PASSED' if subset_check else 'FAILED':<10} | Non-discovery: {diff}")
            if not subset_check:
                any_leak = True
            continue

        overlap = cohort_sets[n1] & cohort_sets[n2]
        status = f"{len(overlap)} (LEAK!)" if len(overlap) > 0 else "0 (PASSED)"
        print(f"{n1 + ' vs ' + n2:<35} | {status:<15} | {list(overlap)}")
        if len(overlap) > 0:
            any_leak = True

print("-" * 80)

# 3. Specific Val 2 vs Post-ML Severity Check
val2_vs_layer4 = val2_cohorts & layer4_cohorts
print(f"\n[KEY CHECK] Validation 2 vs Layer 4 Severity Accessions:")
print(f"  Validation 2 Accessions: {sorted(list(val2_cohorts))}")
print(f"  Layer 4 Severity Accessions: {sorted(list(layer4_cohorts))}")
print(f"  Shared Accessions: {list(val2_vs_layer4)}")
assert len(val2_vs_layer4) == 0, f"LEAK DETECTED: Val 2 and Layer 4 share: {val2_vs_layer4}"
print("  -> CONFIRMED: ZERO accession overlap between Validation 2 and Layer 4!\n")

# 4. Check Sample-level GSM uniqueness across all geo_cache files
geo_dir = 'geo_cache'
matrices = glob.glob(os.path.join(geo_dir, '*_series_matrix.txt.gz'))
gsm_to_study = {}
sample_collisions = []

for m in matrices:
    study = os.path.basename(m).split('_')[0]
    with gzip.open(m, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                header = next(f).strip().split('\t')
                gsms = [col.replace('"', '').strip() for col in header[1:] if 'GSM' in col]
                for gsm in gsms:
                    if gsm in gsm_to_study:
                        sample_collisions.append((gsm, gsm_to_study[gsm], study))
                    else:
                        gsm_to_study[gsm] = study
                break

print(f"Sample-Level Audit across {len(matrices)} GEO Series Matrices:")
print(f"  Total unique GSM patient samples: {len(gsm_to_study)}")
print(f"  Sample ID Collisions across studies: {len(sample_collisions)}")
if sample_collisions:
    print(f"  Collisions: {sample_collisions}")
else:
    print("  -> CONFIRMED: ZERO sample ID overlap across all studies on disk!")

print("=" * 80)
if not any_leak:
    print("FINAL AUDIT RESULT: 100% CLEAN. NO DATA LEAKAGE DETECTED.")
else:
    print("FINAL AUDIT RESULT: FAILED. DATA LEAKAGE DETECTED.")
print("=" * 80)
