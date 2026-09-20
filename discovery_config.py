# -*- coding: utf-8 -*-
"""
Locked Cohort Configuration and Automated Isolation Guard for CSIR Pan-Fibrotic Analysis.

Enforces strict mathematical isolation between Discovery, Validation 1, Validation 2,
and Machine Learning training matrices.
"""

# Explicit 3-Tier Accession Sets per Organ
COHORT_TIERS = {
    "Kidney": {
        "DISCOVERY": {"GSE104066", "GSE66494", "GSE104948", "GSE104954"},
        "VALIDATION_1": {"GSE200818"},
        "VALIDATION_2": {"GSE30529"},
    },
    "Liver": {
        "DISCOVERY": {"GSE77627", "GSE89377", "GSE164760"},
        "VALIDATION_1": {"GSE162694"},
        "VALIDATION_2": {"GSE14323"},
    },
    "Lungs": {
        "DISCOVERY": {"GSE110147", "GSE32537", "GSE53845", "GSE10667"},
        "VALIDATION_1": {"GSE24206"},
        "VALIDATION_2": {"GSE83717"},
    },
    "Skin": {
        "DISCOVERY": {"GSE130955", "GSE95065", "GSE181549"},
        "VALIDATION_1": {"GSE58095"},
        "VALIDATION_2": {"GSE125362"},
    },
}

# Aggregate all Validation 1 and Validation 2 accessions across all organs
VAL1_ALL = set().union(*[tiers["VALIDATION_1"] for tiers in COHORT_TIERS.values()])
VAL2_ALL = set().union(*[tiers["VALIDATION_2"] for tiers in COHORT_TIERS.values()])
DISCOVERY_ALL = set().union(*[tiers["DISCOVERY"] for tiers in COHORT_TIERS.values()])
NON_DISCOVERY_ALL = VAL1_ALL | VAL2_ALL

# Automatically derive exclusions from Discovery
EXCLUDE_FROM_DISCOVERY = {
    organ: tiers["VALIDATION_1"] | tiers["VALIDATION_2"]
    for organ, tiers in COHORT_TIERS.items()
}


def verify_cohort_isolation(ml_studies=None):
    """
    Automated assertion-based verification of cohort isolation across all 3 tiers
    and ML training data.
    Raises AssertionError immediately if any overlap exists.
    """
    header = "=" * 80
    print(header)
    print("[GUARD] AUTOMATED COHORT ISOLATION & TIER VERIFICATION GUARD")
    print(header)
    
    total_discovery = 0
    total_val1 = 0
    total_val2 = 0

    print(f"{'Organ':<10} | {'Discovery Cohorts':<35} | {'Validation 1':<15} | {'Validation 2':<15}")
    print("-" * 80)

    for organ, tiers in COHORT_TIERS.items():
        disc = tiers["DISCOVERY"]
        val1 = tiers["VALIDATION_1"]
        val2 = tiers["VALIDATION_2"]

        # Strict Mathematical Assertions: Zero Overlap Allowed
        disc_val1_overlap = disc & val1
        disc_val2_overlap = disc & val2
        val1_val2_overlap = val1 & val2

        assert len(disc_val1_overlap) == 0, (
            f"FATAL: Cohort overlap detected in {organ} between Discovery and Validation 1: {disc_val1_overlap}"
        )
        assert len(disc_val2_overlap) == 0, (
            f"FATAL: Cohort overlap detected in {organ} between Discovery and Validation 2: {disc_val2_overlap}"
        )
        assert len(val1_val2_overlap) == 0, (
            f"FATAL: Cohort overlap detected in {organ} between Validation 1 and Validation 2: {val1_val2_overlap}"
        )

        total_discovery += len(disc)
        total_val1 += len(val1)
        total_val2 += len(val2)

        disc_str = ", ".join(sorted(disc))
        val1_str = ", ".join(sorted(val1))
        val2_str = ", ".join(sorted(val2))
        print(f"{organ:<10} | {disc_str:<35} | {val1_str:<15} | {val2_str:<15}")

    print("-" * 80)
    print(f"TOTALS     | {total_discovery} Discovery Cohorts                 | {total_val1} Val 1 Cohorts  | {total_val2} Val 2 Cohorts")
    print("STATUS     | [PASSED] Zero data leakage across Discovery, Val 1, and Val 2 tiers.")
    
    # ML Training Isolation Guard (auto-check on-disk matrices if ml_studies is None)
    if ml_studies is None:
        import os
        import pandas as pd
        base_dir = os.path.dirname(os.path.abspath(__file__))
        ml_paths = [
            os.path.join(base_dir, "results", "real_human_patient_ml_training_matrix.csv"),
            os.path.join(base_dir, "results", "real_human_patient_combat_corrected_matrix.csv"),
        ]
        found_studies = set()
        for p in ml_paths:
            if os.path.exists(p):
                df_ml = pd.read_csv(p, usecols=["study"])
                found_studies.update(df_ml["study"].unique())
        if found_studies:
            ml_studies = sorted(list(found_studies))

    if ml_studies is not None:
        print("-" * 80)
        print("[ML GUARD] VERIFYING ML TRAINING DATA STRICT ISOLATION...")
        ml_set = set(ml_studies)
        ml_val1_leak = ml_set & VAL1_ALL
        ml_val2_leak = ml_set & VAL2_ALL
        
        print(f"ML Training Studies: {sorted(list(ml_set))}")
        print(f"Validation 1 Restricted: {sorted(list(VAL1_ALL))}")
        print(f"Validation 2 Restricted: {sorted(list(VAL2_ALL))}")
        
        assert len(ml_val1_leak) == 0, (
            f"FATAL: ML training matrix contains samples from Validation 1 cohorts: {ml_val1_leak}"
        )
        assert len(ml_val2_leak) == 0, (
            f"FATAL: ML training matrix contains samples from Validation 2 cohorts: {ml_val2_leak}"
        )
        print("ML STATUS  | [PASSED] ZERO samples from Validation 1 or Validation 2 present in ML training.")
    print(header + "\n")


# Self-test on import
if __name__ == "__main__":
    verify_cohort_isolation()
