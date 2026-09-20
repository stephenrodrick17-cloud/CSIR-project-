# Pre-Review Master Audit Report: Full Repository Verification

**Audit Date**: September 20, 2026  
**Auditor**: Antigravity Automated Verification Agent  
**Scope**: Complete repository top-to-bottom audit covering all 9 integrity categories prior to mentor review, plus resolution of 5 specific user-identified contradictions.

---

## Executive Summary of Audit Status

| Category | Description | Pre-Audit Status | Action Taken / Fixed | Post-Audit Status |
| :--- | :--- | :---: | :--- | :---: |
| **Cat 1** | Data Authenticity & Stand-in Check | **FLAGGED** (Historical Stand-in) | Audited all `.py` files; identified simulated kidney %TIF stand-in in superseded `fix_real_severity_pipeline.py`; quarantined script to `archive/`. Confirmed 0 synthetic values in final master tables. | **PASSED** |
| **Cat 2** | Cohort Isolation Guard (Discovery/Val1/Val2/ML) | **ATTENTION** (ML Auto-Check Missing) | Live guard passed; extended `verify_cohort_isolation()` in `discovery_config.py` to auto-inspect on-disk ML matrices. Confirmed 0 samples from Val1/Val2 in ML matrices. | **PASSED** |
| **Cat 3** | Hallucinated / Unverified Accessions | **VERIFIED** | Audited all 43 unique `GSE` accessions across all repo text/code. Confirmed 100% of the 16 active study cohorts exist on disk. Documented rejected screening candidates. | **PASSED** |
| **Cat 4** | Funnel Logic & Gene Count Consistency | **ATTENTION** (Cohort Table Count) | Confirmed gene funnel strictly narrows or stays flat ($86 \to 24 \to 18 \to 4$). Synchronized Discovery cohort count ($N=14$) across text and table. Excluded `TNXB` from core panel. | **PASSED** |
| **Cat 5** | One Final, Consistent Gene List | **VERIFIED** | Confirmed only ONE canonical 24-gene master evidence table exists in `README.md` and matches `results/final_master_evidence_table.csv`. | **PASSED** |
| **Cat 6** | Stale Files & Root Artifact Cleanup | **CLEANED** | Quarantined 12 unreferenced/duplicate CSV/PNG files and 4 superseded `.py` scripts from repo root to `archive/` with full inventory documentation. | **PASSED** |
| **Cat 7** | Image-Data Consistency | **REGENERATED** | Created `generate_final_figures.py`; regenerated `plots/study_design_funnel_corrected.png` and `plots/ml_4model_consensus_hub_biomarkers.png` directly from canonical CSV data. | **PASSED** |
| **Cat 8** | Literature Citations for Biological Claims | **UPDATED** | Added citations for VWF capillary rarefaction (Ebina 2004), AEBP1 collinearity grouping (Zou & Hastie 2005), and SSc reverse MR (Lopez-Isac 2019). | **PASSED** |
| **Cat 9** | Internal README Consistency | **RESOLVED** | Harmonized Section 1 ($N=14$), Section 2 table ($N=14$), and Section 3 sample matrix ($N=1,069$ across 9 complete series matrices). Zero numerical contradictions. | **PASSED** |

---

## Detailed Resolution of the Five Key Contradictions

### 1. ML Training Matrix Growth ($799 \to 1,069$ Samples)
* **The Contradiction**: The user asked why the ML matrix count grew after removing the `GSE58095` data leak.
* **Exact Evidence & Arithmetic Reconciliation**:
  * In the earlier 799-sample matrix (`commit 62f6ccd`), 8 cohorts were used: `GSE32537` (217), `GSE164760` (170), `GSE89377` (107), **`GSE58095` (102)**, `GSE66494` (61), `GSE110147` (48), `GSE53845` (48), `GSE10667` (46).
  * `GSE58095` ($n=102$) was an accidental inclusion of the Skin Validation 1 cohort.
  * To fix this data leak, `GSE58095` ($n=102$) was **permanently removed**, and the true Discovery cohorts for Skin were added from raw series matrices:
    * `GSE181549` (Agilent Skin Discovery): **+339 samples** (44 controls, 295 cases)
    * `GSE95065` (Affymetrix Skin Discovery): **+33 samples** (15 controls, 18 cases)
  * **Exact Formula**:
    $$799 - 102 (\text{removed GSE58095}) + 339 (\text{added GSE181549}) + 33 (\text{added GSE95065}) = 1,069 \text{ samples}$$
  * **Live Confirmation**: In `results/real_human_patient_ml_training_matrix.csv` and `results/real_human_patient_combat_corrected_matrix.csv`:
    * `GSE58095` is **100% genuinely absent** (`'GSE58095' in df['study'] == False`).
    * Zero samples from any Validation 1 or Validation 2 cohort are present.

### 2. MDK Reversal ("4/4 Unanimous Votes" $\to$ "Not Supported")
* **The Contradiction**: Why did MDK drop from unanimous ML consensus to "Not Supported" (2/5 stability)?
* **Exact Evidence & Methodological Cause**:
  1. **Dataset Contamination in Early Run**: In `commit 62f6ccd`, ML models were trained on the uncorrected 799-sample matrix containing `GSE58095` (Val1), where MDK had high array signal.
  2. **Pipeline Rigor Upgrades**:
     * ComBat empirical Bayes batch correction was applied across all 9 discovery studies.
     * Training data was balanced across organs (18 controls and 42 cases per organ).
     * Stability was evaluated over **5 random seeds** (42, 123, 456, 789, 999), in which MDK achieved $\ge 3$ votes in only 2 of 5 seeds (stability score: **2/5**).
  3. **Held-Out Validation 2 Failure**: Crucially, in independent held-out Validation 2 (`results/validation2_ecm_core_results.csv`), **MDK only replicated in 1 of 4 organs** (Liver adj. $p = 0.000107$; Kidney $p = 0.225$, Lung $p = 0.907$, Skin $p = 0.282$).
  * **Conclusion**: Per the project's strict multi-organ replication rule ($\ge 2/4$ organs in held-out validation), MDK failed validation and is properly classified as **Not Supported**.

### 3. Discovery Cohort Table Configuration & N=14
* **The Contradiction**: Section 2 table omitted `GSE104066`, `GSE104948`, `GSE104954`, and `GSE77627`, listing only 10 cohorts while text stated $N=14$.
* **Reconciliation**:
  * Reverted Section 2 table in `README.md` to the verified locked 14-cohort configuration:
    * Kidney (4): `GSE66494`, `GSE104066`, `GSE104948`, `GSE104954`
    * Liver (3): `GSE164760`, `GSE89377`, `GSE77627`
    * Lungs (4): `GSE10667`, `GSE110147`, `GSE32537`, `GSE53845`
    * Skin (3): `GSE130955`, `GSE181549`, `GSE95065`
  * Confirmed that `discovery_config.verify_cohort_isolation()` passes with $N=14$ Discovery cohorts, 4 Val 1 cohorts, and 4 Val 2 cohorts.
  * Clarified in Section 3 that 9 of these 14 cohorts had full series matrices ($N=1,069$ samples) for sample-level ComBat ML training.

### 4. Discovery Threshold: $|\log_2\text{FC}| \ge 0.585$ (1.5-fold change)
* **The Contradiction**: Discrepancy between mentions of 0.585 (1.5-fold) and 1.0 (2-fold).
* **Exact Code Grounding**:
  * In `preprocess_build_deg_csvs.py` (line 213) and `run_master_validation_pipeline.py` (line 24):
    `fc_cut = 0.585`
    `p_cut = 0.05`
  * The actual mathematical threshold used to produce the 86 core DEGs and 24 clean ECM genes was:
    $$|\log_2\text{FC}| \ge 0.585 \quad (\ge 1.5\text{-fold change}) \quad \text{and} \quad \text{FDR } p < 0.05$$
  * All mentions in `README.md`, `generate_final_figures.py`, and audit documentation have been synchronized to this exact threshold.

### 5. Removal of `TNXB` from Final Tier 1/Tier 2 Core Evidence Table
* **The Contradiction**: `TNXB` was included in the core master evidence table despite failing Discovery.
* **Exact Resolution**:
  * `TNXB` achieved statistical significance in only 3 of 4 organs in Discovery (failed Kidney, $p = 0.058$).
  * Per the project's strict funnel rule, `TNXB` **cannot be promoted into Tier 1 or Tier 2**.
  * `TNXB` has been **completely removed** from `results/final_master_evidence_table.csv`, the Section 4 table in `README.md`, and the consensus ranking plot.
  * The core panel consists strictly of the **24 Clean Core ECM Genes**.
  * `TNXB` is documented separately in Section 5 of `README.md` as an exploratory non-core finding (noting its independent Val2 concordance and reverse-MR causal signal in SSc, $p = 5.39 \times 10^{-7}$).

---

## Verification of Discovery Differential Expression Re-Run

* **Live Execution**: Re-ran `python preprocess_build_deg_csvs.py` fresh at **20-09-2026 23:18:24**.
* **Fresh Per-Organ DEG Counts** ($|\log_2\text{FC}| \ge 0.585$, FDR $p < 0.05$):
  * **Kidney** (4 cohorts): **11,758** significant DEGs (23,631 total genes)
  * **Liver** (3 cohorts): **2,268** significant DEGs (24,557 total genes)
  * **Lungs** (4 cohorts): **7,803** significant DEGs (22,965 total genes)
  * **Skin** (3 cohorts): **2,979** significant DEGs (26,457 total genes)
* **Downstream Propagation**:
  * Re-ran `python run_master_validation_pipeline.py`:
    * Intersection across all 4 organs: exactly **86 Conserved Pan-Fibrotic Core DEGs**.
    * Matrisome annotation: exactly **24 Clean Core ECM Genes**.
    * Validation 1 replication: **24/24 concordant (100%)**, 12/24 significant in $\ge 2/4$ organs.
    * Validation 2 held-out replication: **18/24 significant in $\ge 2/4$ organs (75%)**.
  * Re-ran `python run_master_consolidation.py`:
    * Produced fresh `results/final_master_evidence_table.csv` (24 rows).
    * Re-generated `plots/ml_4model_consensus_hub_biomarkers.png` and `plots/study_design_funnel_corrected.png`.

---

## Post-ML Severity Validation Isolation (Val 2 vs. Layer 4 Non-Overlap)

* **Audit Objective**: Ensure zero cohort overlap between Validation 2 (held-out validation) and Layer 4 (post-ML independent clinical severity replication cohorts), and verify zero data leakage into ML training.
* **Validation 2 Accessions ($N=4$)**:
  * Kidney: `GSE30529` (Affymetrix Microarray)
  * Liver: `GSE14323` (Affymetrix Microarray)
  * Lung: `GSE83717` (Illumina RNA-seq)
  * Skin: `GSE125362` (Agilent Microarray)
* **Layer 4 Severity Replication Accessions ($N=5$ new cohorts)**:
  * Liver: `GSE84044` (Microarray) + `GSE135251` (RNA-seq)
  * Lung: `GSE38958` (Microarray) + `GSE213001` (RNA-seq)
  * Skin: `GSE9285` (Microarray)
  * Kidney: Confirmed unavailable ($N=0$)
* **Exact Overlap Audit**:
  * `Val 2` $\cap$ `Layer 4`: **$\emptyset$ (0 cohorts, PASSED)**
  * `Val 1` $\cap$ `Layer 4`: **$\emptyset$ (0 cohorts, PASSED)**
  * `Discovery` $\cap$ `Layer 4`: **$\emptyset$ (0 cohorts, PASSED)**
  * `ML Training (N=1,069)` $\cap$ `Val 2`: **$\emptyset$ (0 cohorts, PASSED)**
  * `ML Training (N=1,069)` $\cap$ `Layer 4`: **$\emptyset$ (0 cohorts, PASSED)**
* **Sample-Level Verification**:
  * Scanned **3,133 unique patient GSM IDs** across 32 cached GEO series matrices.
  * Exact cross-study collisions: **0**.
* **Automated Guard Enforced**:
  * Formally encoded into `discovery_config.py` (`LAYER4_SEVERITY_ALL`) with active assertions halting execution if any contamination occurs.

