# Pan-Fibrotic Core Transcriptomic Program & Multi-Model Ensemble Machine Learning Biomarkers Across Human Kidney, Liver, Lung, and Skin

## 1. Locked Multi-Organ Study Design & Architecture

This repository implements a multi-organ transcriptomic discovery and validation framework investigating core extracellular matrix (ECM) remodeling across four human fibrotic diseases:
- **Kidney**: Chronic Kidney Disease / Renal Fibrosis
- **Liver**: Liver Cirrhosis / NASH Fibrosis
- **Lungs**: Idiopathic Pulmonary Fibrosis (IPF)
- **Skin**: Systemic Sclerosis (SSc)

All cohorts are organized into three strictly partitioned tiers:
1. **Tier 1 (Discovery Cohorts, $N=14$)**: Identifies conserved pan-fibrotic differentially expressed genes (DEGs) across all 4 organs ($\ge 2$-fold change, adjusted $p < 0.05$).
2. **Tier 2 (Validation 1 Cohorts, $N=4$)**: First independent cross-platform replication across all 4 organs.
3. **Tier 3 (Validation 2 Held-Out Cohorts, $N=4$)**: Completely isolated, blinded held-out validation cohorts testing pan-fibrotic universality.

---

## 2. Automated Isolation Guard Guarantee

The repository implements an automated architectural guard (`discovery_config.verify_cohort_isolation()`) that runs at the start of every analysis script. It enforces zero cohort overlap across all 12 experimental tiers and asserts that zero samples from Validation 1 or Validation 2 are used in ML training:

| Organ | Discovery Cohorts (Tier 1, $N=14$) | Validation 1 (Tier 2) | Validation 2 Held-Out (Tier 3) |
| :--- | :--- | :--- | :--- |
| **Kidney** | `GSE66494`, `GSE104066`, `GSE104948`, `GSE104954` | `GSE200818` | `GSE30529` |
| **Liver** | `GSE164760`, `GSE89377`, `GSE77627` | `GSE162694` | `GSE14323` |
| **Lungs** | `GSE10667`, `GSE110147`, `GSE32537`, `GSE53845` | `GSE24206` | `GSE83717` |
| **Skin** | `GSE130955`, `GSE181549`, `GSE95065` | `GSE58095` | `GSE125362` |

---

## 3. Pure Discovery Patient Sample Training Matrix ($N=1,069$)

To eliminate data leakage, all Validation 1 (`GSE58095`, `GSE200818`, `GSE162694`, `GSE24206`) and Validation 2 samples (`GSE30529`, `GSE14323`, `GSE83717`, `GSE125362`) are strictly quarantined. While all 14 Discovery cohorts contributed to the initial DEG boundary analyses, 9 cohorts possessed complete, unaggregated per-sample series matrices suitable for ComBat batch correction and sample-level multi-model machine learning training ($N=1,069$ total genuine human biopsy samples):

| Dataset Accession | Organ | Total Samples | Controls | Fibrosis Cases | Platform Type |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **GSE66494** | Kidney | 61 | 8 | 53 | Affymetrix Human Gene 1.0 ST |
| **GSE89377** | Liver | 107 | 13 | 94 | Illumina HumanHT-12 V4.0 |
| **GSE164760** | Liver | 170 | 6 | 164 | Agilent Whole Human Genome Microarray |
| **GSE10667** | Lungs | 46 | 15 | 31 | Affymetrix Human Genome U133 Plus 2.0 |
| **GSE110147** | Lungs | 48 | 11 | 37 | RNA-seq (Illumina HiSeq 2000) |
| **GSE32537** | Lungs | 217 | 50 | 167 | Illumina HumanRef-8 v3.0 |
| **GSE53845** | Lungs | 48 | 8 | 40 | Agilent-014850 Whole Human Genome |
| **GSE95065** | Skin | 33 | 15 | 18 | Affymetrix Human Genome U133A 2.0 |
| **GSE181549** | Skin | 339 | 44 | 295 | Agilent Whole Human Genome 4x44K V2 |
| **TOTAL** | **4 Organs** | **1,069** | **170** | **899** | **100% Pure Discovery Cohort** |

---

## 4. Master Consolidated Evidence Table (24 Clean Core ECM Genes + TNXB)

Uniting Discovery, Validation 1, Validation 2, disease-only histological severity correlation, 5-seed ML stability check, and average diagnostic ROC AUC:

| gene | discovery_status | val1_concordant_organs | val1_significant_organs | val2_concordant_organs | val2_significant_organs | severity_correlation_organs_significant | ml_stability_score | ml_diagnostic_auc | final_evidence_tier |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **COL15A1** | Pass (4/4 Organs) | 4/4 | 2/4 | 4/4 | 4/4 | 0 (None) | 5/5 | 0.843 | **Tier 1 (Full Spectrum)** |
| **COL1A1** | Pass (4/4 Organs) | 4/4 | 2/4 | 2/4 | 2/4 | 1 (Liver) | 5/5 | 0.7904 | **Tier 1 (Full Spectrum)** |
| **SERPINE2** | Pass (4/4 Organs) | 4/4 | 2/4 | 4/4 | 3/4 | 0 (None) | 4/5 | 0.7864 | **Tier 1 (Full Spectrum)** |
| **SERPINF2** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 2/4 | 0 (None) | 5/5 | 0.7635 | **Tier 1 (Full Spectrum)** |
| **COL3A1** | Pass (4/4 Organs) | 4/4 | 2/4 | 4/4 | 3/4 | 0 (None) | 3/5 | 0.7815 | **Tier 2 (Validation-Only)** |
| **COL1A2** | Pass (4/4 Organs) | 4/4 | 2/4 | 4/4 | 3/4 | 1 (Liver) | 0/5 | 0.687 | **Tier 2 (Validation-Only)** |
| **LTBP2** | Pass (4/4 Organs) | 4/4 | 2/4 | 2/4 | 2/4 | 0 (None) | 0/5 | 0.6827 | **Tier 2 (Validation-Only)** |
| **LAMC3** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 2/4 | 0 (None) | 2/5 | 0.6801 | **Tier 2 (Validation-Only)** |
| **CLEC2D** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 2/4 | 0 (None) | 1/5 | 0.6717 | **Tier 2 (Validation-Only)** |
| **PDGFD** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 2/4 | 0 (None) | 0/5 | 0.6161 | **Tier 2 (Validation-Only)** |
| **CCL2** | Pass (4/4 Organs) | 4/4 | 1/4 | 4/4 | 3/4 | 0 (None) | 1/5 | 0.6134 | **Tier 2 (Validation-Only)** |
| **SERPINH1** | Pass (4/4 Organs) | 4/4 | 1/4 | 4/4 | 3/4 | 0 (None) | 0/5 | 0.6057 | **Tier 2 (Validation-Only)** |
| **AEBP1** | Pass (4/4 Organs) | 4/4 | 1/4 | 4/4 | 4/4 | 1 (Liver) | 0/5 | 0.5968 | **Tier 2 (Validation-Only)** |
| **VWF** | Pass (4/4 Organs) | 4/4 | 2/4 | 4/4 | 3/4 | 1 (Liver) | 0/5 | 0.5706 | **Tier 2 (Validation-Only)** |
| **CCL5** | Pass (4/4 Organs) | 4/4 | 2/4 | 2/4 | 2/4 | 0 (None) | 0/5 | 0.5688 | **Tier 2 (Validation-Only)** |
| **CCL19** | Pass (4/4 Organs) | 4/4 | 1/4 | 4/4 | 3/4 | 0 (None) | 0/5 | 0.5687 | **Tier 2 (Validation-Only)** |
| **SVEP1** | Pass (4/4 Organs) | 4/4 | 2/4 | 2/4 | 2/4 | 0 (None) | 1/5 | 0.5396 | **Tier 2 (Validation-Only)** |
| **MFAP4** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 2/4 | 0 (None) | 2/5 | 0.5276 | **Tier 2 (Validation-Only)** |
| **TNXB** | Non-Core (3/4 Organs) | 2/4 | 1/4 | 4/4 | 2/4 | 0 (None) | Not tested (data unavailable) | N/A | **Tier 2 (Validation-Only)** |
| **MDK** | Pass (4/4 Organs) | 4/4 | 2/4 | 2/4 | 1/4 | 0 (None) | 2/5 | 0.7953 | **Not Supported** |
| **FGF14** | Pass (4/4 Organs) | 4/4 | 2/4 | 2/4 | 1/4 | 0 (None) | 0/5 | 0.6736 | **Not Supported** |
| **SPARCL1** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 1/4 | 0 (None) | 0/5 | 0.6156 | **Not Supported** |
| **BMP1** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 1/4 | 0 (None) | 0/5 | 0.5913 | **Not Supported** |
| **CCL21** | Pass (4/4 Organs) | 4/4 | 1/4 | 2/4 | 1/4 | 0 (None) | 0/5 | 0.5675 | **Not Supported** |
| **COLEC11** | Pass (4/4 Organs) | 4/4 | 2/4 | 2/4 | 1/4 | 0 (None) | 0/5 | 0.5196 | **Not Supported** |

---

## 5. Explicit Headline Findings

1. **Tier 1 (Full Spectrum) Qualification**: Exactly four genes qualify for Tier 1 status (`COL15A1`, `COL1A1`, `SERPINE2`, and `SERPINF2`) by demonstrating multi-organ replication in independent Validation 2 cohorts ($\ge 2/4$ organs) AND robust machine learning consensus ($\ge 4/5$ seed stability score with mean AUC $> 0.76$).
2. **`TNXB` Qualification**: `TNXB` demonstrates **4/4 universal Validation 2 concordance** (100% concordance across Kidney, Liver, Lung, and Skin) with significant replication in 2/4 organs (Liver adj_p = 0.000322, Lung adj_p = 0.0242), and a verified reverse-MR causal signal in Systemic Sclerosis ($p = 5.39 \times 10^{-7}$, lead SNP rs6926894; Lopez-Isac et al., *Nat Commun*, 2019). With ML stability unassessed in the discovery microarrays, `TNXB` qualifies for **Tier 2 (Validation-Only)** status.
3. **`COL15A1`**: Emerges as the top pan-fibrotic matrix biomarker with 4/4 universal Validation 2 concordance, 5/5 ML stability score, and the highest individual diagnostic ROC AUC (**0.8430**).
4. **`VWF` Diagnostic**: While universally upregulated in Kidney, Liver, and Skin fibrosis, `VWF` is significantly down-regulated in fibrotic Lung tissue ($\Delta = -3.82$, $p = 5.24 \times 10^{-17}$), consistent with severe pulmonary capillary loss and vascular rarefaction in end-stage IPF (Ebina et al., *Am J Respir Crit Care Med*, 2004), which accounts for its lower pooled cross-organ linear ML performance.
5. **`AEBP1` Empirical Finding**: `AEBP1` demonstrates 4/4 universal Validation 2 replication and significant Liver histological severity correlation ($p = 0.000053$). Its zero selection count in L1-regularized linear classifiers (LASSO/SVM-RFE) reflects the known grouping effect where sparse models select one representative from highly collinear feature sets (Zou & Hastie, *J R Stat Soc B*, 2005) — in this case, primary structural collagens (`COL15A1`, `COL1A1`).

---

## 6. Key Visualizations

### Master Evidence Ranking Across Clean Core ECM Genes
![Master Evidence Ranking](plots/ml_4model_consensus_hub_biomarkers.png)

### Study Design Funnel: Locked Discovery & Multi-Layer Validation
![Study Design Funnel](plots/study_design_funnel_corrected.png)
