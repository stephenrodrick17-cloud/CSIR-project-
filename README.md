# Pan-Fibrotic Core Transcriptomic Program & Multi-Model Ensemble Machine Learning Biomarkers Across Human Kidney, Liver, Lung, and Skin

## 1. Locked Multi-Organ Study Design & Architecture

This repository hosts a multi-organ transcriptomic discovery and validation framework investigating core extracellular matrix (ECM) remodeling across four major human fibrotic diseases:
- **Kidney**: Chronic Kidney Disease / Renal Fibrosis
- **Liver**: Liver Cirrhosis / NASH Fibrosis
- **Lungs**: Idiopathic Pulmonary Fibrosis (IPF)
- **Skin**: Systemic Sclerosis (SSc)

All cohorts are organized into three strictly partitioned tiers:
1. **Tier 1 (Discovery Cohorts, $N=14$)**: Identifies conserved pan-fibrotic differentially expressed genes (DEGs) across all 4 organs ($\ge 2$-fold change, adjusted $p < 0.05$).
2. **Tier 2 (Validation 1 Cohorts, $N=4$)**: First independent replication across all 4 organs.
3. **Tier 3 (Validation 2 Held-Out Cohorts, $N=4$)**: Completely isolated, blinded validation cohort testing final pan-fibrotic universality.

---

## 2. Automated Isolation Guard Guarantee

The repository implements an automated architectural guard (`discovery_config.verify_cohort_isolation()`) that runs at the start of every analysis script. It enforces zero cohort overlap across all 12 experimental tiers:

| Organ | Discovery Cohorts (Tier 1) | Validation 1 (Tier 2) | Validation 2 Held-Out (Tier 3) |
| :--- | :--- | :--- | :--- |
| **Kidney** | `GSE66494` | `GSE200818` | `GSE30529` |
| **Liver** | `GSE164760`, `GSE89377` | `GSE77627` | `GSE14323` |
| **Lungs** | `GSE10667`, `GSE110147`, `GSE32537`, `GSE53845` | `GSE24206` | `GSE83717` |
| **Skin** | `GSE130955`, `GSE181549`, `GSE95065` | `GSE58095` | `GSE125362` |

---

## 3. Executive Summary of Discovery & Multi-Layer Validation

1. **Discovery Analysis**: Identified **86 conserved core genes** significantly dysregulated in 4/4 organ systems, including a dedicated **24 Clean Core ECM Gene Program** (`ecm_clean_genes.csv`).
2. **Validation 1 Replication**: 24/24 (100.0%) Clean Core ECM genes replicated in $\ge 1$ organ; 12/24 (50.0%) in $\ge 2$ organs.
3. **Validation 2 Held-Out Replication**: 24/24 (100.0%) replicated in $\ge 1$ organ; 18/24 (75.0%) replicated in $\ge 2$ organs; `COL15A1` and `AEBP1` replicated in 4/4 organs (100% universal concordance).

---

## 4. Real Patient Sample Training Matrix ($N=799$) & ComBat Batch Correction

To perform rigorous feature selection, individual patient-level expression profiles were extracted from 8 GEO series matrix files across Discovery and Validation 1 cohorts (total $N=799$ real human tissue samples; 147 Healthy Controls, 652 Fibrosis Cases):

| Dataset Accession | Organ | Total Samples | Controls | Fibrosis Cases | Platform Type |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **GSE66494** | Kidney | 61 | 8 | 53 | Affymetrix Human Gene 1.0 ST |
| **GSE89377** | Liver | 107 | 13 | 94 | Illumina HumanHT-12 V4.0 |
| **GSE164760** | Liver | 170 | 6 | 164 | Agilent Whole Human Genome Microarray |
| **GSE10667** | Lungs | 46 | 15 | 31 | Affymetrix Human Genome U133 Plus 2.0 |
| **GSE110147** | Lungs | 48 | 11 | 37 | RNA-seq (Illumina HiSeq 2000) |
| **GSE32537** | Lungs | 217 | 50 | 167 | Illumina HumanRef-8 v3.0 |
| **GSE53845** | Lungs | 48 | 8 | 40 | Agilent-014850 Whole Human Genome |
| **GSE58095** | Skin | 102 | 36 | 66 | Illumina HumanHT-12 V4.0 |
| **TOTAL** | **4 Organs** | **799** | **147** | **652** | **Multi-Platform Human Cohort** |

### Batch Correction & Balancing Protocol
1. **ComBat Batch Correction (`pycombat`)**: Applied empirical Bayes batch correction using `study` as the batch identifier to eliminate cross-platform baseline shifts (e.g., Agilent vs. Affymetrix vs. Illumina vs. RNA-seq).
2. **Organ & Class Balancing**: Constructed a 4-organ stratified cohort (Kidney $N=50$, Liver $N=60$, Lung $N=60$, Skin $N=60$; total $N=230$) with class-weighted estimators (`class_weight='balanced'`) to prevent lung/liver sample size dominance.

---

## 5. Master 4-Model Ensemble ML Feature Selection Results

Four complementary machine learning algorithms were trained on the batch-corrected, organ-balanced matrix:
1. **LASSO Logistic Regression** ($L_1$ penalty, 10-fold Cross-Validation)
2. **SVM-RFE** (Support Vector Machine - Recursive Feature Elimination, linear kernel)
3. **Random Forest** (500 estimators, class-weighted, feature importance threshold $> 1/24$)
4. **XGBoost** (500 boosted gradient trees, max depth 3, learning rate 0.05)

### Final Consensus Ranking Table (24 Clean Core ECM Genes)

| Gene | Total Votes | LASSO | SVM-RFE | Random Forest | XGBoost | Diagnostic ROC AUC | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **COL15A1** | **4 / 4** | Yes ($eta=+0.311$) | Yes (Rank 1) | Yes (Imp=0.1208) | Yes (Imp=0.2128) | **0.8177** | **Unanimous Core Hub** |
| **MDK** | **4 / 4** | Yes ($eta=+0.963$) | Yes (Rank 1) | Yes (Imp=0.0875) | Yes (Imp=0.0652) | **0.8051** | **Unanimous Core Hub** |
| **COL1A1** | **4 / 4** | Yes ($eta=+1.217$) | Yes (Rank 1) | Yes (Imp=0.0777) | Yes (Imp=0.0444) | **0.7935** | **Unanimous Core Hub** |
| **LTBP2** | **4 / 4** | Yes ($eta=+0.926$) | Yes (Rank 1) | Yes (Imp=0.0511) | Yes (Imp=0.0535) | **0.7155** | **Unanimous Core Hub** |
| **SERPINE2** | **3 / 4** | Yes ($eta=+0.019$) | No (Rank 15) | Yes (Imp=0.1050) | Yes (Imp=0.0744) | **0.7811** | **Consensus Hub** |
| **LAMC3** | **3 / 4** | Yes ($eta=-0.414$) | Yes (Rank 1) | Yes (Imp=0.0424) | No (Imp=0.0410) | **0.6086** | **Consensus Hub** |
| **COL3A1** | 2 / 4 | No ($eta=0.000$) | No (Rank 10) | Yes (Imp=0.0839) | Yes (Imp=0.1242) | 0.7723 | Candidate |
| **SERPINF2** | 2 / 4 | Yes ($eta=-0.304$) | No (Rank 8) | Yes (Imp=0.0474) | No (Imp=0.0406) | 0.7052 | Candidate |
| **CLEC2D** | 2 / 4 | Yes ($eta=+0.654$) | Yes (Rank 1) | No (Imp=0.0396) | No (Imp=0.0230) | 0.6930 | Candidate |
| **COL1A2** | 2 / 4 | Yes ($eta=-0.361$) | Yes (Rank 1) | No (Imp=0.0239) | No (Imp=0.0065) | 0.6910 | Candidate |
| **PDGFD** | 2 / 4 | Yes ($eta=+0.198$) | No (Rank 5) | Yes (Imp=0.0476) | No (Imp=0.0372) | 0.6706 | Candidate |
| **SERPINH1** | 2 / 4 | Yes ($eta=+0.372$) | Yes (Rank 1) | No (Imp=0.0198) | No (Imp=0.0323) | 0.6117 | Candidate |
| **CCL5** | 2 / 4 | Yes ($eta=-0.200$) | Yes (Rank 1) | No (Imp=0.0177) | No (Imp=0.0227) | 0.5625 | Candidate |
| **COLEC11** | 2 / 4 | Yes ($eta=-0.360$) | Yes (Rank 1) | No (Imp=0.0318) | No (Imp=0.0300) | 0.5615 | Candidate |
| **FGF14** | 1 / 4 | Yes ($eta=-0.094$) | No (Rank 4) | No (Imp=0.0228) | No (Imp=0.0095) | 0.6517 | Dropped Out |
| **CCL2** | 1 / 4 | Yes ($eta=-0.201$) | No (Rank 11) | No (Imp=0.0249) | No (Imp=0.0211) | 0.6320 | Candidate |
| **SVEP1** | 1 / 4 | Yes ($eta=-0.390$) | No (Rank 6) | No (Imp=0.0225) | No (Imp=0.0377) | 0.5344 | Candidate |
| **AEBP1** | 0 / 4 | No ($eta=0.000$) | No (Rank 3) | No (Imp=0.0178) | No (Imp=0.0215) | 0.6135 | Empirical Observation |
| **BMP1** | 0 / 4 | No ($eta=0.000$) | No (Rank 14) | No (Imp=0.0289) | No (Imp=0.0266) | 0.6049 | Candidate |
| **MFAP4** | 0 / 4 | No ($eta=0.000$) | No (Rank 2) | No (Imp=0.0234) | No (Imp=0.0145) | 0.6035 | Candidate |
| **VWF** | 0 / 4 | No ($eta=0.000$) | No (Rank 9) | No (Imp=0.0170) | No (Imp=0.0096) | 0.5930 | Candidate |
| **SPARCL1** | 0 / 4 | No ($eta=0.000$) | No (Rank 7) | No (Imp=0.0180) | No (Imp=0.0058) | 0.5915 | Candidate |
| **CCL21** | 0 / 4 | No ($eta=0.000$) | No (Rank 13) | No (Imp=0.0154) | No (Imp=0.0370) | 0.5825 | Candidate |
| **CCL19** | 0 / 4 | No ($eta=0.000$) | No (Rank 12) | No (Imp=0.0133) | No (Imp=0.0091) | 0.5731 | Candidate |

---

## 6. Audit Findings: `FGF14`, `MDK`, and `AEBP1`

### 1. Forensic Audit of `FGF14` vs. `MDK`
- **Why `FGF14` scored high before ComBat**: In the raw un-normalized matrix, `GSE164760` (Liver) had 170 samples with 96.5% disease cases (164 disease, 6 controls) and an extreme baseline platform mean of **20.67** (variance 34.46), compared to **0.14** in `GSE66494`. Un-batch-corrected models learned to predict `GSE164760` study membership rather than fibrotic biology. Following ComBat batch correction and organ balancing, `FGF14` **dropped out** from SVM-RFE, Random Forest, and XGBoost (falling to 1 vote).
- **Why `MDK` survived with 4/4 unanimous votes**: Midkine (`MDK`) demonstrated genuine, highly statistically significant upregulation across **all 4 organs individually** post-ComBat (Kidney: diff +4.19, $p = 9.75 \times 10^{-6}$; Liver: diff +1.47, $p = 3.57 \times 10^{-5}$; Lungs: diff +4.08, $p = 5.00 \times 10^{-24}$; Skin: diff +1.48, $p = 4.82 \times 10^{-6}$). Its 4/4 vote count is biologically authentic.

### 2. Empirical Reporting of `AEBP1`
- `AEBP1` (Adipocyte Enhancer-Binding Protein 1 / ACLP) is a matricellular, collagen-binding ECM protein that achieved 100% universal validation in univariate differential expression across all 4 independent Validation 2 cohorts.
- In multivariate machine learning feature selection, sparse regularized linear models (LASSO) and tree ensembles prioritize higher-margin collinear structural collagens (`COL15A1`, `COL1A1`) and growth factors (`MDK`, `LTBP2`, `SERPINE2`), resulting in 0 ML votes for `AEBP1`.
- We report this strictly as an empirical feature selection property of regularized multivariate classifiers, with no unsupported mechanistic claims regarding upstream transcriptional regulation.

---

## 7. Key Visualizations

### 4-Model Feature Consensus on Real Patient Samples
![4-Model Consensus Hub Biomarkers](plots/ml_4model_consensus_hub_biomarkers.png)

### Study Design Funnel: Locked Discovery & Multi-Layer Validation
![Study Design Funnel](plots/study_design_funnel_corrected.png)
