# Graphical Abstract Blueprint: Pan-Fibrotic Core ECM Program & Multi-Model Ensemble Machine Learning Biomarkers

This document provides the complete structural blueprint, exact numbers, dataset accessions, gene counts, and layout instructions for creating the **Graphical Abstract** (for BioRender, Adobe Illustrator, Inkscape, or PowerPoint).

---

## 1. High-Level Visual Flowchart & Architecture

```
┌───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                           GRAPHICAL ABSTRACT ARCHITECTURE                                          │
└───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

 [PANEL A: DISCOVERY]                      [PANEL B: BEFORE ML VALIDATION]               [PANEL C: DURING ML]
 14 GEO Cohorts (4 Organs)                 Validation 1 (Internal, N=4)                  1,069 Pure Discovery Biopsies
 1,069 Patient Biopsies                    4/4 Organs Direction Concordant (24/24)       ComBat Batch Harmonized
 Moderated eBayes |log2FC| >= 0.585                   │                                            │
 FDR p < 0.05                                         ▼                                            ▼
           │                               Validation 2 (Held-Out, N=4)                  4 ML Models x 5 Random Seeds
           ▼                               GSE30529, GSE14323, GSE83717, GSE125362       Random Forest, XGBoost,
 86 Conserved Core DEGs                    Non-Parametric Mann-Whitney U                 LASSO, SVM-RFE
           │                                          │                                            │
           ▼ (Matrisome Database)                     ▼                                            ▼
 24 Clean Core ECM Genes ────────────────► 18 / 24 Genes Survived                        4 High-Stability Hubs
                                           (Replicated in >= 2/4 Organs)                 (>= 4/5 Seeds, AUC > 0.76)
                                                      │                                            │
                                                      └──────────────────────┬─────────────────────┘
                                                                             │
                                                                             ▼
                                                           [PANEL D: AFTER ML EVIDENCE SYNTHESIS]
                                                           ══════════════════════════════════════
                                                           ★ TIER 1: FULL-SPECTRUM HUBS (N = 4)
                                                             • COL15A1 (#1 Hub, AUC 0.843, 5/5 Seeds)
                                                             • COL1A1  (Fibrillar Standard, AUC 0.790)
                                                             • SERPINE2 (Upregulated Antiprotease)
                                                             • SERPINF2 (Downregulated Antiprotease)
                                                           ──────────────────────────────────────
                                                           ★ TIER 2: VALIDATION-ONLY (N = 14)
                                                             (Replicated >= 2/4 Val2, ML Stability < 4/5)
                                                           ★ NOT SUPPORTED (N = 6)
                                                             (Replicated in only 1/4 Val2)
                                                                             │
                                                                             ▼
                                                           [PANEL E: SUPPLEMENTARY EXTERNAL SEVERITY]
                                                           4 External Cohorts (N = 474 Biopsies)
                                                           • Liver: GSE84044 (Scheuer S0-S4, N=124)
                                                           • Liver: GSE135251 (Kleiner F0-F4, N=216)
                                                           • Lung: GSE38958 (% Pred. DLCO, N=60)
                                                           • Skin: GSE9285 (mRSS Score, N=74)
```

---

## 2. Panel Breakdown & Exact Data Manifest

### Panel A: Cross-Organ Multi-Cohort Discovery
* **Goal**: Illustrate the un-biased extraction of conserved molecular drivers across distinct fibrotic organs.
* **Tissues Represented**: Human Kidney, Liver, Lung, and Skin.
* **Datasets ($N = 14$)**:
  - **Kidney** ($N=4$): `GSE104066` ($n=73$), `GSE66494` ($n=61$), `GSE104948` ($n=21$), `GSE104954` ($n=28$)
  - **Liver** ($N=3$): `GSE77627` ($n=58$), `GSE89377` ($n=107$), `GSE164760` ($n=170$)
  - **Lung** ($N=4$): `GSE10667` ($n=46$), `GSE110147` ($n=48$), `GSE32537` ($n=217$), `GSE53845` ($n=48$)
  - **Skin** ($N=3$): `GSE130955` ($n=24$), `GSE95065` ($n=33$), `GSE181549` ($n=339$)
* **Total Discovery Samples**: **1,069 human biopsies** (170 non-fibrotic controls, 899 fibrotic cases).
* **Funnel Counts**:
  1. Genome-Wide Transcripts: $\sim 25,000$ per platform.
  2. 4-Organ Conserved Core: **86 genes** (significant in all 4 organs, $|\log_2\text{FC}| \ge 0.585$, FDR $p < 0.05$).
  3. Clean Core Matrisome: **24 genes** (mapping to Collagens, ECM Regulators, Secreted Factors, Glycoproteins).

---

### Panel B: Before ML Validation (Validation 1 & Validation 2)
* **Goal**: Demonstrate rigorous multi-center biological replication before applying machine learning.
* **Validation 1 (Internal Replication, $N = 4$ Cohorts, $n = 384$ Biopsies)**:
  - Kidney: `GSE200818` ($n=64$)
  - Liver: `GSE162694` ($n=138$)
  - Lung: `GSE24206` ($n=23$)
  - Skin: `GSE58095` ($n=159$)
  - **Outcome**: **24/24 genes (100%)** had unanimous direction concordance with fibrosis across all 4 organs; 12/24 were statistically significant in $\ge 2/4$ organs.
* **Validation 2 (Held-Out Multi-Center Replication, $N = 4$ Cohorts, $n = 105$ Biopsies)**:
  - Kidney: `GSE30529` (Affymetrix Microarray, 10 DKD vs. 12 Controls)
  - Liver: `GSE14323` (Affymetrix Microarray, 41 Cirrhosis vs. 19 Controls)
  - Lung: `GSE83717` (Illumina RNA-seq, 6 IPF vs. 5 Controls)
  - Skin: `GSE125362` (Agilent Microarray, 8 SSc vs. 4 Controls)
* **Statistical Method**: Two-sided Mann-Whitney U test (non-parametric).
* **Before ML Survival Count**:
  - **18 Genes Survived** (Statistically significant in $\ge 2/4$ completely held-out organ cohorts):
    - *4/4 Organs*: `COL15A1`, `AEBP1` (2 genes)
    - *3/4 Organs*: `SERPINE2`, `COL3A1`, `COL1A2`, `CCL2`, `SERPINH1`, `VWF`, `CCL19` (7 genes)
    - *2/4 Organs*: `COL1A1`, `SERPINF2`, `LTBP2`, `LAMC3`, `CLEC2D`, `PDGFD`, `CCL5`, `SVEP1`, `MFAP4` (9 genes)
  - **6 Genes Failed** (Replicated in only 1/4 organ):
    - `MDK`, `FGF14`, `SPARCL1`, `BMP1`, `CCL21`, `COLEC11`.

---

### Panel C: Machine Learning Ensemble Training (During ML)
* **Goal**: Rank biomarkers by diagnostic power and multi-seed stability on real human patient biopsies.
* **Training Dataset**: Unified 9-cohort Discovery matrix of **1,069 genuine human patient biopsy samples** (ComBat batch-harmonized across platforms).
* **Isolation Guarantee**: Zero samples from Validation 1, Validation 2, or External datasets used in training.
* **4 ML Architectures**:
  1. Random Forest (Non-linear bagging ensemble)
  2. XGBoost (Gradient boosted decision trees)
  3. LASSO (L1-penalized sparse logistic regression)
  4. SVM-RFE (Support Vector Machine Recursive Feature Elimination)
* **Protocol**: Evaluated across **5 independent random seeds** (42, 123, 456, 789, 2024).
* **Stability Rule**: High stability required feature selection in $\ge 4/5$ seeds with mean diagnostic AUC $> 0.76$.

---

### Panel D: After ML Evidence Synthesis & Final Biomarker Tiers
* **Goal**: Present the final evidence classification uniting held-out replication and ML stability.
* **Total Clean Core Genes Evaluated**: 24.
* **Tier Stratification**:
  1. **Tier 1: Full-Spectrum Universal Biomarkers ($N = 4$)**
     - Passed Validation 2 ($\ge 2/4$ organs) **AND** High ML Stability ($\ge 4/5$ seeds):
     - **`COL15A1`**: #1 ranked universal matrix marker across all layers (Val 2: 4/4 organs, ML Stability: 5/5 seeds, Mean AUC: **0.8430**).
     - **`COL1A1`**: Canonical fibrillar collagen standard (Val 2: 2/4 organs, ML Stability: 5/5 seeds, Mean AUC: **0.7904**).
     - **`SERPINE2`**: Upregulated serine protease inhibitor / antiprotease axis (Val 2: 3/4 organs, ML Stability: 4/5 seeds, Mean AUC: **0.7864**).
     - **`SERPINF2`**: Downregulated alpha-2-antiplasmin / antiprotease shutdown (Val 2: 2/4 organs, ML Stability: 5/5 seeds, Mean AUC: **0.7635**).
  2. **Tier 2: Biological Validation-Only ($N = 14$)**
     - Passed Validation 2 ($\ge 2/4$ organs), but ML Stability $< 4/5$ seeds:
     - `AEBP1`, `COL3A1`, `COL1A2`, `VWF`, `CCL2`, `SERPINH1`, `CCL19`, `LTBP2`, `CCL5`, `SVEP1`, `LAMC3`, `CLEC2D`, `PDGFD`, `MFAP4`.
  3. **Not Supported ($N = 6$)**
     - Failed multi-organ Validation 2 (replicated in only 1/4 organ):
     - `MDK`, `FGF14`, `COLEC11`, `SPARCL1`, `BMP1`, `CCL21`.

---

### Panel E: Supplementary External Severity & Dose-Response Validation (After ML)
* **Goal**: Test whether priority Tier 1 hubs track progressive clinical disease severity and physiological decline.
* **4 External GEO Cohorts ($N = 474$ Human Biopsies)**:
  - **Liver (Microarray)**: `GSE84044` ($n = 124$) $\rightarrow$ Scheuer Fibrosis Staging (S0 to S4)
  - **Liver (RNA-seq)**: `GSE135251` ($n = 216$) $\rightarrow$ Kleiner Fibrosis Staging (F0 to F4) & NAS Score (0 to 8)
  - **Lung (Microarray)**: `GSE38958` ($n = 60$) $\rightarrow$ % Predicted DLCO & % Predicted FVC
  - **Skin (Microarray)**: `GSE9285` ($n = 74$) $\rightarrow$ Modified Rodnan Skin Score (mRSS 0 to 51)
* **Key Findings**:
  - `COL1A1` & `SERPINE2`: Cross-platform, multi-organ monotonic increase with advancing fibrosis stage.
  - `COL15A1`: Profound inverse correlation with pulmonary diffusing capacity (% Predicted DLCO, $\rho = -0.581, p = 5.73 \times 10^{-6}$).
  - `SERPINF2`: Progressive negative correlation tracking parenchymal functional loss.

---

## 3. Visual Styling & Color Recommendations for Graphical Abstract

| Element | Recommended Color / Hex | Symbol / Icon |
| :--- | :--- | :--- |
| **Kidney** | Crimson / Deep Red (`#D32F2F`) | Kidney silhouette |
| **Liver** | Amber / Burnt Orange (`#F57C00`) | Liver silhouette |
| **Lung** | Cerulean / Teal Blue (`#0288D1`) | Lungs silhouette |
| **Skin** | Purple / Violet (`#7B1FA2`) | Epidermis/Dermis icon |
| **Tier 1 Biomarkers** | Emerald Green / Gold Border (`#2E7D32`) | Star / Diamond badge |
| **Tier 2 Biomarkers** | Slate Blue (`#455A64`) | Circle / Checkmark badge |
| **Excluded Genes** | Muted Grey (`#9E9E9E`) | X / Dash badge |
| **Flowchart Arrows** | Dark Charcoal (`#263238`) | Thick curved vector arrow |

---

## 4. Accompanying Data Files in This Folder

1. [`pipeline_stages_and_datasets.csv`](pipeline_stages_and_datasets.csv): Complete list of all 27 datasets categorized by stage, platform, sample count, and usage.
2. [`gene_funnel_counts.csv`](gene_funnel_counts.csv): Step-by-step numbers, inputs, outputs, attrition criteria, and exact surviving gene lists.
3. [`all_24_genes_full_trajectory.csv`](all_24_genes_full_trajectory.csv): The complete 24-gene matrix tracking Discovery, Val 1, Val 2 (Before ML), ML Stability (During ML), and Final Tier (After ML).
