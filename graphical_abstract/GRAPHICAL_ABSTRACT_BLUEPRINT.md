# Graphical Abstract Blueprint: Pan-Fibrotic Core ECM Program & Multi-Model Consensus Machine Learning Hub Biomarkers

This document provides the complete structural blueprint, exact numbers, dataset accessions, gene counts, and layout instructions for creating the **Graphical Abstract** (for BioRender, Adobe Illustrator, Inkscape, or PowerPoint), aligned with the audited **3/2/4 Hub Gene Hierarchy**.

---

## 1. High-Level Visual Flowchart & Architecture

```
┌───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                           GRAPHICAL ABSTRACT ARCHITECTURE                                          │
└───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

 [PANEL A: DISCOVERY]                      [PANEL B: PRE-ML REPLICATION]                 [PANEL C: CONSENSUS ML & WGCNA]
 14 GEO Cohorts (4 Organs)                 Validation 1 (Internal, N=4)                  1,069 Discovery Biopsies (9 Cohorts)
 1,069 Patient Biopsies                    4/4 Organs Direction Concordant (24/24)       ComBat Batch Harmonized
 Moderated eBayes |log2FC| >= 0.585                   │                                            │
 FDR p < 0.05                                         ▼                                            ▼
           │                               Validation 2 (Held-Out, N=4)                  5 Feature Selection Methods
           ▼                               GSE30529, GSE14323, GSE83717, GSE125362       LASSO, SVM-RFE, RF, XGBoost,
 86 Conserved Core DEGs                    Non-Parametric Mann-Whitney U                 TOM Co-expression / WGCNA
           │                                          │                                            │
           ▼ (Matrisome Database)                     ▼                                            ▼
 24 Clean Core ECM Genes ────────────────► 18 / 24 Genes Survived                        9 Consensus Hub Genes (>=3/5)
                                           (Replicated in >= 2/4 Organs)                 5 Unanimous Hubs (5/5 Votes)
                                                      │                                            │
                                                      └──────────────────────┬─────────────────────┘
                                                                             │
                                                                             ▼
                                                           [PANEL D: STAGE 3 DUAL-PLATFORM VALIDATION]
                                                           ══════════════════════════════════════════
                                                           ★ 3 STRICT CROSS-PLATFORM VALIDATED HUBS
                                                             (p < 0.05 & Direction-Concordant on BOTH Microarray & RNA-seq)
                                                             • COL15A1 (Basement membrane organizer, AUC 0.843)
                                                             • COL3A1  (Early repair / fibrillar collagen, AUC 0.852)
                                                             • SERPINE2 (Extracellular antiprotease regulator, AUC 0.786)
                                                           ──────────────────────────────────────────
                                                           ★ 2 DIRECTION-CONCORDANT UNANIMOUS HUBS
                                                             (5/5 ML Votes, 100% Direction-Concordant, RNA-seq p >= 0.05)
                                                             • COL1A1  (Primary fibrillar collagen, RNA-seq p = 0.101)
                                                             • SERPINF2 (Alpha-2-antiplasmin, RNA-seq p = 0.788)
                                                           ──────────────────────────────────────────
                                                           ★ 4 STAGE 2 CANDIDATE HUBS (Failed Stage 3)
                                                             • LAMC3, LTBP2, SVEP1 (Direction discordant)
                                                             • MDK (RNA-seq non-significant, p = 0.760)
                                                                             │
                                                                             ▼
                                                           [PANEL E: DOWNSTREAM PUBLICATION SUITES (MODS 1–6)]
                                                           • Mod 1: PPI Network (STRING v12 Collagen Triad)
                                                           • Mod 2: Early Detection (GSE84044 Scheuer Stages)
                                                           • Mod 3: Diagnostic Nomogram & DCA (AUC 0.900 / 0.924)
                                                           • Mod 4: Immune Infiltration (Stromal/Myofibroblast)
                                                           • Mod 5: HPA IHC Baseline & Patient Biopsy RNA
                                                           • Mod 6: Functional Enrichment (GO/KEGG/DO via Enrichr)
```

---

## 2. Panel Breakdown & Exact Data Manifest

### Panel A: Cross-Organ Multi-Cohort Discovery
* **Goal**: Illustrate the unbiased extraction of conserved molecular drivers across distinct fibrotic organs.
* **Tissues Represented**: Human Kidney, Liver, Lung, and Skin.
* **Datasets ($N = 14$)**:
  - **Kidney** ($N=4$): `GSE104066` ($n=73$), `GSE66494` ($n=61$), `GSE104948` ($n=21$), `GSE104954` ($n=28$)
  - **Liver** ($N=3$): `GSE77627` ($n=58$), `GSE89377` ($n=107$), `GSE164760` ($n=170$)
  - **Lung** ($N=4$): `GSE10667` ($n=46$), `GSE110147` ($n=48$), `GSE32537` ($n=217$), `GSE53845` ($n=48$)
  - **Skin** ($N=3$): `GSE130955` ($n=24$), `GSE95065` ($n=33$), `GSE181549` ($n=339$)
* **Total Discovery Samples**: **1,069 human biopsies** (170 non-fibrotic controls, 899 fibrotic cases across 9 cohorts with paired controls).
* **Funnel Counts**:
  1. Genome-Wide Transcripts: $\sim 25,000$ per platform.
  2. 4-Organ Conserved Core: **86 genes** (significant in all 4 organs, $|\log_2\text{FC}| \ge 0.585$, FDR $p < 0.05$).
  3. Clean Core Matrisome: **24 genes** (mapping to Collagens, ECM Regulators, Secreted Factors, Glycoproteins).

---

### Panel B: Pre-ML Replication (Validation 1 & Validation 2)
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
* **Pre-ML Survival Count**:
  - **18 Genes Survived** (Statistically significant in $\ge 2/4$ held-out organ cohorts):
    - *4/4 Organs*: `COL15A1`, `AEBP1` (2 genes)
    - *3/4 Organs*: `SERPINE2`, `COL3A1`, `COL1A2`, `CCL2`, `SERPINH1`, `VWF`, `CCL19` (7 genes)
    - *2/4 Organs*: `COL1A1`, `SERPINF2`, `LTBP2`, `LAMC3`, `CLEC2D`, `PDGFD`, `CCL5`, `SVEP1`, `MFAP4` (9 genes)
  - **6 Genes Failed** (Replicated in only 1/4 organ):
    - `MDK`, `FGF14`, `SPARCL1`, `BMP1`, `CCL21`, `COLEC11`.

---

### Panel C: Consensus Machine Learning & WGCNA (During ML)
* **Goal**: Rank biomarkers by diagnostic power and multi-seed stability on real human patient biopsies.
* **Training Dataset**: Unified 9-cohort Discovery matrix of **1,069 genuine human patient biopsy samples** (ComBat batch-harmonized across platforms).
* **Isolation Guarantee**: Zero samples from Validation 1, Validation 2, or External datasets used in training.
* **5 Feature Selection Architectures**:
  1. LASSO (L1-penalized sparse logistic regression)
  2. SVM-RFE (Support Vector Machine Recursive Feature Elimination)
  3. Random Forest (Class-weighted bagging ensemble)
  4. XGBoost (Gradient boosted decision trees)
  5. TOM Co-expression Clustering / WGCNA (Topological Overlap Matrix soft-thresholding)
* **Protocol**: Evaluated across **5 independent random seeds** (42, 123, 456, 789, 999).
* **Consensus Selection**:
  - **9 Candidate Hub Genes** ($\ge 3/5$ votes): `COL15A1`, `COL1A1`, `COL3A1`, `SERPINE2`, `SERPINF2`, `LAMC3`, `LTBP2`, `MDK`, `SVEP1`.
  - **5 Unanimous Hubs** (5/5 votes): `COL15A1`, `COL1A1`, `COL3A1`, `SERPINE2`, `SERPINF2`.

---

### Panel D: Stage 3 Cross-Platform Validation & Final Hub Hierarchy
* **Goal**: Establish the locked 3-tier biomarker classification across platforms.
* **Tiers**:
  1. **3 Strict Cross-Platform Validated Hubs**:
     - `COL15A1`, `COL3A1`, `SERPINE2` ($p < 0.05$ & concordant direction on both Microarray and RNA-seq).
  2. **2 Direction-Concordant Unanimous ML Hubs**:
     - `COL1A1`, `SERPINF2` (5/5 ML votes, 100% direction-concordant; RNA-seq $p = 0.101$ / $p = 0.788$).
  3. **4 Stage 2 Candidate Hubs**:
     - `LAMC3`, `LTBP2`, `MDK`, `SVEP1` (failed Stage 3 dual-platform replication).

---

### Panel E: Stage 4 Downstream Publication Suites (Modules 1–6)
* **Module 1: PPI Network (STRING v12)** — High-confidence ($\ge 0.700$) direct interactions forming the core Collagen Triad (`COL1A1`–`COL3A1`–`COL15A1`).
* **Module 2: Early Detection Triptych (GSE84044)** — Scheuer S0 vs. S1–S2 subclinical ROC (Ensemble AUC 0.716 / 0.719) + trajectory switch dynamics + DCA.
* **Module 3: Diagnostic Nomogram & DCA (N=1,069)** — Multivariable scoring ruler (5-Gene AUC: 0.9002; 9-Gene AUC: 0.9238, evaluated in-sample).
* **Module 4: Immune & Stromal Infiltration (GSE84044)** — Correlations with activated myofibroblasts (`ACTA2`), matrix stroma (`POSTN`), macrophages (`CD68`, `CD163`), and endothelial cells (`PECAM1`).
* **Module 5: HPA IHC Baseline & Patient Biopsy RNA Characterization** — Normal baseline IHC across 4 human organs + empirical patient biopsy $\log_2\text{FC}$.
* **Module 6: Functional Enrichment (GO / KEGG / DO via Enrichr)** — Live Enrichr query: ECM Organization, Collagen Fibril Organization, ECM-receptor interaction.

---

## 3. Visual Styling & Color Recommendations for Graphical Abstract

| Element | Recommended Color / Hex | Symbol / Icon |
| :--- | :--- | :--- |
| **Kidney** | Crimson / Deep Red (`#D32F2F`) | Kidney silhouette |
| **Liver** | Amber / Burnt Orange (`#F57C00`) | Liver silhouette |
| **Lung** | Cerulean / Teal Blue (`#0288D1`) | Lungs silhouette |
| **Skin** | Purple / Violet (`#7B1FA2`) | Epidermis/Dermis icon |
| **3 Strict Validated Hubs** | Royal Blue / Gold Badge (`#1E40AF`) | Star / Shield badge |
| **2 Direction-Concordant Hubs** | Jade Green (`#059669`) | Diamond badge |
| **4 Stage 2 Candidate Hubs** | Slate Gray (`#64748B`) | Circle badge |
| **Excluded Genes** | Muted Grey (`#9E9E9E`) | Dash / X badge |
| **Flowchart Arrows** | Dark Charcoal (`#263238`) | Thick curved vector arrow |
