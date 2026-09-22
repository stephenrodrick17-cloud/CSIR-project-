# Canonical Publication Figures Directory (`plots/`)

This directory contains the 15 final, publication-grade figures for the **CSIR Pan-Fibrotic Multi-Organ Biomarker Study**. Every figure has been generated directly from locked, audited patient cohorts and public database APIs.

---

### Core Hierarchy & Locked Parameters
- **Discovery Cohorts**: $N = 14$ independent cohorts across 4 human organs (Kidney: 4, Liver: 3, Lung: 4, Skin: 3).
- **Discovery Sample Size**: $N = 1,069$ human patient tissue samples (170 non-fibrotic controls, 899 fibrotic cases across 9 cohorts with paired controls).
- **Core Significance Filter**: eBayes moderated $t$-test $|\log_2\text{FC}| \ge 0.585$ (1.5-fold), Benjamini-Hochberg FDR $p < 0.05$.
- **Pan-Fibrotic Conserved Core**: **86 DEGs** shared across all 4 organs.
- **Clean Core ECM Program**: **24 genes** overlapping the Human Matrisome Masterlist ($N=1,027$).
- **Machine Learning & WGCNA Hub Discovery**: 5-method consensus on $N=1,069$ pure discovery samples (LASSO, SVM-RFE, RF, XGBoost, TOM Co-expression Clustering / Python WGCNA) $\to$ **9 Multi-Model Consensus Hub Genes** ($\ge 3/5$ votes): `COL15A1`, `COL1A1`, `COL3A1`, `SERPINE2`, `SERPINF2`, `LAMC3`, `LTBP2`, `MDK`, `SVEP1`.
- **Unanimous Consensus Hubs (5/5 Votes)**: `COL15A1`, `COL1A1`, `COL3A1`, `SERPINE2`, `SERPINF2` (selected unanimously by all 5 feature selection architectures).
- **Cross-Platform Validation (Microarray vs. RNA-seq)**:
  - **3 Strict Cross-Platform Validated Hubs**: `COL15A1`, `COL3A1`, `SERPINE2` ($p < 0.05$, AUC $> 0.76$, concordant direction on both Microarray and RNA-seq).
  - **2 Direction-Concordant Unanimous ML Hubs**: `COL1A1`, `SERPINF2` (100% direction-concordant across platforms; RNA-seq $p = 0.101$ / $p = 0.788$).
  - **4 Stage 2 Candidate Hubs**: `LAMC3`, `LTBP2`, `MDK`, `SVEP1` (failed Stage 3 dual-platform validation).

---

### Upstream Discovery & Conserved Core Figures (3 Figures)

| File | Pipeline Stage | Description | Source Script |
| :--- | :--- | :--- | :--- |
| `venn_4organ_manual_ellipses.png` | **Stage 1: Discovery Intersection** | 4-Organ Venn diagram showing DEGs for Kidney (11,758), Liver (2,268), Lung (7,803), Skin (2,979), and the conserved 86-gene core intersection. | `venn_4organs.py` |
| `upset_plot_4organs_corrected.png` | **Stage 1: Discovery Intersection** | UpSet plot showing all 15 subset intersections across the 4 organs; mathematically and cardinality identical to the 4-organ Venn. | `generate_upset_plot.py` |
| `venn_organ_vs_ecm_compendium.png` | **Stage 1: ECM Filtering** | Multi-panel Venn compendium: 2x2 grid of organ DEGs vs. ECM Masterlist (1,027 genes), 3-way ECM intersection, and 4-way pan-fibrotic core (86 DEGs $\to$ 24 ECM genes). | `venn_organ_vs_ecm.py` |

---

### Stage 4 Parallel Downstream Publication Suites: 5-Gene vs. 9-Gene (12 Figures across 6 Modules)

Both suites are generated separately and preserved side-by-side using 100% real human patient cohorts and public databases.

| Module | 5-Gene Unanimous Consensus Suite | 9-Gene Multi-Model Consensus Suite | Key Methodological Details |
| :--- | :--- | :--- | :--- |
| **1. PPI Network (STRING v12)** | `hub_genes_ppi_network_5genes.png`<br>`results/hub_genes_ppi_network_nodes_5genes.csv`<br>`results/hub_genes_ppi_network_edges_5genes.csv` | `hub_genes_ppi_network_9genes.png`<br>`results/hub_genes_ppi_network_nodes_9genes.csv`<br>`results/hub_genes_ppi_network_edges_9genes.csv` | **5-Gene**: 17 nodes, 68 edges ($p < 10^{-16}$), highlighting core collagens and serpins with top functional interactors.<br>**9-Gene**: 9 pure consensus hub nodes with direct high-confidence interactions ($\ge 0.700$, collagen triad `COL1A1`–`COL3A1`–`COL15A1`). |
| **2. Early Detection Triptych** | `hub_genes_early_detection_triptych_5genes.png`<br>`results/hub_genes_early_stage_roc_metrics_5genes.csv` | `hub_genes_early_detection_triptych_9genes.png`<br>`results/hub_genes_early_stage_roc_metrics_9genes.csv` | **5-Gene**: Early ROC (S0 vs S1–S2, $N=96$) ensemble $\text{AUC} = 0.716$; Net benefit $+0.396$ at $p_t=25\%$.<br>**9-Gene**: Early ROC ensemble $\text{AUC} = 0.719$; Net benefit $+0.383$ at $p_t=25\%$. Evaluated on histological Scheuer stages (`GSE84044`). |
| **3. Diagnostic Nomogram & DCA** | `hub_genes_nomogram_and_dca_5genes.png`<br>`results/hub_genes_nomogram_parameters_5genes.csv` | `hub_genes_nomogram_and_dca_9genes.png`<br>`results/hub_genes_nomogram_parameters_9genes.csv` | **5-Gene**: Multivariable logistic model across $N=1,069$ biopsies: Diagnostic $\text{AUC} = 0.9002$ [95% CI: 0.8747–0.9246], $\text{Brier} = 0.0797$, 0–300 point ruler.<br>**9-Gene**: Diagnostic $\text{AUC} = 0.9238$ [95% CI: 0.9023–0.9438], $\text{Brier} = 0.0705$, 0–450 point ruler. *(In-sample evaluation on discovery matrix)*. |
| **4. Immune Infiltration** | `hub_genes_immune_infiltration_5genes.png`<br>`results/hub_genes_immune_correlations_5genes.csv` | `hub_genes_immune_infiltration_9genes.png`<br>`results/hub_genes_immune_correlations_9genes.csv` | **5-Gene**: $10 \times 5$ Spearman correlation heatmap with BH-FDR + 6 representative scatter plots across $N=124$ liver biopsies (`GSE84044`).<br>**9-Gene**: $10 \times 9$ Spearman correlation heatmap + full $3 \times 3$ grid covering all 9 hub genes vs. validated stromal/immune markers. |
| **5. HPA IHC & RNA Characterization** | `hub_genes_hpa_ihc_summary_5genes.png`<br>`results/hub_genes_hpa_ihc_validation_5genes.csv` | `hub_genes_hpa_ihc_summary_9genes.png`<br>`results/hub_genes_hpa_ihc_validation_9genes.csv` | **5-Gene**: Panel A displays genuine HPA v23 normal-tissue baseline protein IHC across 4 organs (`SERPINE2` explicitly noted as "No IHC Data in HPA"); Panel B displays empirical patient biopsy $\log_2\text{FC}$ across discovery cohorts.<br>**9-Gene**: Panel A displays HPA v23 normal IHC across all 9 hubs; Panel B displays patient transcript-level $\log_2\text{FC}$ across all 4 organs. |
| **6. Functional Enrichment (GO / KEGG / DO)** | `enrichment_5genes.png`<br>`results/go_*_5genes.csv`<br>`results/kegg_5genes.csv`<br>`results/do_5genes.csv` | `enrichment_9genes.png`<br>`results/go_*_9genes.csv`<br>`results/kegg_9genes.csv`<br>`results/do_9genes.csv` | Querying live Enrichr REST API (`gseapy v1.3.1`) across GO BP/MF/CC (2023), KEGG (2021), and DisGeNET Disease Ontology. Top pathways: Extracellular Matrix Organization, Collagen Fibril Organization, ECM-receptor interaction. |
