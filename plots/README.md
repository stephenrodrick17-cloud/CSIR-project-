# Canonical Publication Figures Directory (`plots/`)

This directory contains the final, publication-grade figures for the **CSIR Pan-Fibrotic Multi-Organ Biomarker Study**. Every figure here has been verified against the locked dataset and source files.

---

### Core Numbers Enforced Across All Figures
- **Discovery Cohorts**: $N = 14$ independent cohorts across 4 human organs (Kidney: 4, Liver: 3, Lung: 4, Skin: 3).
- **Discovery Sample Size**: $N = 1,069$ human patient tissue samples (170 non-fibrotic controls, 899 fibrotic cases).
- **Core Significance Filter**: eBayes moderated $t$-test $|\log_2\text{FC}| \ge 0.585$ (1.5-fold), Benjamini-Hochberg FDR $p < 0.05$.
- **Pan-Fibrotic Conserved Core**: **86 DEGs** shared across all 4 organs.
- **Clean Core ECM Program**: **24 genes** overlapping the Human Matrisome Masterlist ($N=1,027$).
- **Machine Learning & WGCNA Hub Discovery**: 5-method consensus on $N=1,069$ pure discovery samples (LASSO, SVM-RFE, RF, XGBoost, TOM Co-expression Clustering / Python WGCNA) -> **9 Multi-Model Consensus Hub Genes** ($\ge 3/5$ votes): `COL15A1`, `COL1A1`, `COL3A1`, `SERPINE2`, `SERPINF2`, `LAMC3`, `LTBP2`, `MDK`, `SVEP1`.
- **Unanimous Consensus Hubs (5/5 Votes)**: `COL15A1`, `COL1A1`, `COL3A1`, `SERPINE2`, `SERPINF2` (selected unanimously by all 5 feature selection architectures).
- **Cross-Platform Validation (Microarray vs. RNA-seq)**: Under strict dual-platform criteria ($p < 0.05$, AUC $> 0.76$, concordant direction on both platforms), exactly **3 Hub Genes strictly pass**: `COL15A1`, `COL3A1`, and `SERPINE2`. `COL1A1` and `SERPINF2` are 100% direction-concordant across platforms but non-significant in RNA-seq ($p = 0.216$ and $p = 0.788$).

---

### Canonical Figure Index: Upstream Discovery & Conserved Core

| File | Pipeline Stage | Description | Source Script / Input Data |
| :--- | :--- | :--- | :--- |
| `venn_4organ_manual_ellipses.png` | **Stage 1: Discovery Intersection** | 4-Organ Venn diagram showing DEGs for Kidney (11,758), Liver (2,268), Lung (7,803), Skin (2,979), and the conserved 86-gene core intersection. | `venn_4organs.py` |
| `upset_plot_4organs_corrected.png` | **Stage 1: Discovery Intersection** | UpSet plot showing all 15 subset intersections across the 4 organs; mathematically and cardinality identical to the 4-organ Venn. | `generate_upset_plot.py` |
| `venn_organ_vs_ecm_compendium.png` | **Stage 1: ECM Filtering** | Multi-panel Venn compendium: 2x2 grid of organ DEGs vs. ECM Masterlist (1,027 genes), 3-way ECM intersection, and 4-way pan-fibrotic core (86 DEGs -> 24 ECM genes). | `venn_organ_vs_ecm.py` |

---

### Stage 4 Parallel Downstream Publication Suites: 5-Gene Unanimous Hubs vs. 9-Gene Multi-Model Hubs

Both suites are fully generated, verified, and preserved side-by-side using 100% real human patient cohorts and public databases.

| Stage 4 Analysis | 5-Gene Unanimous Consensus Suite | 9-Gene Multi-Model Consensus Suite | Key Differences & Characteristics |
| :--- | :--- | :--- | :--- |
| **1. PPI Network** | `hub_genes_ppi_network_5genes.png`<br>`results/hub_genes_ppi_network_nodes_5genes.csv`<br>`results/hub_genes_ppi_network_edges_5genes.csv` | `hub_genes_ppi_network_9genes.png`<br>`results/hub_genes_ppi_network_nodes_9genes.csv`<br>`results/hub_genes_ppi_network_edges_9genes.csv` | **5-Gene**: 17 nodes, 68 edges ($p < 10^{-16}$), focused on core collagens (`COL1A1`, `COL3A1`, `COL15A1`) and serpins (`SERPINE2`, `SERPINF2`).<br>**9-Gene**: Exactly 9 nodes (100% pure consensus hubs: `COL15A1`, `COL1A1`, `COL3A1`, `SERPINE2`, `SERPINF2`, `LAMC3`, `LTBP2`, `MDK`, `SVEP1`) with 3 high-confidence direct STRING interactions ($\ge 0.700$, the core collagen triad: `COL1A1`, `COL3A1`, `COL15A1`); zero bridging proteins or threshold dilution, displaying the 6 non-collagen hubs as isolated nodes as-is. |
| **2. Early Detection Triptych** | `hub_genes_early_detection_triptych_5genes.png`<br>`results/hub_genes_early_stage_roc_metrics_5genes.csv` | `hub_genes_early_detection_triptych_9genes.png`<br>`results/hub_genes_early_stage_roc_metrics_9genes.csv` | **5-Gene**: Early ROC (S0 vs S1-S2, $N=96$) ensemble $\text{AUC} = 0.716$; Net benefit $+0.396$ at $p_t=25\%$.<br>**9-Gene**: Early ROC ensemble $\text{AUC} = 0.719$; Net benefit $+0.383$ at $p_t=25\%$. |
| **3. Diagnostic Nomogram & DCA** | `hub_genes_nomogram_and_dca_5genes.png`<br>`results/hub_genes_nomogram_parameters_5genes.csv` | `hub_genes_nomogram_and_dca_9genes.png`<br>`results/hub_genes_nomogram_parameters_9genes.csv` | **5-Gene**: Multivariable logistic model across 1,069 biopsies has Diagnostic $\text{AUC} = 0.9002$ [95% CI: 0.8747–0.9246], $\text{Brier} = 0.0797$, 0–300 point ruler.<br>**9-Gene**: Diagnostic $\text{AUC} = 0.9238$ [95% CI: 0.9023–0.9438], $\text{Brier} = 0.0705$, 0–450 point ruler. *(Note: Evaluated in-sample on N=1,069; 9-gene gain expected from adding parameters).* |
| **4. Immune Infiltration** | `hub_genes_immune_infiltration_5genes.png`<br>`results/hub_genes_immune_correlations_5genes.csv` | `hub_genes_immune_infiltration_9genes.png`<br>`results/hub_genes_immune_correlations_9genes.csv` | **5-Gene**: $10 \times 5$ Spearman correlation heatmap with Benjamini-Hochberg FDR + 6 regressions across $N=124$ liver biopsies (`GSE84044`).<br>**9-Gene**: $10 \times 9$ Spearman correlation heatmap + 6 representative regressions across $N=124$ biopsies. |
| **5. HPA IHC Baseline & Patient RNA Characterization** | `hub_genes_hpa_ihc_summary_5genes.png`<br>`results/hub_genes_hpa_ihc_validation_5genes.csv` | `hub_genes_hpa_ihc_summary_9genes.png`<br>`results/hub_genes_hpa_ihc_validation_9genes.csv` | **5-Gene**: Panel A displays genuine HPA v23 normal-tissue baseline protein IHC across 4 organs (`SERPINE2` marked "No IHC Data in HPA"); Panel B displays empirical patient biopsy transcript-level $\log_2\text{FC}$ across $N=1,069$ discovery samples *(Note: $\log_2\text{FC}$ values reflect Stage 1 organ-specific discovery biopsies, distinct from Stage 3 pooled cross-platform validation)*.<br>**9-Gene**: Panel A displays genuine HPA v23 normal-tissue baseline protein IHC across all 9 consensus hubs using verified monospecific antibodies; Panel B displays empirical patient biopsy transcript-level $\log_2\text{FC}$ across all 4 organs. |


