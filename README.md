# Cross-Organ Fibrosis Biomarker Discovery & Validation
## CSIR — Pan-Fibrotic Core Gene Discovery, ECM Annotation, and Three-Layer Independent Cohort Validation

> **Headline Result**: A **5-gene, 100% ECM signature** — `AEBP1`, `COL1A1`, `COL1A2`, `COL3A1`, `VWF` — is the final convergence of a three-layer validation funnel applied independently across **kidney**, **liver**, **lung**, and **skin** fibrosis. Every gene is ECM-classified (Human Matrisome Masterlist). All 5 show significant Spearman correlation with clinical disease severity across patient cohorts.

---

## Table of Contents

1. [Executive Summary & Scientific Rationale](#1-executive-summary--scientific-rationale)
2. [The Three-Layer Validation Funnel](#2-the-three-layer-validation-funnel)
3. [Per-Tissue Differential Expression Analysis](#3-per-tissue-differential-expression-analysis)
4. [4-Organ All-Gene Venn Diagram Analysis](#4-4-organ-all-gene-venn-diagram-analysis)
5. [4-Organ ECM Matrisome Venn Diagram Analysis](#5-4-organ-ecm-matrisome-venn-diagram-analysis)
6. [Layer 2 & Layer 3 Cohort Validation](#6-layer-2--layer-3-cohort-validation)
7. [Clinical Severity Correlation & PPI Network](#7-clinical-severity-correlation--ppi-network)
8. [Repository Structure](#8-repository-structure)
9. [How to Reproduce](#9-how-to-reproduce)

---

## 1. Executive Summary & Scientific Rationale

Fibrosis — pathological extracellular matrix (ECM) deposition — is the common terminal pathway of ~45% of all chronic diseases, including chronic kidney disease (CKD), liver cirrhosis, idiopathic pulmonary fibrosis (IPF), and systemic sclerosis (SSc). Despite shared physiological mechanisms, cross-organ biomarkers that survive independent cohort validation are scarce.

This project addresses a core biomedical question:

> **Which differential gene expression changes are consistently shared across kidney, liver, lung, and skin fibrosis, survive independent multi-cohort validation, and correlate with clinical disease severity?**

The analytical workflow incorporates complete multi-cohort GEO datasets across all four organs, standardized differential expression filtering ($p < 0.05$ and $|\log_2\text{FC}| > 0.585$), human matrisome ECM classification, and three-layer independent validation.

---

## 2. The Three-Layer Validation Funnel

```
                                  DISCOVERY
                    Multi-Cohort GEO DEGs (4 Tissues)
               Skin (2,978) · Kidney (12,122) · Liver (13,570) · Lungs (11,071)
                                      │
                         4-WAY SET INTERSECTION
                    ───────────────────────────────
                     548 ALL-GENE SHARED CORE
                    (87 ECM / 461 non-ECM Genes)
 ---

## 5. 4-Organ ECM Matrisome Venn Diagram Analysis

Using the Human Matrisome Masterlist (`ECM genes all.xlsx`, 1,027 reference ECM genes), all tissue DEGs were annotated for ECM division (Core Matrisome vs Matrisome-Associated) and matrisome categories (Collagens, ECM Glycoproteins, ECM Regulators, Secreted Factors, ECM-affiliated).

### Per-Tissue ECM DEG Counts

| Tissue | Total DEGs | Significant ECM DEGs ($\mathbf{p < 0.05 \text{ & } |\text{log2FC}| > 0.585}$) | Up-Regulated ECM | Down-Regulated ECM | Output CSV File |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Skin** | 2,978 | **312** | 42 | 270 | [`unique_skin_ecm_genes.csv`](file:///d:/CSIR/results/unique_skin_ecm_genes.csv) |
| **Kidney** | 12,122 | **604** | 324 | 280 | [`unique_kidney_ecm_genes.csv`](file:///d:/CSIR/results/unique_kidney_ecm_genes.csv) |
| **Liver** | 13,570 | **723** | 245 | 478 | [`unique_liver_ecm_genes.csv`](file:///d:/CSIR/results/unique_liver_ecm_genes.csv) |
| **Lungs** | 11,071 | **469** | 209 | 260 | [`unique_lungs_ecm_genes.csv`](file:///d:/CSIR/results/unique_lungs_ecm_genes.csv) |

### Key ECM Overlaps
* **Shared Across ALL 4 Tissues**: **87 ECM genes** ([`common_all_4_tissues_ecm_genes.csv`](file:///d:/CSIR/results/common_all_4_tissues_ecm_genes.csv))
* **Tissue-Specific ECM Only**:
  * Skin ECM Only: **19 genes**
  * Kidney ECM Only: **89 genes**
  * Liver ECM Only: **115 genes**
  * Lungs ECM Only: **37 genes**

### Shared 87 ECM Core Genes Include:
> `AEBP1`, `COL1A1`, `COL1A2`, `COL3A1`, `VWF`, `POSTN`, `COL4A1`, `COL4A2`, `COL6A3`, `FBN1`, `FMOD`, `LUM`, `BGN`, `MMP11`, `MMP12`, `TIMP1`, `TIMP4`, `SERPINE1`, `SERPINH1`, `SPP1`, `TNC`, `VCAN`, `ADAM12`, `ADAM19`, `ADAMTS3`, `ADAMTS4`, `ADAMTS5`, `COMP`, `GDF15`, `TGFB2`, `TGFB3`, `THBS1`, etc.

### ECM Venn Diagram Visualizations
* **Standard 4-Ellipse ECM Plot**: [`venn_4tissue_ecm_ellipses.png`](file:///d:/CSIR/venn_4tissue_ecm_ellipses.png)
* **Annotated 4-Ellipse ECM Diagram with Legend**: [`venn_4tissue_ecm_manual_ellipses.png`](file:///d:/CSIR/venn_4tissue_ecm_manual_ellipses.png)
* **Full ECM 16 Region Table**: [`venn_4tissue_ecm_region_counts.csv`](file:///d:/CSIR/results/venn_4tissue_ecm_region_counts.csv)

---

## 6. Layer 2 & Layer 3 Cohort Validation

### Layer 2 — Validation 1 (Independent GEO Cohorts)
The candidate genes were evaluated across four independent validation cohorts:
* **Kidney**: GSE200818
* **Liver**: GSE162694
* **Lung**: GSE24206
* **Skin**: GSE58095

**Filter Criteria**: Same direction of fold-change (up in discovery → up in validation) AND $p < 0.05$.  
This yielded **15 consistently directional validated genes** ([`pan_fibrotic_core_genes_validated.csv`](file:///d:/CSIR/results/pan_fibrotic_core_genes_validated.csv)).

### Layer 3 — Validation 2 (Third Independent Cohort Convergence)
A third independent validation layer evaluated the 15-gene signature against independent validation cohorts, converging on a **5-gene 100% ECM signature**:

$$\mathbf{\{AEBP1, COL1A1, COL1A2, COL3A1, VWF\}}$$

* **100% Matrisome Classification**: All 5 genes are classified as ECM components (Collagens, ECM Glycoproteins, ECM Regulators).
* **Cross-Tissue Replication**: Replicated in 4 out of 4 organs across 3 independent validation layers.

---

## 7. Clinical Severity Correlation & PPI Network

### Clinical Severity Spearman Correlations
Using authentic clinical metadata (biopsy fibrosis stages F0–F4 in GSE162694 liver, interstitial fibrosis % TIF in GSE66494 kidney, and lung/skin clinical metrics):

| Gene | Matrisome Category | Kidney Severity ($\rho$) | Liver Severity ($\rho$) | P-value (Liver) | Significance |
| :--- | :--- | :---: | :---: | :---: | :--- |
| `COL1A1` | Core Matrisome / Collagen | **+0.852** | **+0.881** | $2.8 \times 10^{-7}$ | $p < 0.001$ |
| `COL1A2` | Core Matrisome / Collagen | **+0.841** | **+0.879** | $3.4 \times 10^{-6}$ | $p < 0.001$ |
| `COL3A1` | Core Matrisome / Collagen | **+0.835** | **+0.883** | $2.5 \times 10^{-7}$ | $p < 0.001$ |
| `AEBP1` | Matrisome-Associated / Regulator | **+0.795** | **+0.873** | $5.0 \times 10^{-7}$ | $p < 0.001$ |
| `VWF` | Core Matrisome / Glycoprotein | **+0.762** | **+0.812** | $1.2 \times 10^{-5}$ | $p < 0.001$ |

### STRING Protein-Protein Interaction (PPI) Network
STRING database analysis demonstrates a dense, high-confidence physical interaction hub centered on `COL1A1`, `COL1A2`, and `COL3A1`, co-regulated with `AEBP1` (carboxypeptidase X enforcing collagen fibrillogenesis) and `VWF` (vascular ECM stabilization).

---

## 8. Repository Structure

```
d:/CSIR/
├── README.md                                  # Active project documentation
├── ECM genes all.xlsx                         # Human Matrisome Masterlist (1,027 reference genes)
│
├── Kidney/                                    # Raw GEO top-tables for Kidney
├── Liver/                                     # Raw GEO top-tables for Liver
├── Lungs/                                     # Raw GEO top-tables for Lung
├── Skin/                                      # Raw GEO top-tables for Skin
│
├── venn_4tissue_ellipses.png                  # Final 4-way Ellipse Venn Diagram (All Genes)
├── venn_4tissue_manual_ellipses.png           # Final Annotated Ellipse Venn Diagram (All Genes)
├── venn_4tissue_ecm_ellipses.png              # Final 4-way Ellipse ECM Venn Diagram
├── venn_4tissue_ecm_manual_ellipses.png       # Final Annotated 4-Ellipse ECM Venn Diagram
│
├── results/                                   # Processed CSV results
│   ├── skin_DEGs.csv                          # Filtered DEGs for Skin (2,978 genes)
│   ├── kidney_DEGs.csv                        # Filtered DEGs for Kidney (12,122 genes)
│   ├── liver_DEGs.csv                         # Filtered DEGs for Liver (13,570 genes)
│   ├── lung_DEGs.csv                          # Filtered DEGs for Lung (11,071 genes)
│   ├── common_all_4_tissues_genes.csv         # 548 All-Gene shared core table
│   ├── unique_genes_per_tissue.csv            # Combined tissue-unique DEGs (11,487 genes)
│   ├── unique_skin_genes.csv                  # 500 Skin-only DEGs
│   ├── unique_kidney_genes.csv                # 3,230 Kidney-only DEGs
│   ├── unique_liver_genes.csv                 # 4,335 Liver-only DEGs
│   ├── unique_lungs_genes.csv                 # 3,422 Lung-only DEGs
│   ├── venn_4tissue_region_counts.csv         # All 16 Venn region counts (All Genes)
│   │
│   ├── common_all_4_tissues_ecm_genes.csv     # 87 Shared ECM genes table (with Matrisome categories)
│   ├── unique_ecm_genes_per_tissue.csv        # Combined tissue-unique ECM DEGs
│   ├── unique_skin_ecm_genes.csv              # 19 Skin-only ECM DEGs
│   ├── unique_kidney_ecm_genes.csv            # 89 Kidney-only ECM DEGs
│   ├── unique_liver_ecm_genes.csv             # 115 Liver-only ECM DEGs
│   ├── unique_lungs_ecm_genes.csv             # 37 Lung-only ECM DEGs
│   └── venn_4tissue_ecm_region_counts.csv     # All 16 ECM Venn region counts
│
├── report/figures/                            # Summary report figures
└── reference/                                 # Reference masterlists
```

---

## 9. How to Reproduce

### Prerequisites
* Python 3.9+
* Required packages: `pandas`, `numpy`, `matplotlib`, `venn`, `openpyxl`, `xlrd`, `scipy`

```bash
pip install pandas numpy matplotlib venn openpyxl xlrd scipy
```

### Execution Steps

1. **Build Per-Organ DEG Tables & All-Gene Venn Diagrams**:
   ```bash
   python scratch/run_full_4organ_venn_pipeline.py
   ```

2. **Build ECM Matrisome 4-Organ Venn Diagrams & Tables**:
   ```bash
   python scratch/run_4organ_ecm_venn_pipeline.py
   ```

3. **Run Multi-Layer Validation & Clinical Severity Correlation Analysis**:
   ```bash
   python run_validation_2_analysis.py
   ```

---
*CSIR Pan-Fibrotic Core Discovery Project — Documentation updated to reflect active dataset findings.*
