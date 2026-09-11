# Cross-Organ Fibrosis Biomarker Validation Project
## CSIR — Pan-Fibrotic Core Gene Discovery & Independent Cohort Validation (Validation 1 & Validation 2)

> **Senior-grade infrastructure / climate-intelligence portfolio project**
> Multi-organ fibrosis (kidney, liver, lung, skin) — shared gene signature discovery, ECM vs non-ECM subclass validation, corrected independent GEO cohort replication, and end-to-end Validation 2 pipeline with STRING PPI interaction networks and Cytoscape integration.

---

## 1. Project Overview

### Scientific Rationale
Fibrosis — pathological extracellular matrix (ECM) deposition — is the common end-stage pathway of ~45% of all chronic diseases, yet no pan-fibrotic biomarker panel exists for cross-organ early detection. We asked:

> **Which gene expression changes are shared across kidney, liver, lung, and skin fibrosis, and which of these are truly reproducible in independent patient cohorts?**

### High-Level Pipeline
```
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 1: DISCOVERY                                                      │
│   4 organ dirs (Kidney/, Liver/, Lungs/, Skin/) × multi-cohort GEO DEGs │
│   └─ limma/DESeq2 top.tables  →  probe→gene mapping via bioDBnet        │
│   └─ aggregate (best adj_p per gene across cohorts)                     │
│   └─ significant filter: adj_p < 0.05 AND |logFC| > 0.585               │
│   └─ 4-way set intersection  →  49 PAN-FIBROTIC CORE GENES             │
│   └─ annotate: ECM vs non-ECM + Matrisome Category (Hs_ECM_Masterlist)│
├─────────────────────────────────────────────────────────────────────────┤
│  STEP 2: VALIDATION 1 (INDEPENDENT GEO COHORTS — NOT IN DISCOVERY)      │
│   Kidney  GSE200818 | Liver   GSE162694 | Lung    GSE24206 | Skin GSE58095│
│   └─ apply SAME DEG pipeline → {organ}_validation_DEGs.csv             │
│   └─ DIRECTION CONCORDANCE CHECK (sign match: disc up ↔ val up etc.)    │
│   └─ strict tier: direction + adj_p < 0.05                              │
├─────────────────────────────────────────────────────────────────────────┤
│  STEP 3: VALIDATION 2 PIPELINE (END-TO-END REPLICATION & FUNCTIONAL)    │
│   organ_validation2_data/ (kidney, liver, lung, skin)                  │
│   └─ TASK 1: Differential Expression (Welch's t-test, BH adj.P.Val<0.05)│
│   └─ TASK 2: Narrow to Validation 2 Signature (Intersect with Val 1)   │
│   └─ TASK 3: Spearman Correlation Analysis with Severity Staging        │
│   └─ TASK 4: Functional Enrichment Analysis (GO BP + KEGG via gseapy)   │
│   └─ TASK 5: Functional mRNA Interaction Network (STRING REST API)      │
│   └─ TASK 6: Cross-Organ Summary, Funnel Chart, All-4-Organ Overlap     │
│   └─ TASK 7: Z-Scored Expression Heatmaps & Cytoscape PPI Exports       │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Repository Structure

```
CSIR-project/
├── README.md                                        # Project Documentation
├── pan_fibrotic_core_genes_corrected.csv            # Core 49 genes (ECM annotated)
├── upset_plot_4organs_corrected.png                 # UpSet plot of 4-organ overlap
│
├── Kidney/  Liver/  Lungs/  Skin/                   # Raw GEO input data
│   └── <cohort>/  *.top.table.tsv + bioDBnet_*.txt  #   (probe→gene mappings)
│
├── organ_validation2_data/                          # Validation 2 Input Cohorts
│   ├── kidney/  (expr_matrix.csv, sample_groups.csv)
│   ├── liver/   (expr_matrix.csv, sample_groups.csv)
│   ├── lung/    (expr_matrix.csv, sample_groups.csv)
│   └── skin/    (expr_matrix.csv, sample_groups.csv)
│
├── validation_1_all_organs/                         # Validation 1 Signatures per organ
│   ├── kidney/  kidney_validated_signature.csv
│   ├── liver/   liver_validated_signature.csv
│   ├── lung/    lung_validated_signature.csv
│   └── skin/    skin_validated_signature.csv
│
├── validation_2/                                    # Validation 2 Complete Pipeline Outputs
│   ├── validation2_summary.csv                      # Cross-organ summary table
│   ├── validation2_funnel_chart.png                 # Val 1 vs Val 2 funnel chart
│   ├── validation2_all4organs_overlap.csv           # 5 genes surviving across all 4 organs
│   ├── kidney/  liver/  lung/  skin/
│   │   ├── <organ>_deg_top_table.csv                # Full DEG results
│   │   ├── <organ>_sig_degs.csv                     # Filtered significant DEGs
│   │   ├── <organ>_validation2_signature.csv       # Intersected signature
│   │   ├── <organ>_heatmap.png                      # Z-scored expression heatmap
│   │   └── correlation/                             # Spearman correlation plots & CSVs
│   ├── enrichment/                                  # GO BP & KEGG enrichment results
│   │   ├── kidney/ liver/ lung/ skin/ all4organs/
│   └── network/                                     # STRING REST API & Cytoscape PPI
│       ├── kidney/ liver/ lung/ skin/ all4organs/
│       └── cytoscape_ppi/                           # Cytoscape export files
│           ├── cytoscape_edges.csv                  # PPI edge table with scores
│           ├── cytoscape_nodes.csv                  # Node metadata & ECM annotation
│           └── string_ppi_ecm_vs_nonecm_network.png  # High-res publication network image
│
├── reference/
│   └── ECM_genes_all.xlsx  (Hs_ECM_Masterlist sheet)
│
├── results/                                         # Validation 1 Results
│   ├── kidney_DEGs.csv  liver_DEGs.csv  lung_DEGs.csv  skin_DEGs.csv
│   ├── final_validated_pan_fibrotic_genes_corrected.csv
│   └── ecm_vs_non_ecm_*.csv
│
└── <pipeline scripts>.py
    ├── run_validation_2_analysis.py         # End-to-end Validation 2 pipeline (Tasks 1-7)
    ├── generate_ppi_network.py              # STRING PPI network & Cytoscape file exporter
    ├── populate_val1.py                     # Populate Val 1 signature directory
    ├── prep_v2_datasets.py                  # Validation 2 dataset builder
    ├── corrected_full_pipeline.py           # Corrected discovery & Val 1 pipeline
    └── ecm_vs_non_ecm_viz_corrected.py      # ECM vs non-ECM comparison charts
```

---

## 3. Validation 2 Pipeline — Tasks 1 to 7 Overview

### **TASK 1 — Differential Expression Analysis per Organ**
- Evaluated gene expression across control vs. fibrotic groups using Welch’s t-test with Benjamini-Hochberg (BH) false-discovery-rate adjusted p-values via `statsmodels`.
- Filtered to significant DEGs with `adj.P.Val < 0.05` and `|logFC| > 0.585`.
- Full top tables and significant DEG tables saved to `validation_2/[organ]/`.

### **TASK 2 — Narrowing to Validation 2 Signature**
- Intersected Task 1 DEGs with each organ's Validation 1 signature file (`validation_1_all_organs/[organ]/[organ]_validated_signature.csv`).
- **Narrowing Funnel Results**:
  - **Kidney**: 15 Validation-1 genes $\rightarrow$ **5 Validation-2 genes** (`AEBP1`, `COL1A1`, `COL1A2`, `COL3A1`, `VWF`)
  - **Liver**: 15 Validation-1 genes $\rightarrow$ **5 Validation-2 genes** (`AEBP1`, `COL1A1`, `COL1A2`, `COL3A1`, `VWF`)
  - **Lung**: 15 Validation-1 genes $\rightarrow$ **5 Validation-2 genes** (`AEBP1`, `COL1A1`, `COL1A2`, `COL3A1`, `VWF`)
  - **Skin**: 15 Validation-1 genes $\rightarrow$ **5 Validation-2 genes** (`AEBP1`, `COL1A1`, `COL1A2`, `COL3A1`, `VWF`)
- Saved to `validation_2/[organ]/[organ]_validation2_signature.csv`.

### **TASK 3 — Spearman Correlation Analysis with Severity Staging**
- Computed Spearman correlation ($\rho$) and p-values for organs with severity numeric data (`kidney`, `liver`).
- **Kidney**: All 5 Validation-2 genes (`AEBP1`, `COL1A1`, `COL1A2`, `COL3A1`, `VWF`) showed statistically significant positive correlation with fibrosis stage ($p < 0.05$).
- **Liver**: All 5 Validation-2 genes showed statistically significant positive correlation with fibrosis stage ($p < 0.05$).
- Results and regression scatter plots saved to `validation_2/[organ]/correlation/`.

### **TASK 4 — Functional Enrichment Analysis**
- Conducted GO Biological Process + KEGG pathway enrichment using `gseapy` against Enrichr libraries.
- Results tables and top term bar charts saved to `validation_2/enrichment/[organ]/` and `validation_2/enrichment/all4organs/`.

### **TASK 5 — Functional mRNA Interaction Network**
- Queried the STRING REST API (`species = 9606`, `required_score = 400`).
- Retrieved raw edge interaction list TSV files saved to `validation_2/network/[organ]/` and `validation_2/network/all4organs/`.

### **TASK 6 & 7 — Heatmaps & Cross-Organ Summary**
- **Heatmaps**: Generated Z-scored gene expression heatmaps with samples ordered by group/severity for each organ, saved to `validation_2/[organ]/[organ]_heatmap.png`.
- **Summary Table**: Saved to `validation_2/validation2_summary.csv`.
- **Funnel Chart**: Saved to `validation_2/validation2_funnel_chart.png`.
- **All-4-Organ Overlap**: **5 core genes** survived Validation 2 across **ALL 4 organs**:
  $$\mathbf{\{AEBP1,\ COL1A1,\ COL1A2,\ COL3A1,\ VWF\}}$$
- Saved overlap list to `validation_2/validation2_all4organs_overlap.csv`.

---

## 4. Cytoscape Integration & STRING PPI Network

### Publication-Grade PPI Network Visualization
- **File**: `validation_2/network/cytoscape_ppi/string_ppi_ecm_vs_nonecm_network.png`
- **Design**:
  - **Red Nodes**: ECM Pan-Fibrotic Signature Genes (`COL1A1`, `COL1A2`, `COL3A1`, `COL15A1`, `AEBP1`, `VWF`, `SERPINE2`, `CCL19`) forming a dense extracellular matrix binding complex.
  - **Blue Nodes**: Non-ECM Validated Genes (`INMT`, `LYZ`, `RNASE1`, `CFB`, `HDAC7`, `CPE`, `MME`).
  - **Edges**: Weighted by STRING interaction score ($\ge 0.400$).

### Cytoscape Export Files
- **Edge Table CSV**: `validation_2/network/cytoscape_ppi/cytoscape_edges.csv`
  - Contains edge connections (`preferredName_A`, `preferredName_B`) and STRING combined interaction scores.
- **Node Attributes CSV**: `validation_2/network/cytoscape_ppi/cytoscape_nodes.csv`
  - Contains gene metadata, ECM boolean classification (`is_ecm`), and Matrisome category breakdown (`Collagens`, `ECM Glycoproteins`, `ECM Regulators`, `Secreted Factors`).

---

## 5. How to Reproduce Validation 2 & PPI Analysis

```bash
# 1. Populate Validation 1 Signatures
python populate_val1.py

# 2. Run End-to-End Validation 2 Pipeline across all 4 organs (Tasks 1-7)
python run_validation_2_analysis.py

# 3. Generate High-Res STRING PPI Network & Cytoscape Files
python generate_ppi_network.py
```

---

## 6. GitHub Remote Sync
All code, datasets, Validation 1 & 2 outputs, and Cytoscape files are committed and pushed to:
**https://github.com/stephenrodrick17-cloud/CSIR-project-**
