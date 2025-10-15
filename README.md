# Single-cell and Molecular Insights into the Immunopathogenic Mechanisms of Oral Mucous Membrane Pemphigoid

This repository contains the code for the for the manuscript:

> **Single-cell and Molecular Insights into the Immunopathogenic Mechanisms of Oral Mucous Membrane Pemphigoid**

---
##  📖 Description
This repository includes scripts for processing single-cell RNA sequencing data, performing various analyses such as clustering, differential expression analysis, and gene set enrichment analysis, as well as visualizing the results through heatmaps and dot plots.

---
## 🧊 Key Components
- Fig1.R: Single-cell profiling and unbiased clustering of cells in buccal tissues of MMP samples
- Fig2.R: Activation of the WNT pathway in MMP samples
- Fig3.R: Leukocyte migration activation triggered by T cell subsets in the MMP Environment
- Fig4.R: Subpopulations of myeloid cells play a pro-inflammatory and pro-fibrotic role in the buccal mucosa of MMP patients
- Fig5.R: Integrates peripheral blood and tissue samples, and use the scRepertoire and UpSetR packages for T-cell receptor (TCR) analysis.

---
## ⚙️ Requirements
R Version
- R (version 4.2 or higher)
- Required R packages:
- Seurat - For single-cell RNA sequencing analysis.
- ggplot2 - For data visualization.
- dplyr - For data manipulation.
- Harmony - For integrating multiple single-cell datasets.
- GSVA - For Gene Set Variation Analysis.
- ComplexHeatmap - For generating heatmaps with annotations.
- circlize - For circular visualizations and complex data visualization.
---

## ⬇️ Requirements
To install the required packages, run the following R code:

- install.packages(c("ggplot2", "dplyr", "circlize", "ComplexHeatmap"))
- install.packages("Seurat")
- devtools::install_github("jokergoo/ComplexHeatmap")
- install.packages("GSVA")
- install.packages("harmony")
- install.packages("CellChat")
---

