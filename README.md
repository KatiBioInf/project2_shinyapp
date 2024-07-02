# project2_shinyapp
This repository contains the code and data for the shiny app that presents the results of the differential gene expression and gene set enrichment analysis of 6 samples (3 normal and 3 tumor) from the [TCGA-KIRC](https://portal.gdc.cancer.gov/projects/TCGA-KIRC) (Kidney Renal Clear Cell Carcinoma) project.

The app is available here.

The repository includes the following files:

1. data_for_shiny_app.RData: all data necessary to run the app
2. app.R: code for shiny app
3. gsea_tcga_kirc.R: underlying gene set enrichment analysis using [clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html) and [enrichplot](https://www.bioconductor.org/packages/release/bioc/html/enrichplot.html)
