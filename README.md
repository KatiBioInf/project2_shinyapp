# project2_shinyapp
This repository contains the code and data for the shiny app that presents the main results of the differential gene expression and gene set enrichment analysis of 6 samples (3 normal and 3 tumor) from the [TCGA-KIRC](https://portal.gdc.cancer.gov/projects/TCGA-KIRC) (Kidney Renal Clear Cell Carcinoma) project.

The app is available [here](https://katibioinf.shinyapps.io/DGE_and_GSEA_shiny_app/).

The repository includes the following files:

1. data_for_shiny_app.RData: all data necessary to run the app
2. app.R: code for shiny app
3. dge_tcga_kirc_results.rds: output data from the dge analysis = input for the gsea
4. gsea_tcga_kirc.R: underlying gene set enrichment analysis using [clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html) and [enrichplot](https://www.bioconductor.org/packages/release/bioc/html/enrichplot.html)

The details of the differential gene expression analysis is available in [my other repository](https://github.com/KatiBioInf/project1_DESeq2).
