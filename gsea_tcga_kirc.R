library("clusterProfiler")
library("pathview")
library("enrichplot")

######### STEP 1: importaing annotation table

# home sapiens annotation database

organism = "org.Hs.eg.db"

# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
 
############ STEP 2: prepare gene list

df_tcga_kirc <- readRDS("dge_tcga_kirc_results.rds")

# all genes are included, even the non-significant ones

# we want the log2 fold change 
original_gene_list <- df_tcga_kirc$log2FoldChange

# name the vector with gene names
names(original_gene_list) <- df_tcga_kirc$ensid

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
# we get then the genes (gene names) and log2Fold changes sorted from highest (positive values) to lowest (negative values)
gene_list = sort(gene_list, decreasing = TRUE)

########## STEP 3: exploring the annotation file

# getting to know the annotation file

keytypes(org.Hs.eg.db)

keys <- keys(org.Hs.eg.db)

select(org.Hs.eg.db, keys=keys, columns=c("ENSEMBL", "SYMBOL","GENETYPE"))

########## STEP 4: performing gene set enrichment analysis

# creating gse object

## documentation

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

gseGO_results <- gse@result

################################################################################
# dotplot
################################################################################

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

################################################################################
# enrichment map
################################################################################

pairwise_termism_gse <- pairwise_termsim(gse, showCategory = 100)

emapplot(pairwise_termism_gse)

################################################################################
# netplot
################################################################################

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 5, layout="star")
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3, layout="fr")

################################################################################
# ridgeplot
################################################################################

ridgeplot(gse)

################################################################################
# GSEA / running enrichment plot
################################################################################

gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1)

gseaplot(gse, by = "all", title = gse@result$Description[5], geneSetID = 5)

gseaplot(gse, by = "all", title = gse@result$Description[12], geneSetID = 12)

gseaplot(gse, by = "all", title = gse@result$Description[40], geneSetID = 40)

gse@result$Description[1:50]
