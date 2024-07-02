# Libraries----
library(shiny)
library(shinydashboard)
library(shinyjs)
library(tidyverse)
library(scales)
library(ggplot2)
library(dashboardthemes)
library(shinythemes)
library(clusterProfiler)
library(pathview)
library(enrichplot)

# load data

load("data_for_shinyapp.RData")

# data cleaning

# only top 200 genes in dge

dge_top_200 <- dge_model_results %>% 
  arrange(padj) %>% 
  top_n(-200) %>% 
  dplyr::mutate(keep="yes") %>% 
  select(ensid, keep)

df_dot_chart <- df_dot_chart %>% 
  left_join(dge_top_200, by="ensid") %>% 
  dplyr::filter(keep=="yes") %>% 
  select(-keep)

dge_model_results <- dge_model_results %>% 
  left_join(dge_top_200, by="ensid") %>% 
  dplyr::filter(keep=="yes") %>% 
  select(-keep)
  

# the app

ui <- shinyUI(fluidPage(
  theme = shinytheme("sandstone"),
  
  titlePanel(h1(strong("TCGA-KIRC - DGE and GSEA"))),
  
  sidebarLayout(
    sidebarPanel(
      #conditionalPanel(condition = "input.tab_selected==1", h4("Introduction to app features")),
      conditionalPanel(style = "font-size: 16px", condition = "input.tab_selected==2", selectInput(inputId = "ind_gene", label = "Select gene symbol", choices = sort(unique(dge_model_results$gene_name)), selected = NULL)),
      conditionalPanel(style = "font-size: 16px", condition = "input.tab_selected==3", selectInput(inputId = "GO_id", label = "Select GO ID", choices = sort(unique(gseGO_results$ID)), selected = NULL),
                       ),
    ),
    
    mainPanel(
      
      tabsetPanel(
        type="tabs",
        id="tab_selected",
        selected=1,
        tabPanel(title="Documentation", value = 1,
                 h2("Introduction"),
                 p("This app presents the results of the gene expression and gene set enrichment analysis based on six samples (3 normal and 3 tumor) from the", a(em("TCGA-KIRC"), href="https://portal.gdc.cancer.gov/projects/TCGA-KIRC"), "project.",
                   style = "font-size: 110%"),
                 p("The GitHub repository for this app is available", a(em("here."), href="https://github.com/KatiBioInf/project2_shinyapp"),
                   style = "font-size: 110%"),
                 h2("Tabs"),
                 p("The tab 'differential gene expression' presents the results from the", a(em("DESeq2"), href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html"), "analysis. For details visit", a(em("my other GitHub repository."), href="https://github.com/KatiBioInf/project1_DESeq2"),
                   "The plot shows the differences in the normalized expression counts per sample type for the gene symbol selected by the user in the left-side panel (scroll-down menu). 
                   The differential gene expression results (incl. log2-fold changes and adjusted p-values) for the same selected gene is presented in the table below the plot. Only the top 200 genes based on adjusted p-value are included in this app.",
                   style = "font-size: 110%"),
                 p("The tab 'gene set enrichment' presents the results of the gene set enrichment analysis using the GO database via the", a(em("clusterProfiler"), href="https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html"), "R/Bioconductor package. We present gene sets with p-value of 0.05 or lower.
                   The user can choose the GO id to inspect using the scroll-down menu on the left-side.",
                   style = "font-size: 110%"),
                 icon = icon("question-circle")),
        tabPanel(title="Differential gene expression", 
                 fluidRow(br(),
                          br(),
                          p("Normalized count per sample type", style = "font-size: 130%"),
                          br(),
                          plotOutput("plot_gene", width = 550, height = 450),
                          br(),
                          br(),
                          p("DESeq2 results", style = "font-size: 130%"),
                          br(),
                          tableOutput("table_DGE")), value=2),
        tabPanel(title="Gene set enrichment", 
                 fluidRow(
                   br(),
                   br(),
                   p("Barcode plot", style = "font-size: 130%"),
                   br(),
          plotOutput("plot_GO", width = 500, height = 600), 
          br(),
          br(),
          p("GO analysis", style = "font-size: 130%"),
          br(),
          tableOutput("table_gse")), value=3)
      )
    )
  )
)
)


shinyServer(server <- function(input, output){
  
  # data cleaning
  prepped_data1 <- reactive({
    df_single_gene <- df_dot_chart %>% 
      dplyr::select(c(2:8)) %>% 
      dplyr::filter(gene_name %in% input$ind_gene) %>% 
      dplyr::select(-gene_name) %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column(var="sample") 
    
    
    colnames(df_single_gene) <- c("sample", "norm_count")
    
    df_single_gene$sample <- gsub("\\.", "-", df_single_gene$sample)
    
    df_single_gene <- df_single_gene %>% 
      left_join(sampleinfo, by="sample")
    
  })
  
  
  prepped_data2 <- reactive({
    data2 <- dge_model_results %>% 
      dplyr::filter(gene_name %in% input$ind_gene)
  })
  
  prepped_data3 <- reactive({
    data3 <- gseGO_results %>%
      filter( ID %in% input$GO_id) 
  })
  
  index <- reactive({
    which(gseGO_results$ID %in% input$GO_id)
    
  })
  
  output$plot_gene <- renderPlot({
    ggplot(prepped_data1(), aes(x=types, y=norm_count, fill=types)) + 
      geom_dotplot(binaxis='y', stackdir='center')+
      scale_fill_manual(values=c("#9FF3BD", "#F43911"))+
      theme_bw()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 16))+
      ggtitle(req(input$ind_gene))+
      xlab("sample type")+
      ylab("normalized count")
  })
  
  output$table_DGE <- renderTable({
    prepped_data2()

  })
  
  output$plot_GO <- renderPlot({
    gseaplot(gse, by = "all", title = req(input$GO_id), geneSetID = index())
  })
  
  output$table_gse <- renderTable({
    prepped_data3()
    
  })
  
  

})

shinyApp(ui, server)