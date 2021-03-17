
# library prep ------------------------------------------------------------
library(shiny)
library(shinycssloaders)
library(cowplot)
library(Seurat)
library(ggplot2)

# data prep ---------------------------------------------------------------
# file_name <- "~/OneDrive - 숭실대학교 - Soongsil University/BI/KU/Yale/HSC_scRNA_seq/Version1/20210225_c04mer/20210302_cxcr4_HR.rds"
cxcr4 <- readRDS('./R/20210302_cxcr4_HR.rds')

cxcr4_meta <- cxcr4
cxcr4_exp_data <- cxcr4@assays$MAGIC@data

# Shiny codes
ui <- fluidPage(
  # title
  titlePanel("Cell Browser"),

  # line break
  hr(),

  sidebarLayout(
    sidebarPanel(
      # select the condition
      radioButtons(
        inputId = "input_condition",
        label = "cells:",
        choices = c("of all conditions" = "all",
                    "split be biological conditions" = "split"),
      ),
      # select the gene
      selectInput(
        inputId = "input_gene",
        label = "Gene:",
        choices = c(Choose = "", rownames(cxcr4_exp_data)),
        selected = "Cxcl12",
        selectize = TRUE
      )
    ),

    # mainPanel for visualize plots
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "UMAP plot",
          fluidRow(
            column(
              width = 6,
              withSpinner(plotOutput("UMAP_feature_plot"))
            ),
            column(
              width = 6,
              withSpinner(plotOutput("UMAP_dimplot"))
            )
          )
        ),
        tabPanel(
          title = "Violin Plot",
          withSpinner(plotOutput("Violin_plot"))
        )
      )
    )
  )
)
server <- function(input, output){
  # generate featureplot
  featurePlot_reactive <- reactive({
    req(input$input_gene)
    feature_plot <- FeaturePlot(cxcr4, features = input$input_gene, cols = c("light gray", "red"))+
      labs(title =paste0("Gene Expression of ", input$input_gene))+
      theme(aspect.ratio = 1)
    if(input$input_condition=="split"){
      feature_plot <-FeaturePlot(cxcr4, features = input$input_gene, cols = c("light gray", "red"), split.by = "orig.ident")+
        theme(aspect.ratio = 1)
    }
    feature_plot
  })
  dimplot_reactive <- reactive({
    DimPlot(cxcr4) +
      theme(aspect.ratio = 1)
  })
  vlnplot_reactive <- reactive({
    req(input$input_gene)
    vlnplot <- VlnPlot(cxcr4, features = input$input_gene)
    if(input$input_condition=="split"){
      vlnplot <-VlnPlot(cxcr4, features = input$input_gene, split.by = "orig.ident", split.plot = T)+
        theme(aspect.ratio = 1)
    }
    vlnplot
  })
  output$UMAP_feature_plot <- renderPlot(featurePlot_reactive())
  output$UMAP_dimplot <- renderPlot(dimplot_reactive())
  output$Violin_plot <- renderPlot(vlnplot_reactive())
}
shinyApp(ui=ui,server=server)
