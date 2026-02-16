# COVID-19 Lung Atlas Interactive Explorer
# Based on Melms et al. 2021 "A molecular single-cell lung atlas of lethal COVID-19"
# Data: GSE171524

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(dplyr)
library(ggplot2)
library(viridis)

# ============================================================================
# Load Data
# ============================================================================

# UMAP and metadata
umap_data <- read.csv("data/shiny_app_data_umap_sampled.csv", row.names = 1)

# Gene expression
gene_expr <- read.csv("data/shiny_app_data_gene_expression.csv", row.names = 1)

# DE results
de_global <- read.csv("data/shiny_app_data_de_results.csv")
de_celltype <- read.csv("data/shiny_app_data_de_celltype.csv")

# Cell-cell interactions
interactions <- read.csv("data/shiny_app_data_interactions.csv")

# Define colors
cell_type_colors <- c(
  "Epithelial" = "#2ECC71",
  "Myeloid" = "#E74C3C",
  "Fibroblast" = "#9B59B6",
  "Endothelial" = "#F39C12",
  "T_cell" = "#3498DB",
  "B_cell" = "#1ABC9C",
  "NK_cell" = "#E91E63",
  "Plasma_cell" = "#00BCD4",
  "Mast_cell" = "#795548"
)

condition_colors <- c("COVID" = "#E74C3C", "Control" = "#3498DB")

# Get available genes
available_genes <- colnames(gene_expr)

# ============================================================================
# UI
# ============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "COVID-19 Lung Atlas"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("UMAP Explorer", tabName = "umap", icon = icon("project-diagram")),
      menuItem("Gene Expression", tabName = "expression", icon = icon("dna")),
      menuItem("Differential Expression", tabName = "de", icon = icon("chart-bar")),
      menuItem("Cell Proportions", tabName = "proportions", icon = icon("chart-pie")),
      menuItem("Cell-Cell Interactions", tabName = "interactions", icon = icon("exchange-alt")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f4f4f4; }
        .box { border-top: 3px solid #3c8dbc; }
      "))
    ),
    
    tabItems(
      # ====== UMAP Tab ======
      tabItem(
        tabName = "umap",
        fluidRow(
          box(
            title = "UMAP Visualization", status = "primary", solidHeader = TRUE,
            width = 8, height = "600px",
            plotlyOutput("umap_plot", height = "550px")
          ),
          box(
            title = "Options", status = "info", solidHeader = TRUE,
            width = 4,
            selectInput("color_by", "Color by:",
                        choices = c("Cell Type" = "cell_type",
                                    "Condition" = "condition",
                                    "Sample" = "sample_id",
                                    "Leiden Cluster" = "leiden"),
                        selected = "cell_type"),
            selectInput("split_by", "Split by:",
                        choices = c("None" = "none",
                                    "Condition" = "condition"),
                        selected = "none"),
            hr(),
            h4("Summary Statistics"),
            verbatimTextOutput("umap_summary")
          )
        )
      ),
      
      # ====== Gene Expression Tab ======
      tabItem(
        tabName = "expression",
        fluidRow(
          box(
            title = "Gene Expression on UMAP", status = "primary", solidHeader = TRUE,
            width = 8, height = "600px",
            plotlyOutput("gene_umap", height = "550px")
          ),
          box(
            title = "Gene Selection", status = "info", solidHeader = TRUE,
            width = 4,
            selectizeInput("gene_select", "Select Gene:",
                           choices = available_genes,
                           selected = "NEAT1"),
            hr(),
            h4("Expression by Condition"),
            plotOutput("gene_violin", height = "200px"),
            hr(),
            h4("Expression by Cell Type"),
            plotOutput("gene_dotplot", height = "200px")
          )
        )
      ),
      
      # ====== Differential Expression Tab ======
      tabItem(
        tabName = "de",
        fluidRow(
          box(
            title = "Volcano Plot", status = "primary", solidHeader = TRUE,
            width = 8, height = "500px",
            plotlyOutput("volcano_plot", height = "450px")
          ),
          box(
            title = "Options", status = "info", solidHeader = TRUE,
            width = 4,
            selectInput("de_celltype", "Cell Type:",
                        choices = c("Global (All Cells)" = "global",
                                    unique(de_celltype$cell_type)),
                        selected = "global"),
            sliderInput("log2fc_thresh", "Log2FC Threshold:",
                        min = 0, max = 2, value = 0.5, step = 0.1),
            sliderInput("padj_thresh", "-Log10(padj) Threshold:",
                        min = 0, max = 10, value = 1.3, step = 0.1),
            hr(),
            h4("DEG Summary"),
            verbatimTextOutput("de_summary")
          )
        ),
        fluidRow(
          box(
            title = "Differentially Expressed Genes", status = "primary", solidHeader = TRUE,
            width = 12,
            DTOutput("de_table")
          )
        )
      ),
      
      # ====== Cell Proportions Tab ======
      tabItem(
        tabName = "proportions",
        fluidRow(
          box(
            title = "Cell Type Proportions by Condition", status = "primary", solidHeader = TRUE,
            width = 6, height = "500px",
            plotlyOutput("proportion_bar", height = "450px")
          ),
          box(
            title = "Cell Type Distribution", status = "primary", solidHeader = TRUE,
            width = 6, height = "500px",
            plotlyOutput("proportion_pie", height = "450px")
          )
        ),
        fluidRow(
          box(
            title = "Proportion Comparison (COVID vs Control)", status = "info", solidHeader = TRUE,
            width = 12,
            DTOutput("proportion_table")
          )
        )
      ),
      
      # ====== Cell-Cell Interactions Tab ======
      tabItem(
        tabName = "interactions",
        fluidRow(
          box(
            title = "Ligand-Receptor Interactions", status = "primary", solidHeader = TRUE,
            width = 8, height = "500px",
            plotlyOutput("interaction_plot", height = "450px")
          ),
          box(
            title = "Options", status = "info", solidHeader = TRUE,
            width = 4,
            selectInput("pathway_select", "Pathway:",
                        choices = unique(interactions$pathway),
                        selected = unique(interactions$pathway)[1]),
            sliderInput("interaction_thresh", "Min |Log2FC|:",
                        min = 0, max = 3, value = 1, step = 0.5),
            hr(),
            h4("Top Changed Interactions"),
            verbatimTextOutput("interaction_summary")
          )
        ),
        fluidRow(
          box(
            title = "Interaction Details", status = "info", solidHeader = TRUE,
            width = 12,
            DTOutput("interaction_table")
          )
        )
      ),
      
      # ====== About Tab ======
      tabItem(
        tabName = "about",
        fluidRow(
          box(
            title = "About This App", status = "primary", solidHeader = TRUE,
            width = 12,
            h3("COVID-19 Lung Atlas Explorer"),
            p("This interactive application allows exploration of single-cell RNA sequencing data from COVID-19 lung tissue."),
            hr(),
            h4("Data Source"),
            p("Melms et al. (2021) 'A molecular single-cell lung atlas of lethal COVID-19' Nature"),
            p("GEO Accession: GSE171524"),
            hr(),
            h4("Dataset Summary"),
            tags$ul(
              tags$li("94,027 cells from 27 samples"),
              tags$li("19 COVID-19 patients, 7 controls"),
              tags$li("9 major cell types identified"),
              tags$li("Single-nucleus RNA sequencing")
            ),
            hr(),
            h4("Key Findings"),
            tags$ul(
              tags$li("Myeloid cells significantly increased in COVID-19"),
              tags$li("Epithelial cells (AT1/AT2) decreased"),
              tags$li("Pathological fibroblasts expanded (CTHRC1+)"),
              tags$li("NEAT1/MALAT1 upregulated in macrophages"),
              tags$li("Impaired efferocytosis (AXL/MERTK down)")
            ),
            hr(),
            h4("Analysis Pipeline"),
            p("Scanpy, scVI integration, Leiden clustering, Wilcoxon DE, gseapy enrichment")
          )
        )
      )
    )
  )
)

# ============================================================================
# Server
# ============================================================================

server <- function(input, output, session) {
  
  # ====== UMAP Tab ======
  output$umap_plot <- renderPlotly({
    color_var <- input$color_by
    
    if (input$split_by == "none") {
      p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = .data[[color_var]],
                                  text = paste("Cell:", rownames(umap_data),
                                               "<br>Type:", cell_type,
                                               "<br>Condition:", condition))) +
        geom_point(size = 0.5, alpha = 0.6) +
        theme_minimal() +
        labs(title = paste("UMAP colored by", color_var))
      
      if (color_var == "cell_type") {
        p <- p + scale_color_manual(values = cell_type_colors)
      } else if (color_var == "condition") {
        p <- p + scale_color_manual(values = condition_colors)
      }
    } else {
      p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = .data[[color_var]])) +
        geom_point(size = 0.5, alpha = 0.6) +
        facet_wrap(~ condition) +
        theme_minimal() +
        labs(title = paste("UMAP colored by", color_var, "split by condition"))
      
      if (color_var == "cell_type") {
        p <- p + scale_color_manual(values = cell_type_colors)
      }
    }
    
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(itemsizing = "constant"))
  })
  
  output$umap_summary <- renderPrint({
    cat("Total cells:", nrow(umap_data), "\n")
    cat("\nCells by condition:\n")
    print(table(umap_data$condition))
    cat("\nCells by cell type:\n")
    print(table(umap_data$cell_type))
  })
  
  # ====== Gene Expression Tab ======
  output$gene_umap <- renderPlotly({
    gene <- input$gene_select
    
    if (gene %in% colnames(gene_expr)) {
      plot_data <- umap_data
      plot_data$expression <- gene_expr[rownames(plot_data), gene]
      
      p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = expression,
                                  text = paste("Expression:", round(expression, 2)))) +
        geom_point(size = 0.5, alpha = 0.7) +
        scale_color_viridis_c(option = "viridis") +
        theme_minimal() +
        labs(title = paste(gene, "Expression"), color = "Expr")
      
      ggplotly(p, tooltip = "text")
    }
  })
  
  output$gene_violin <- renderPlot({
    gene <- input$gene_select
    
    if (gene %in% colnames(gene_expr)) {
      plot_data <- umap_data
      plot_data$expression <- gene_expr[rownames(plot_data), gene]
      
      ggplot(plot_data, aes(x = condition, y = expression, fill = condition)) +
        geom_violin(scale = "width") +
        geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
        scale_fill_manual(values = condition_colors) +
        theme_minimal() +
        labs(x = "", y = "Expression") +
        theme(legend.position = "none")
    }
  })
  
  output$gene_dotplot <- renderPlot({
    gene <- input$gene_select
    
    if (gene %in% colnames(gene_expr)) {
      plot_data <- umap_data
      plot_data$expression <- gene_expr[rownames(plot_data), gene]
      
      summary_data <- plot_data %>%
        group_by(cell_type) %>%
        summarise(
          mean_expr = mean(expression),
          pct_expr = mean(expression > 0) * 100
        )
      
      ggplot(summary_data, aes(x = cell_type, y = 1, size = pct_expr, color = mean_expr)) +
        geom_point() +
        scale_color_viridis_c() +
        scale_size_continuous(range = c(2, 10)) +
        theme_minimal() +
        labs(x = "", y = "", size = "% Expr", color = "Mean") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })
  
  # ====== Differential Expression Tab ======
  de_data <- reactive({
    if (input$de_celltype == "global") {
      de_global
    } else {
      de_celltype %>% filter(cell_type == input$de_celltype)
    }
  })
  
  output$volcano_plot <- renderPlotly({
    df <- de_data()
    df$neg_log_padj <- -log10(df$padj + 1e-300)
    df$significance <- "NS"
    df$significance[df$log2FC > input$log2fc_thresh & df$neg_log_padj > input$padj_thresh] <- "Up"
    df$significance[df$log2FC < -input$log2fc_thresh & df$neg_log_padj > input$padj_thresh] <- "Down"
    
    colors <- c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "#AAAAAA")
    
    p <- ggplot(df, aes(x = log2FC, y = neg_log_padj, color = significance,
                        text = paste("Gene:", gene, "<br>Log2FC:", round(log2FC, 2),
                                     "<br>padj:", signif(padj, 3)))) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = input$padj_thresh, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = c(-input$log2fc_thresh, input$log2fc_thresh), 
                 linetype = "dashed", color = "gray") +
      theme_minimal() +
      labs(x = "Log2 Fold Change", y = "-Log10(padj)", 
           title = "Volcano Plot: COVID vs Control")
    
    ggplotly(p, tooltip = "text")
  })
  
  output$de_summary <- renderPrint({
    df <- de_data()
    df$neg_log_padj <- -log10(df$padj + 1e-300)
    
    up <- sum(df$log2FC > input$log2fc_thresh & df$neg_log_padj > input$padj_thresh)
    down <- sum(df$log2FC < -input$log2fc_thresh & df$neg_log_padj > input$padj_thresh)
    
    cat("Upregulated:", up, "\n")
    cat("Downregulated:", down, "\n")
    cat("Total DEGs:", up + down, "\n")
  })
  
  output$de_table <- renderDT({
    df <- de_data()
    df <- df %>%
      filter(abs(log2FC) > input$log2fc_thresh & padj < 10^(-input$padj_thresh)) %>%
      arrange(padj) %>%
      select(gene, log2FC, pvals, padj) %>%
      head(100)
    
    datatable(df, options = list(pageLength = 10)) %>%
      formatSignif(columns = c("log2FC", "pvals", "padj"), digits = 3)
  })
  
  # ====== Cell Proportions Tab ======
  output$proportion_bar <- renderPlotly({
    prop_data <- umap_data %>%
      group_by(condition, cell_type) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(condition) %>%
      mutate(prop = n / sum(n) * 100)
    
    p <- ggplot(prop_data, aes(x = cell_type, y = prop, fill = condition)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = condition_colors) +
      theme_minimal() +
      labs(x = "", y = "Proportion (%)", title = "Cell Type Proportions") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
  })
  
  output$proportion_pie <- renderPlotly({
    prop_data <- umap_data %>%
      group_by(cell_type) %>%
      summarise(n = n()) %>%
      mutate(prop = n / sum(n) * 100)
    
    plot_ly(prop_data, labels = ~cell_type, values = ~n, type = "pie",
            marker = list(colors = cell_type_colors[prop_data$cell_type])) %>%
      layout(title = "Overall Cell Type Distribution")
  })
  
  output$proportion_table <- renderDT({
    prop_data <- umap_data %>%
      group_by(condition, cell_type) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(condition) %>%
      mutate(prop = n / sum(n) * 100) %>%
      select(-n) %>%
      tidyr::pivot_wider(names_from = condition, values_from = prop) %>%
      mutate(Difference = COVID - Control) %>%
      arrange(desc(abs(Difference)))
    
    datatable(prop_data, options = list(pageLength = 10)) %>%
      formatRound(columns = c("COVID", "Control", "Difference"), digits = 2)
  })
  
  # ====== Cell-Cell Interactions Tab ======
  output$interaction_plot <- renderPlotly({
    df <- interactions %>%
      filter(pathway == input$pathway_select & abs(log2FC) > input$interaction_thresh) %>%
      arrange(desc(abs(log2FC))) %>%
      head(20)
    
    df$label <- paste(df$sender, "->", df$receiver, "\n", df$ligand, "-", df$receptor)
    df$direction <- ifelse(df$log2FC > 0, "Increased", "Decreased")
    
    p <- ggplot(df, aes(x = reorder(label, log2FC), y = log2FC, fill = direction)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Increased" = "#E74C3C", "Decreased" = "#3498DB")) +
      coord_flip() +
      theme_minimal() +
      labs(x = "", y = "Log2 Fold Change (COVID/Control)", 
           title = paste(input$pathway_select, "Interactions"))
    
    ggplotly(p)
  })
  
  output$interaction_summary <- renderPrint({
    df <- interactions %>%
      filter(pathway == input$pathway_select & abs(log2FC) > input$interaction_thresh)
    
    cat("Pathway:", input$pathway_select, "\n")
    cat("Interactions shown:", nrow(df), "\n")
    cat("\nTop increased:\n")
    top_up <- df %>% arrange(desc(log2FC)) %>% head(3)
    for (i in 1:nrow(top_up)) {
      cat(paste0("  ", top_up$sender[i], " -> ", top_up$receiver[i], 
                 " (", top_up$ligand[i], "-", top_up$receptor[i], ")\n"))
    }
  })
  
  output$interaction_table <- renderDT({
    df <- interactions %>%
      filter(pathway == input$pathway_select) %>%
      arrange(desc(abs(log2FC))) %>%
      select(sender, receiver, ligand, receptor, COVID, Control, log2FC)
    
    datatable(df, options = list(pageLength = 10)) %>%
      formatRound(columns = c("COVID", "Control", "log2FC"), digits = 3)
  })
}

# ============================================================================
# Run App
# ============================================================================

shinyApp(ui = ui, server = server)
