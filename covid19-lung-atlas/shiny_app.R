# COVID-19 Lung Atlas Interactive Explorer
# Based on Melms et al. 2021 "A molecular single-cell lung atlas of lethal COVID-19"
# Data: GSE171524
library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(data.table)

# =============================================================================
# Load Data (with error handling) — uses fread() for speed
# =============================================================================
load_data <- function(path, desc, row_names = FALSE) {
  if (!file.exists(path)) {
    warning(paste("Missing data file:", path))
    return(NULL)
  }
  tryCatch({
    dt <- fread(path, showProgress = FALSE)
    df <- as.data.frame(dt)
    if (row_names && ncol(df) > 1) {
      rownames(df) <- df[[1]]
      df <- df[, -1, drop = FALSE]
    }
    df
  }, error = function(e) {
    warning(paste("Error loading", desc, ":", e$message))
    NULL
  })
}

message("Loading data...")
t0 <- Sys.time()
umap_data    <- load_data("data/shiny_app_data_umap_sampled.csv",
                           "UMAP data", row_names = TRUE)
gene_expr    <- load_data("data/shiny_app_data_gene_expression.csv",
                           "Gene expression", row_names = TRUE)
de_global    <- load_data("data/shiny_app_data_de_results.csv",
                           "Global DE")
de_celltype  <- load_data("data/shiny_app_data_de_celltype.csv",
                           "Cell type DE")
interactions <- load_data("data/shiny_app_data_interactions.csv",
                           "Interactions")
message(sprintf("Data loaded in %.1f seconds", difftime(Sys.time(), t0, units = "secs")))

# Ensure 'leiden' column exists (may not be present in all datasets)
if (!is.null(umap_data) && !"leiden" %in% colnames(umap_data)) {
  umap_data$leiden <- "0"
}

# Define colors
cell_type_colors <- c(
  "Epithelial" = "#2ECC71", "Myeloid" = "#E74C3C", "Fibroblast" = "#9B59B6",
  "Endothelial" = "#F39C12", "T_cell" = "#3498DB", "B_cell" = "#1ABC9C",
  "NK_cell" = "#E91E63", "Plasma_cell" = "#00BCD4", "Mast_cell" = "#795548"
)
condition_colors <- c("COVID" = "#E74C3C", "Control" = "#3498DB")

# Available genes for selection
available_genes <- if (!is.null(gene_expr)) sort(colnames(gene_expr)) else character(0)

# DE cell types
de_ct_choices <- if (!is.null(de_celltype) && "cell_type" %in% colnames(de_celltype)) {
  c("Global (All Cells)" = "global", setNames(unique(de_celltype$cell_type), unique(de_celltype$cell_type)))
} else {
  c("Global (All Cells)" = "global")
}

# Interaction pathways
pathway_choices <- if (!is.null(interactions) && "pathway" %in% colnames(interactions)) {
  unique(interactions$pathway)
} else {
  "None"
}

# COVID-19 gene signatures for the signature tab
signatures <- list(
  "Myeloid Dysfunction (lncRNA)" = c("NEAT1", "MALAT1"),
  "Efferocytosis" = c("AXL", "MERTK", "GAS6"),
  "DATP (Damage-Associated)" = c("KRT8", "CLDN4", "CDKN1A", "KRT17", "SOX4"),
  "Fibrosis" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "CTHRC1"),
  "T Cell Exhaustion" = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4"),
  "Interferon Response" = c("ISG15", "IFI6", "IFI27", "MX1", "MX2", "IFIT1"),
  "Inflammatory Cytokines" = c("IL1B", "IL6", "TNF", "CXCL8", "CCL2"),
  "Epithelial AT1" = c("AGER", "PDPN", "CLIC5", "HOPX"),
  "Epithelial AT2" = c("SFTPC", "SFTPA1", "SFTPB", "ABCA3"),
  "Alveolar Macrophage" = c("MARCO", "FABP4", "MCEMP1", "PPARG"),
  "Monocyte-Derived Macrophage" = c("FCN1", "S100A8", "S100A9", "S100A12")
)

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  title = "COVID-19 Lung Atlas",
  skin = "blue",
  dashboardHeader(
    title = tags$span(
      tags$i(class = "fas fa-lungs", style = "margin-right: 6px;"),
      "COVID-19 Lung Atlas"
    ),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "sidebar",
      menuItem("UMAP Explorer",          tabName = "umap",         icon = icon("diagram-project")),
      menuItem("Gene Expression",        tabName = "expression",   icon = icon("dna")),
      menuItem("Gene Signatures",        tabName = "signatures",   icon = icon("layer-group")),
      menuItem("Differential Expression",tabName = "de",           icon = icon("chart-bar")),
      menuItem("Cell Proportions",       tabName = "proportions",  icon = icon("chart-pie")),
      menuItem("Sample Explorer",        tabName = "samples",      icon = icon("flask")),
      menuItem("Cell-Cell Interactions", tabName = "interactions", icon = icon("arrows-alt-h")),
      menuItem("About",                  tabName = "about",        icon = icon("info-circle"))
    ),
    hr(),
    div(style = "padding: 10px; color: #b8c7ce; font-size: 11px;",
        p(icon("database"), paste("Cells:", format(nrow(umap_data), big.mark = ","))),
        p(icon("flask"),    paste("Samples:", length(unique(umap_data$sample_id)))),
        p(icon("dna"),      paste("Genes:", length(available_genes)))
    )
  ),

  dashboardBody(
    tags$head(
      # FontAwesome 6 (includes FA5 icon names as aliases)
      tags$link(
        rel  = "stylesheet",
        href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.1/css/all.min.css"
      ),
      tags$style(HTML("
      /* Layout */
      .content-wrapper { background-color: #f0f2f5; }

      /* Boxes */
      .box {
        border-radius: 6px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.10);
        border-top: none;
      }
      .box-header {
        border-bottom: 1px solid #f0f0f0;
        border-radius: 6px 6px 0 0;
      }
      .box.box-solid.box-primary > .box-header { background: #3c8dbc; }
      .box.box-solid.box-info    > .box-header { background: #00c0ef; }

      /* Sidebar */
      .main-sidebar { box-shadow: 2px 0 6px rgba(0,0,0,0.15); }
      .sidebar-menu > li > a { font-size: 13px; padding: 10px 15px; }
      .sidebar-menu > li.active > a { border-left: 3px solid #fff; }

      /* Header title icon */
      .logo .fas { vertical-align: middle; }

      /* Tables */
      .dataTables_wrapper { font-size: 13px; }

      /* Inputs */
      .selectize-input { font-size: 13px; }
      .btn-download { margin-top: 5px; }

      /* Info boxes */
      .info-box { min-height: 80px; border-radius: 6px; }
      .small-box { border-radius: 6px; }
      .small-box .icon { font-size: 60px; }

      /* Value boxes in About tab */
      .small-box h3 { font-size: 28px; }

      /* Download buttons */
      .btn-sm { font-size: 12px; }
    "))
    ),

    tabItems(
      # ====== UMAP Tab ======
      tabItem(tabName = "umap",
        fluidRow(
          box(title = "UMAP Visualization", status = "primary", solidHeader = TRUE,
              width = 8, height = "620px",
              plotlyOutput("umap_plot", height = "560px")
          ),
          box(title = "Options", status = "info", solidHeader = TRUE, width = 4,
              selectInput("color_by", "Color by:",
                          choices = c("Cell Type" = "cell_type", "Condition" = "condition",
                                      "Sample" = "sample_id", "Leiden Cluster" = "leiden"),
                          selected = "cell_type"),
              selectInput("split_by", "Split by:",
                          choices = c("None" = "none", "Condition" = "condition"),
                          selected = "none"),
              sliderInput("point_size", "Point size:", min = 0.1, max = 3, value = 0.8, step = 0.1),
              sliderInput("point_alpha", "Opacity:", min = 0.1, max = 1, value = 0.6, step = 0.1),
              hr(),
              h4("Dataset Summary"),
              verbatimTextOutput("umap_summary"),
              downloadButton("download_umap", "Download Plot Data", class = "btn-sm btn-download")
          )
        )
      ),

      # ====== Gene Expression Tab ======
      tabItem(tabName = "expression",
        fluidRow(
          box(title = "Gene Expression on UMAP", status = "primary", solidHeader = TRUE,
              width = 8, height = "620px",
              plotlyOutput("gene_umap", height = "560px")
          ),
          box(title = "Gene Selection", status = "info", solidHeader = TRUE, width = 4,
              selectizeInput("gene_select", "Select Gene:",
                             choices = NULL,  # Server-side update for performance
                             selected = NULL,
                             options = list(placeholder = "Type gene name...")),
              hr(),
              h4("Expression by Condition"),
              plotOutput("gene_violin", height = "200px"),
              hr(),
              h4("Expression by Cell Type"),
              plotOutput("gene_dotplot", height = "250px")
          )
        )
      ),

      # ====== Gene Signatures Tab ======
      tabItem(tabName = "signatures",
        fluidRow(
          box(title = "COVID-19 Gene Signature Explorer", status = "primary", solidHeader = TRUE,
              width = 12,
              selectInput("sig_select", "Select Signature:",
                          choices = names(signatures),
                          selected = names(signatures)[1]),
              p(textOutput("sig_genes_text"), style = "color: #666; font-style: italic;")
          )
        ),
        fluidRow(
          box(title = "Signature Expression by Cell Type", status = "primary", solidHeader = TRUE,
              width = 6, height = "450px",
              plotOutput("sig_heatmap", height = "400px")
          ),
          box(title = "Signature Expression: COVID vs Control", status = "primary", solidHeader = TRUE,
              width = 6, height = "450px",
              plotOutput("sig_violin", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Individual Gene Expression", status = "info", solidHeader = TRUE,
              width = 12, height = "400px",
              plotlyOutput("sig_umap_grid", height = "350px")
          )
        )
      ),

      # ====== Differential Expression Tab ======
      tabItem(tabName = "de",
        fluidRow(
          box(title = "Volcano Plot", status = "primary", solidHeader = TRUE,
              width = 8, height = "520px",
              plotlyOutput("volcano_plot", height = "470px")
          ),
          box(title = "Options", status = "info", solidHeader = TRUE, width = 4,
              selectInput("de_celltype", "Cell Type:", choices = de_ct_choices, selected = "global"),
              sliderInput("log2fc_thresh", "Log2FC Threshold:", min = 0, max = 3, value = 0.5, step = 0.1),
              sliderInput("padj_thresh", "-Log10(padj) Threshold:", min = 0, max = 20, value = 1.3, step = 0.1),
              hr(),
              h4("DEG Summary"),
              verbatimTextOutput("de_summary"),
              downloadButton("download_de", "Download DEG Table", class = "btn-sm btn-download")
          )
        ),
        fluidRow(
          box(title = "Differentially Expressed Genes", status = "primary", solidHeader = TRUE,
              width = 12,
              DTOutput("de_table")
          )
        )
      ),

      # ====== Cell Proportions Tab ======
      tabItem(tabName = "proportions",
        fluidRow(
          box(title = "Cell Type Proportions by Condition", status = "primary", solidHeader = TRUE,
              width = 6, height = "520px",
              plotlyOutput("proportion_bar", height = "470px")
          ),
          box(title = "Cell Type Distribution", status = "primary", solidHeader = TRUE,
              width = 6, height = "520px",
              plotlyOutput("proportion_pie", height = "470px")
          )
        ),
        fluidRow(
          box(title = "Proportion Comparison (COVID vs Control)", status = "info", solidHeader = TRUE,
              width = 12,
              DTOutput("proportion_table")
          )
        )
      ),

      # ====== Sample Explorer Tab ======
      tabItem(tabName = "samples",
        fluidRow(
          box(title = "Cells per Sample", status = "primary", solidHeader = TRUE,
              width = 12, height = "400px",
              plotlyOutput("sample_bar", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Sample Details", status = "info", solidHeader = TRUE, width = 6,
              DTOutput("sample_table")
          ),
          box(title = "Highlight Sample on UMAP", status = "primary", solidHeader = TRUE, width = 6,
              selectInput("highlight_sample", "Select Sample:", choices = NULL),
              plotlyOutput("sample_umap", height = "350px")
          )
        )
      ),

      # ====== Cell-Cell Interactions Tab ======
      tabItem(tabName = "interactions",
        fluidRow(
          box(title = "Ligand-Receptor Interactions", status = "primary", solidHeader = TRUE,
              width = 8, height = "520px",
              plotlyOutput("interaction_plot", height = "470px")
          ),
          box(title = "Options", status = "info", solidHeader = TRUE, width = 4,
              selectInput("pathway_select", "Pathway:", choices = pathway_choices,
                          selected = pathway_choices[1]),
              sliderInput("interaction_thresh", "Min |Log2FC|:", min = 0, max = 5, value = 1, step = 0.25),
              hr(),
              h4("Top Changed Interactions"),
              verbatimTextOutput("interaction_summary")
          )
        ),
        fluidRow(
          box(title = "Interaction Details", status = "info", solidHeader = TRUE, width = 12,
              DTOutput("interaction_table")
          )
        )
      ),

      # ====== About Tab ======
      tabItem(tabName = "about",
        fluidRow(
          box(title = "About This Application", status = "primary", solidHeader = TRUE, width = 8,
              h3("COVID-19 Lung Atlas Explorer"),
              p("Interactive exploration of single-nucleus RNA sequencing data from COVID-19 and healthy control lung tissue."),
              hr(),
              h4("Data Source"),
              p(strong("Melms et al. (2021)"), "'A molecular single-cell lung atlas of lethal COVID-19'", em("Nature"), "590, 635-641."),
              p("GEO Accession:", tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524",
                                         "GSE171524", target = "_blank")),
              hr(),
              h4("Key Findings"),
              tags$ol(
                tags$li(strong("Myeloid expansion:"), "Significant increase in myeloid cells (macrophages, monocytes) in COVID-19 lungs"),
                tags$li(strong("Epithelial loss:"), "Reduced AT1/AT2 alveolar epithelial cells"),
                tags$li(strong("Pathological fibroblasts:"), "CTHRC1+ fibroblasts expanded ~2.5-fold, driving fibrosis"),
                tags$li(strong("Macrophage dysfunction:"), "NEAT1/MALAT1 upregulated; AXL/MERTK downregulated (impaired efferocytosis)"),
                tags$li(strong("DATP cells:"), "Damage-associated transient progenitors (KRT8+/CLDN4+) expanded in COVID-19")
              ),
              hr(),
              h4("Analysis Pipeline"),
              tags$ul(
                tags$li("Data loading and QC: Scanpy, Scrublet (doublet detection)"),
                tags$li("Integration: scVI (variational inference for batch correction)"),
                tags$li("Clustering: Leiden algorithm on scVI latent space"),
                tags$li("Annotation: Marker gene scoring + manual curation"),
                tags$li("Differential expression: Wilcoxon rank-sum test + pseudobulk"),
                tags$li("Visualization: This R Shiny application")
              )
          ),
          box(title = "Dataset Summary", status = "info", solidHeader = TRUE, width = 4,
              valueBox(format(nrow(umap_data), big.mark = ","), "Total Cells (sampled)", icon = icon("dot-circle"),   color = "blue",   width = 12),
              valueBox(length(unique(umap_data$sample_id)), "Samples",         icon = icon("flask"),       color = "green",  width = 12),
              valueBox(length(unique(umap_data$cell_type)), "Cell Types",      icon = icon("layer-group"), color = "purple", width = 12),
              valueBox(length(available_genes),             "Genes Available", icon = icon("dna"),         color = "orange", width = 12)
          )
        )
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  # Server-side gene selectize (much faster for many genes)
  updateSelectizeInput(session, "gene_select",
                       choices = available_genes,
                       selected = if ("NEAT1" %in% available_genes) "NEAT1" else available_genes[1],
                       server = TRUE)

  # Update sample choices
  updateSelectInput(session, "highlight_sample",
                    choices = sort(unique(umap_data$sample_id)),
                    selected = sort(unique(umap_data$sample_id))[1])

  # ====== UMAP Tab ======
  output$umap_plot <- renderPlotly({
    req(umap_data)
    color_var <- input$color_by

    if (input$split_by == "none") {
      p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = .data[[color_var]],
                                  text = paste("Type:", cell_type,
                                               "<br>Condition:", condition,
                                               "<br>Sample:", sample_id))) +
        geom_point(size = input$point_size, alpha = input$point_alpha) +
        theme_minimal(base_size = 12) +
        labs(title = paste("UMAP -", gsub("_", " ", color_var)), color = "")

      if (color_var == "cell_type") {
        p <- p + scale_color_manual(values = cell_type_colors)
      } else if (color_var == "condition") {
        p <- p + scale_color_manual(values = condition_colors)
      }
    } else {
      p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = .data[[color_var]])) +
        geom_point(size = input$point_size, alpha = input$point_alpha) +
        facet_wrap(~ condition) +
        theme_minimal(base_size = 12) +
        labs(title = paste("UMAP -", gsub("_", " ", color_var), "(split by condition)"), color = "")

      if (color_var == "cell_type") p <- p + scale_color_manual(values = cell_type_colors)
    }

    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(itemsizing = "constant", font = list(size = 10)))
  })

  output$umap_summary <- renderPrint({
    req(umap_data)
    cat("Cells displayed:", format(nrow(umap_data), big.mark = ","), "\n\n")
    cat("By condition:\n")
    cond_tbl <- table(umap_data$condition)
    for (nm in names(cond_tbl)) cat(sprintf("  %-10s %s\n", nm, format(cond_tbl[nm], big.mark = ",")))
    cat("\nBy cell type:\n")
    ct_tbl <- sort(table(umap_data$cell_type), decreasing = TRUE)
    for (nm in names(ct_tbl)) cat(sprintf("  %-15s %s\n", nm, format(ct_tbl[nm], big.mark = ",")))
  })

  output$download_umap <- downloadHandler(
    filename = function() paste0("umap_data_", Sys.Date(), ".csv"),
    content = function(file) write.csv(umap_data, file, row.names = TRUE)
  )

  # ====== Gene Expression Tab ======
  output$gene_umap <- renderPlotly({
    req(input$gene_select, gene_expr, umap_data)
    gene <- input$gene_select
    if (!(gene %in% colnames(gene_expr))) return(NULL)

    plot_data <- umap_data
    matched_rows <- intersect(rownames(plot_data), rownames(gene_expr))
    plot_data <- plot_data[matched_rows, ]
    plot_data$expression <- gene_expr[matched_rows, gene]

    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = expression,
                                text = paste("Expr:", round(expression, 2),
                                             "<br>Type:", cell_type,
                                             "<br>Cond:", condition))) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_viridis_c(option = "magma", direction = -1) +
      theme_minimal(base_size = 12) +
      labs(title = paste(gene, "Expression"), color = "Expr")

    ggplotly(p, tooltip = "text")
  })

  output$gene_violin <- renderPlot({
    req(input$gene_select, gene_expr, umap_data)
    gene <- input$gene_select
    if (!(gene %in% colnames(gene_expr))) return(NULL)

    plot_data <- umap_data
    matched_rows <- intersect(rownames(plot_data), rownames(gene_expr))
    plot_data <- plot_data[matched_rows, ]
    plot_data$expression <- gene_expr[matched_rows, gene]

    ggplot(plot_data, aes(x = condition, y = expression, fill = condition)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3, alpha = 0.8) +
      scale_fill_manual(values = condition_colors) +
      theme_minimal(base_size = 11) +
      labs(x = "", y = "Expression", title = paste(gene, "- COVID vs Control")) +
      theme(legend.position = "none")
  })

  output$gene_dotplot <- renderPlot({
    req(input$gene_select, gene_expr, umap_data)
    gene <- input$gene_select
    if (!(gene %in% colnames(gene_expr))) return(NULL)

    plot_data <- umap_data
    matched_rows <- intersect(rownames(plot_data), rownames(gene_expr))
    plot_data <- plot_data[matched_rows, ]
    plot_data$expression <- gene_expr[matched_rows, gene]

    summary_data <- plot_data %>%
      group_by(cell_type, condition) %>%
      summarise(mean_expr = mean(expression), pct_expr = mean(expression > 0) * 100, .groups = "drop")

    ggplot(summary_data, aes(x = cell_type, y = condition, size = pct_expr, color = mean_expr)) +
      geom_point() +
      scale_color_viridis_c(option = "magma", direction = -1) +
      scale_size_continuous(range = c(1, 10), name = "% Expressing") +
      theme_minimal(base_size = 11) +
      labs(x = "", y = "", color = "Mean Expr", title = paste(gene, "by Cell Type")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # ====== Gene Signatures Tab ======
  output$sig_genes_text <- renderText({
    sig_genes <- signatures[[input$sig_select]]
    avail <- sig_genes[sig_genes %in% colnames(gene_expr)]
    paste("Genes:", paste(avail, collapse = ", "),
          paste0("(", length(avail), "/", length(sig_genes), " available)"))
  })

  output$sig_heatmap <- renderPlot({
    req(gene_expr, umap_data)
    sig_genes <- signatures[[input$sig_select]]
    avail <- sig_genes[sig_genes %in% colnames(gene_expr)]
    if (length(avail) == 0) return(NULL)

    matched_rows <- intersect(rownames(umap_data), rownames(gene_expr))
    plot_data <- umap_data[matched_rows, ]
    plot_data[avail] <- gene_expr[matched_rows, avail]

    # Mean expression per cell type
    heat_data <- plot_data %>%
      group_by(cell_type) %>%
      summarise(across(all_of(avail), mean), .groups = "drop") %>%
      pivot_longer(cols = all_of(avail), names_to = "gene", values_to = "expression") %>%
      group_by(gene) %>%
      mutate(scaled = scale(expression)[, 1]) %>%
      ungroup()

    ggplot(heat_data, aes(x = gene, y = cell_type, fill = scaled)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0, name = "Z-score") +
      theme_minimal(base_size = 12) +
      labs(x = "", y = "", title = paste(input$sig_select, "- Mean Expression")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
  })

  output$sig_violin <- renderPlot({
    req(gene_expr, umap_data)
    sig_genes <- signatures[[input$sig_select]]
    avail <- sig_genes[sig_genes %in% colnames(gene_expr)]
    if (length(avail) == 0) return(NULL)

    matched_rows <- intersect(rownames(umap_data), rownames(gene_expr))
    plot_data <- umap_data[matched_rows, ]

    # Compute mean signature score per cell
    expr_mat <- as.matrix(gene_expr[matched_rows, avail, drop = FALSE])
    plot_data$sig_score <- rowMeans(expr_mat)

    ggplot(plot_data, aes(x = cell_type, y = sig_score, fill = condition)) +
      geom_violin(scale = "width", position = position_dodge(width = 0.8), alpha = 0.7) +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), outlier.size = 0.2, alpha = 0.8) +
      scale_fill_manual(values = condition_colors) +
      theme_minimal(base_size = 12) +
      labs(x = "", y = "Mean Signature Score", title = paste(input$sig_select, "Score"), fill = "Condition") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$sig_umap_grid <- renderPlotly({
    req(gene_expr, umap_data)
    sig_genes <- signatures[[input$sig_select]]
    avail <- sig_genes[sig_genes %in% colnames(gene_expr)]
    if (length(avail) == 0) return(NULL)

    matched_rows <- intersect(rownames(umap_data), rownames(gene_expr))
    plot_data <- umap_data[matched_rows, ]
    expr_mat <- as.matrix(gene_expr[matched_rows, avail, drop = FALSE])
    plot_data$sig_score <- rowMeans(expr_mat)

    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = sig_score,
                                text = paste("Score:", round(sig_score, 3),
                                             "<br>Type:", cell_type))) +
      geom_point(size = 0.3, alpha = 0.6) +
      scale_color_viridis_c(option = "inferno") +
      theme_minimal(base_size = 11) +
      labs(title = paste(input$sig_select, "Score on UMAP"), color = "Score")

    ggplotly(p, tooltip = "text")
  })

  # ====== Differential Expression Tab ======
  de_data <- reactive({
    req(de_global)
    if (input$de_celltype == "global") {
      de_global
    } else if (!is.null(de_celltype)) {
      de_celltype %>% filter(cell_type == input$de_celltype)
    } else {
      de_global
    }
  })

  output$volcano_plot <- renderPlotly({
    df <- de_data()
    req(nrow(df) > 0)
    df$neg_log_padj <- -log10(pmax(df$padj, 1e-300))
    df$significance <- "NS"
    df$significance[df$log2FC > input$log2fc_thresh & df$neg_log_padj > input$padj_thresh] <- "Up"
    df$significance[df$log2FC < -input$log2fc_thresh & df$neg_log_padj > input$padj_thresh] <- "Down"

    colors <- c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "#CCCCCC")

    p <- ggplot(df, aes(x = log2FC, y = neg_log_padj, color = significance,
                        text = paste("Gene:", gene, "<br>Log2FC:", round(log2FC, 3),
                                     "<br>padj:", signif(padj, 3)))) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = input$padj_thresh, linetype = "dashed", color = "gray60") +
      geom_vline(xintercept = c(-input$log2fc_thresh, input$log2fc_thresh),
                 linetype = "dashed", color = "gray60") +
      theme_minimal(base_size = 12) +
      labs(x = "Log2 Fold Change (COVID vs Control)", y = "-Log10(padj)",
           title = ifelse(input$de_celltype == "global", "Volcano: All Cells",
                          paste("Volcano:", input$de_celltype)),
           color = "")

    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", y = -0.1))
  })

  output$de_summary <- renderPrint({
    df <- de_data()
    df$neg_log_padj <- -log10(pmax(df$padj, 1e-300))
    up <- sum(df$log2FC > input$log2fc_thresh & df$neg_log_padj > input$padj_thresh, na.rm = TRUE)
    down <- sum(df$log2FC < -input$log2fc_thresh & df$neg_log_padj > input$padj_thresh, na.rm = TRUE)
    cat(sprintf("Upregulated:   %d\n", up))
    cat(sprintf("Downregulated: %d\n", down))
    cat(sprintf("Total DEGs:    %d\n", up + down))
  })

  output$de_table <- renderDT({
    df <- de_data()
    # Select columns that exist
    display_cols <- intersect(c("gene", "log2FC", "pvals", "padj", "scores"), colnames(df))
    df <- df %>%
      filter(abs(log2FC) > input$log2fc_thresh & padj < 10^(-input$padj_thresh)) %>%
      arrange(padj) %>%
      select(all_of(display_cols)) %>%
      head(200)

    num_cols <- intersect(c("log2FC", "pvals", "padj", "scores"), display_cols)
    datatable(df, options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE, class = "compact stripe") %>%
      formatSignif(columns = num_cols, digits = 3)
  })

  output$download_de <- downloadHandler(
    filename = function() paste0("DEGs_", input$de_celltype, "_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- de_data() %>%
        filter(abs(log2FC) > input$log2fc_thresh & padj < 10^(-input$padj_thresh)) %>%
        arrange(padj)
      write.csv(df, file, row.names = FALSE)
    }
  )

  # ====== Cell Proportions Tab ======
  output$proportion_bar <- renderPlotly({
    req(umap_data)
    prop_data <- umap_data %>%
      group_by(condition, cell_type) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(condition) %>%
      mutate(prop = n / sum(n) * 100)

    p <- ggplot(prop_data, aes(x = reorder(cell_type, -prop), y = prop, fill = condition,
                                text = paste(cell_type, "\n", condition, ":", round(prop, 1), "%"))) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      scale_fill_manual(values = condition_colors) +
      theme_minimal(base_size = 12) +
      labs(x = "", y = "Proportion (%)", title = "Cell Type Proportions", fill = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p, tooltip = "text")
  })

  output$proportion_pie <- renderPlotly({
    req(umap_data)
    prop_data <- umap_data %>%
      group_by(cell_type) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(prop = round(n / sum(n) * 100, 1))

    plot_ly(prop_data, labels = ~cell_type, values = ~n, type = "pie",
            textinfo = "label+percent",
            marker = list(colors = unname(cell_type_colors[prop_data$cell_type])),
            hovertemplate = "%{label}<br>%{value} cells (%{percent})<extra></extra>") %>%
      layout(title = list(text = "Overall Cell Type Distribution"),
             showlegend = TRUE, legend = list(font = list(size = 10)))
  })

  output$proportion_table <- renderDT({
    req(umap_data)
    prop_data <- umap_data %>%
      group_by(condition, cell_type) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(condition) %>%
      mutate(prop = n / sum(n) * 100) %>%
      select(-n) %>%
      pivot_wider(names_from = condition, values_from = prop, values_fill = 0) %>%
      mutate(Difference = COVID - Control,
             Direction = ifelse(Difference > 0, "Increased in COVID", "Decreased in COVID")) %>%
      arrange(desc(abs(Difference)))

    datatable(prop_data, options = list(pageLength = 10), rownames = FALSE, class = "compact stripe") %>%
      formatRound(columns = c("COVID", "Control", "Difference"), digits = 2) %>%
      formatStyle("Difference",
                  color = styleInterval(0, c("#3498DB", "#E74C3C")),
                  fontWeight = "bold")
  })

  # ====== Sample Explorer Tab ======
  output$sample_bar <- renderPlotly({
    req(umap_data)
    sample_data <- umap_data %>%
      group_by(sample_id, condition) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(condition, sample_id)

    p <- ggplot(sample_data, aes(x = reorder(sample_id, -n), y = n, fill = condition,
                                  text = paste(sample_id, "\n", condition, "\n", n, "cells"))) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = condition_colors) +
      theme_minimal(base_size = 12) +
      labs(x = "Sample", y = "Number of Cells", title = "Cells per Sample", fill = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

    ggplotly(p, tooltip = "text")
  })

  output$sample_table <- renderDT({
    req(umap_data)
    sample_summary <- umap_data %>%
      group_by(sample_id, condition) %>%
      summarise(
        n_cells = n(),
        n_cell_types = n_distinct(cell_type),
        top_cell_type = names(sort(table(cell_type), decreasing = TRUE))[1],
        .groups = "drop"
      ) %>%
      arrange(condition, sample_id)

    datatable(sample_summary, options = list(pageLength = 15), rownames = FALSE, class = "compact stripe")
  })

  output$sample_umap <- renderPlotly({
    req(umap_data, input$highlight_sample)
    plot_data <- umap_data
    plot_data$highlight <- ifelse(plot_data$sample_id == input$highlight_sample, input$highlight_sample, "Other")

    cond <- unique(plot_data$condition[plot_data$sample_id == input$highlight_sample])
    highlight_color <- if (length(cond) > 0 && cond[1] == "COVID") "#E74C3C" else "#3498DB"

    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = highlight, alpha = highlight,
                                text = paste("Sample:", sample_id, "<br>Type:", cell_type))) +
      geom_point(size = 0.4) +
      scale_color_manual(values = c(setNames(highlight_color, input$highlight_sample), "Other" = "#DDDDDD")) +
      scale_alpha_manual(values = c(setNames(0.8, input$highlight_sample), "Other" = 0.15)) +
      theme_minimal(base_size = 11) +
      labs(title = paste("Sample:", input$highlight_sample), color = "", alpha = "") +
      guides(alpha = "none")

    ggplotly(p, tooltip = "text") %>%
      layout(showlegend = FALSE)
  })

  # ====== Cell-Cell Interactions Tab ======
  output$interaction_plot <- renderPlotly({
    req(interactions)
    df <- interactions %>%
      filter(pathway == input$pathway_select & abs(log2FC) > input$interaction_thresh) %>%
      arrange(desc(abs(log2FC))) %>%
      head(20)

    if (nrow(df) == 0) {
      return(plotly_empty() %>% layout(title = "No interactions match current filters"))
    }

    # Short y-axis label: "Ligand > Receptor" and full detail in tooltip
    df$label <- paste0(df$ligand, " \u2192 ", df$receptor)
    # Deduplicate labels by appending sender/receiver initials when needed
    if (anyDuplicated(df$label)) {
      df$label <- paste0(df$ligand, " \u2192 ", df$receptor, "  [", df$sender, " > ", df$receiver, "]")
    }
    df$direction <- ifelse(df$log2FC > 0, "Increased", "Decreased")

    p <- ggplot(df, aes(x = reorder(label, log2FC), y = log2FC, fill = direction,
                        text = paste0("<b>", ligand, " \u2192 ", receptor, "</b>",
                                      "<br>Sender: ", sender,
                                      "<br>Receiver: ", receiver,
                                      "<br>Log2FC: ", round(log2FC, 2)))) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("Increased" = "#E74C3C", "Decreased" = "#3498DB")) +
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_text(size = 10)) +
      labs(x = "", y = "Log2 Fold Change (COVID / Control)",
           title = paste(input$pathway_select, "Pathway"), fill = "")

    ggplotly(p, tooltip = "text") %>%
      layout(margin = list(l = 150))
  })

  output$interaction_summary <- renderPrint({
    req(interactions)
    df <- interactions %>%
      filter(pathway == input$pathway_select & abs(log2FC) > input$interaction_thresh)

    cat("Pathway:", input$pathway_select, "\n")
    cat("Interactions:", nrow(df), "\n")

    if (nrow(df) > 0) {
      cat("\nTop increased:\n")
      top_up <- df %>% arrange(desc(log2FC)) %>% head(3)
      for (i in seq_len(nrow(top_up))) {
        cat(sprintf("  %s -> %s (%s-%s) FC=%.2f\n",
                    top_up$sender[i], top_up$receiver[i],
                    top_up$ligand[i], top_up$receptor[i], top_up$log2FC[i]))
      }
      cat("\nTop decreased:\n")
      top_down <- df %>% arrange(log2FC) %>% head(3)
      for (i in seq_len(nrow(top_down))) {
        cat(sprintf("  %s -> %s (%s-%s) FC=%.2f\n",
                    top_down$sender[i], top_down$receiver[i],
                    top_down$ligand[i], top_down$receptor[i], top_down$log2FC[i]))
      }
    } else {
      cat("\nNo interactions above threshold.\nTry lowering the Min |Log2FC| slider.")
    }
  })

  output$interaction_table <- renderDT({
    req(interactions)
    df <- interactions %>%
      filter(pathway == input$pathway_select) %>%
      arrange(desc(abs(log2FC)))

    display_cols <- intersect(c("sender", "receiver", "ligand", "receptor", "COVID", "Control", "log2FC"), colnames(df))
    df <- df %>% select(all_of(display_cols))

    num_cols <- intersect(c("COVID", "Control", "log2FC"), display_cols)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE, class = "compact stripe") %>%
      formatRound(columns = num_cols, digits = 3) %>%
      formatStyle("log2FC",
                  color = styleInterval(0, c("#3498DB", "#E74C3C")),
                  fontWeight = "bold")
  })
}

# Run App
shinyApp(ui = ui, server = server)
