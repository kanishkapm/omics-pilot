# OmicsPilot — Single-cell RNA-seq analysis with AI assistant
# GitHub: https://github.com/kanishkapm/omics-pilot
#
# Setup Instructions:
# 1. Create .Renviron in your home directory with: OPENAI_API_KEY="sk-proj-YOUR_KEY"
# 2. Run: shiny::runApp()
#
# DO NOT commit .Renviron to GitHub - it's protected by .gitignore

options(shiny.maxRequestSize = 1000 * 1000 * 1000) # 1 GB

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rmarkdown)
library(DT)
library(EnhancedVolcano)
library(pheatmap)
library(stringdist)
library(httr2)
suppressWarnings(suppressMessages(require(shinycssloaders)))

custom_faq_path <- file.path(tempdir(), "omics_custom_faqs.csv")

load_custom_faq <- function() {
  if (file.exists(custom_faq_path)) {
    df <- tryCatch(read.csv(custom_faq_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df) && all(c("question", "answer") %in% colnames(df))) {
      return(df %>% distinct(question, .keep_all = TRUE))
    }
  }
  data.frame(question = character(0), answer = character(0), stringsAsFactors = FALSE)
}

save_custom_faq <- function(df) {
  tryCatch({
    write.csv(df %>% distinct(question, .keep_all = TRUE), custom_faq_path, row.names = FALSE)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

# Load OpenAI API key from environment
load_openai_key <- function() {
  key <- Sys.getenv("OPENAI_API_KEY")
  if (key == "") {
    return(NULL)
  }
  key
}

# Query OpenAI API
query_openai <- function(query, api_key) {
  if (is.null(api_key)) return(NULL)
  
  tryCatch({
    resp <- request("https://api.openai.com/v1/chat/completions") %>%
      req_headers("Authorization" = paste("Bearer", api_key),
                  "Content-Type" = "application/json") %>%
      req_body_json(list(
        model = "gpt-4o-mini",
        messages = list(
          list(role = "system", content = "You are OmicsPilot, a bioinformatics assistant for single-cell RNA-seq analysis. Be concise and technical. Keep responses under 200 words."),
          list(role = "user", content = query)
        ),
        temperature = 0.7,
        max_tokens = 300
      )) %>%
      req_perform()
    
    result <- resp_body_json(resp)
    if (!is.null(result$choices) && length(result$choices) > 0) {
      return(result$choices[[1]]$message$content)
    }
    NULL
  }, error = function(e) {
    message("OpenAI API error: ", e$message)
    return(NULL)
  })
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background: #f7fafc; color: #1f2937; }
      .sidebar .well { background: #ffffff; border-radius: 8px; border: 1px solid #e2e8f0; padding: 12px; }
      .main-panel { padding-left: 18px; }
      .app-title { font-weight: 700; color: #0f172a; }
      .bio-badge { background: linear-gradient(90deg,#0ea5a2,#06b6d4); color: white; padding:4px 8px; border-radius:6px; font-weight:600; }
      .ai-msg-user { background: #e6fffa; padding:10px; border-radius:8px; margin-bottom:8px; }
      .ai-msg-assistant { background: #f1f5f9; padding:10px; border-radius:8px; margin-bottom:8px; }
      .small-note { font-size: 12px; color:#64748b; }
      .panel-box { background: white; padding: 12px; border-radius: 8px; border: 1px solid #e6eef4; margin-bottom: 12px; }
      .footer-note { font-size: 12px; color: #94a3b8; margin-top: 8px; }
      .api-status { font-size: 11px; padding: 4px 8px; border-radius: 4px; display: inline-block; }
      .api-active { background: #d1fae5; color: #065f46; }
      .api-inactive { background: #fee2e2; color: #7f1d1d; }
    "))
  ),
  
  div(class="app-header", fluidRow(
    column(9, h2("OmicsPilot", class="app-title"), span(class="bio-badge", "Bioinformatics + AI")),
    column(3, tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0b/DNA_icon.svg/1200px-DNA_icon.svg.png", height="48px", style="float:right;"))
  )),
  hr(),
  
  sidebarLayout(
    sidebarPanel(class = "sidebar",
                 div(class="panel-box",
                     h4("Input (10X files)"),
                     fileInput("matrix_file", "Matrix (.mtx or .mtx.gz)", accept = c(".mtx", ".mtx.gz", ".gz")),
                     fileInput("features_file", "Features (.tsv or .tsv.gz)", accept = c(".tsv", ".tsv.gz", ".gz")),
                     fileInput("barcodes_file", "Barcodes (.tsv or .tsv.gz)", accept = c(".tsv", ".tsv.gz", ".gz")),
                     helpText("Each file must be < 1 GB. Browser and server validated."),
                     actionButton("run_analysis", "Run Analysis", class = "btn-primary")
                 ),
                 
                 div(class="panel-box",
                     h4("AI Assistant"),
                     uiOutput("api_status_ui"),
                     br(),
                     textAreaInput("ai_query", "Ask the assistant", rows = 3, placeholder = "E.g. 'Summarize QC metrics' or 'Top markers for cluster 1'"),
                     fluidRow(column(6, actionButton("ai_send", "Send", class="btn-success")), column(6, actionButton("ai_clear", "Clear Chat"))),
                     br(),
                     fluidRow(column(4, actionButton("example_upload", "How to upload", class="btn-xs")),
                              column(4, actionButton("example_summary", "Summarize QC", class="btn-xs")),
                              column(4, actionButton("example_markers", "Top markers C0", class="btn-xs"))),
                     br(),
                     p(class="small-note", "Hybrid: Local FAQ + OpenAI for complex questions")
                 ),
                 
                 div(class="panel-box",
                     h4("Teach assistant"),
                     textInput("teach_question", "Question (short)"),
                     textAreaInput("teach_answer", "Answer", rows = 3),
                     actionButton("teach_add", "Add Q&A", class = "btn-success"),
                     br(), br(),
                     fileInput("upload_faq_csv", "Upload FAQ CSV (question,answer)", accept = ".csv"),
                     downloadButton("download_faqs", "Download current FAQs")
                 ),
                 
                 div(class="panel-box",
                     h4("Report"),
                     helpText("HTML report with plots and interpretations."),
                     downloadButton("download_report", "Download HTML Report"),
                     br(),
                     p(class="footer-note", "Requires report.Rmd in app directory.")
                 )
    ),
    
    mainPanel(class="main-panel",
              tabsetPanel(
                tabPanel("Chat",
                         br(),
                         uiOutput("ai_status"),
                         br(),
                         uiOutput("ai_chat_ui"),
                         hr(),
                         h4("Custom FAQs"),
                         DTOutput("faqs_table")
                ),
                tabPanel("Data",
                         fluidRow(
                           column(12, div(class="panel-box", h4("Uploaded files"), textOutput("matrix_info"), textOutput("features_info"), textOutput("barcodes_info")))
                         ),
                         fluidRow(
                           column(12, div(class="panel-box", h4("QC table"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(DTOutput("qc_table")) else DTOutput("qc_table")))
                         )
                ),
                tabPanel("QC Plots",
                         fluidRow(column(12, div(class="panel-box", h4("Violin plots"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(plotOutput("qc_violin")) else plotOutput("qc_violin")))),
                         fluidRow(column(6, div(class="panel-box", h4("nCount vs percent.mt"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(plotOutput("feature_scatter_mt")) else plotOutput("feature_scatter_mt"))),
                                  column(6, div(class="panel-box", h4("nCount vs nFeature"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(plotOutput("feature_scatter_features")) else plotOutput("feature_scatter_features"))))
                ),
                tabPanel("UMAP & Clusters",
                         fluidRow(column(12, div(class="panel-box", h4("UMAP"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(plotOutput("umap_plot")) else plotOutput("umap_plot")))),
                         fluidRow(column(12, div(class="panel-box", h4("UMAP summary (quick)"), verbatimTextOutput("umap_summary"))))
                ),
                tabPanel("Markers & GO",
                         fluidRow(column(6, div(class="panel-box", h4("Volcano (cluster 0)"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(plotOutput("volcano_plot")) else plotOutput("volcano_plot"))),
                                  column(6, div(class="panel-box", h4("Heatmap (top markers)"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(plotOutput("heatmap_plot")) else plotOutput("heatmap_plot")))),
                         fluidRow(column(12, div(class="panel-box", h4("Top GO terms sample"), if (requireNamespace("shinycssloaders", quietly=TRUE)) withSpinner(plotOutput("go_cluster_0")) else plotOutput("go_cluster_0"))))
                )
              )
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    seurat_obj = NULL,
    go_results = NULL,
    markers = NULL,
    chat = list(),
    custom_faq = load_custom_faq(),
    openai_key = load_openai_key(),
    api_enabled = !is.null(load_openai_key())
  )
  
  # API Status UI
  output$api_status_ui <- renderUI({
    status_class <- if (rv$api_enabled) "api-active" else "api-inactive"
    status_text <- if (rv$api_enabled) "✓ OpenAI Ready" else "✗ No API Key"
    tags$span(class = paste("api-status", status_class), status_text)
  })
  
  # ---------- Upload info ----------
  output$matrix_info <- renderText({
    if (!is.null(input$matrix_file)) paste("Matrix:", input$matrix_file$name, "(", format(input$matrix_file$size, big.mark=","), "bytes)") else "Matrix: none"
  })
  output$features_info <- renderText({
    if (!is.null(input$features_file)) paste("Features:", input$features_file$name, "(", format(input$features_file$size, big.mark=","), "bytes)") else "Features: none"
  })
  output$barcodes_info <- renderText({
    if (!is.null(input$barcodes_file)) paste("Barcodes:", input$barcodes_file$name, "(", format(input$barcodes_file$size, big.mark=","), "bytes)") else "Barcodes: none"
  })
  
  validate_input_file <- function(f, expected_type = c("matrix","features","barcodes")) {
    if (is.null(f)) return("missing")
    if (!is.null(f$size) && f$size >= 1000*1000*1000) {
      return(paste0("file_too_large: ", f$name, " (", f$size, " bytes)"))
    }
    name <- tolower(f$name)
    if (expected_type == "matrix") {
      if (!grepl("\\.mtx(\\.gz)?$", name)) return("invalid_type_matrix")
    } else {
      if (!grepl("\\.tsv(\\.gz)?$", name)) return("invalid_type_tsv")
    }
    return("ok")
  }
  
  # ---------- Analysis pipeline ----------
  observeEvent(input$run_analysis, {
    if (is.null(input$matrix_file) || is.null(input$features_file) || is.null(input$barcodes_file)) {
      showModal(modalDialog(title = "Missing files", "Please upload matrix, features and barcodes files before running analysis.", easyClose = TRUE))
      return()
    }
    v1 <- validate_input_file(input$matrix_file, "matrix")
    v2 <- validate_input_file(input$features_file, "features")
    v3 <- validate_input_file(input$barcodes_file, "barcodes")
    if (v1 != "ok" || v2 != "ok" || v3 != "ok") {
      msgs <- c()
      if (v1 != "ok") msgs <- c(msgs, paste("Matrix:", v1))
      if (v2 != "ok") msgs <- c(msgs, paste("Features:", v2))
      if (v3 != "ok") msgs <- c(msgs, paste("Barcodes:", v3))
      showModal(modalDialog(title = "File validation failed", paste(msgs, collapse = "\n"), easyClose = TRUE))
      return()
    }
    
    withProgress(message = 'Running Analysis', value = 0, {
      temp_dir <- tempdir()
      
      m_name <- input$matrix_file$name
      m_gz <- grepl("\\.gz$", m_name, ignore.case = TRUE)
      dest_matrix <- file.path(temp_dir, ifelse(m_gz, "matrix.mtx.gz", "matrix.mtx"))
      file.copy(input$matrix_file$datapath, dest_matrix, overwrite = TRUE)
      
      f_name <- input$features_file$name
      f_gz <- grepl("\\.gz$", f_name, ignore.case = TRUE)
      dest_features <- file.path(temp_dir, ifelse(f_gz, "features.tsv.gz", "features.tsv"))
      file.copy(input$features_file$datapath, dest_features, overwrite = TRUE)
      
      b_name <- input$barcodes_file$name
      b_gz <- grepl("\\.gz$", b_name, ignore.case = TRUE)
      dest_barcodes <- file.path(temp_dir, ifelse(b_gz, "barcodes.tsv.gz", "barcodes.tsv"))
      file.copy(input$barcodes_file$datapath, dest_barcodes, overwrite = TRUE)
      
      incProgress(0.1, detail = "Loading data")
      data <- tryCatch({
        Read10X(data.dir = temp_dir)
      }, error = function(e) {
        showModal(modalDialog(title = "Error loading data", paste("Read10X failed:", e$message), easyClose = TRUE))
        return(NULL)
      })
      req(data)
      
      seurat_obj <- CreateSeuratObject(counts = data, project = "SingleCellAnalysis", min.cells = 3, min.features = 200)
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      incProgress(0.2, detail = "Filtering cells")
      
      seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
      incProgress(0.3, detail = "Normalizing data")
      
      seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
      incProgress(0.4, detail = "Finding variable features")
      
      seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
      incProgress(0.5, detail = "Scaling data")
      
      seurat_obj <- ScaleData(seurat_obj)
      incProgress(0.6, detail = "Running PCA")
      
      seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
      incProgress(0.7, detail = "Running UMAP")
      
      seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
      incProgress(0.8, detail = "Clustering")
      
      seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
      seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
      incProgress(0.9, detail = "Finding DEGs")
      
      markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      
      markers_entrez <- tryCatch({
        bitr(unique(markers$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      }, error = function(e) NULL)
      if (!is.null(markers_entrez)) {
        markers <- left_join(markers, markers_entrez, by = c("gene" = "SYMBOL"))
      } else {
        markers$ENTREZID <- NA
      }
      
      go_results <- list()
      for (cluster in unique(markers$cluster)) {
        cluster_genes <- markers %>% filter(cluster == !!cluster) %>% pull(ENTREZID) %>% na.omit()
        if (length(cluster_genes) > 0) {
          go_res <- tryCatch({
            enrichGO(gene = cluster_genes, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
          }, error = function(e) NULL)
          go_results[[paste0("cluster_", cluster)]] <- go_res
        }
      }
      
      rv$seurat_obj <- seurat_obj
      rv$go_results <- go_results
      rv$markers <- markers
      incProgress(1, detail = "Analysis complete")
      
      showModal(modalDialog(title = "Analysis Complete", "Single-cell analysis finished. You can view plots and download the HTML report.", easyClose = TRUE, footer = modalButton("OK")))
    })
  })
  
  # ---------- Plots & outputs ----------
  output$qc_violin <- renderPlot({ req(rv$seurat_obj); VlnPlot(rv$seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) })
  output$feature_scatter_mt <- renderPlot({ req(rv$seurat_obj); FeatureScatter(rv$seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") })
  output$feature_scatter_features <- renderPlot({ req(rv$seurat_obj); FeatureScatter(rv$seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") })
  output$umap_plot <- renderPlot({ req(rv$seurat_obj); DimPlot(rv$seurat_obj, reduction = "umap", label = TRUE) })
  
  output$go_cluster_0 <- renderPlot({ req(rv$go_results); if (!is.null(rv$go_results[["cluster_0"]])) barplot(rv$go_results[["cluster_0"]], showCategory = 10) })
  output$volcano_plot <- renderPlot({
    req(rv$markers)
    cluster_0_markers <- rv$markers %>% filter(cluster == 0)
    if (nrow(cluster_0_markers) > 0) {
      EnhancedVolcano(cluster_0_markers, lab = cluster_0_markers$gene, x = 'avg_log2FC', y = 'p_val_adj',
                      title = 'Cluster 0 vs Others', pCutoff = 0.05, FCcutoff = 0.25)
    }
  })
  output$heatmap_plot <- renderPlot({
    req(rv$seurat_obj, rv$markers)
    top10 <- rv$markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
    if (nrow(top10) > 0) DoHeatmap(rv$seurat_obj, features = unique(top10$gene)) + NoLegend()
  })
  
  output$qc_table <- renderDT({ req(rv$seurat_obj); datatable(rv$seurat_obj@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]) })
  output$umap_summary <- renderText({
    if (is.null(rv$seurat_obj)) return("No UMAP (run analysis).")
    cl_col <- if ("seurat_clusters" %in% colnames(rv$seurat_obj@meta.data)) "seurat_clusters" else NULL
    if (is.null(cl_col)) return("No cluster labels.")
    tab <- table(rv$seurat_obj@meta.data[[cl_col]])
    paste("Clusters:", paste(names(tab), tab, sep=":", collapse="; "))
  })
  
  # ---------- Assistant ----------
  built_in_faq <- list(
    "how to upload" = "Use the Upload panel: provide matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz), then Run Analysis.",
    "how to run analysis" = "Run Analysis executes filtering, NormalizeData, PCA, UMAP, clustering and marker detection (FindAllMarkers).",
    "what is percent.mt" = "percent.mt is percent mitochondrial counts per cell. High values suggest low-quality cells; common threshold: <5%.",
    "how to download report" = "Use the Download HTML Report button. The report.Rmd must exist in the app directory and will be rendered to HTML."
  )
  
  summarize_qc <- function() {
    if (is.null(rv$seurat_obj)) return("No Seurat object available. Run analysis first.")
    md <- rv$seurat_obj@meta.data
    if (!all(c("nFeature_RNA", "nCount_RNA", "percent.mt") %in% colnames(md))) return("QC fields missing.")
    s <- function(x) paste0("min=", signif(min(x, na.rm=TRUE),3), ", median=", signif(median(x, na.rm=TRUE),3), ", mean=", signif(mean(x, na.rm=TRUE),3), ", max=", signif(max(x, na.rm=TRUE),3))
    n_total <- nrow(md); n_mt5 <- sum(md$percent.mt >= 5, na.rm=TRUE); n_lowfeat <- sum(md$nFeature_RNA <= 200, na.rm=TRUE)
    paste0("QC: ", n_total, " cells. nFeature_RNA: ", s(md$nFeature_RNA), ". nCount_RNA: ", s(md$nCount_RNA), ". percent.mt: ", s(md$percent.mt),
           ". Cells with percent.mt >=5%: ", n_mt5, " (", round(100 * n_mt5 / n_total,2), "%). Cells with nFeature_RNA <=200: ", n_lowfeat, " (", round(100 * n_lowfeat / n_total,2), "%).")
  }
  
  summarize_umap_clusters <- function() {
    if (is.null(rv$seurat_obj)) return("No Seurat object available. Run analysis first.")
    cl_col <- if ("seurat_clusters" %in% colnames(rv$seurat_obj@meta.data)) "seurat_clusters" else NULL
    if (is.null(cl_col)) return("No cluster labels.")
    tab <- table(rv$seurat_obj@meta.data[[cl_col]])
    paste("Clusters:", length(tab), "  Counts per cluster:", paste(names(tab), as.integer(tab), sep=":", collapse=", "))
  }
  
  summarize_top_markers <- function(cluster = 0, n = 10) {
    if (is.null(rv$markers)) return("No marker results available.")
    df <- rv$markers %>% filter(cluster == !!cluster) %>% arrange(p_val_adj, desc(avg_log2FC))
    if (nrow(df) == 0) return(paste0("No markers for cluster ", cluster))
    paste("Top markers:", paste(head(df$gene, n), collapse = ", "))
  }
  
  summarize_go <- function(cluster = 0, n = 5) {
    key <- paste0("cluster_", cluster)
    if (is.null(rv$go_results) || is.null(rv$go_results[[key]])) return(paste0("No GO for cluster ", cluster))
    df <- as.data.frame(rv$go_results[[key]])
    if (nrow(df) == 0) return(paste0("No significant GO terms for cluster ", cluster))
    paste("Top GO:", paste(head(df$Description, n), collapse = "; "))
  }
  
  match_custom_faq <- function(query) {
    df <- rv$custom_faq
    if (is.null(df) || nrow(df) == 0) return(NULL)
    q <- tolower(trimws(query))
    for (i in seq_len(nrow(df))) {
      if (grepl(tolower(df$question[i]), q, fixed = TRUE) || grepl(q, tolower(df$question[i]), fixed = TRUE)) return(df$answer[i])
    }
    keys <- tolower(df$question)
    d <- stringdist::stringdist(q, keys, method = "jw")
    best <- which.min(d)
    if (length(best) == 1 && d[best] < 0.30) return(df$answer[best])
    NULL
  }
  
  is_greeting <- function(q) {
    ql <- tolower(trimws(q))
    grepl("^(hi|hello|hey|greetings)\\b", ql)
  }
  
  find_local_answer <- function(query) {
    q <- tolower(trimws(query))
    if (is_greeting(query)) {
      return("Hi! I'm OmicsPilot. Ask me about QC, UMAP, markers, GO enrichment, or let me learn new FAQs!")
    }
    cf <- tryCatch(match_custom_faq(query), error = function(e) NULL)
    if (!is.null(cf)) return(cf)
    for (k in names(built_in_faq)) if (grepl(k, q, fixed = TRUE)) return(built_in_faq[[k]])
    if (grepl("qc|quality control|qc summary|summarize qc|describe qc", q)) return(summarize_qc())
    if (grepl("umap|clusters|cluster summary|cells per cluster|how many clusters", q)) return(summarize_umap_clusters())
    if (grepl("top markers|markers for cluster|top 5 markers", q) && !is.null(rv$markers)) {
      m <- regmatches(q, regexpr("[0-9]+", q))
      cluster_num <- ifelse(length(m) > 0, as.numeric(m), 0)
      return(summarize_top_markers(cluster = cluster_num, n = 10))
    }
    if (grepl("go terms|go enrichment", q) && !is.null(rv$go_results)) {
      m <- regmatches(q, regexpr("[0-9]+", q))
      cluster_num <- ifelse(length(m) > 0, as.numeric(m), 0)
      return(summarize_go(cluster = cluster_num, n = 5))
    }
    NULL
  }
  
  # ---------- Custom FAQ UI ----------
  output$faqs_table <- renderDT({
    df <- rv$custom_faq
    if (is.null(df) || nrow(df) == 0) {
      datatable(data.frame(Note = "No custom FAQs trained"), options = list(dom='t'))
    } else {
      datatable(df, options = list(pageLength = 5), rownames = FALSE)
    }
  })
  
  observeEvent(input$teach_add, {
    q <- trimws(input$teach_question); a <- trimws(input$teach_answer)
    if (q == "" || a == "") {
      showModal(modalDialog("Both question and answer are required.", easyClose = TRUE)); return()
    }
    df <- rv$custom_faq
    new <- data.frame(question = q, answer = a, stringsAsFactors = FALSE)
    df <- bind_rows(new, df) %>% distinct(question, .keep_all = TRUE)
    rv$custom_faq <- df
    save_custom_faq(df)
    showModal(modalDialog("Added custom Q&A.", easyClose = TRUE))
    updateTextInput(session, "teach_question", value = ""); updateTextAreaInput(session, "teach_answer", value = "")
  })
  
  observeEvent(input$upload_faq_csv, {
    fi <- input$upload_faq_csv; if (is.null(fi)) return()
    df <- tryCatch(read.csv(fi$datapath, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) { showModal(modalDialog("Failed to read CSV.", easyClose = TRUE)); return() }
    names(df) <- tolower(names(df))
    if (!all(c("question","answer") %in% names(df))) { showModal(modalDialog("CSV must contain 'question' and 'answer' columns.", easyClose = TRUE)); return() }
    df <- df %>% dplyr::select(question, answer)
    combined <- bind_rows(df, rv$custom_faq) %>% distinct(question, .keep_all = TRUE)
    rv$custom_faq <- combined; save_custom_faq(combined)
    showModal(modalDialog(paste0("Uploaded and merged ", nrow(df), " FAQs. Total now: ", nrow(combined)), easyClose = TRUE))
  })
  
  output$download_faqs <- downloadHandler(
    filename = function() paste0("omics_custom_faqs_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- rv$custom_faq
      if (is.null(df) || nrow(df) == 0) write.csv(data.frame(question=character(0), answer=character(0)), file, row.names = FALSE)
      else write.csv(df, file, row.names = FALSE)
    }
  )
  
  # ---------- Chat interaction ----------
  output$ai_status <- renderUI({
    nf <- if (!is.null(rv$custom_faq)) nrow(rv$custom_faq) else 0
    api_status <- if (rv$api_enabled) "ENABLED" else "DISABLED"
    tagList(
      tags$div(style = "color:green;", paste0("Hybrid Assistant: ", api_status)), 
      tags$div(class="small-note", paste0("Custom FAQs: ", nf, " | Local + OpenAI"))
    )
  })
  
  observeEvent(input$ai_clear, { rv$chat <- list() })
  observeEvent(input$example_upload, { updateTextAreaInput(session, "ai_query", value = "How do I upload data?") })
  observeEvent(input$example_summary, { updateTextAreaInput(session, "ai_query", value = "Summarize QC metrics for my dataset") })
  observeEvent(input$example_markers, { updateTextAreaInput(session, "ai_query", value = "What are the top markers for cluster 0?") })
  
  output$ai_chat_ui <- renderUI({
    msgs <- rv$chat
    if (length(msgs) == 0) return(tags$p("No messages yet. Use the input box to ask questions."))
    tags$div(
      lapply(seq_along(msgs), function(i) {
        m <- msgs[[i]]
        cls <- if (m$role == "user") "ai-msg-user" else "ai-msg-assistant"
        source_tag <- if (m$role == "assistant" && !is.null(m$source)) {
          tags$span(style = "font-size:9px; color:#94a3b8; margin-left:8px;", paste0("[", m$source, "]"))
        } else { NULL }
        tags$div(class=cls, 
                 tags$strong(ifelse(m$role=="user","You: ","Assistant: ")), 
                 tags$span(m$text), 
                 source_tag,
                 tags$div(style="font-size:10px;color:#64748b;", format(m$time, "%Y-%m-%d %H:%M:%S")))
      })
    )
  })
  
  observeEvent(input$ai_send, {
    req(input$ai_query)
    user_text <- trimws(input$ai_query)
    if (user_text == "") return()
    
    rv$chat[[length(rv$chat) + 1]] <- list(role = "user", text = user_text, time = Sys.time(), source = NULL)
    
    # Try local answers first
    local_ans <- tryCatch(find_local_answer(user_text), error = function(e) NULL)
    if (!is.null(local_ans)) {
      rv$chat[[length(rv$chat) + 1]] <- list(role = "assistant", text = local_ans, time = Sys.time(), source = "local")
      updateTextAreaInput(session, "ai_query", value = "")
      return()
    }
    
    # Fall back to OpenAI if available
    if (rv$api_enabled) {
      openai_ans <- tryCatch(query_openai(user_text, rv$openai_key), error = function(e) NULL)
      if (!is.null(openai_ans)) {
        rv$chat[[length(rv$chat) + 1]] <- list(role = "assistant", text = openai_ans, time = Sys.time(), source = "OpenAI")
        updateTextAreaInput(session, "ai_query", value = "")
        return()
      }
    }
    
    # Final fallback
    rv$chat[[length(rv$chat) + 1]] <- list(role = "assistant", 
                                           text = "I couldn't generate an answer. Try teaching me with the 'Teach assistant' panel.", 
                                           time = Sys.time(), source = "fallback")
    updateTextAreaInput(session, "ai_query", value = "")
  })
  
  # ---------- Report rendering (HTML only) ----------
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("omics_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      req(rv$seurat_obj, rv$go_results, rv$markers)
      
      rmd_path <- "report.Rmd"
      if (!file.exists(rmd_path)) {
        err_html <- "<html><body><h2>report.Rmd not found</h2><p>Please add a 'report.Rmd' file to the app directory. See the app instructions.</p></body></html>"
        writeLines(err_html, file)
        showModal(modalDialog(
          title = "Missing report.Rmd",
          "report.Rmd not found in the app directory. A diagnostic HTML file was provided for download.",
          easyClose = TRUE
        ))
        return()
      }
      
      temp_work_dir <- tempdir()
      out_file <- file.path(temp_work_dir, paste0("omics_report_", Sys.Date(), ".html"))
      
      params <- list(
        seurat_obj = rv$seurat_obj,
        go_results = rv$go_results,
        markers = rv$markers
      )
      
      tryCatch({
        rmarkdown::render(
          input = rmd_path,
          output_format = "html_document",
          output_file = out_file,
          params = params,
          envir = new.env(parent = globalenv()),
          knit_root_dir = temp_work_dir,
          quiet = TRUE
        )
        if (!file.copy(out_file, file, overwrite = TRUE)) {
          stop("Failed to copy rendered report to download location")
        }
        showModal(modalDialog(
          title = "Report Generated",
          "Your HTML report has been generated and is ready to download.",
          easyClose = TRUE
        ))
      }, error = function(e) {
        diag_html <- paste0("<html><body><h2>Report generation failed</h2><pre>", htmltools::htmlEscape(e$message), "</pre><p>See R console for full traceback.</p></body></html>")
        writeLines(diag_html, file)
        showModal(modalDialog(
          title = "Report Generation Failed",
          paste0("Error: ", e$message, ". A diagnostic HTML file was provided for download."),
          easyClose = TRUE
        ))
      })
    }
  )
  
}

shinyApp(ui, server)