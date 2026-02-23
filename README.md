#  OmicsPilot

**Single-cell RNA-seq analysis with AI-powered bioinformatics assistant**

A Shiny web application for analyzing 10X Genomics single-cell RNA-seq data with integrated OpenAI assistance. Process quality control, clustering, marker gene detection, and GO enrichment—all with an intelligent assistant to answer your bioinformatics questions.

![R](https://img.shields.io/badge/R-4.0+-blue)
![Shiny](https://img.shields.io/badge/Shiny-Web-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

**[Try OmicsPilot Online](https://kanishka-2004.shinyapps.io/SCRNA/)** - No installation required!

---

##  Features

- **10X Genomics Support** - Direct import of matrix.mtx, features.tsv, barcodes.tsv files
- **Quality Control** - Violin plots, scatter plots, QC metrics (nFeature_RNA, nCount_RNA, percent.mt)
- **Advanced Analysis** - PCA, UMAP dimensionality reduction, clustering (resolution 0.5)
- **Marker Detection** - FindAllMarkers for differential expression analysis
- **GO Enrichment** - Gene Ontology enrichment analysis per cluster
- **Hybrid AI Assistant** - Local FAQ + OpenAI GPT-4o-mini for bioinformatics Q&A
- **Custom FAQ Training** - Teach the assistant domain-specific answers
- **Auto Report** - Generate professional HTML reports with plots and interpretations
- **Session Management** - Upload/download custom FAQs as CSV

---

## Quick Start

### Try Online (No Installation)

**No setup needed! Access the live app here:**  
**https://kanishka-2004.shinyapps.io/SCRNA/**

Just upload your 10X files and start analyzing!

---

### Run Locally

#### Requirements

- **R 4.0+**
- **RStudio** (recommended)
- **OpenAI API key** (for AI features)

### Installation

1. **Clone the repository:**
```bash
git clone https://github.com/kanishkapm/omics-pilot.git
cd omics-pilot
```

2. **Install R dependencies:**
```r
install.packages(c("shiny", "Seurat", "ggplot2", "dplyr", "DT", "rmarkdown", "httr2", "usethis", "shinycssloaders"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "EnhancedVolcano", "pheatmap"))
```

3. **Create `.Renviron` file with your OpenAI API key:**
```r
usethis::edit_r_environ()
```
Add this line:
```
OPENAI_API_KEY=sk-proj-YOUR_API_KEY_HERE
```
Save and restart R.

4. **Run the app:**
```r
setwd("~/omics-pilot")  # Or your project directory
shiny::runApp()
```

The app will open at `http://localhost:3838`

---

## 🌍 Online Demo

**Live Application:** https://kanishka-2004.shinyapps.io/SCRNA/

Try OmicsPilot without installing anything:
- Upload your 10X Genomics data
- Run analysis in the cloud
- Use AI assistant (with your own API key in `.Renviron`)
- Download reports

**Note:** The online version also requires a `.Renviron` file with your OpenAI API key for full AI features.

---

## Usage Guide

### Step 1: Upload Data

1. Navigate to the **"Input (10X files)"** panel
2. Upload three files:
   - **Matrix**: `matrix.mtx` or `matrix.mtx.gz`
   - **Features**: `features.tsv` or `features.tsv.gz`
   - **Barcodes**: `barcodes.tsv` or `barcodes.tsv.gz`
3. Click **"Run Analysis"**

### Step 2: View Quality Control

- **Data tab**: See uploaded file info and QC table
- **QC Plots tab**: Violin plots, scatter plots (nCount vs percent.mt, nCount vs nFeature)
- Check for:
  - nFeature_RNA > 200 (complexity)
  - percent.mt < 5% (mitochondrial contamination)

### Step 3: Explore Clusters

- **UMAP & Clusters tab**: Visualize cell clusters in 2D space
- Review cluster sizes and distribution

### Step 4: Analyze Markers

- **Markers & GO tab**: 
  - Volcano plot showing differentially expressed genes in Cluster 0
  - Heatmap of top markers across all clusters
  - GO enrichment terms

### Step 5: Talk to the AI Assistant

Ask questions in the **Chat panel**:
- "Summarize QC metrics for my dataset"
- "What are the top markers for cluster 0?"
- "What GO terms are enriched?"
- Any custom questions (uses OpenAI if local FAQ doesn't match)

### Step 6: Download Report

- Click **"Download HTML Report"** to generate a comprehensive report with:
  - QC metrics and visualizations
  - UMAP plot and cluster summary
  - Marker gene analysis
  - GO enrichment results
  - Interpretations

---

## AI Assistant

### Local Mode (No API calls)
- Greetings & basic questions
- Built-in FAQ responses
- QC summarization
- Cluster/marker queries

### OpenAI Mode (Requires API key)
- Complex bioinformatics questions
- Contextual explanations
- Custom interpretations
- Any question not in local FAQ

**Status Indicator:**
- **✓ OpenAI Ready** - API key loaded
- **✗ No API Key** - Using local FAQ only

### Teach the Assistant

1. Go to **"Teach assistant"** panel
2. Enter a question and answer
3. Click **"Add Q&A"**
4. Download FAQs as CSV for backup

---

## File Structure

```
omics-pilot/
├── README.md              # This file
├── app.R                  # Main Shiny application
├── report.Rmd             # HTML report template
├── .gitignore             # Git ignore file (protects .Renviron)
└── .Renviron              # API key (local only, NOT in GitHub)
```

---

## Security

**Important: Never commit `.Renviron` to GitHub**

- `.Renviron` contains your OpenAI API key
- It's automatically ignored by `.gitignore`
- Each user creates their own `.Renviron` locally
- Regenerate keys if accidentally exposed: https://platform.openai.com/api-keys

---

## Pipeline Details

### Analysis Steps

1. **Data Loading** - Read 10X Genomics sparse matrix format
2. **QC Metrics** - Calculate nFeature_RNA, nCount_RNA, percent.mt
3. **Filtering** - Remove cells: nFeature < 200, nFeature > 2500, percent.mt > 5%
4. **Normalization** - LogNormalize (scale factor 10,000)
5. **Variable Features** - VST method, 2000 features
6. **Scaling** - Center and scale data
7. **PCA** - Top 10 PCs
8. **UMAP** - 2D visualization (dims 1:10)
9. **Clustering** - Shared nearest neighbor (resolution 0.5)
10. **Marker Detection** - FindAllMarkers (logfc.threshold 0.25, min.pct 0.25)
11. **GO Enrichment** - BP terms (p-value cutoff 0.05)

---

## Examples

### Example 1: Basic Analysis
```
1. Upload 3 10X files
2. Click "Run Analysis"
3. View plots in each tab
4. Download report
```

### Example 2: Ask AI Questions
```
You: "What is percent.mt?"
Assistant: "percent.mt is the percentage of mitochondrial reads per cell. High values (>5%) suggest low-quality cells or stressed cells."

You: "Summarize QC metrics"
Assistant: "QC: 2500 cells. nFeature_RNA: min=200, median=1800, mean=1950, max=2500..."
```

### Example 3: Custom FAQ
```
Question: "What cell types are in my data?"
Answer: "Based on marker analysis, we see immune cells (CD45+), epithelial cells (EPCAM+), and fibroblasts (COL1A1+)"
```

---

## Troubleshooting

### Problem: API Key not loading
```r
# Check if key is set
Sys.getenv("OPENAI_API_KEY")
# Should return your key, not empty string
```

### Problem: report.Rmd not found
- Ensure `report.Rmd` is in the same directory as `app.R`
- File must be named exactly `report.Rmd`

### Problem: Package installation fails
```r
# Try installing one package at a time
install.packages("Seurat")  # May take 5-10 minutes
```

### Problem: File size too large
- Max file size: 1 GB per file
- Compress with gzip (.gz extension)
- Use sparse matrix format (.mtx)

---

## Resources

- **Seurat Documentation**: https://satijalab.org/seurat/
- **10X Genomics**: https://www.10xgenomics.com/
- **Gene Ontology**: http://geneontology.org/
- **clusterProfiler**: https://bioconductor.org/packages/clusterProfiler/
- **OpenAI API**: https://platform.openai.com/docs

---

## Contributing

Found a bug or have a feature request?

1. Fork the repository
2. Create a branch: `git checkout -b feature/your-feature`
3. Make changes
4. Commit: `git commit -m "Add feature X"`
5. Push: `git push origin feature/your-feature`
6. Open a Pull Request

---

## 📄 License

This project is licensed under the **MIT License** - see the LICENSE file for details.

---

## Authors

**Kanishka P M**
- GitHub: [@kanishkapm](https://github.com/kanishkapm)
  
**Indhushree Reddy V**

**Vahini**

For questions/support: Open an issue on GitHub

---

## Acknowledgments

- Built with [Shiny](https://shiny.posit.co/)
- Analysis powered by [Seurat](https://satijalab.org/seurat/)
- Enrichment analysis via [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)
- AI assistant by [OpenAI](https://openai.com/)

---

## If you find OmicsPilot useful, please star this repository!

---

**Last Updated:** 2025 
**Version:** 1.0.0
