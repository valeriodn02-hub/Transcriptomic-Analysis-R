# RNA-Seq Differential Expression Analysis Pipeline

# Project Overview
This project contains an R-based pipeline for the analysis of bulk RNA-Seq data. The workflow proceeds from raw count importation (Salmon output) to functional enrichment analysis, identifying differentially expressed genes (DEGs) associated with specific experimental conditions.

**Key Objectives:**
* Data preprocessing and quality control.
* Variance stabilizing transformation for dimensionality reduction (PCA).
* Statistical testing for differential expression using generalized linear models.
* Functional enrichment analysis (Gene Ontology & Hallmark pathways).

# Tools & Libraries
* **Data Handling:** `SummarizedExperiment`, `dplyr`, `openxlsx`
* **Normalization & Statistics:** `DESeq2`, `edgeR`
* **Visualization:** `ggplot2`, `factoextra`, `pheatmap`, `RColorBrewer`
* **Enrichment:** `clusterProfiler`

---

# 1. Data Loading & Preprocessing
The pipeline starts by importing Salmon transcript quantification files and converting them to gene-level counts. Data is then structured into a `DGEList` object for edgeR analysis.

```r
library(edgeR)
library(DESeq2)
library(SummarizedExperiment)

# Load gene counts (Example structure)
# gene_counts is a SummarizedExperiment object from Salmon
matrice_conteggi <- assay(gene_counts, "counts")

# Round counts to integers (required for DESeq2/edgeR)
df_conteggi <- floor(as.data.frame(matrice_conteggi))

# Create DGEList object combining counts and gene metadata
DGE_list0 <- DGEList(counts = df_conteggi, genes = df_dati_geni)

# Assign groups to samples
DGE_list0$samples$group <- as.factor(pData$groupsGep)

# Filtering low-expression genes to increase statistical power
keep <- filterByExpr(y = DGE_list0$counts, group = DGE_list0$samples$group)
DGE_list0 <- DGE_list0[keep, , keep.lib.sizes=FALSE]

#2. Dimensionality Reduction (PCA)
To assess sample grouping and identify potential batch effects, I utilized Variance Stabilizing Transformation (VST) followed by Principal Component Analysis (PCA)

library(FactoMineR)
library(factoextra)

# Variance Stabilizing Transformation using DESeq2
dds <- DESeqDataSetFromMatrix(countData = DGE_list0$counts,
                              colData = DGE_list0$samples,
                              design = ~group)
vst_dds <- vst(dds, blind = TRUE)
vst_assay <- assay(vst_dds)

# Selection of top variable genes for PCA
varianza_geni <- apply(vst_assay, 1, var)
top_genes <- names(sort(varianza_geni, decreasing=TRUE))[1:500]
pca_matrix <- t(vst_assay[top_genes, ])

# Compute PCA
res.pca <- PCA(pca_matrix, graph = FALSE)

# Visualization of Individuals (Samples)
fviz_pca_ind(res.pca,
             geom.ind = c("point", "text"),
             repel = TRUE,
             habillage = pData$groupsGep, # Color by group
             addEllipses = FALSE,
             pointsize = 3) +
  theme_minimal()

#3. Differential Expression Analysis (edgeR)
Differential expression was performed using a Generalized Linear Model (GLM) approach (Likelihood Ratio Test)

# Create design matrix
design <- model.matrix(~ 0 + group, data = DGE_list0$samples)

# Estimate Dispersion
dglist <- estimateDisp(DGE_list0, design = design, robust = TRUE)

# Fit GLM
fit <- glmFit(dglist, design = design, robust = TRUE)

# Define Contrasts (e.g., Condition A vs Condition B)
my_contrasts <- makeContrasts(
  Contrast1 = groupHunk_dox_pos - groupVector_dox_pos,
  levels = design
)

# Perform Likelihood Ratio Test
lrt <- glmLRT(fit, contrast = my_contrasts[,"Contrast1"])

# Extract Top Tags (DEGs)
top_tags <- topTags(lrt, n = nrow(fit))
DEGs <- top_tags$table

# Filter significant genes (FDR < 0.05 and |logFC| > 0.58)
sig_genes <- DEGs[DEGs$FDR < 0.05 & abs(DEGs$logFC) > 0.58, ]

#4. Visualization: Volcano Plot & Heatmap
Custom visualization of the results using ggplot2 and pheatmap

library(ggplot2)
library(ggrepel)

# Volcano Plot with custom thresholds
ggplot(DEGs, aes(x = logFC, y = -log10(FDR), color = Threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot: HUNK vs Vector", x = "log2 Fold Change", y = "-log10 FDR") +
  geom_text_repel(data = top_genes_labels, aes(label = gene_name), size = 3)

# Heatmap of Top DEGs
library(pheatmap)
pheatmap(vst_assay[rownames(sig_genes)[1:50], ], 
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         annotation_col = sample_metadata,
         show_rownames = TRUE)

#5. Functional Enrichment Analysis
Performing Over-Representation Analysis (ORA) using clusterProfiler on Gene Ontology (GO) and Hallmark pathways

library(clusterProfiler)

# Prepare gene lists
genes_up <- rownames(sig_genes[sig_genes$logFC > 0, ])
genes_down <- rownames(sig_genes[sig_genes$logFC < 0, ])

# Enrichment Function Loop (for multiple ontologies)
enrich_results <- enricher(
    gene = genes_up,
    universe = all_detected_genes,
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    TERM2GENE = term2gene_category # Loaded from .gmt files
)

# Visualization of Enrichment
barplot(enrich_results, showCategory=10) + ggtitle("Enrichment Analysis - UpRegulated")

Valerio Di Nunzio, Bioinformatics & Data Analysis
