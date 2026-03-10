#Requirements
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")
if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy")
if (!requireNamespace("oligo", quietly = TRUE)) BiocManager::install("oligo")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!requireNamespace("impute", quietly = TRUE)) BiocManager::install("impute")
if (!requireNamespace("preprocessCore", quietly = TRUE)) BiocManager::install("preprocessCore")
if (!requireNamespace("flashClust", quietly = TRUE)) install.packages("flashClust")
if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) install.packages("dynamicTreeCut")

#Data collection
geo_id <- "GSE18520"
gse <- getGEO(geo_id, GSEMatrix = TRUE)

phenoData_full <- pData(phenoData(gse[[1]]))
phenoData <- phenoData_full[, 1, drop = FALSE]
head(phenoData)
pheno <- phenoData
pheno$SampleID <- rownames(pheno)
pheno <- pheno[, c("SampleID", colnames(pheno)[1])]

#Load Packages
library(WGCNA)
library(affy)
library(oligo)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db) # Replace with your appropriate chip annotation package
library(flashClust)
library(RColorBrewer)

# Enable multi-threading for WGCNA
set.seed(12345)
enableWGCNAThreads()

# Directory set and file review
cel_dir <- "F:/minor project/GSE18520_RAW"
cel_files <- list.files(cel_dir, pattern = "\\.CEL$|\\.cel$", full.names = TRUE)

if(length(cel_files) == 0) {
  stop("No CEL files found in the specified directory.")
} else {
  cat("Found", length(cel_files), "CEL files.\n")
}

# Read CEL files
raw_data <- ReadAffy(filenames = cel_files)

# Get and print CDF name
cdf_name <- raw_data@cdfName
cat("CDF Name detected:", cdf_name, "\n")

# Platform check
if(cdf_name == "HG-U133_Plus_2") {
  cat("Detected Affymetrix HG-U133_Plus_2 platform.\n")
  platform <- "affymetrix"
} else {
  cat("Unknown platform, defaulting to Affymetrix processing.\n")
  platform <- "affymetrix"
}

# Read CEL files (Affymetrix HG-U133 Plus 2.0 → affy is correct)
raw_data <- ReadAffy(filenames = cel_files)

# Extract sample names
sample_names <- sampleNames(raw_data)

# Quality control – boxplot of raw intensities
pdf("QC_Raw_Data_Boxplot.pdf", width = 12, height = 8)
boxplot(raw_data,
        main = "Raw Data Intensity",
        cex.axis = 0.7,
        las = 2)
dev.off()

# RMA normalization
cat("Performing RMA normalization...\n")
eset <- affy::rma(raw_data)

# Expression matrix
expr_matrix <- exprs(eset)

# Get probe to gene mapping
cat("Mapping probes to gene symbols...\n")
probe_ids <- rownames(expr_matrix)

#annotation package
if(exists("pkg_used") && pkg_used == "oligo") {
  # For newer arrays, try to get annotation from featureData
  annotations <- annotation(raw_data)
  cat("Array annotation: ", annotations, "\n")
  
# Map probes to gene symbols using hgu133plus2.db
  gene_symbols <- mapIds(hgu133plus2.db, 
                         keys = probe_ids,
                         column = "SYMBOL", 
                         keytype = "PROBEID",
                         multiVals = "first")
} else {
  
# For older Affymetrix arrays
  gene_symbols <- mapIds(hgu133plus2.db, 
                         keys = probe_ids,
                         column = "SYMBOL", 
                         keytype = "PROBEID",
                         multiVals = "first")
}

# Debugging the Row Mismatch 
dim(probe_ids)       # Should return (54675, 1) or (54675,)
dim(gene_symbols)    # Should return (54675, 1) or (54675,)
dim(expr_matrix)     # Should return (44650, X)

#Match expr_matrix to probe_ids and gene_symbols
common_probes <- intersect(probe_ids, rownames(expr_matrix))  # Get shared probe IDs
expr_matrix <- expr_matrix[common_probes, , drop = FALSE]
probe_ids <- probe_ids[common_probes]
gene_symbols <- gene_symbols[common_probes]

# Combine into a data frame
expr_data <- data.frame(ProbeID = probe_ids, 
                        GeneSymbol = gene_symbols,
                        expr_matrix,
                        stringsAsFactors = FALSE)

# Remove probes without gene symbols
expr_data <- expr_data[!is.na(expr_data$GeneSymbol), ]

# Collapse multiple probes to unique genes
cat("Collapsing multiple probes to unique genes...\n")

# Step 1: Compute mean expression per probe
expr_data$MeanExpr <- rowMeans(expr_data[, 3:(ncol(expr_data)-1)], na.rm = TRUE)
head(expr_data)

# Step 2: Sort by GeneSymbol and descending mean expression
# This ensures that for genes with multiple probes, the probe with highest mean is first
expr_data <- expr_data[order(expr_data$GeneSymbol, -expr_data$MeanExpr), ]

# Step 3: Keep only the first occurrence per gene
unique_genes <- !duplicated(expr_data$GeneSymbol)
expr_data <- expr_data[unique_genes, ]

# Remove temporary MeanExpr column
expr_data$MeanExpr <- NULL

# Set rownames to gene symbols
rownames(expr_data) <- expr_data$GeneSymbol

# Prepare final expression matrix
final_expr <- expr_data[, -c(1, 2)] # Remove ProbeID & GeneSymbol columns
colnames(final_expr) <- gsub("\\.CEL$", "", colnames(final_expr)) # Remove ".CEL" from GSM IDS
data <- final_expr

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

# 1. Load phenotype and expression data
pheno <- read.csv("pheno.csv")
data <- read.csv("Processed_Expression_Data.csv")

# 2. Format phenotype data
# Move sample IDs to rownames
samp2 <- pheno[,-1]
rownames(samp2) <- pheno[,1]
pheno <- samp2

# Remove unnecessary phenotype columns
pheno[2:43] <- list(NULL)

# 3. Quality control – gene & sample filtering--------------------------------
# Reload expression data with genes as rownames
data <- read.csv(
  "Processed_Expression_Data.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Identify good genes and samples
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

# 4. Outlier sample detection-------------------------------------------------
# Method 1: Hierarchical clustering
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# Method 2: PCA-based visualization
library(ggrepel)

pca <- prcomp(t(data))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(pca.dat)), size = 3) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %')) +
  theme_bw(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

# Keep phenotype and expression data aligned
colData <- pheno
data.subset <- data

# Verify sample name consistency
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

# Transpose expression matrix (samples as rows)
new.norm <- (t(data.subset))

# 5. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 30, by = 2))

# Define candidate soft-thresholding powers
sft <- pickSoftThreshold(new.norm,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
sft.data <- sft$fitIndices

# Visualize scale-free topology and connectivity
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Convert matrix to numeric
new.norm[] <- sapply(new.norm, as.numeric)

# 7. Module detection using blockwiseModules (TOM)-----------------------------
soft_power <- 10
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(new.norm,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
cor <- temp_cor

# 8. Module Eigengenes and visualization ------------------------------------
module_eigengenes <- bwnet$MEs

# Get number of genes for each module
table(bwnet$colors)

# Plot dendrogram with module colors
nGenes <- length(bwnet$dendrograms[[1]]$order)
print(nGenes)
print(length(bwnet$colors))
print(length(bwnet$unmergedColors))

colorVector <- cbind(
  bwnet$unmergedColors[bwnet$dendrograms[[1]]$order],
  bwnet$colors[bwnet$dendrograms[[1]]$order]
)

plotDendroAndColors(
  bwnet$dendrograms[[1]], 
  colorVector,
  c("unmerged", "merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05
)

# 9. Trait processing and module–trait correlation-----------------------------
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('Normal Ovary', title), 1, 0)) %>% 
  select(2)

# Binarize categorical variables
colData$title <- factor(colData$title, levels = c("Normal Ovary", "Ovarian Tumor"))
severity.out <- binarizeCategoricalColumns(colData$title,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

traits <- cbind(traits, severity.out)

# Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18],
             y = names(heatmap.data)[1:18],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# 10. Extract genes from blue module-----------------------------------
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'blue') %>% 
  rownames()

# 11. Intramodular analysis – hub gene identification --------------------

# Define numbers of genes and samples
nSamples <- nrow(new.norm)
nGenes <- ncol(new.norm)

# Calculate the module membership and the associated p-values
module.membership.measure <- cor(module_eigengenes, new.norm, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:18,1:18]

# Calculate the gene significance and associated p-values
gene.signf.corr <- cor(new.norm, traits$`data.Ovarian Tumor.vs.all`, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(63)


# Load expression matrix
expr <- read.csv("Processed_Expression_Data.csv", row.names = 1, check.names = FALSE)

# Load metadata
metadata <- read.csv("Sample_Info.csv", stringsAsFactors = FALSE)

# Ensure column names in expr exactly match GSM_ID
all(colnames(expr) %in% metadata$Samples)
# Should return TRUE

# Match column order
metadata_ordered <- metadata[match(colnames(expr), metadata$Samples), ]

# Check
all(metadata_ordered$Samples == colnames(expr))  # Should be TRUE

# Make group factor
group <- factor(metadata_ordered$Type)

library(limma)

# Design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit model
fit <- lmFit(expr, design)

# Make contrast
contrast.matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEGs
deg <- topTable(fit2, coef="Tumor_vs_Normal", number=Inf, adjust.method="fdr")

# Add regulation column
deg$Regulation <- "NotSig"
deg$Regulation[deg$logFC > 1 & deg$adj.P.Val < 0.05] <- "Up"
deg$Regulation[deg$logFC < -1 & deg$adj.P.Val < 0.05] <- "Down"

# Summary
table(deg$Regulation)

deg_plot <- deg
deg_plot$logFC <- -deg_plot$logFC   # Normal - Tumor

desired_genes <- c("BCLAF1", "MRFAP1", "MORF4L1", "HNRNPM", "NORAD", "LINC00899")

# MA plot
library(ggplot2)
library(ggrepel)

# Prepare MA plot data
ma_plot <- deg
ma_plot$Significant <- "NotSig"

ma_plot$Significant[ma_plot$logFC > 1 & ma_plot$adj.P.Val < 0.05] <- "Up"
ma_plot$Significant[ma_plot$logFC < -1 & ma_plot$adj.P.Val < 0.05] <- "Down"

# Add flag for desired genes
ma_plot$Highlight <- ifelse(rownames(ma_plot) %in% desired_genes, "Target", "Other")

# MA plot
ggplot(ma_plot, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
  
  # Highlight target genes
  geom_point(
    data = subset(ma_plot, Highlight == "Target"),
    color = "black",
    size = 3
  ) +
  
  # Add labels to target genes
  geom_text_repel(
    data = subset(ma_plot, Highlight == "Target"),
    aes(label = rownames(subset(ma_plot, Highlight == "Target"))),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(values = c(
    "Down" = "red",
    "Up" = "blue",
    "NotSig" = "grey"
  )) +
  
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  
  labs(
    title = "MA Plot Showing Differential Expression (Normal vs Tumor)",
    x = "Average log2 Expression",
    y = "log2 Fold Change (Normal − Tumor)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Volcano plot
library(ggplot2)
library(ggrepel)

# Prepare volcano plot data
vol_plot <- deg
vol_plot$Significant <- "NotSig"

vol_plot$Significant[vol_plot$logFC > 1 & vol_plot$adj.P.Val < 0.05] <- "Up"
vol_plot$Significant[vol_plot$logFC < -1 & vol_plot$adj.P.Val < 0.05] <- "Down"

# Flag desired genes
vol_plot$Highlight <- ifelse(rownames(vol_plot) %in% desired_genes, "Target", "Other")

vol_plot$Gene <- rownames(vol_plot)

# Volcano plot
ggplot(vol_plot, aes(x = logFC, y = -log10(adj.P.Val))) +
  
  # Main points
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
  
  # Highlight target genes
  geom_point(
    data = subset(vol_plot, Highlight == "Target"),
    color = "black",
    size = 3
  ) +
  
  # Labels for highlighted genes
  geom_text_repel(
    data = subset(vol_plot, Highlight == "Target"),
    aes(label = Gene),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  # Threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  # SAME colors as MA plot
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "NotSig" = "grey"
  )) +
  
  labs(
    title = "Volcano Plot Showing Differential Expression (Normal vs Tumor)",
    x = "log2 Fold Change (Normal − Tumor)",
    y = "-log10 Adjusted P-value"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#Heatmap
#Genes with an adjusted P-value (Benjamini–Hochberg FDR) < 0.05 and an absolute log2 fold change > 1 were considered significantly differentially expressed. Among these, the top 50 genes were selected based on the smallest adjusted P-values, representing the most statistically robust differential expression between normal and tumor samples.These genes were visualized using a heatmap with row-wise Z-score normalization (standardized across samples, independently of other genes.)
library(pheatmap)

# 1. Select top 50 DEGs (significant only)
top50_genes <- deg %>%
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
  dplyr::arrange(adj.P.Val) %>%
  head(50)

# 2. Extract expression matrix for these genes
heatmap_mat <- expr[rownames(top50_genes), ]

# 3. Z-score scaling (row-wise)
heatmap_mat <- t(scale(t(heatmap_mat)))

# 4. Create sample annotation
annotation_col <- data.frame(
  Group = metadata_ordered$Type
)
rownames(annotation_col) <- metadata_ordered$Samples

# 5. Define annotation colors
ann_colors <- list(
  Group = c(
    Normal = "#1f78b4",
    Tumor  = "#e31a1c"
  )
)

# 6. Plot heatmap
pheatmap(
  heatmap_mat,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Heatmap of Top 50 Differentially Expressed Genes"
)

1. **KEGG top 10**

\# --- Packages ---

library(dplyr)

library(tidyr)

library(stringr)

library(readxl)   # to read Excel files

library(writexl)  # to save output



\# --- Read KEGG file ---

infile <- "KEGG.xlsx"   # <-- your uploaded file

df <- read\_excel(infile)



\# Clean columns (keep names like "P-value")

names(df) <- trimws(names(df))



df <- df %>%

&nbsp; mutate(

&nbsp;   Term   = str\_trim(as.character(Term)),

&nbsp;   Genes  = str\_trim(as.character(Genes)),

&nbsp;   p\_num  = suppressWarnings(as.numeric(gsub("\[^0-9eE.+-]", "", `P-value`))),

&nbsp;   term\_l = tolower(Term)

&nbsp; )



\# --- Regex patterns for your target KEGG terms ---

patterns <- list(

&nbsp; "Pathways of neurodegeneration - multiple diseases" =

&nbsp;   "pathways.\*neurodegeneration.\*multiple.\*disease",

&nbsp; "MAPK signaling pathway" =

&nbsp;   "mapk.\*signaling",

&nbsp; "Calcium signaling pathway" =

&nbsp;   "calcium.\*signaling",

&nbsp; "Neuroactive ligand-receptor interaction" =

&nbsp;   "neuroactive.\*ligand.\*receptor.\*interaction",

&nbsp; "Chemokine signaling pathway" =

&nbsp;   "chemokine.\*signaling",

&nbsp; "TNF signaling pathway" =

&nbsp;   "tnf.\*signaling",

&nbsp; "Rap1 signaling pathway" =

&nbsp;   "rap1.\*signaling",

&nbsp; "cAMP signaling pathway" =

&nbsp;   "camp.\*signaling",

&nbsp; "PI3K-Akt signaling pathway" =

&nbsp;   "pi3k.\*akt.\*signaling",

&nbsp; "cGMP-PKG signaling pathway" =

&nbsp;   "cgmp.\*pkg.\*signaling"

)



\# --- Pick exactly one row per target term ---

picked <- lapply(names(patterns), function(lbl) {

&nbsp; pat <- patterns\[\[lbl]]

&nbsp; m <- df %>% filter(str\_detect(term\_l, regex(pat, ignore\_case = TRUE)))

&nbsp; if (nrow(m) == 0) return(NULL)

&nbsp; m %>% arrange(p\_num) %>% slice(1)

})



picked <- picked\[!sapply(picked, is.null)]

sel <- bind\_rows(picked)



if (nrow(sel) == 0) stop("No KEGG terms matched. Check the Term column in KEGG.xlsx")



\# --- Expand Genes row-wise ---

final <- sel %>%

&nbsp; mutate(Gene = str\_split(Genes, ";")) %>%

&nbsp; unnest(Gene) %>%

&nbsp; mutate(Gene = str\_trim(Gene)) %>%

&nbsp; select(Term, Gene, `P-value`)



\# --- Save output ---

outfile <- "KEGG\_top10\_from\_screenshot\_expanded.xlsx"

write\_xlsx(final, outfile)



message("Wrote: ", outfile,

&nbsp;       " | rows: ", nrow(final),

&nbsp;       " | KEGG terms requested: 10 | KEGG terms matched: ", dplyr::n\_distinct(sel$Term))



\# Optional: print any missing terms

missing <- setdiff(names(patterns), unique(sel$Term))

if (length(missing) > 0) message("Not matched (by pattern): ", paste(missing, collapse = "; "))






2. **GO top 10**

\# --- Packages ---

library(dplyr)

library(tidyr)

library(stringr)



\# --- Read file (keep original column names like "P-value") ---

infile <- "GO\_Biological\_Process\_2025\_table.csv"  # <- change path if needed

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)



\# Clean columns and values

names(df) <- trimws(names(df))

df <- df %>%

&nbsp; mutate(

&nbsp;   Term   = str\_trim(as.character(Term)),

&nbsp;   Genes  = str\_trim(as.character(Genes)),

&nbsp;   # parse numeric p-values robustly

&nbsp;   p\_num  = suppressWarnings(as.numeric(gsub("\[^0-9eE.+-]", "", `P-value`))),

&nbsp;   term\_l = tolower(Term)

&nbsp; )



\# --- Regex patterns for the 10 terms (case-insensitive) ---

patterns <- list(

&nbsp; "Sensory Perception of Smell" =

&nbsp;   "sensory\\\\s+perception.\*smell",

&nbsp; "Cell-Cell Adhesion" =

&nbsp;   "cell\\\\s\*-?\\\\s\*cell.\*adhesion",

&nbsp; "Inflammatory Response" =

&nbsp;   "inflamm.\*response",

&nbsp; "Chemical Synaptic Transmission" =

&nbsp;   "chemical\\\\s+synaptic\\\\s+transmission",

&nbsp; "Positive Regulation of Transcription" =

&nbsp;   "positive.\*regulation.\*transcription",

&nbsp; "Adenylate Cyclase-Inhibiting G-Protein Coupled Receptor Signaling Pathway" =

&nbsp;   "adenylate\\\\s\*cyclase.\*g\\\\s\*-?\\\\s\*protein.\*coupled.\*receptor.\*signaling",

&nbsp; "Sensory Perception of Chemical Stimulus" =

&nbsp;   "sensory\\\\s+perception.\*chemical.\*stimulus",

&nbsp; "Anterograde Trans-Synaptic Signaling" =

&nbsp;   "anterograde.\*synaptic.\*signaling",

&nbsp; "Phospholipase C-Activating G-Protein Coupled Receptor Signaling Pathway" =

&nbsp;   "phospholipase\\\\s\*c.\*g\\\\s\*-?\\\\s\*protein.\*coupled.\*receptor.\*signaling",

&nbsp; "Generation of Neurons" =

&nbsp;   "generation.\*neurons"

)



\# --- Pick exactly one row per target term (smallest p-value if multiple) ---

picked <- lapply(names(patterns), function(lbl) {

&nbsp; pat <- patterns\[\[lbl]]

&nbsp; m <- df %>% filter(str\_detect(term\_l, regex(pat, ignore\_case = TRUE)))

&nbsp; if (nrow(m) == 0) return(NULL)

&nbsp; m %>% arrange(p\_num) %>% slice(1)  # keep the most significant match

})



picked <- picked\[!sapply(picked, is.null)]

sel <- bind\_rows(picked)



if (nrow(sel) == 0) stop("No matching terms found. Check the 'Term' values in your CSV.")



\# --- Expand Genes row-wise ---

final <- sel %>%

&nbsp; mutate(Gene = str\_split(Genes, ";")) %>%

&nbsp; unnest(Gene) %>%

&nbsp; mutate(Gene = str\_trim(Gene)) %>%

&nbsp; select(Term, Gene, `P-value`)



\# --- Save ---

outfile <- "GO\_BP\_top10\_from\_screenshot\_expanded.csv"

write.csv(final, outfile, row.names = FALSE)



message("Wrote: ", outfile,

&nbsp;       " | rows: ", nrow(final),

&nbsp;       " | terms requested: 10 | terms matched: ", dplyr::n\_distinct(sel$Term))



\# Optional: print any terms that didn't match

missing <- setdiff(names(patterns), unique(names(picked)))

if (length(missing) > 0) message("Not matched (by pattern): ", paste(missing, collapse = "; "))



3. **Get exp data only for lncrnas n mrnas each**

library(readr)
library(dplyr)

expr <- read_csv("Processed_Expression_Data.csv")

# Ensure first column is gene names
colnames(expr)[1] <- "Gene"

# Standardize gene symbols
expr$Gene <- toupper(expr$Gene)

rna_map <- read_csv("RNAs only blue.csv")

lnc_genes  <- rna_map$lncRNAs %>% na.omit() %>% toupper()
mrna_genes <- rna_map$mRNAs   %>% na.omit() %>% toupper()

lnc_expr <- expr %>% filter(Gene %in% lnc_genes)
mrna_expr <- expr %>% filter(Gene %in% mrna_genes)

write_csv(lnc_expr,  "Blue_module_lncRNA_expression.csv")
write_csv(mrna_expr, "Blue_module_mRNA_expression.csv")


4. **If u want pvals for lncrna-mrna other than corr vals (getting in coexp.py)**
# Load libraries

library(Hmisc)   # for correlation with p-values

library(dplyr)   # for data wrangling



\# Step 1: Read input CSV files

lnc\_exp <- read.csv("lncRNA\_expression.csv", row.names = 1)   # rows=lncRNAs, cols=samples

mrna\_exp <- read.csv("mRNA\_expression.csv", row.names = 1)    # rows=mRNAs, cols=samples



\# Step 2: Make sure samples (columns) match in order

common\_samples <- intersect(colnames(lnc\_exp), colnames(mrna\_exp))

lnc\_exp <- lnc\_exp\[, common\_samples]

mrna\_exp <- mrna\_exp\[, common\_samples]



\# Step 3: Compute correlations

res <- rcorr(t(as.matrix(lnc\_exp)), t(as.matrix(mrna\_exp)), type = "pearson")



\# Extract correlation and p-value matrices

cor\_matrix <- res$r\[1:nrow(lnc\_exp), (nrow(lnc\_exp)+1):nrow(res$r)]

p\_matrix   <- res$P\[1:nrow(lnc\_exp), (nrow(lnc\_exp)+1):nrow(res$P)]



\# Step 4: Convert to edge list

edges <- data.frame(

&nbsp; lncRNA   = rep(rownames(lnc\_exp), each = nrow(mrna\_exp)),

&nbsp; mRNA     = rep(rownames(mrna\_exp), times = nrow(lnc\_exp)),

&nbsp; Correlation = as.vector(cor\_matrix),

&nbsp; Pvalue      = as.vector(p\_matrix)

)



\# Step 5: Filter significant correlations

edges\_filtered <- edges %>%

&nbsp; filter(abs(Correlation) >= 0.7 \& Pvalue <= 0.05)   # you can adjust thresholds



\# Step 6: Save results

write.csv(edges\_filtered, "lncRNA\_mRNA.csv", row.names = FALSE)



\# Optional: also save full correlation results

write.csv(edges, "lncRNA\_mRNA\_all\_correlations.csv", row.names = FALSE)




5. **lncrnas top 10**

library(readr)
library(dplyr)

# Load edge list generated for Cytoscape
edges <- read_csv("mRNA_lncRNA_network_for_cytoscape.csv")

# Check structure (optional but useful)
# colnames(edges) should be: Source, Target, Correlation

# Count how many mRNAs each lncRNA connects with
lnc_degree <- edges %>% 
  group_by(Target) %>%
  summarise(
    Degree = n(),
    MeanCorrelation = mean(abs(Correlation), na.rm = TRUE)
  ) %>%
  arrange(desc(Degree))

# Top 10 lncRNAs by number of connections
top10_lnc <- head(lnc_degree, 10)

print(top10_lnc)

# Save results
write.csv(top10_lnc, "lncRNAs_mRNA_Top10.csv", row.names = FALSE)




6. **mirnas top 10**

# ------------------ Packages ------------------
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(janitor)

# ------------------ Input ------------------
infile <- "miRTarBase_2017_table.csv"
df <- read_csv(infile)

# ------------------ Clean column names ------------------
df <- df %>% clean_names()
# Columns now:
# term | overlap | p_value | adjusted_p_value | old_p_value |
# old_adjusted_p_value | odds_ratio | combined_score | genes

# ------------------ miRNAs of interest ------------------
mirnas <- c(
  "hsa-miR-615-3p",
  "hsa-miR-16-5p",
  "hsa-miR-92a-3p",
  "hsa-miR-484",
  "hsa-let-7b-5p",
  "hsa-miR-193b-3p",
  "hsa-miR-186-5p",
  "hsa-miR-320a",
  "hsa-miR-877-3p",
  "hsa-miR-30a-5p"
)

# ------------------ Pick best row per miRNA ------------------
picked <- lapply(mirnas, function(mir) {
  m <- df %>%
    filter(str_detect(tolower(term), tolower(mir)))

  if (nrow(m) == 0) return(NULL)

  m %>% arrange(adjusted_p_value) %>% slice(1)
})

picked <- picked[!sapply(picked, is.null)]
sel <- bind_rows(picked)

# ------------------ Expand genes ------------------
cyto_edges <- sel %>%
  separate_rows(genes, sep = ";|,") %>%
  mutate(genes = str_trim(genes)) %>%
  transmute(
    source = term,                # miRNA
    target = genes,               # mRNA
    p_value = adjusted_p_value
  )

# ------------------ Save CSV ------------------
outfile <- "miRNA_10_Cytoscape_edges.csv"
write_csv(cyto_edges, outfile)

# ------------------ Report ------------------
message("Saved: ", outfile)
message("Total edges: ", nrow(cyto_edges))
message("miRNAs matched: ", n_distinct(cyto_edges$source))

# ------------------ Missing miRNAs ------------------
missing <- setdiff(mirnas, unique(sel$term))
if (length(missing) > 0)
  message("Not matched: ", paste(missing, collapse = ", "))


7. **Scores assign to each module**

library(readr)

library(readxl)

library(dplyr)

library(stringr)



\# === Load data ===

gene\_sign\_corr <- read\_csv("gene sign corr.csv")       

trans\_mm <- read\_csv("trans\_MM.csv")                   

blue\_module <- read\_excel("Blue module.xlsx")          



\# === Prepare GS dataframe ===

gs\_df <- gene\_sign\_corr %>%

&nbsp; rename(Gene = 1, GS = 2) %>% 

&nbsp; mutate(Gene = str\_trim(toupper(Gene)))  # standardize



\# === Prepare MM dataframe ===

mm\_df <- trans\_mm %>%

&nbsp; rename(Gene = 1) %>%

&nbsp; mutate(Gene = str\_trim(toupper(Gene)))  # standardize



\# === Prepare Blue module gene list ===

blue\_genes <- blue\_module %>%

&nbsp; rename(Gene = 1) %>%

&nbsp; mutate(Gene = str\_trim(toupper(Gene)))  # standardize



\# === Merge GS and MM ===

merged\_df <- inner\_join(gs\_df, mm\_df, by = "Gene")



\# === Extract Blue module genes (no threshold) ===

blue\_all <- merged\_df %>%

&nbsp; filter(Gene %in% blue\_genes$Gene) %>%

&nbsp; select(Gene, GS, MM\_blue = MEblue) %>%

&nbsp; mutate(Module = "Blue")



\# === Save results ===

write\_csv(blue\_all, "blue\_module\_all\_genes.csv")



\# === Print preview ===

print(head(blue\_all))

cat("Blue module total genes:", nrow(blue\_all), "\\n")



8. **Module threshold set to find real hub genes**

library(dplyr)

library(readr)



\# ------------------------------

\# 1. Load your data

\# ------------------------------

blue <- read\_csv("blue\_module\_all\_genes.csv")

turq <- read\_csv("turquoise\_module\_all\_genes.csv")



\# Detect column names (in case they differ)

detect\_col <- function(df, pattern) {

&nbsp; col <- grep(pattern, names(df), ignore.case = TRUE, value = TRUE)

&nbsp; if (length(col) == 0) stop(paste("No column matching", pattern, "found"))

&nbsp; return(col\[1])

}



gene\_col\_b <- detect\_col(blue, "gene")

mm\_col\_b   <- detect\_col(blue, "MM|ModuleMembership")

gs\_col\_b   <- detect\_col(blue, "GS|GeneSignificance")



gene\_col\_t <- detect\_col(turq, "gene")

mm\_col\_t   <- detect\_col(turq, "MM|ModuleMembership")

gs\_col\_t   <- detect\_col(turq, "GS|GeneSignificance")



\# ------------------------------

\# 2. Apply thresholds

\# ------------------------------

blue\_filt <- blue %>%

&nbsp; filter(abs(.data\[\[mm\_col\_b]]) >= 0.6,

&nbsp;        abs(.data\[\[gs\_col\_b]]) >= 0.24)



turq\_filt <- turq %>%

&nbsp; filter(abs(.data\[\[mm\_col\_t]]) >= 0.6,

&nbsp;        abs(.data\[\[gs\_col\_t]]) >= 0.24)



\# ------------------------------

\# 3. Find overlap \& merge in desired format

\# ------------------------------

overlap\_genes <- intersect(blue\_filt\[\[gene\_col\_b]], turq\_filt\[\[gene\_col\_t]])



merged <- data.frame(Gene = overlap\_genes) %>%

&nbsp; left\_join(blue %>%

&nbsp;             select(Gene = all\_of(gene\_col\_b),

&nbsp;                    Blue\_MM = all\_of(mm\_col\_b),

&nbsp;                    Blue\_GS = all\_of(gs\_col\_b)),

&nbsp;           by = "Gene") %>%

&nbsp; left\_join(turq %>%

&nbsp;             select(Gene = all\_of(gene\_col\_t),

&nbsp;                    Turquoise\_MM = all\_of(mm\_col\_t),

&nbsp;                    Turquoise\_GS = all\_of(gs\_col\_t)),

&nbsp;           by = "Gene") %>%

&nbsp; arrange(Gene)



\# ------------------------------

\# 4. Save output

\# ------------------------------

write\_csv(merged, "Real\_Hub\_Genes\_MM0.6\_GS0.24.csv")



cat("Found", nrow(merged), "overlapping hub genes.\\n")





9. **Extract thes rest of mrnas from module after getting lncrnas from lncsea**

# Example: list of your 29 lncRNAs
lncRNAs <- c("LINC001", "LINC002", "LINC003", ..., "LINC029")  # replace with actual names

# Read full blue module file
blue <- read.csv("blue_module_all_genes.csv", stringsAsFactors = FALSE)

# Separate lncRNAs
blue_lnc <- blue[blue$Gene %in% lncRNAs, ]

# Separate mRNAs
blue_mrna <- blue[!blue$Gene %in% lncRNAs, ]

# Optional: save them
write.csv(blue_lnc, "Blue_module_lncRNAs.csv", row.names = FALSE)
write.csv(blue_mrna, "Blue_module_mRNAs.csv", row.names = FALSE)



10. **FDR**
blue <- read.csv("blue_module_all_genes.csv", stringsAsFactors = FALSE)

head(blue)

nSamples <- nrow(data)
blue$GS_pvalue <- corPvalueStudent(blue$GS, nSamples)
blue$GS_FDR <- p.adjust(blue$GS_pvalue, method = "BH")
hub_candidates <- subset(
  blue,
  GS_FDR < 0.05 & abs(MM_blue) > 0.6
)
blue$GS_pvalue[blue$GS_pvalue == 0] <- .Machine$double.xmin
blue$GS_FDR[blue$GS_FDR == 0] <- .Machine$double.xmin






11. **miRNA-mRNA overlap (mirTARBase)**
library(readr)
library(dplyr)

# Step 1: Read the files
mirna_targets <- read_csv("mirTARBase.csv")
blue_mrnas <- read_csv("Blue_module_mRNAs.csv")

# Step 2: Inspect columns
cat("Columns in miRTarBase file:\n")
print(colnames(mirna_targets))
cat("\nColumns in Blue module mRNAs file:\n")
print(colnames(blue_mrnas))

# Step 3: Find overlap
# Keep only miRNA, Target Gene, Support Type
common_mrnas <- mirna_targets %>%
  filter(`Target Gene` %in% blue_mrnas$Gene) %>%
  select(miRNA, `Target Gene`, `Support Type`)

# Step 4: Summary
cat("Number of overlapping interactions:", nrow(common_mrnas), "\n")
cat("Number of unique mRNAs:", length(unique(common_mrnas$`Target Gene`)), "\n")

# Step 5: Save overlap for Cytoscape
write_csv(common_mrnas, "miRTarBase_BlueModule_overlap.csv")
cat("✅ Overlap file saved as 'miRTarBase_BlueModule_overlap.csv'\n")





12. **mirna-mrna (Targetscan)**
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Step 1: Read files
targetscan <- read_csv("TargetScan_microRNA_2017_table.csv")
blue_mrnas <- read_csv("Blue_module_mRNAs.csv")

# Step 2: Inspect columns (already done, but kept for safety)
cat("TargetScan columns:\n")
print(colnames(targetscan))

cat("\nBlue module mRNA columns:\n")
print(colnames(blue_mrnas))

# Step 3: Separate comma-separated genes into rows
targetscan_long <- targetscan %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = str_trim(Genes))  # remove extra spaces

mirna_mrna_edges <- targetscan_long %>%
  filter(Genes %in% blue_mrnas$Gene) %>%   # overlap only
  select(
    miRNA = mirna,
    mRNA  = Genes,
    P_value = `P-value`,
    Combined_Score = `Combined Score`
  )

# Step 6: Optional filtering (recommended)
# keep only strong miRNA signals
mirna_mrna_edges <- mirna_mrna_edges %>%
  filter(Adj_Pvalue < 0.05)

# Step 7: Summary
cat("Total miRNA–mRNA interactions:", nrow(mirna_mrna_edges), "\n")
cat("Unique miRNAs:", length(unique(mirna_mrna_edges$miRNA)), "\n")
cat("Unique mRNAs:", length(unique(mirna_mrna_edges$mRNA)), "\n")

# Step 8: Save for Cytoscape
write_csv(mirna_mrna_edges, "TargetScan_BlueModule_miRNA_mRNA_edges.csv")

cat("✅ File saved: TargetScan_BlueModule_miRNA_mRNA_edges.csv\n")




13. **Triplets**
library(dplyr)
library(readr)

# ===============================
# 1. Define the 18 hub genes
# ===============================
hub_genes <- c(
  "HNRNPM","KHDRBS1","HNRNPDL","MDH1","ASNSD1","BCLAF1",
  "DDX5","RABGGTB","SRSF2","MORF4L1","UNC50","MRFAP1",
  "ZC3H11A","ASH1L","CCT4","TOMM20","EEF1B2","EIF3E"
)

# ===============================
# 2. Load input files
# ===============================
lnc_mrna <- read_csv("mRNA_lncRNA_network_for_cytoscape.csv")
mirna_mrna <- read_csv("miRTarBase_BlueModule_miRNA_mRNA_edges.csv")

# ===============================
# 3. Standardize column names
# ===============================
lnc_mrna <- lnc_mrna %>%
  rename(
    lncRNA = lncRNA,
    mRNA   = mRNA
  )

mirna_mrna <- mirna_mrna %>%
  rename(
    miRNA = miRNA,
    mRNA  = mRNA
  )

# ===============================
# 4. Standardize gene symbols
# ===============================
lnc_mrna$mRNA <- toupper(trimws(lnc_mrna$mRNA))
lnc_mrna$lncRNA <- toupper(trimws(lnc_mrna$lncRNA))

mirna_mrna$mRNA <- toupper(trimws(mirna_mrna$mRNA))
mirna_mrna$miRNA <- toupper(trimws(mirna_mrna$miRNA))

hub_genes <- toupper(hub_genes)

# ===============================
# 5. Identify common mRNAs
# ===============================
common_mRNAs <- intersect(
  unique(lnc_mrna$mRNA),
  unique(mirna_mrna$mRNA)
)

# ===============================
# 6. Generate ALL triplets
# ===============================
triplets_all <- lnc_mrna %>%
  filter(mRNA %in% common_mRNAs) %>%
  inner_join(
    mirna_mrna %>% filter(mRNA %in% common_mRNAs),
    by = "mRNA"
  ) %>%
  select(lncRNA, mRNA, miRNA)

# ===============================
# 7. Extract triplets for the 18 hub genes
# ===============================
hub_triplets <- triplets_all %>%
  filter(mRNA %in% hub_genes)

# ===============================
# 8. Export files
# ===============================
write_csv(triplets_all,
          "ALL_mRNA_lncRNA_miRNA_triplets_miRTarBase.csv")

write_csv(hub_triplets,
          "Hub18_mRNA_lncRNA_miRNA_triplets_miRTarBase.csv")

# ===============================
# 9. Summary statistics
# ===============================
cat("===== TRIPLET SUMMARY =====\n")
cat("Total triplets formed:", nrow(triplets_all), "\n")
cat("Total hub triplet interactions:", nrow(hub_triplets), "\n\n")

cat("Unique hub genes forming triplets:\n")
print(unique(hub_triplets$mRNA))

cat("\nlncRNAs involved in hub triplets:",
    length(unique(hub_triplets$lncRNA)), "\n")

cat("miRNAs involved in hub triplets:",
    length(unique(hub_triplets$miRNA)), "\n")




14a. **=COUNTA(UNIQUE(B:B)) to find unique ones within a column**
14b. **=FILTER(A:A, ISNUMBER(MATCH(A:A, B:B, 0))) to find out the intersect btw 2 columns**




15. **Total unique mrnas in my cerna (old)**

# read csv
df <- read.csv("Book1.csv", stringsAsFactors = FALSE)

# combine both columns
all_genes <- c(df$Gene, df$mrna)

# count occurrences across both columns
gene_counts <- table(all_genes)

# extract genes appearing exactly once
unique_genes <- names(gene_counts[gene_counts == 1])

# print result
cat("Unique genes across Gene and mrna columns:\n")
cat(unique_genes, sep = "\n")
cat("Number of unique genes:", length(unique_genes), "\n")




16. **Candidate genes**
library(dplyr)
library(readr)

# ------------------------------
# 1. Load data
# ------------------------------
blue <- read_csv("FDR hub candidates blue.csv")
turq <- read_csv("FDR hub candidates turq.csv")

# ------------------------------
# 2. Apply thresholds SEPARATELY
# ------------------------------
blue_hubs <- blue %>%
  filter(abs(MM_blue) >= 0.6,
         abs(GS) >= 0.24) %>%
  transmute(
    Gene   = Gene,
    GS     = GS,
    MM     = MM_blue,
    Module = "Blue"
  )

turq_hubs <- turq %>%
  filter(abs(MM_turquoise) >= 0.6,
         abs(GS) >= 0.24) %>%
  transmute(
    Gene   = Gene,
    GS     = GS,
    MM     = MM_turquoise,
    Module = "Turquoise"
  )

# ------------------------------
# 3. Combine results (NO intersection)
# ------------------------------
final_hubs <- bind_rows(blue_hubs, turq_hubs) %>%
  arrange(Module, desc(abs(MM)), desc(abs(GS)))

# ------------------------------
# 4. Save output
# ------------------------------
write_csv(final_hubs, "Module_Hub_Genes_MM0.6_GS0.24.csv")

cat("Blue hubs:", nrow(blue_hubs), "\n")
cat("Turquoise hubs:", nrow(turq_hubs), "\n")
cat("Total hubs:", nrow(final_hubs), "\n")



17.**out of 115 which ones form triplets**
library(dplyr)
library(readr)
library(tidyr)

# ====== 1. Load hub genes ======
hub_genes <- c("AASDHPPT","ACADM","AIDA","ANAPC13","ANXA2","ARL5A","ARL6IP5",
               "ASH1L","ASNSD1","ATP5F1A","BCLAF1","BIRC2","BTF3","BZW1","CAPZA2",
               "CCT4","CD46","COX4I1","CSDE1","CUL5","CYFIP1","DDX1","DDX17","DDX5",
               "DDX50","DSTN","EEF1B2","EFCAB14","EIF2A","EIF2S3","EIF3E","EIF3F",
               "EIF4G2","EML4","FAM120A","GLOD4","GOLM2","GPBP1","H3-3B","HAPSTR1",
               "HBP1","HECTD1","HNRNPDL","HNRNPM","HSD17B12","IREB2","ITGAV","KHDRBS1",
               "LIN7C","MDH1","MICU2","MORF4L1","MRFAP1","NACA","NDFIP2","NPTN","NUDCD2",
               "OSBPL9","PCMT1","PDS5A","PIGY","PKN2","PMPCB","PPP1CC","PPP2CB","PPP2R5C",
               "PRKAR1A","PSMA1","PSMA2","PSMC1","PSMD10","PTPN12","RABGGTB","RBMX","RPF1",
               "RSL24D1","RTRAF","RWDD4","SAR1A","SARAF","SCP2","SEC11A","SEH1L","SELENOF",
               "SEPTIN7","SLMAP","SLTM","SMAD2","SMIM15","SPCS2","SRP72","SRSF2","STXBP3",
               "SUCLA2","SYPL1","TGOLN2","TM9SF2","TMEM14A","TMEM263","TMEM50A","TOMM20",
               "TRMT5","TSN","UBE2G1","UNC50","UQCRC2","USP33","VPS26A","VPS54","WASL",
               "ZC3H11A","ZFAND6","ZFR","ZMYND11")

# ====== 2. Load networks ======
lnc_mrna <- read_csv("mRNA_lncRNA_network_for_cytoscape.csv") # mRNA, lncRNA, Correlation
mirna_mrna <- read_csv("miRNA_mRNA_overlap.csv")               # miRNA, Target Gene, Support Type

# Standardize names
lnc_mrna$mRNA <- toupper(trimws(lnc_mrna$mRNA))
lnc_mrna$lncRNA <- toupper(trimws(lnc_mrna$lncRNA))
mirna_mrna$`Target Gene` <- toupper(trimws(mirna_mrna$`Target Gene`))
mirna_mrna$miRNA <- toupper(trimws(mirna_mrna$miRNA))
hub_genes <- toupper(trimws(hub_genes))

# ====== 3. Filter for mRNAs present in both networks ======
mRNAs_in_both <- intersect(unique(lnc_mrna$mRNA), unique(mirna_mrna$`Target Gene`))

lnc_mrna_filt <- lnc_mrna %>% filter(mRNA %in% mRNAs_in_both)
mirna_mrna_filt <- mirna_mrna %>% filter(`Target Gene` %in% mRNAs_in_both)

# ====== 4. Generate triplets ======
# Cross join lncRNAs and miRNAs per mRNA to form full triplets
triplets <- lnc_mrna_filt %>%
  inner_join(mirna_mrna_filt, by=c("mRNA"="Target Gene")) %>%
  select(lncRNA, mRNA, miRNA, Correlation, `Support Type`) %>%
  mutate(Is_Hub = ifelse(mRNA %in% hub_genes, "Yes", "No"))

# ====== 5. Summary ======
total_triplets <- nrow(triplets)
cat("Total triplets (lncRNA-mRNA-miRNA):", total_triplets, "\n")

cat("Unique lncRNAs:", length(unique(triplets$lncRNA)), "\n")
cat("Unique miRNAs:", length(unique(triplets$miRNA)), "\n")
cat("Unique mRNAs:", length(unique(triplets$mRNA)), "\n")
cat("Hub genes forming triplets:", length(unique(triplets$mRNA[triplets$Is_Hub=="Yes"])), "\n")
hub_gene_names <- unique(triplets$mRNA[triplets$Is_Hub == "Yes"])
cat("Hub genes forming triplets:\n")
print(hub_gene_names)

# ====== 6. Separate hub triplets ======
hub_triplets <- triplets %>% filter(Is_Hub == "Yes")
