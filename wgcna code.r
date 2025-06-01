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
#Phenodata fetch
phenoData_full <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1)]

#Load Packages
library(WGCNA)
library(affy)
library(oligo)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db) # Replace with your appropriate chip annotation package
library(flashClust)
library(RColorBrewer)
set.seed(12345)
enableWGCNAThreads()
library(affy)  # or library(oligo)

# Directory set and files review
cel_dir <- "C:/Users/ASUS/Downloads/GSE18520/"
cel_files <- list.files(cel_dir, pattern = "\\.CEL$|\\.cel$", full.names = TRUE)
if(length(cel_files) == 0) {
  stop("No CEL files found in the specified directory.")
} else {
  cat("Found", length(cel_files), "CEL files.\n")
}

# Install if needed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GEOquery")

# Check if files exist
cel_files <- list.files("C:/Users/ASUS/Downloads/GSE18520/", pattern = "*.CEL", full.names = TRUE)
raw_data <- ReadAffy(filenames = cel_files)

# Get the CDF name from your raw data
raw_data@cdfName

#check files are in affy n hg-u133-plus2 platform detected
if(raw_data@cdfName == "HG-U133_Plus_2") {
  cat("Detected Affymetrix HG-U133_Plus_2 platform.\n")
  platform <- "affymetrix"
} else {
  cat("Unknown platform, defaulting to Affymetrix processing.\n")
  platform <- "affymetrix"
}

# Read CEL files using affy
if(platform == "affymetrix") 
  # Read CEL files using affy or oligo package
  # Try first with oligo (for newer arrays)
  tryCatch({
    cat("Attempting to read with oligo package (for newer arrays)...\n")
    raw_data <- read.celfiles(cel_files)
    pkg_used <- "affy"
  }, error = function(e) {
    cat("Oligo failed, trying with affy package (for older arrays)...\n")
    raw_data <- ReadAffy(filenames = cel_files)
    pkg_used <- "affy"
  })

# Extract sample names
if(exists("pkg_used") && pkg_used == "affy") {
  sample_names <- sampleNames(raw_data)
} else {
  sample_names <- sampleNames(raw_data)
}

# RMA normalization
cat("Performing RMA normalization...\n")
if(exists("pkg_used") && pkg_used == "affy") {
  eset <- rma(raw_data) # Background correction, normalization, and summarization
} else {
  eset <- rma(raw_data) # Background correction, normalization, and summarization
}

# Convert to expression matrix
expr_matrix <- exprs(eset)

# Get probe to gene mapping
cat("Mapping probes to gene symbols...\n")
probe_ids <- rownames(expr_matrix)

#annotation package
if(exists("pkg_used") && pkg_used == "oligo") {
  # For newer arrays, try to get annotation from featureData
  annotations <- annotation(raw_data)
  cat("Array annotation: ", annotations, "\n")
  
  # For demonstration - using hgu133plus2.db
  # Replace with the correct annotation package for your array
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
expr_data <- data.frame(ProbeID = probe_ids, 
                        GeneSymbol = gene_symbols,
                        expr_matrix,
                        stringsAsFactors = FALSE)
expr_data <- expr_data[!is.na(expr_data$GeneSymbol), ]
cat("Collapsing multiple probes to unique genes...\n")
expr_data$MeanExpr <- rowMeans(expr_data[, 3:(ncol(expr_data)-1)], na.rm = TRUE)
head(expr_data)
expr_data <- expr_data[!is.na(expr_data$GeneSymbol), ]
cat("Collapsing multiple probes to unique genes...\n")

# Calculate mean expression for each probe
expr_data$MeanExpr <- rowMeans(expr_data[, 3:(ncol(expr_data)-1)], na.rm = TRUE)

# Sort by gene symbol and mean expression
expr_data <- expr_data[order(expr_data$GeneSymbol, -expr_data$MeanExpr), ]
unique_genes <- !duplicated(expr_data$GeneSymbol)
expr_data <- expr_data[unique_genes, ]
expr_data$MeanExpr <- NULL
rownames(expr_data) <- expr_data$GeneSymbol
final_expr <- expr_data[, -c(1, 2)]

# Transpose for WGCNA (genes as columns, samples as rows)
datExpr <- t(final_expr)
rownames(datExpr) <- colnames(final_expr)
write.csv(final_expr, file = "Processed_Expression_Data.csv")

#column to row names
samp2 <- pheno[,-1]
rownames(samp2) <- pheno[,1]
pheno <- samp2
pheno[2:43] <- list(NULL)

# 2. QC - outlier detection ------------------------------------------------
# Detect outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# pca - method 2
pca <- prcomp(t(data))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

#Samples to be excluded
#samples.to.be.excluded <- c('GSM4622647', 'GSM462650', 'GSM462652', 'GSM462651', 'GSM462648', 'GSM462649')
#data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

#colData <- pheno %>% 
  #filter(!row.names(.) %in% samples.to.be.excluded)
colData <- pheno
data.subset <- data

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))
new.norm <- (t(data.subset))

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(new.norm,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices
# visualization to pick power

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

# convert matrix to numeric
new.norm[] <- sapply(new.norm, as.numeric)
soft_power <- 10
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(new.norm,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
# Check dimensions
nGenes <- length(bwnet$dendrograms[[1]]$order)
print(nGenes)
print(length(bwnet$colors))
print(length(bwnet$unmergedColors))

# Make sure color vectors match the order length
colorVector <- cbind(
  bwnet$unmergedColors[bwnet$dendrograms[[1]]$order],
  bwnet$colors[bwnet$dendrograms[[1]]$order]
)

# Now plot
plotDendroAndColors(
  bwnet$dendrograms[[1]], 9
  colorVector,
  c("unmerged", "merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05
)

# create traits file - binarize categorical variables
# incase error: write.csv(coldata, "ww.csv"), remove spaces from excel file, coldata <- read.csv("ww.csv")
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('Normal Ovary', title), 1, 0)) %>% 
  select(2)

# binarize categorical variables

colData$title <- factor(colData$title, levels = c("Normal Ovary", "Ovarian Tumor"))

severity.out <- binarizeCategoricalColumns(colData$title,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)


traits <- cbind(traits, severity.out)
# visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18],
             y = names(heatmap.data)[1:18],
             col = c("blue1", "skyblue", "white", "pink", "red"))
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'grey') %>% 
  rownames()

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Define numbers of genes and samples
nSamples <- nrow(new.norm)
nGenes <- ncol(new.norm)

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, new.norm, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:18,1:18]

# Calculate the gene significance and associated p-values
write.csv(gene.signf.corr, "gene_sig.csv")
gene.signf.corr <- cor(new.norm, traits$`data.Ovarian Tumor.vs.all`, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(63)
