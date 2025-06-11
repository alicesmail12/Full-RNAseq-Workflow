# Full RNAseq Workflow in R
Basic analysis of publically available PPP RNAseq dataset.

**Publication for Dataset:** Baum et al. (2022) *Pustular psoriasis: Molecular pathways and effects of spesolimab in generalized pustular psoriasis.* The Journal of Allergy and Clinical Immunology. 149(4):1402-1412.

## Overview of Steps
### Loading and Checking Data
1. Read in a **count matrix** (read counts for each sample) and a **meta table** (containing information about each sample). For example, if using a GEO dataset use:
```
# Packages
library(GEOquery)

# Read count and meta data
metaData <- getGEO(filename="~/Downloads/GSE000000_series_matrix.txt.gz")
countMatrix <- read.table('GSE000000_read_counts.txt', fill=TRUE, header=TRUE) 
```
In the example script I use a meta table and count table from www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11144?query=GPP%20baum%20schmid. 
### Variance Partitioning
2. Next I can perform **variance partitioning** to see how each meta variable contributes to the variation in the counts of each gene. First the data needs to be normalised using **DESeq2**, and then fitExtractVarPartModel() can be used to model the effect of each variable on each gene across all the samples.
```
library(variancePartition)
library(DESeq2)

# DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countMatrix, colData=metaData, design=~1)

# Pre-processing suggested by variancePartition package
dds <- estimateSizeFactors(dds)

# Identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds) > 1) >= 0.5 * ncol(dds)

# Compute log2 Fragments Per Million
quantLog <- log2(fpm(dds)[isexpr, ] + 1)

# Define formula using each variable
form <- ~ (1 | group) + (1 | ageBracket) + (1 | sex) 

# Run variancePartition analysis
varPart <- fitExtractVarPartModel(quantLog, form, metaData)

# Visualise
vp <- sortCols(varPart)
plotPercentBars(vp[1:20, ])
plotVarPart(vp)
```
### Differential Expression Analysis
3. Next I can use DESeq2 to normalise the counts and perform **DE analysis**.
```
# Check sample names match between meta and count data
all(colnames(countMatrix) %in% rownames(metaData))
all(colnames(countMatrix) == rownames(metaData))

# Create deseq object
dds <- DESeqDataSetFromMatrix(countData=countMatrix, colData=metaData, design=~group+sex)
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast=c("group","Case","Control")) %>% as.data.frame() %>% rownames_to_column('ensembl_gene_id')
resultsNames(dds)

# Get results with shrinkage applied
resLFC <- lfcShrink(dds, coef="group_Control_vs_Case", type="apeglm")

# Get gene symbols
gene_symbol <- select(EnsDb.Hsapiens.v86, 
                      keys=as.character(res$ensembl_gene_id) ,
                      keytype="GENEID",
                      columns=c("GENEID", "GENENAME"))

# Add gene symbols
colnames(gene_symbol) <- c("ensembl_gene_id", "geneName")
resLFC <- left_join(resLFC, gene_symbol, by="ensembl_gene_id") %>% dplyr::filter(!is.na(geneName))
```
4. Now I can also make lots of plots!
```
# Plot one gene
plotCounts(dds, gene='ENSG00000017427', intgroup="group")

# Plot PCA of normalised gene counts
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("group"))
```





























