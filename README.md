# RNA Sequencing Analysis in R
Basic analysis of publically available PPP RNA-seq dataset.

**Publication for Dataset:** Baum et al. (2022) *Pustular psoriasis: Molecular pathways and effects of spesolimab in generalized pustular psoriasis.* The Journal of Allergy and Clinical Immunology. 149(4):1402-1412.

## Overview
### Loading and Checking Data
1. Read in a **count matrix** (read counts for each sample) and a **meta table** (containing information about each sample). For example, if using a GEO dataset:
```
# Packages
library(GEOquery)
library(EnsDb.Hsapiens.v86)

# Read count and meta data
metaData <- getGEO(filename="~/Downloads/GSE000000_series_matrix.txt.gz")
countMatrix <- read.table('GSE000000_read_counts.txt', fill=TRUE, header=TRUE)

# Get gene symbols from biomart
gene_symbol <- select(EnsDb.Hsapiens.v86, 
                      keys=as.character(countMatrix$geneID) ,
                      keytype="GENEID",
                      columns=c("GENEID", "GENENAME"))

# Add gene symbols
countMatrixSym <- merge(countMatrix, gene_symbol, by.x='geneID', by.y='GENEID')

# Remove duplicate genes
countMatrixSym <- countMatrixSym[!duplicated(countMatrixSym[c('GENENAME')]), ]
```
2. Make sure the column names in the **count matrix** match the rownames in the **meta table**.
```
# Create column and count data
metaMatrix <- metaFilt %>% remove_rownames %>% column_to_rownames(var="Sample")
countMatrix <- countMatrixSym %>% remove_rownames %>% column_to_rownames(var="GENENAME")

# Remove gene column
countMatrix <- countMatrix %>% dplyr::select(-c(geneID))

# Make sure column names match row names
all(colnames(countMatrix) %in% rownames(metaMatrix))
all(colnames(countMatrix) == rownames(metaMatrix))
```
### Variance Partitioning
3. Next I can perform **variance partitioning** to see how each meta variable contributes to the variation in the counts of each gene. First the data needs to be normalised using **DESeq2**, and then fitExtractVarPartModel() can be used to model the effect of each variable on each gene across all the samples.
```
# Packages
library(variancePartition)
library(DESeq2)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = metaMatrix,
                              design = ~1)

# Median of ratios method for normalisation
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
```
### Differential Expression Analysis
4. Next I can use **DESeq2** to normalise the counts and perform **DE analysis**.
```
# Create deseq object
dds <- DESeqDataSetFromMatrix(countData=countMatrix, colData=metaData, design=~group+sex)
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast=c("group","Case","Control")) %>% as.data.frame() %>% rownames_to_column('ensembl_gene_id')
resultsNames(dds)

# Get results with shrinkage applied
resLFC <- lfcShrink(dds, coef="group_Control_vs_Case", type="apeglm")

# Top genes
resLFC %>% arrange(padj)
```
### Gene Set Enrichment Analysis (GSEA)
5. Finally, I like using **fgsea** to see if any pathways are up or downregulated.
```
# Packages
library(fgsea)

# Create gene list
resLFCList <- resLFC$log2FoldChange
names(resLFCList) <- resLFC$Gene

# Clean
resLFCList <- na.omit(resLFCList)
resLFCList <- sort(resLFCList, decreasing=T)
resLFCList <- resLFCList[!duplicated(names(resLFCList))]

# Get GO BP gmt file
GO_file <- "GOBP.gmt"
GO=fgsea::gmtPathways(GO_file)

# Run FGSEA
set.seed(1)
resGSEA <- fgsea::fgsea(pathways=GO,
                        stats=resLFCList,
                        minSize=10, 
                        maxSize=200)

# Get significant pathways
resGSEAFiltPADJ <- resGSEA %>%
  dplyr::filter(padj<=0.05) %>% 
  arrange(pval)
```





























