# Open libraries
library(tidyverse)
library(biomaRt)
library(DESeq2)
library(extrafont)
library(pheatmap)
library(EnsDb.Hsapiens.v86)
library(variancePartition)

# Set directory
workingDir = {wd}
setwd(workingDir)

## FUNCTIONS ###################################################################
makeVolcanoPlot <- function(resLFC){
  
  # Edit colours and labels
  resLFC$Label <- ifelse((resLFC$padj<0.05)&(abs(resLFC$log2FoldChange)>6), resLFC$Gene, NA)
  resLFC$Colour <- ifelse((resLFC$padj<0.05)&((resLFC$log2FoldChange)>0.5), 'Increased', 'None')
  resLFC$Colour <- ifelse((resLFC$padj<0.05)&((resLFC$log2FoldChange)<(-0.5)), 'Decreased', resLFC$Colour)
  resLFCPlot <- resLFC %>% dplyr::filter(!is.na(padj))
  
  # Plot
  plot <- ggplot(resLFCPlot, aes(x=log2FoldChange, y=-log10(padj), label=Label, colour=Colour, size=abs(log2FoldChange)*(-log10(padj))))+
    geom_hline(yintercept=-log10(0.05), colour='grey')+
    geom_vline(xintercept=c(0.5,-0.5), colour='grey')+
    geom_point(alpha=0.75)+
    ggrepel::geom_text_repel(family='Franklin Gothic Book', colour='black', size=4)+
    theme_classic()+
    theme(text=element_text(family='Franklin Gothic Book'))+
    scale_colour_manual(values=c('#ba5346', '#75a450', '#ebebeb'))+
    guides(colour='none', size='none')
  
  # Return
  return(plot)
}

# Plot DESEQ normalised counts
NormCountsPlot <- function(gene, dds=dds, metaData=metaData){
  
  # Extract normalised counts
  normalizedCounts <- counts(dds, normalized=TRUE) %>% as.data.frame() %>% rownames_to_column('Gene')
  
  # Filter for genes
  PlotCounts <- normalizedCounts %>% dplyr::filter(Gene %in% gene) 
  PlotCounts <- PlotCounts[-1] %>% t() %>% as.data.frame() %>% rownames_to_column('Patient')
  colnames(PlotCounts) <- c('Patient', 'Gene')

  # Add metadata
  PlotCounts <- merge(PlotCounts, metaData %>% dplyr::rename('Patient'='Source.Name'))
  
  # Plot
  plot <- ggplot(PlotCounts, aes(x=Condition, y=Gene, fill=Condition, colour=Condition))+
    geom_violin(alpha=0.5)+
    geom_point(size=4)+
    theme_classic()+
    theme(text=element_text(family='Franklin Gothic Book', colour='black', size=14),
          title=element_text(family='Franklin Gothic Book', colour='black', size=16))+
    scale_fill_manual(values=c('#90bdcf', '#cfc963', '#ba5346', '#75a450'))+
    scale_colour_manual(values=c('#90bdcf', '#cfc963', '#ba5346', '#75a450'))+
    labs(fill='Group', colour='Group', x='', y='NormalisedExpression')+
    ggtitle(paste(gene))+
    guides(colour='none', fill='none')
  
  # Return
  return(plot)
  
}

## DATA HANDLING ###############################################################
# Get meta data
meta <- read.table(paste0(workingDir, 'E-MTAB-11144.sdrf.txt'), sep='\t', fill=TRUE, header=TRUE) %>%
  dplyr::rename("Condition"="Characteristics.disease.",
                "Drug"="Factor.Value.compound.")

# Get relevant samples (L-PPP)
metaFilt <- meta %>%
  dplyr::filter(Condition %in% c('generalized pustular psoriasis', 'normal')) %>%
  dplyr::filter(Characteristics.organism.part. == 'skin') %>%
  dplyr::filter(Characteristics.sampling.site. %in% c('lesion', 'sole', 'palm'))

# Get sample names
sampleNames <- metaFilt$Source.Name

# Get both count files
countMatrixPPPHV <- read.table(paste0(workingDir, 'count_matrix_PPP_HV_allGenes.txt'), sep='\t', fill=TRUE, header=TRUE)
countMatrixGPP <- read.table(paste0(workingDir, 'counts_matrix_GPP_allGenes.txt'), sep='\t', fill=TRUE, header=TRUE) %>% dplyr::rename('geneID'='X')
countMatrix <- merge(countMatrixPPPHV, countMatrixGPP, by='geneID')

# Get only GPP and HV samples
countMatrixFilt <- countMatrix %>% dplyr::select('geneID', all_of(sampleNames))

# Get biomart
gene_symbol <- select(EnsDb.Hsapiens.v86, 
                      keys=as.character(countMatrix$geneID) ,
                      keytype="GENEID",
                      columns=c("GENEID", "GENENAME"))

# Add gene symbols
countMatrixFiltSym <- merge(countMatrixFilt, gene_symbol, by.x='geneID', by.y='GENEID')

# Remove duplicate genes
countMatrixFiltSym <- countMatrixFiltSym[!duplicated(countMatrixFiltSym[c('GENENAME')]), ]

# DESEQ ########################################################################
# Create column data
metaMatrix <- metaFilt %>%
  remove_rownames %>% 
  column_to_rownames(var="Source.Name")

# Create count data
countMatrix <- countMatrixFiltSym %>% 
  remove_rownames %>% 
  column_to_rownames(var="GENENAME")

# Remove gene column
countMatrix <- countMatrix %>% dplyr::select(-c(geneID))

# Make sure column names match row names
all(colnames(countMatrix) %in% rownames(metaMatrix))
all(colnames(countMatrix) == rownames(metaMatrix))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = metaMatrix,
                              design = ~Condition)

# Filter
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Median of ratios method for normalisation
dds <- estimateSizeFactors(dds)
normalisedCounts <- counts(dds, normalized=TRUE) %>% as.data.frame()

# vsd transformed data
vsd <- vst(dds, blind=F)

# Get PCA coords
pcaData <- plotPCA(vsd, intgroup=c("Condition", "Drug"), returnData=TRUE)

# Get variance explained
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA plot
plot <- ggplot(pcaData %>% rownames_to_column('Patient'), aes(PC1, PC2, color=Condition, shape=Drug)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_point(size=6, alpha=1)+
  theme_classic()+
  theme(text=element_text(family='Franklin Gothic Book', colour='black', size=14))+
  scale_color_manual(values=c('#90bdcf', '#ba5346', '#ba5346', '#75a450'))+
  scale_y_continuous(expand=c(0.1,0.1))+
  labs(colour='Condition', shape='Drug')

# Save
ggsave('GPP-PCA-Plot.png', plot, dpi=300)

# VARIANCE PARTITIONING ########################################################
# Identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds) > 1) >= 0 * ncol(dds)

# Compute log2 Fragments Per Million
quantLog <- log2(fpm(dds)[isexpr, ] + 1)

# Define formula using each variable
form <- ~ (Condition) + (Drug)

# Run variancePartition analysis
varPart <- fitExtractVarPartModel(exprObj=quantLog, formula=form, data=metaFilt)

# Visualise
vp <- sortCols(varPart)

# Violin plot
plot <- plotVarPart(vp) + 
  theme_classic()+
  theme(text=element_text(family='Franklin Gothic Book', colour='black', size=14))+
  scale_fill_manual(values=c('#90bdcf', '#cfc963', '#ba5346', '#ebebeb'))+
  guides(fill='none')

# Save
ggsave('GPP-VarPart-Plot.png', plot, dpi=300)

# DE ANALYSIS ##################################################################
# Run DESEQ
dds <- DESeq(dds)
dds$Condition <- relevel(dds$Condition, ref="normal")

# Results
resultsNames(dds)
res <- results(dds, contrast=c("Condition","normal","generalized pustular psoriasis")) %>% as.data.frame() %>% rownames_to_column('Gene')

# LFC shrinkage
resLFC <- lfcShrink(dds, coef="Condition_normal_vs_generalized.pustular.psoriasis", type="apeglm") %>% as.data.frame() %>% rownames_to_column('Gene')

# Plot 1 gene
plotCounts(dds, gene='SPINK9', intgroup="Condition")

# Top genes
resLFC %>% arrange(padj)

# Flip L2FC
resLFC$log2FoldChange <- resLFC$log2FoldChange*-1

# Volcano plot
makeVolcanoPlot(resLFC) 

# Normalised counts
plot <- NormCountsPlot('TMEM45B', dds=dds, metaData=metaFilt)

# Save
ggsave('GPP-DESEQ-TMEM45B.png', plot, dpi=300)

# GSEA #########################################################################
# Create gene list
resLFCList <- resLFC$log2FoldChange
names(resLFCList) <- resLFC$Gene

# Clean
resLFCList <- na.omit(resLFCList)
resLFCList <- sort(resLFCList, decreasing=T)
resLFCList <- resLFCList[!duplicated(names(resLFCList))]

# Get GO BP gmt file
GO_file <- "GOMF.txt"
GO=fgsea::gmtPathways(GO_file)

# Run FGSEA
set.seed(1)
resGSEA <- fgsea::fgsea(pathways=GO,
                        stats=resLFCList,
                        minSize=10, 
                        maxSize=200)

# Get significant pathways
resGSEAPADJ <- resGSEA %>%
  dplyr::filter(padj<=0.05) %>% 
  arrange(pval)

# Enrichment labels
resGSEA$Enrichment=ifelse((resGSEA$NES > 0)&(resGSEA$padj<0.05), "Upregulated", 'None')
resGSEA$Enrichment=ifelse((resGSEA$NES < 0)&(resGSEA$padj<0.05), "Downregulated", resGSEA$Enrichment)

# Edit GOBP term names
resGSEA$pathway <- gsub(x=resGSEA$pathway, pattern="GOMF_", replacement="") 
resGSEA$pathway <- gsub(x=resGSEA$pathway, pattern="Reactome ", replacement="") 
resGSEA$pathway <- gsub(x=resGSEA$pathway, pattern="_", replacement=" ")
resGSEA$pathway <- str_to_sentence(resGSEA$pathway)
resGSEA$pathway <- gsub(x=resGSEA$pathway, pattern="Fcgr", replacement="FCGR") 
resGSEA$pathway <- gsub(x=resGSEA$pathway, pattern="Cell cell", replacement="Cell-cell") 
resGSEA$pathway <- gsub(x=resGSEA$pathway, pattern="Mhc", replacement="MHC") 
resGSEA$pathway <- gsub(x=resGSEA$pathway, pattern="Ccr", replacement="CCR") 

# GSEA volcano plot
plot <- ggplot(resGSEA, aes(x=NES, y=-log10(padj), size=size)) +
  geom_hline(yintercept=-log10(0.05), color="#ebebeb", linewidth=0.5)+
  geom_vline(xintercept=1, color="#ebebeb", linewidth=0.5)+
  geom_vline(xintercept=-1, color="#ebebeb", linewidth=0.5) +
  geom_point(aes(colour=Enrichment), alpha=0.6) +
  scale_size_continuous(range=c(0,15))+
  theme_bw() +
  ggrepel::geom_text_repel(aes(label=stringr::str_wrap(pathway, 40)), family='Franklin Gothic Book', 
                           data=subset(resGSEA, NES>2.25 & padj<0.05),
                           nudge_y=-0.05, min.segment.length = 4,
                           lineheight=0.8,
                           size=4) +
  xlab("NES") + ylab("-Log10 adjusted p-value") +
  theme(text=element_text(size=16, family='Franklin Gothic Book'),
        strip.background=element_rect(colour=NA, fill='white'),
        panel.grid=element_blank()) +
  scale_color_manual(values=c("#35264b",'#ebebeb', "#ba5346"))+
  scale_x_continuous(limits=c(-5,5))+
  guides(size='none', colour='none')+
  ggtitle('GO Molecular Functions Dysregulated in Generalized Pustular Psoriasis')

# Save
ggsave('GSEA-GOMF-GPP.png', plot, dpi=300)

