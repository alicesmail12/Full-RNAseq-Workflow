# Full RNAseq Workflow in R
## Loading and Checking Data
1. Read in a **count matrix** (read counts for each sample) and a **meta table** (containing information about each sample). For example, if using a GEO dataset use:
```
# Read count and meta data
metaData <- getGEO(filename="~/Downloads/GSE000000_series_matrix.txt.gz")
countMatrix <- read.table('GSE000000_read_counts.txt', fill=TRUE, header=TRUE) 
```
2. Next perform a **PCA** on the count matrix.
```
# PCA using prcomp
results <- prcomp(countMatrix, scale=F)$rotation %>% as.data.frame() %>% rownames_to_column('Sample')

# Plot
ggplot(results, aes(x=PC1, y=PC2, colour=Sample, label=Sample))+
  geom_point()+
  geom_text(colour='black')
```
