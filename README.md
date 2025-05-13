# Full RNAseq Workflow in R
## Loading Data
1. Read in a **count matrix** (read counts for each sample) and a **meta table** (containing information about each sample).

If using a GEO dataset use:
```
metaData <- getGEO(filename="~/Downloads/GSE000000_series_matrix.txt.gz")
```
