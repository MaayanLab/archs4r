# archs4r
R package to load and query ARCHS4 data in R


## Installation

Install directly from GitHub:

```R
install.packages("devtools")  # or "remotes"
library("devtools")
install_github("MaayanLab/archs4r")
```

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
```

## Usage

```R
h5file = "human_gene_v2.latest.h5"

# Search metadata for a pattern (e.g. "liver")
df_meta <- a4.data.meta(h5file, "liver")

# Randomly select 5 samples
df_rand <- a4.data.rand(h5file, 5, seed = 123)

# Select samples from a given series
df_series <- a4.data.series(h5file, "GSE64016")

# Select specific samples by their geo_accession IDs
df_samples <- a4.data.samples(h5file, c("GSM1158284","GSM1482938","GSM1562817"))
```
