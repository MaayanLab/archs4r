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
# Search metadata for a pattern (e.g. "liver")
df_meta <- a4.data.meta("mydata.h5", "liver")

# Randomly select 5 samples
df_rand <- a4.data.rand("mydata.h5", 5, seed = 123)

# Select samples from a given series
df_series <- a4.data.series("mydata.h5", "GSE12345")

# Select specific samples by their geo_accession IDs
df_samples <- a4.data.samples("mydata.h5", c("GSM111", "GSM222"))
```
