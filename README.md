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

### Metadata
```R
h5file = "human_gene_v2.latest.h5"

# Search for samples whose metadata fields (among a subset) contain “liver”
df_meta <- a4.meta.meta(h5file, "liver")

# Get all metadata for field
samples <- a4.meta.field(h5file, "geo_accession")
genes <- a4.meta.field(h5file, "symbol")

# Select samples from a given series
df_series <- a4.meta.series(h5file, "GSE64016")

# Extract metadata for specific samples:
df_samples <- a4.meta.samples(h5file, c("GSM12345", "GSM67890"))

```

### Data
```R
h5file = "human_gene_v2.latest.h5"

# Search metadata for a pattern (e.g. "liver")
df_meta <- a4.data.meta(h5file, "liver")

# Randomly select 5 samples
df_rand <- a4.data.rand(h5file, 20, seed = 123)

# Select samples from a given series
df_series <- a4.data.series(h5file, "GSE64016")

# Select specific samples by their geo_accession IDs
df_samples <- a4.data.samples(h5file, c("GSM1158284","GSM1482938","GSM1562817"))
```
