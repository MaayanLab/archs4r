
![archs4r](https://github.com/user-attachments/assets/ad4111c0-9ad4-42a2-87fa-bd03d332948a)

# archs4r - Official R Package for ARCHS4 Data

**archs4r** is an R package designed to streamline the loading and querying of ARCHS4 data directly within the R environment. ARCHS4 (All RNA-seq and ChIP-seq Sample and Signature Search) is a comprehensive resource offering access to a vast collection of gene expression data derived from RNA-seq experiments. This package empowers users to efficiently access, manipulate, and analyze ARCHS4 data for a wide range of bioinformatics applications.

## Installation

To get started with `archs4r`, you need R installed on your system. The package depends on the `rhdf5` package from Bioconductor to handle HDF5 files. Follow these steps to install `archs4r`:

1. **Install the `devtools` package** (if not already installed):
```R
 install.packages("devtools")
 ```
2. **Install the `archs4` package**
```R
library("devtools")
install_github("MaayanLab/archs4r")
```
3. **Manual install of `Bioconductor` dependencies**
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
BiocManager::install("prerpocessCore")
```

### Obtaining ARCHS4 Data

The archs4r package requires ARCHS4 data files in HDF5 format, which contain gene expression data and associated metadata. **ARCHS4r** can interact with gene and transcript-level files. Download the latest versions from the <a href="https://archs4.org/download" target="_blank">ARCHS4 download page</a>. After downloading, note the file path (e.g., "path/to/human_gene_v2.latest.h5") as it will be needed for package functions.

<img src="https://github.com/user-attachments/assets/f2de5f94-a790-42fb-b99e-cd77c74642eb" alt="image" width="700">


## Usage

The archs4r package is organized into three main modules: Metadata, Data, and Utilities. Below are detailed descriptions and examples for each.

### Metadata
### Metadata

The **Metadata** module in `archs4r` provides tools to extract and query metadata associated with genes, transcripts, and samples from ARCHS4 HDF5 files. These functions allow users to filter and retrieve metadata based on search terms, specific fields, sample IDs, or GEO series IDs.

#### `a4.meta.meta`

Searches for samples whose metadata matches a given pattern across specified fields.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file (e.g., `"human_gene_v2.latest.h5"`).
  - `search_term`: A string to search for in metadata (e.g., `"liver"`).
  - `meta_fields`: A character vector of metadata fields to search. Defaults to `c("characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")`.
  - `remove_sc`: Logical; if `TRUE`, filters out samples with a single-cell probability ≥ 0.5. Defaults to `FALSE`.
  - `silent`: Logical; if `FALSE`, prints the number of matches per field. Defaults to `FALSE`.

- **Returns**: A transposed data frame where rows are samples matching the `search_term` and columns are the metadata fields.

- **Example**:
  ```R
  library(archs4r)
  h5file <- "human_gene_v2.latest.h5"
  df_meta <- a4.meta.meta(h5file, "liver", meta_fields = c("title", "source_name_ch1"), remove_sc = TRUE)
  ```

#### `a4.meta.field`

Retrieves all values for a specific metadata field from samples, genes, or transcripts.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file.
  - `field`: A string specifying the metadata field (e.g., `"geo_accession"`, `"symbol"`).

- **Returns**: A character vector of values for the specified field, converted to uppercase. Raises an error if the field doesn’t exist in `/meta/samples`, `/meta/genes`, or `/meta/transcripts`.

- **Example**:
  ```R
  samples <- a4.meta.field(h5file, "geo_accession")
  genes <- a4.meta.field(h5file, "symbol")
  ```

#### `a4.meta.samples`

Extracts metadata for a specific set of samples identified by their GEO accession IDs.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file.
  - `samples`: A character vector of GEO accession IDs (e.g., `c("GSM12345", "GSM67890")`).
  - `meta_fields`: A character vector of metadata fields to retrieve. Defaults to `c("geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")`.
  - `silent`: Logical; if `FALSE`, prints extraction details for each field. Defaults to `FALSE`.

- **Returns**: A transposed data frame where rows are the specified samples and columns are the requested metadata fields.

- **Example**:
  ```R
  df_samples <- a4.meta.samples(h5file, c("GSM12345", "GSM67890"), meta_fields = c("title", "series_id"))
  ```

#### `a4.meta.series`

Fetches metadata for all samples within a specified GEO series.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file.
  - `series`: A string specifying the GEO series ID (e.g., `"GSE64016"`).
  - `meta_fields`: A character vector of metadata fields to retrieve. Defaults to `c("geo_accession", "series_id", "characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")`.
  - `silent`: Logical; if `FALSE`, prints extraction details for each field. Defaults to `FALSE`.

- **Returns**: A transposed data frame where rows are samples from the specified series and columns are the requested metadata fields.

- **Example**:
  ```R
  df_series <- a4.meta.series(h5file, "GSE64016", meta_fields = c("geo_accession", "characteristics_ch1"))
  ```

### Data

The data module has endpoints similar to those of the metadata module but instead returns gene/transcript expression matrices.

```R
library("archs4r")

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

### Utilities

The utilities module contains some helpful functions such as data normalization and gene filtering. Gene symbols are not unique and by aggregating the gene rows are deduplicated by summing up all counts belonging to the same gene symbol. For normalization there is quantile normalization, TMM, and CPM.

```R
library("archs4r")

h5file = "human_gene_v2.latest.h5"

# List H5 file structure and fields
a4.ls(h5file)

exp = a4.data.rand(h5file, 100)
normalized_exp = a4.normalize((exp, method = "log_quantile") # method options: log_quantile, cpm, tmm, quantile

# filter genes with low expression
fexp = a4.filter_genes(exp, readThreshold = 20, sampleThreshold = 0.02, deterministic = TRUE, aggregate = TRUE)

# Merge counts when ensembl ids point to the sample gene symbol. Counts are added.
dexp = a4.aggregate_duplicate_genes(exp)
```
