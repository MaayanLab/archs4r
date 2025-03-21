
![archs4r](https://github.com/user-attachments/assets/ad4111c0-9ad4-42a2-87fa-bd03d332948a)

# archs4r - Official R Package for ARCHS4 Data

**archs4r** is an R package designed to streamline the loading and querying of ARCHS4 data directly within the R environment. ARCHS4 (All RNA-seq and ChIP-seq Sample and Signature Search) is a comprehensive resource offering access to a vast collection of gene expression data derived from RNA-seq experiments. This package empowers users to efficiently access, manipulate, and analyze ARCHS4 data for a wide range of bioinformatics applications.

## Table of Contents
- [Installation](#installation)
- [Obtaining ARCHS4 Data](#obtaining-archs4-data)
- [Usage](#usage)
    - [🚀 Quick Start Guide](#usage)
    - [Metadata](#metadata)
        - [a4.meta.meta](#a4.meta.meta) - Search metadata across fields
        - [a4.meta.field](#a4.meta.field) - Retrieve values for a specific metadata field
        - [a4.meta.samples](#a4.meta.samples) - Extract metadata for specific samples
        - [a4.meta.series](#a4.meta.series) - Fetch metadata for a GEO series
    - [Data](#data)
        - [a4.data.meta](#a4.data.meta) - Fetch expression data by metadata search
        - [a4.data.rand](#a4.data.rand) - Randomly sample expression data
        - [a4.data.series](#a4.data.series) - Get expression data for a GEO series
        - [a4.data.samples](#a4.data.samples) - Extract expression data for specific samples
        - [a4.data.index](#a4.data.index) - Helper function for index-based data extraction
    - [Utilities](#utilities)
        - [a4.ls](#a4.ls) - List ARCHS4 HDF5 file structure
        - [a4.normalize](#a4.normalize) - Normalize expression count matrix
        - [a4.aggregate_duplicate_genes](#a4.aggregate_duplicate_genes) - Aggregate duplicate gene entries
        - [a4.filter_genes](#a4.filter_genes) - Filter genes by read threshold

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

### 🚀 Quick Start Guide

This workflow provides a foundation for analyzing myoblast-related RNA-seq data from ARCHS4. Adjust parameters (e.g., search_term, readThreshold) as needed for your research question.

```R
# Load library
library(archs4r)

# Set file path
h5file <- "path/to/human_gene_v2.latest.h5"

# List file contents
structure <- a4.ls(h5file)
print(structure)

# Extract myoblast expression
myoblast_exp <- a4.data.meta(h5file, "myoblast", c("title", "source_name_ch1"), remove_sc = TRUE)

# Get metadata for samples
myoblast_samples <- colnames(myoblast_exp)
myoblast_meta <- a4.meta.samples(h5file, myoblast_samples, c("geo_accession", "title", "source_name_ch1", "series_id"))

# Filter genes and aggregate duplicates
filtered_exp <- a4.filter_genes(myoblast_exp, readThreshold = 20, sampleThreshold = 0.02, aggregate = TRUE)

# Log quantile normalization
normalized_exp <- a4.normalize(filtered_exp, method = "log_quantile")

# Inspect results
dim(normalized_exp)
print(head(normalized_exp[, 1:5]))
```


### Metadata

The **Metadata** module in `archs4r` provides tools to extract and query metadata associated with genes, transcripts, and samples from ARCHS4 HDF5 files. These functions allow users to filter and retrieve metadata based on search terms, specific fields, sample IDs, or GEO series IDs.

#### Supported metadata fields

Metadata is grouped into **gene/transcripts** (depending on file), **info**, and **samples**. List file content using **a4.ls(h5file)**.

| Group          | Field                  | Type    |
|----------------|------------------------|---------|
| **meta/genes** | biotype                | str     |
|                | ensembl_gene           | str     |
|                | symbol                 | str     |
| **meta/info**  | author                 | str     |
|                | contact                | str     |
|                | creation-date          | str     |
|                | laboratory             | str     |
|                | version                | str     |
| **meta/samples** | alignedreads         | float64 |
|                | channel_count          | str     |
|                | characteristics_ch1    | str     |
|                | contact_address        | str     |
|                | contact_city           | str     |
|                | contact_country        | str     |
|                | contact_institute      | str     |
|                | contact_name           | str     |
|                | contact_zip            | str     |
|                | data_processing        | str     |
|                | extract_protocol_ch1   | str     |
|                | geo_accession          | str     |
|                | instrument_model       | str     |
|                | last_update_date       | str     |
|                | library_selection      | str     |
|                | library_source         | str     |
|                | library_strategy       | str     |
|                | molecule_ch1           | str     |
|                | organism_ch1           | str     |
|                | platform_id            | str     |
|                | relation               | str     |
|                | sample                 | str     |
|                | series_id              | str     |
|                | singlecellprobability  | float64 |
|                | source_name_ch1        | str     |
|                | status                 | str     |
|                | submission_date        | str     |
|                | taxid_ch1              | str     |
|                | title                  | str     |
|                | type                   | str     |

#### `a4.meta.meta`

Searches for samples whose metadata matches a given pattern across specified fields.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file (e.g., `"human_gene_v2.latest.h5"`).
  - `search_term`: A string to search for in metadata (e.g., `"liver"`).
  - `meta_fields`: A character vector of metadata fields to search. Defaults to `c("characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")`.
  - `remove_sc`: Logical; if `TRUE`, filters out samples with a single-cell probability ≥ 0.5. Defaults to `FALSE`.
  - `silent`: Logical; if `FALSE`, prints the number of matches per field. Defaults to `FALSE`.

- **Returns**: A data frame where rows are samples matching the `search_term` and columns are the metadata fields.

- **Example**:
  ```R
  library(archs4r)
  h5file <- "human_gene_v2.latest.h5"
  df_meta <- a4.meta.meta(h5file, "liver", meta_fields = c("title", "source_name_ch1"), remove_sc = TRUE)
  ```

<a id="a4.meta.field"></a>
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

The **Data** module in `archs4r` provides functions to retrieve gene expression matrices from ARCHS4 HDF5 files based on metadata searches, random sampling, GEO series, or specific sample IDs. These functions rely on the underlying `a4.data.index` helper function to extract expression data.

#### `a4.data.meta`

Fetches expression data for samples whose metadata matches a given search term across specified fields.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file (e.g., `"human_gene_v2.latest.h5"`).
  - `search_term`: A string to search for in metadata using regular expressions (e.g., `"liver"`).
  - `meta_fields`: A character vector of metadata fields to search. Defaults to `c("characteristics_ch1", "source_name_ch1", "title")`.
  - `remove_sc`: Logical; if `TRUE`, excludes samples with single-cell probability ≥ 0.5. Defaults to `FALSE`.
  - `silent`: Logical; if `FALSE`, prints search progress. Defaults to `FALSE`.

- **Returns**: A matrix of expression counts where rows are genes and columns are matching samples.

- **Example**:
  ```R
  library(archs4r)
  h5file <- "human_gene_v2.latest.h5"
  df_meta <- a4.data.meta(h5file, "liver", meta_fields = c("title", "source_name_ch1"), remove_sc = TRUE)
  ```

#### `a4.data.rand`

Randomly selects a specified number of samples and returns their expression data.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file.
  - `number`: Integer; the number of samples to randomly select.
  - `seed`: Integer; seed for reproducibility of random sampling. Defaults to `1`.
  - `remove_sc`: Logical; if `TRUE`, excludes single-cell samples (probability ≥ 0.5). Defaults to `FALSE`.
  - `silent`: Logical; if `FALSE`, passed to `a4.data.index` for verbosity. Defaults to `FALSE`.

- **Returns**: A matrix of expression counts where rows are genes and columns are randomly selected samples.

- **Example**:
  ```R
  df_rand <- a4.data.rand(h5file, 20, seed = 123, remove_sc = TRUE)
  ```

#### `a4.data.series`

Retrieves expression data for all samples in a specified GEO series.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file.
  - `series_id`: A string specifying the GEO series ID (e.g., `"GSE64016"`).
  - `silent`: Logical; if `FALSE`, passed to `a4.data.index` for verbosity. Defaults to `FALSE`.

- **Returns**: A matrix of expression counts where rows are genes and columns are samples from the series. Returns `NULL` if no samples match the series ID.

- **Example**:
  ```R
  df_series <- a4.data.series(h5file, "GSE64016")
  ```

#### `a4.data.samples`

Extracts expression data for a specific set of samples identified by their GEO accession IDs.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file.
  - `sample_ids`: A character vector of GEO accession IDs (e.g., `c("GSM12345", "GSM67890")`).
  - `silent`: Logical; if `FALSE`, passed to `a4.data.index` for verbosity. Defaults to `FALSE`.

- **Returns**: A matrix of expression counts where rows are genes and columns are the specified samples. Returns `NULL` if no samples match the IDs.

- **Example**:
  ```R
  df_samples <- a4.data.samples(h5file, c("GSM1158284", "GSM1482938"))
  ```

#### `a4.data.index`

A helper function that extracts expression data for specified sample and gene indices.

- **Parameters**:
  - `file`: Path to the ARCHS4 HDF5 file.
  - `sample_idx`: Integer vector of sample indices to extract.
  - `gene_idx`: Integer vector of gene indices to extract. Defaults to `integer(0)` (all genes).
  - `silent`: Logical; if `FALSE`, allows verbosity (though not explicitly used here). Defaults to `FALSE`.

- **Returns**: A transposed matrix of expression counts where rows are genes and columns are samples, with row names as gene symbols and column names as GEO accession IDs.

- **Note**: This is typically called internally by other `a4.data.*` functions but can be used directly for custom index-based queries.

- **Example**:
  ```R
  idx <- c(1, 2, 3)  # Example sample indices
  exp_data <- a4.data.index(h5file, idx)
  ```

### Utilities

The **Utilities** module in `archs4r` offers a set of helper functions for preprocessing and analyzing ARCHS4 expression data. These include tools for inspecting file structure, normalizing count matrices, filtering genes, and aggregating duplicate gene entries. Supported normalization methods include quantile normalization, log quantile normalization, counts per million (CPM), and trimmed mean of M-values (TMM).

#### `a4.ls`

Lists the structure and fields of an ARCHS4 HDF5 file, including datasets under `/data`, `/meta/genes`, `/meta/info`, and `/meta/samples`.

- **Parameters**:
  - `h5file`: Path to the ARCHS4 HDF5 file (e.g., `"human_gene_v2.latest.h5"`).

- **Returns**: A data frame with columns `group`, `name`, `otype`, `dclass`, and `dim`, detailing the file’s structure.

- **Example**:
  ```R
  library(archs4r)
  h5file <- "human_gene_v2.latest.h5"
  structure <- a4.ls(h5file)
  ```

Here’s an expanded version of the Markdown documentation for the `a4.normalize` function, with detailed explanations of each normalization option to help users understand their purpose and application:

#### `a4.normalize`

Normalizes an expression count matrix using a specified method to adjust for sequencing depth, library size differences, or other technical biases in RNA-seq data. This function supports four normalization techniques: quantile normalization (`"quantile"`), log quantile normalization (`"log_quantile"`), counts per million (`"cpm"`), and trimmed mean of M-values (`"tmm"`), each suited to different analysis scenarios.

- **Parameters**:
  - `counts`: A matrix or data frame of raw expression counts, where rows represent genes and columns represent samples.
  - `method`: A string specifying the normalization method:
    - `"quantile"`: Applies quantile normalization to force all samples to have the same distribution of counts, useful for removing technical variation while preserving biological differences.
    - `"log_quantile"`: Applies quantile normalization to log2-transformed counts (log2(counts + 1)), reducing the impact of extreme values and stabilizing variance across samples. This is the default method.
    - `"cpm"`: Normalizes counts to counts per million, scaling each sample’s counts by its total library size divided by 1 million, providing a simple measure of relative expression.
    - `"tmm"`: Uses the trimmed mean of M-values method on log2-transformed counts (log2(counts + 1)) to compute normalization factors that adjust for library size and composition biases, trimming outliers to improve robustness.
    Defaults to `"log_quantile"`.
  - `tmm_outlier`: Numeric; the percentage (e.g., `0.05` for 5%) of extreme values to trim from each end when calculating TMM normalization factors. Only used with `method = "tmm"`. Defaults to `0.05`.

- **Returns**: A normalized matrix with the same row names (gene symbols) and column names (sample IDs) as the input, containing adjusted expression values.

- **Normalization Options Explained**:
  - **`quantile`**: This method assumes that the underlying distribution of gene expression should be similar across samples. It ranks the counts within each sample, then replaces them with values from a reference distribution (typically the average quantile across all samples). It’s effective for removing technical noise but may obscure subtle biological differences if sample conditions vary widely.
  - **`log_quantile`**: Similar to `"quantile"`, but first transforms the counts to a log2 scale (adding 1 to avoid log(0)) before normalization. The log transformation reduces the influence of highly expressed genes, making it suitable for datasets with skewed distributions or when downstream analyses (e.g., differential expression) assume log-normal data. This is the default due to its balance of robustness and interpretability.
  - **`cpm`**: Converts raw counts to counts per million by dividing each count by the sample’s total count and multiplying by 1 million. This simple scaling accounts for differences in sequencing depth but does not adjust for gene length or composition biases, making it a quick option for exploratory analysis or when absolute expression levels are less critical.
  - **`tmm`**: The trimmed mean of M-values method calculates a normalization factor for each sample based on the log2-transformed counts, trimming a specified percentage of extreme values (controlled by `tmm_outlier`) to reduce the impact of outliers (e.g., highly expressed or rare genes). It assumes most genes are not differentially expressed and adjusts for both library size and RNA composition differences, making it robust for differential expression studies.

- **Example**:
  ```R
  library(archs4r)
  h5file <- "human_gene_v2.latest.h5"
  exp <- a4.data.rand(h5file, 100)  # Get 100 random samples
  # Normalize with log quantile (default)
  normalized_exp <- a4.normalize(exp, method = "log_quantile")
  # Normalize with CPM
  cpm_exp <- a4.normalize(exp, method = "cpm")
  # Normalize with TMM, trimming 10% outliers
  tmm_exp <- a4.normalize(exp, method = "tmm", tmm_outlier = 0.10)
  ```


#### `a4.aggregate_duplicate_genes`

Aggregates duplicate gene entries by summing their counts.

- **Parameters**:
  - `exp`: A matrix or data frame of expression counts with gene names as row names.

- **Returns**: A data frame with deduplicated gene rows, where counts for duplicate gene symbols are summed.

- **Example**:
  ```R
  dexp <- a4.aggregate_duplicate_genes(exp)
  ```

#### `a4.filter_genes`

Filters genes based on a minimum read threshold across a fraction of samples.

- **Parameters**:
  - `exp`: A matrix or data frame of expression counts.
  - `readThreshold`: Numeric; minimum count threshold for a gene to be retained. Defaults to `20`.
  - `sampleThreshold`: Numeric; minimum fraction of samples where the count exceeds `readThreshold`. Defaults to `0.02`.
  - `deterministic`: Logical; if `TRUE`, sets a fixed seed (`42`) for reproducibility. Defaults to `TRUE`.
  - `aggregate`: Logical; if `TRUE`, aggregates duplicate genes before filtering. Defaults to `TRUE`.

- **Returns**: A filtered matrix or data frame with only genes meeting the criteria.

- **Example**:
  ```R
  fexp <- a4.filter_genes(exp, readThreshold = 20, sampleThreshold = 0.02)
  ```
