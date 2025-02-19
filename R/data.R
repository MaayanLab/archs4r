
# Clean a string by removing underscores, dashes, quotes, slashes, spaces, and dots,
# then making it uppercase.
clean_string <- function(x) {
  toupper(gsub("(_|-|'|/| |\\.)", "", x))
}

# meta: searches the metadata of an HDF5 file for a given search term.
a4.data.meta <- function(file, search_term,
                 meta_fields = c("geo_accession", "series_id", "characteristics_ch1",
                                 "extract_protocol_ch1", "source_name_ch1", "title"),
                 remove_sc = FALSE, silent = FALSE) {
  # Clean the search term
  search_term_clean <- clean_string(search_term)
  if (!silent)
    cat("Searching for any occurrence of", search_term_clean, "as regular expression\n")
  meta_local(file, search_term_clean, meta_fields, remove_sc, silent)
}

# meta_local: For each metadata field present in the HDF5 file, search for match(es)
# with the prepared search term; optionally remove samples with singlecellprobability ≥ 0.5.
meta_local <- function(file, search_term,
                       meta_fields = c("geo_accession", "series_id", "characteristics_ch1",
                                       "extract_protocol_ch1", "source_name_ch1", "title"),
                       remove_sc = FALSE, silent = FALSE) {
  # List available entries under the meta/samples group
  fid <- H5Fopen(h5file)
  gid <- H5Gopen(fid, "/meta/samples")
  meta_fields <- h5ls(gid, recursive = FALSE)$name
  H5Gclose(gid)
  H5Fclose(fid)
  idx <- integer(0)
  for (field in meta_fields) {
    if (field %in% samples_fields) {
      # Read the metadata values
      meta_values <- h5read(file, paste0("meta/samples/", field))
      meta_clean <- clean_string(meta_values)
      # Find indices where the search term is present using regular expression matching
      matches <- which(grepl(search_term, meta_clean, perl = TRUE))
      idx <- union(idx, matches)
    }
  }
  if (remove_sc) {
    # Read the single-cell probability values and select those below 0.5
    singleprob <- h5read(file, "meta/samples/singlecellprobability")
    idx <- intersect(idx, which(singleprob < 0.5))
  }
  idx <- sort(unique(idx))
  counts <- a4.data.index(file, idx, silent = silent)
  return(counts)
}

# rand: randomly select a given number of samples.
a4.data.rand <- function(file, number, seed = 1, remove_sc = FALSE, silent = FALSE) {
  set.seed(seed)
  rand_local(file, number, remove_sc, silent)
}

rand_local <- function(file, number, remove_sc, silent = FALSE) {
  gsm_ids <- h5read(file, "meta/samples/geo_accession")
  if (remove_sc) {
    singleprob <- h5read(file, "meta/samples/singlecellprobability")
    eligible <- which(singleprob < 0.5)
    idx <- sort(sample(eligible, number))
  } else {
    total <- length(gsm_ids)
    idx <- sort(sample(seq_len(total), number))
  }
  return(a4.data.index(file, idx, silent = silent))
}

# series: retrieve samples belonging to a specific series_id.
a4.data.series <- function(file, series_id, silent = FALSE) {
  series_local(file, series_id, silent)
}

series_local <- function(file, series_id, silent = FALSE) {
  series_vec <- h5read(file, "meta/samples/series_id")
  idx <- which(series_vec == series_id)
  if (length(idx) > 0)
    return(a4.data.index(file, idx, silent = silent))
  else
    return(NULL)
}

# samples: retrieve samples by a list of geo_accession identifiers.
a4.data.samples <- function(file, sample_ids, silent = FALSE) {
  samples_local(file, sample_ids, silent)
}

samples_local <- function(file, sample_ids, silent = FALSE) {
  sample_ids <- as.character(sample_ids)
  gsm_ids <- h5read(file, "meta/samples/geo_accession")
  idx <- which(gsm_ids %in% sample_ids)
  if (length(idx) > 0)
    return(a4.data.index(file, idx, silent = silent))
  else
    return(NULL)
}

# index: load gene expression data for given sample indices (and an optional set of genes).
# It determines the gene names via get_encoding(), reads the expression dataset,
# subsets the data, and returns a data.frame.
a4.data.index <- function(file, sample_idx, gene_idx = integer(0), silent = FALSE) {
  sample_idx <- sort(sample_idx)
  gene_idx <- sort(gene_idx)
  genes <- h5read(file, "meta/genes/symbol")
  if (length(gene_idx) == 0)
      gene_idx <- seq_along(genes)
  gsm_ids <- h5read(file, "meta/samples/geo_accession")
  gsm_ids <- gsm_ids[sample_idx]

  exp_data <- h5read(file, "data/expression", index = list(gene_idx, sample_idx))
  rownames(exp_data) <- genes[gene_idx]
  colnames(exp_data) <- gsm_ids
  return(exp_data)
}

# get_encoding: Determines which metadata field contains gene (or transcript) identifiers.
get_encoding <- function(file) {
  ls_meta <- h5ls(file, recursive = TRUE)
  # Check if the 'meta/genes' group exists and whether it has a gene_symbol or symbol field.
  genes_group <- ls_meta$name[ls_meta$group == "/meta/genes"]
  if ("gene_symbol" %in% genes_group)
    return("meta/genes/gene_symbol")
  else if ("symbol" %in% genes_group)
    return("meta/genes/symbol")
  
  # Otherwise, look for a transcripts group with ensembl_id.
  transcripts_group <- ls_meta$name[ls_meta$group == "/meta/transcripts"]
  if ("ensembl_id" %in% transcripts_group)
    return("meta/transcripts/ensembl_id")
  
  stop("error in gene/transcript meta data")
}
