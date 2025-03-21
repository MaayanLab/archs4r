a4.data.meta <- function(file, search_term,
                 meta_fields = c("characteristics_ch1", "source_name_ch1", "title"),
                 remove_sc = FALSE, silent = FALSE) {

  if (!silent)
    cat("Searching for any occurrence of", search_term, "as regular expression\n")
  
  idx <- integer(0)
  for (field in meta_fields) {
      # Read the metadata values
      meta_values <- h5read(file, paste0("meta/samples/", field))
      # Find indices where the search term is present using regular expression matching
      matches <- which(grepl(search_term, meta_values, perl = TRUE, ignore.case = TRUE))
      idx <- union(idx, matches)
      h5closeAll()
  }
  if (remove_sc) {
    # Read the single-cell probability values and select those below 0.5
    singleprob <- h5read(file, "meta/samples/singlecellprobability")
    idx <- intersect(idx, which(singleprob < 0.5))
    h5closeAll()
  }
  idx <- sort(unique(idx))
  counts <- a4.data.index(file, idx, silent = silent)
  return(counts)
}


a4.data.rand <- function(file, number, seed = 1, remove_sc = FALSE, silent = FALSE) {
  set.seed(seed)
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

a4.data.series <- function(file, series_id, silent = FALSE) {
    series_vec <- h5read(file, "meta/samples/series_id")
    idx <- which(series_vec == series_id)
    if (length(idx) > 0)
      return(a4.data.index(file, idx, silent = silent))
    else
      return(NULL)
}

a4.data.samples <- function(file, sample_ids, silent = FALSE) {
    sample_ids <- as.character(sample_ids)
    gsm_ids <- h5read(file, "meta/samples/geo_accession")
    idx <- which(gsm_ids %in% sample_ids)
    if (length(idx) > 0)
      return(a4.data.index(file, idx, silent = silent))
    else
      return(NULL)
}

a4.data.index <- function(file, sample_idx, gene_idx = integer(0), silent = FALSE) {
  sample_idx <- sort(sample_idx)
  gene_idx <- sort(gene_idx)
  genes <- h5read(file, "meta/genes/symbol")
  if (length(gene_idx) == 0)
      gene_idx <- seq_along(genes)
  gsm_ids <- h5read(file, "meta/samples/geo_accession")
  gsm_ids <- gsm_ids[sample_idx]

  exp_data <- t(h5read(file, "data/expression", index = list(sample_idx, gene_idx)))
  h5closeAll()
  rownames(exp_data) <- genes[gene_idx]
  colnames(exp_data) <- gsm_ids
  return(exp_data)
}

get_encoding <- function(h5file) {
  # Open the HDF5 file
  fid <- H5Fopen(h5file)
  
  # Open the /meta/ group
  ggid <- H5Gopen(fid, "/meta/")
  
  # List contents of /meta/
  meta_contents <- h5ls(ggid, recursive = FALSE)
  
  # Clean up
  H5Gclose(ggid)
  H5Fclose(fid)
  
  # Get the names of datasets/groups in /meta/
  names <- meta_contents$name
  
  # Check for gene-related or transcript-related entries
  has_genes <- any(grepl("gene[s|_ids]?", names, ignore.case = TRUE))
  has_transcripts <- any(grepl("transcript[s|_ids]?", names, ignore.case = TRUE))
  
  # Return based on presence
  if (has_genes && !has_transcripts) {
    return("/meta/genes/symbol")
  } else if (has_transcripts && !has_genes) {
    return("/meta/transcripts/ensembl_id")
  } else {
    return("unknown")  # If both or neither are present, return "unknown"
  }
}