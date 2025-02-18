library("rhdf5")

# Helper: convert a vector of strings to uppercase.
to_upper <- function(x) {
  toupper(as.character(x))
}

# 1. meta()
# Search for samples based on a search term in selected metadata fields.
a4.meta.meta <- function(file, search_term,
                 meta_fields = c("characteristics_ch1", "extract_protocol_ch1",
                                 "source_name_ch1", "title"),
                 remove_sc = FALSE, silent = FALSE) {
  # Convert the search term to uppercase for case–insensitive matching.
  search_term <- toupper(search_term)
  
  # Read geo_accession values from the samples group and convert to uppercase.
  geo_acc <- to_upper(h5read(file, "meta/samples/geo_accession"))
  
  # List available fields under meta/samples.
  samples_list <- h5ls(file, "/meta/samples")$name
  
  meta_list <- list()  # To store each meta field’s values.
  used_fields <- c()   # Keep track of fields we successfully read.
  
  # Loop over each field provided.
  for (field in meta_fields) {
    if (field %in% samples_list) {
      # Try reading and converting to uppercase.
      field_data <- to_upper(h5read(file, paste0("meta/samples/", field)))
      meta_list[[field]] <- field_data
      used_fields <- c(used_fields, field)
    }
  }
  
  # Create a data.frame with rows = meta_fields and columns = samples (geo_accession)
  # (Each row is one meta_field; columns are ordered by geo_accession.)
  meta_df <- as.data.frame(do.call(rbind, meta_list), stringsAsFactors = FALSE)
  rownames(meta_df) <- used_fields
  colnames(meta_df) <- geo_acc
  
  # Collect sample indices where any meta field contains the search term.
  idx <- integer(0)
  for (i in seq_len(nrow(meta_df))) {
    # Use grepl to search the (already uppercase) strings.
    matches <- which(grepl(search_term, meta_df[i, ], perl = TRUE))
    idx <- union(idx, matches)
    if (!silent) {
      # Optionally print progress for each field.
      cat(sprintf("Field '%s': found %d matches\n", rownames(meta_df)[i], length(matches)))
    }
  }
  
  if (remove_sc) {
    # Read singlecellprobability and select samples with probability < 0.5.
    singleprob <- as.numeric(h5read(file, "meta/samples/singlecellprobability"))
    idx <- intersect(idx, which(singleprob < 0.5))
  }
  idx <- sort(unique(idx))
  
  # Return the selected metadata with samples as rows.
  # (Transpose so rows = samples and columns = meta fields.)
  return(t(meta_df[, idx, drop = FALSE]))
}

# 2. field()
# Retrieve a specified field from one of the meta groups.
a4.meta.field <- function(file, field) {
  # List keys available in the samples, genes, and transcripts groups.
  sample_fields <- h5ls(file, "/meta/samples")$name
  gene_fields <- tryCatch(h5ls(file, "/meta/genes")$name,
                          error = function(e) NULL)
  transcript_fields <- tryCatch(h5ls(file, "/meta/transcripts")$name,
                                error = function(e) NULL)
  
  if (field %in% sample_fields) {
    out <- h5read(file, paste0("meta/samples/", field))
    return(to_upper(out))
  } else if (!is.null(gene_fields) && (field %in% gene_fields)) {
    out <- h5read(file, paste0("meta/genes/", field))
    return(to_upper(out))
  } else if (!is.null(transcript_fields) && (field %in% transcript_fields)) {
    out <- h5read(file, paste0("meta/transcripts/", field))
    return(to_upper(out))
  } else {
    stop("Specified field does not exist. Choose from supported sample meta fields or gene meta fields.")
  }
}

# 3. samples()
# Extract metadata for specified samples.
a4.meta.samples <- function(file, samples,
                    meta_fields = c("geo_accession", "series_id",
                                    "characteristics_ch1", "extract_protocol_ch1",
                                    "source_name_ch1", "title"),
                    silent = FALSE) {
  # Convert samples to uppercase.
  samples <- to_upper(samples)
  geo_acc <- to_upper(h5read(file, "meta/samples/geo_accession"))
  
  # Determine indices of samples that are in the specified set.
  idx <- which(geo_acc %in% samples)
  
  meta_list <- list()
  used_fields <- c()
  # Loop over the meta fields and extract only for the selected sample indices.
  for (field in meta_fields) {
    sample_fields <- h5ls(file, "/meta/samples")$name
    if (field %in% sample_fields) {
      field_data <- to_upper(h5read(file, paste0("meta/samples/", field)))
      meta_list[[field]] <- field_data[idx]
      used_fields <- c(used_fields, field)
      if (!silent) {
        cat(sprintf("Extracting field '%s' for %d samples\n", field, length(idx)))
      }
    }
  }
  
  meta_df <- as.data.frame(do.call(rbind, meta_list), stringsAsFactors = FALSE)
  rownames(meta_df) <- used_fields
  colnames(meta_df) <- geo_acc[idx]
  # Return the transposed data.frame so that samples are rows.
  inter <- intersect(colnames(meta_df), samples)
  return(t(meta_df[, inter, drop = FALSE]))
}

# 4. series()
# Extract metadata for samples in a specified series.
a4.meta.series <- function(file, series,
                   meta_fields = c("geo_accession", "series_id",
                                   "characteristics_ch1", "extract_protocol_ch1",
                                   "source_name_ch1", "title"),
                   silent = FALSE) {
  series <- toupper(series)
  series_vec <- to_upper(h5read(file, "meta/samples/series_id"))
  
  idx <- which(series_vec == series)
  meta_list <- list()
  used_fields <- c()
  
  for (field in meta_fields) {
    sample_fields <- h5ls(file, "/meta/samples")$name
    if (field %in% sample_fields) {
      field_data <- to_upper(h5read(file, paste0("meta/samples/", field)))
      meta_list[[field]] <- field_data[idx]
      used_fields <- c(used_fields, field)
      if (!silent) {
        cat(sprintf("Series '%s': field '%s' has %d entries\n", series, field, length(idx)))
      }
    }
  }
  
  meta_df <- as.data.frame(do.call(rbind, meta_list), stringsAsFactors = FALSE)
  rownames(meta_df) <- used_fields
  colnames(meta_df) <- to_upper(h5read(file, "meta/samples/geo_accession")[idx])
  return(t(meta_df))
}

# 5. get_meta()
# Return all sample metadata as a named list.
get_meta <- function(file) {
  meta_list <- list()
  sample_fields <- h5ls(file, "/meta/samples")$name
  for (field in sample_fields) {
    data <- h5read(file, paste0("meta/samples/", field))
    if (!is.character(data)) {
      data <- as.character(data)
    }
    meta_list[[field]] <- data
  }
  return(meta_list)
}

# 6. get_meta_sample_field()
# Return one specific field from meta/samples.
get_meta_sample_field <- function(file, field) {
  data <- h5read(file, paste0("meta/samples/", field))
  if (!is.character(data)) {
    data <- as.character(data)
  }
  return(data)
}

# 7. get_meta_gene_field()
# Return one specific field from meta/genes.
get_meta_gene_field <- function(file, field) {
  data <- h5read(file, paste0("meta/genes/", field))
  if (!is.character(data)) {
    data <- as.character(data)
  }
  return(data)
}