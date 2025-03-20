
# Helper: convert a vector of strings to uppercase.
to_upper <- function(x) {
  toupper(as.character(x))
}

a4.meta.meta <- function(file, search_term,
                 meta_fields = c("characteristics_ch1", "extract_protocol_ch1",
                                 "source_name_ch1", "title"),
                 remove_sc = FALSE, silent = FALSE) {

  geo_acc <- h5read(file, "meta/samples/geo_accession")
  fid <- H5Fopen(h5file)
  gid <- H5Gopen(fid, "/meta/samples")
  samples_list = h5ls(gid, recursive = FALSE)$name
  H5Gclose(gid)
  H5Fclose(fid)
  h5closeAll()
  
  meta_list <- list()
  used_fields <- c()
  
  for (field in meta_fields) {
    if (field %in% samples_list) {
      field_data <- h5read(file, paste0("meta/samples/", field))
      meta_list[[field]] <- field_data
      used_fields <- c(used_fields, field)
    }
  }
  
  meta_df <- as.data.frame(do.call(rbind, meta_list), stringsAsFactors = FALSE)
  rownames(meta_df) <- used_fields
  colnames(meta_df) <- geo_acc
  
  idx <- integer(0)
  for (i in seq_len(nrow(meta_df))) {
    matches <- which(grepl(search_term, meta_df[i, ], perl = TRUE, ignore.case = TRUE))
    idx <- union(idx, matches)
    if (!silent) {
      cat(sprintf("Field '%s': found %d matches\n", rownames(meta_df)[i], length(matches)))
    }
  }
  
  if (remove_sc) {
    singleprob <- as.numeric(h5read(file, "meta/samples/singlecellprobability"))
    idx <- intersect(idx, which(singleprob < 0.5))
  }
  idx <- sort(unique(idx))
  return(t(meta_df[, idx, drop = FALSE]))
}

get_group_names <- function(fid, group_path) {
  # Try to open the group
  gid <- try(H5Gopen(fid, group_path), silent = TRUE)
  
  # If successful, get names and close group
  if (!inherits(gid, "try-error")) {
    names <- h5ls(gid, recursive = FALSE)$name
    H5Gclose(gid)
    return(names)
  } else {
    # If group doesn't exist, return empty character vector
    return(c())
  }
}



a4.meta.field <- function(file, field) {

  # Open the file once and reuse the file ID
  fid <- H5Fopen(h5file)

  # Get names from each group
  sample_fields <- get_group_names(fid, "/meta/samples")
  gene_fields <- get_group_names(fid, "/meta/genes")
  transcript_fields <- get_group_names(fid, "/meta/transcripts")

  # Close the file and all remaining open handles
  H5Fclose(fid)
  h5closeAll()

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

a4.meta.samples <- function(file, samples,
                    meta_fields = c("geo_accession", "series_id",
                                    "characteristics_ch1", "extract_protocol_ch1",
                                    "source_name_ch1", "title"),
                    silent = FALSE) {

  geo_acc <- h5read(file, "meta/samples/geo_accession")
  
  # Determine indices of samples that are in the specified set.
  idx <- which(geo_acc %in% samples)
  
  meta_list <- list()
  used_fields <- c()
  # Loop over the meta fields and extract only for the selected sample indices.
  for (field in meta_fields) {
      # Open the file once and reuse the file ID
      fid <- H5Fopen(h5file)
      sample_fields <- get_group_names(fid, "/meta/samples")
      H5Fclose(fid)
      h5closeAll()

    if (field %in% sample_fields) {
      field_data <- h5read(file, paste0("meta/samples/", field))
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

a4.meta.series <- function(file, series,
                    meta_fields = c("geo_accession", "series_id",
                                    "characteristics_ch1", "extract_protocol_ch1",
                                    "source_name_ch1", "title"),
                    silent = FALSE) {

  geo_acc <- h5read(file, "meta/samples/geo_accession")
  series_id <- h5read(file, "meta/samples/series_id")
  
  # Determine indices of samples that are in the specified set.
  idx <- which(series_id == series)
  
  meta_list <- list()
  used_fields <- c()
  # Loop over the meta fields and extract only for the selected sample indices.
  for (field in meta_fields) {
      # Open the file once and reuse the file ID
      fid <- H5Fopen(h5file)
      sample_fields <- get_group_names(fid, "/meta/samples")
      H5Fclose(fid)
      h5closeAll()

    if (field %in% sample_fields) {
      field_data <- h5read(file, paste0("meta/samples/", field))
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
  inter <- intersect(colnames(meta_df), geo_acc[idx])
  return(t(meta_df[, inter, drop = FALSE]))
}
