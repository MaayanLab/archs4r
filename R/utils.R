
# 1. get_config()
# Read a JSON configuration file from the relative data folder.
a4.get_config <- function() {
  # Here we assume the config file (config.json) is in a folder called "data"
  # located in the same directory as the current working directory.
  config_url <- file.path("data", "config.json")
  data <- fromJSON(config_url)
  return(data)
}

# 2. versions()
# Return the names of available versions of human gene counts.
a4.versions <- function() {
  config <- get_config()
  vers <- names(config[["GENE_COUNTS"]][["HUMAN"]])  # Returns the keys
  return(vers)
}


# Normalize a count matrix (data.frame or matrix) according to the chosen method.
a4.normalize <- function(counts, method = "log_quantile", tmm_outlier = 0.05) {
  norm_exp <- NULL
  if (method == "quantile") {
    # preprocessCore::normalize.quantiles expects a matrix of numeric values.
    norm_exp <- normalize.quantiles(as.matrix(counts))
  } else if (method == "log_quantile") {
    norm_exp <- normalize.quantiles(log2(1 + as.matrix(counts)))
  } else if (method == "cpm") {
    norm_exp <- cpm_normalization(counts)
  } else if (method == "tmm") {
    norm_exp <- tmm_norm(counts, percentage = tmm_outlier)
  } else {
    stop(paste("Unsupported normalization method:", method))
  }
  
  # Keep as matrix instead of converting to data.frame
  # Assign row and column names directly
  rownames(norm_exp) <- rownames(counts)
  colnames(norm_exp) <- colnames(counts)
  
  # Debugging output (optional)
  print("Normalization completed")
  
  # Return the matrix
  return(norm_exp)
}


# 4. tmm_norm()
# Apply a trimmed mean normalization on the log2(counts + 1) values.
tmm_norm <- function(exp, percentage = 0.05) {
  # Compute log-transformed expression (log2(1 + exp))
  lexp <- log2(1 + as.matrix(exp))
  # Compute a TMM factor for each sample (column) using trimmed_mean()
  tmm <- trimmed_mean(lexp, percentage)
  # Replicate tmm values as rows (each column gets its normalization factor)
  nf <- matrix(rep(tmm, each = nrow(lexp)), nrow = nrow(lexp))
  temp <- lexp / nf
  return(temp)
}

# 5. trimmed_mean()
# For each column of a matrix, compute the trimmed mean of positive values.
trimmed_mean <- function(matrix_in, percentage) {
  matrix_in <- as.matrix(matrix_in)
  n_cols <- ncol(matrix_in)
  trimmed_means <- numeric(n_cols)
  
  for (col in seq_len(n_cols)) {
    data <- matrix_in[, col]
    # Use only positive values
    data <- data[data > 0]
    if (length(data) == 0) {
      trimmed_means[col] <- NA
      next
    }
    n_trim <- floor(length(data) * percentage)
    sorted_values <- sort(data)
    if (length(sorted_values) > 2 * n_trim) {
      trimmed_values <- sorted_values[(n_trim + 1):(length(sorted_values) - n_trim)]
    } else {
      trimmed_values <- sorted_values
    }
    trimmed_means[col] <- mean(trimmed_values)
  }
  return(trimmed_means)
}

# 6. cpm_normalization()
# Normalize counts to counts per million.
cpm_normalization <- function(df) {
  df <- as.matrix(df)
  sample_sum <- colSums(df)
  scaling_factor <- sample_sum / 1e6
  normalized_df <- sweep(df, 2, scaling_factor, FUN = "/")
  return(normalized_df)
}

get_dataset_dims <- function(h5file, dataset_path) {
  fid <- H5Fopen(h5file)
  on.exit(H5Fclose(fid))
  did <- H5Dopen(fid, dataset_path)
  on.exit(H5Dclose(did), add = TRUE)
  space <- H5Dget_space(did)
  on.exit(H5Sclose(space), add = TRUE)
  dims <- H5Sget_simple_extent_dims(space)$size
  return(rev(dims))  # Reverse dimensions as per your example
}


a4.ls <- function(h5file) {
  fid <- H5Fopen(h5file)
  gid <- H5Gopen(fid, "/meta/samples")
  ggid <- H5Gopen(fid, "/meta/genes")
  iid <- H5Gopen(fid, "/meta/info")
  d1 = h5ls(iid, recursive = FALSE)
  d2 = h5ls(gid, recursive = FALSE)
  d3 = h5ls(ggid, recursive = FALSE)
  # Clean up by closing the group and file handles
  H5Gclose(ggid)
  H5Gclose(gid)
  H5Fclose(fid)

  exp_dims <- get_dataset_dims(h5file, "data/expression")
  exp_type <- 'INTEGER'

  # Create a dataframe for "data/expression"
  exp_df <- data.frame(
    group = "/data",
    name = "expression",
    otype = "H5I_DATASET",
    dclass = exp_type,
    dim = paste(exp_dims, collapse = " x "),
    stringsAsFactors = FALSE
  )

  d1 <- d1 %>% mutate(group = "/meta/info")
  d2 <- d2 %>% mutate(group = "/meta/samples")
  d3 <- d3 %>% mutate(group = "/meta/genes")

  combined_df <- bind_rows(exp_df, d1, d3, d2)
  combined_df
}

a4.aggregate_duplicate_genes <- function(exp) {
  # Use rowsum: rownames(exp) is used as a grouping factor.
  # Convert data.frame to matrix if needed.
  exp_mat <- as.matrix(exp)
  aggregated <- rowsum(exp_mat, group = rownames(exp_mat))
  return(as.data.frame(aggregated, stringsAsFactors = FALSE))
}

# 9. filter_genes()
# Filter genes based on a minimum read threshold in a fraction of samples.
a4.filter_genes <- function(exp, readThreshold = 20, sampleThreshold = 0.02,
                         deterministic = TRUE, aggregate = TRUE) {
  if (deterministic) {
    set.seed(42)
  }
  if (aggregate) {
    exp <- a4.aggregate_duplicate_genes(exp)
  }
  # Count for each gene (row): number of samples with count > readThreshold.
  kk <- rowSums(exp > readThreshold)
  # Identify rows in which at least sampleThreshold fraction of samples exceed the threshold.
  valid_rows <- which(kk >= ncol(exp) * sampleThreshold)
  return(exp[valid_rows, , drop = FALSE])
}