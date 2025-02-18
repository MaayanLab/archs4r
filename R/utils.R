library("jsonlite")
library("preprocessCore")

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

# 3. normalize()
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
  # Return the result as a data.frame with the same row and column names as counts.
  norm_exp <- as.data.frame(norm_exp, stringsAsFactors = FALSE)
  rownames(norm_exp) <- rownames(counts)
  colnames(norm_exp) <- colnames(counts)
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

# 7. ls()
# List all meta data groups and meta data fields in an H5 file.
a4.ls <- function(file) {
  # List objects in the HDF5 file using h5ls()
  ls_result <- h5ls(file)
  # A helper function to print each object similar to the Python code.
  print_data <- function(name, obj, prefix="") {
    if (obj$class == "H5D") {  # Dataset
      # Try to mimic data type reporting and shape.
      data_type <- as.character(obj$dclass)
      shape <- paste(obj$dim, collapse = " x ")
      cat(sprintf("%s%-20s  %-6s | %s\n", prefix, name, data_type, shape))
    } else {
      cat(sprintf("%s%-26s\n", prefix, name))
    }
    # Attributes (if any) may be obtained by h5readAttributes
    # (Not all objects have attributes; skip them if none)
    attr_list <- tryCatch(h5readAttributes(file, obj$fullname), error = function(e) list())
    if (length(attr_list) > 0) {
      for (key in names(attr_list)) {
        cat(sprintf("%s   %-11s : %s\n", prefix, key, as.character(attr_list[[key]])))
      }
    }
  }
  # Print top-level objects
  top_level <- unique(ls_result$group)
  for (grp in unique(ls_result$group)) {
    # Get objects in this group.
    sub_objs <- ls_result[ls_result$group == grp, ]
    for (i in seq_len(nrow(sub_objs))) {
      print_data(sub_objs$name[i], sub_objs[i, ])
    }
  }
}

# 8. aggregate_duplicate_genes()
# Sum values for duplicate genes (rows with identical row names).
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