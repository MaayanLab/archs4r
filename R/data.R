library("rhdf5")
library("stringr")


a4.data.meta <- function(file, search_term, meta_fields=c("characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")){
  search_term <- toupper(str_replace_all(search_term, "_|-|'|/| |\\.", ""))
  idx = c()
  for(mfield in meta_fields){
    meta = h5read(file, paste0("meta/samples/", mfield))
    cleaned_meta <- toupper(gsub("_|-|'|/|\\.", "", meta))
    idx <- append(idx, grep(search_term, cleaned_meta))
  }
  idx = sort(unique(idx))
  return(a4.data.index(file, idx))
}

a4.data.series <- function(file, series){
  h5_series = h5read(file, "meta/samples/series_id")
  samples = h5read(file, "meta/samples/geo_accession")
  genes = h5read(file, "meta/genes/symbol")

  series_locations = which(h5_series == series)

  expression = t(h5read(file, "data/expression", index=list(series_locations, 1:length(genes))))
  rownames(expression) = genes
  colnames(expression) = samples[series_locations]

  return(expression)
}

a4.data.samples <- function(file, samples) {
  samples <- unique(samples)

  h5_samples = h5read(file, "meta/samples/geo_accession")
  genes = h5read(file, "meta/genes/symbol")

  sample_locations = which(h5_samples %in% samples)

  expression = t(h5read(file, "data/expression", index=list(sample_locations, 1:length(genes))))
  rownames(expression) = genes
  colnames(expression) = h5_samples[sample_locations]

  return(expression)
}

a4.data.rand <- function(file, N) {
  samples = h5read(file, "meta/samples/geo_accession")
  random_idx <- sort(sample(0:length(samples), N, replace = FALSE))
  return(a4.data.index(file, random_idx))
}

a4.data.index <- function(file, idx) {
  samples <- unique(samples)

  samples = h5read(file, "meta/samples/geo_accession")[idx]
  genes = h5read(file, "meta/genes/symbol")

  expression = t(h5read(file, "data/expression", index=list(idx, 1:length(genes))))
  rownames(expression) = genes
  colnames(expression) = samples

  return(expression)
}
