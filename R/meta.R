a4.meta.meta <- function(file, search_term, meta_fields=c("characteristics_ch1", "extract_protocol_ch1", "source_name_ch1", "title")){
  search_term <- toupper(str_replace_all(search_term, "_|-|'|/| |\\.", ""))
  idx = c()
  for(mfield in meta_fields){
    meta = h5read(file, paste0("meta/samples/", mfield))
    cleaned_meta <- toupper(gsub("_|-|'|/|\\.", "", meta))
    idx <- append(idx, grep(search_term, cleaned_meta))
  }
  idx = sort(unique(idx))

  res = data.frame(h5read(destination_file, "meta/samples", bit64conversion = "bit64"))[idx,]
  rownames(res)=res$geo_accession
  return(res)
}

a4.meta.series <- function(file, series_id){
  series = h5read(file, "meta/samples/series_id")
  idx = which(series == series_id)
  res = data.frame(h5read(destination_file, "meta/samples", bit64conversion = "bit64"))[idx,]
  rownames(res)=res$geo_accession
  return(res)
}

a4.meta.samples <- function(file, samples){
  samples <- unique(samples)
  h5_samples = h5read(file, "meta/samples/geo_accession")
  genes = h5read(file, "meta/genes/symbol")
  idx = which(h5_samples %in% samples)
  res = data.frame(h5read(destination_file, "meta/samples", bit64conversion = "bit64"))[idx,]
  rownames(res)=res$geo_accession
  return(res)
}
