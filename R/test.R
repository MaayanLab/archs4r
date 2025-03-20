
h5file = "/Users/maayanlab/Downloads/human_gene_v2.6.h5"

samples = h5read(h5file, "meta/samples/geo_accession")
genes = h5read(h5file, "meta/genes/symbol")

exp = a4.data.samples(h5file, c("GSM1158284","GSM1482938","GSM1562817"))

#exp = a4.data.series(h5file, series="GSE64016")


#exp = a4.data.meta(h5file, 'myoblast')


#res = a4.meta.series(h5file, "GSE40859")

