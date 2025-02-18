
destination_file = "/Users/maayanlab/Documents/human_gene_v2.6.h5"

samples = h5read(destination_file, "meta/samples/geo_accession")
genes = h5read(destination_file, "meta/genes/symbol")


#exp = a4.data.samples('/Users/maayanlab/Documents/human_gene_v2.2.h5', samples=c("GSM1158284","GSM1482938","GSM1562817"))

#exp = a4.data.series('/Users/maayanlab/Documents/human_gene_v2.2.h5', series="GSE64016")


#a4.data.meta('/Users/maayanlab/Documents/human_gene_v2.2.h5', 'myoblast')


#res = a4.meta.series(destination_file, "GSE40859")
