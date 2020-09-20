rtracklayer::import.bed('BioData/eclip/')
hepg2 <- read_tsv('BioData/eclip/HepG2.tsv', col_names = FALSE)
k562 <- read_tsv('BioData/eclip/K562.tsv', col_names = FALSE)
merged <- rbind(hepg2, k562)
merged <- merged[merged$X7 > 1 & merged$X8 > 1, ]

gencode <- rtracklayer::import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- GenomicFeatures::makeTxDbFromGRanges(gencode)
genes <- gencode[gencode$type == 'gene']

unique(merged$X4)
gene_name <- stringr::str_split(merged$X4, '_', simplify = TRUE)[, 1]
merged$X4 <- gene_name
sum(merged$X2 < merged$X3)
colnames(merged) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'x1', 'x2', 'x3', 'x4')
merged <- GenomicRanges::makeGRangesFromDataFrame(
  merged, keep.extra.columns = TRUE 
)
rtracklayer::export.bed(merged, 'vkr/merged.bed')

grlist_merged <- split(merged, merged$name)

library(parallel)
library(GenomicFeatures)
res <- mclapply(
  names(grlist_merged),
  function(gene_name) {
    pgene <- genes[genes$gene_name == gene_name]
    subsetByOverlaps(grlist_merged[[gene_name]], pgene)
  },
  mc.cores = 32
)
names(res) <- names(grlist_merged)
res <- GRangesList(res)
res <- unlist(res)
rtracklayer::export.bed(res, 'vkr/self_peaks.bed')

unique(res$name)

# self-rbp genes with uORFs
ut <- rtracklayer::import.bed('vkr/ut_uorfs_genomic.bed')
tis <- rtracklayer::import.bed('vkr/tisdb_uorfs_genomic.bed')
un <- c(ut, tis)
un <- un[!duplicated(un)]

genes_with_uorfs <- subsetByOverlaps(genes, un)
genes_with_uorfs <- genes_with_uorfs$gene_name
uorf_prediction <- readxl::read_excel('Supp/Supplemental_Data_Tables_.xlsx', sheet = 4, skip = 2)
unique(uorf_prediction$gene)
genes_with_uorfs <- c(genes_with_uorfs, unique(uorf_prediction$gene))
genes_with_uorfs <- unique(genes_with_uorfs)
gene_with_uORF_and_rbp <- unique(res$name[res$name %in% genes_with_uorfs])
