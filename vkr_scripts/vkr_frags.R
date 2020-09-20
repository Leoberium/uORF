library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(parallel)
cores <- detectCores()

# importing uorfs
ut_tx <- import.bed('vkr/ut_uorfs_tx.bed')
ut_gene <- import.bed('vkr/ut_uorfs_genomic.bed')
tisdb_tx <- import.bed('vkr/tisdb_uorfs_tx.bed')
tisdb_gene <- import.bed('vkr/tisdb_uorfs_genomic.bed')

# importing gencode
gencode <- import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- makeTxDbFromGRanges(gencode)

# uniting genomic uORFs
united_gene <- c(ut_gene, tisdb_gene)
sum(duplicated(united_gene)) # 274 duplicates
united_gene <- united_gene[!duplicated(united_gene)] # removing them

# loading transcript ranges
full_tx <- exonsBy(x = gencode_db, by = 'tx', use.names = TRUE)
five_tx <- fiveUTRsByTranscript(x = gencode_db, use.names = TRUE)
cds_tx <- cdsBy(x = gencode_db, by = 'tx', use.names = TRUE)

# filtering tx to have 5'UTR
full_tx <- full_tx[names(full_tx) %in% names(five_tx)]
cds_tx <- cds_tx[names(cds_tx) %in% names(five_tx)]

# mapping uORF to transcripts
united_tx <- c(ut_tx, tisdb_tx)

# filtering out duplicates
united_tx <- united_tx[!duplicated(united_tx)]

# leaving only one uORF per tx
library(magrittr)
uorfs <- GenomicRanges::makeGRangesFromDataFrame(
  df = data.frame(united_tx) %>% 
    dplyr::group_by(seqnames) %>% 
    dplyr::top_n(n = 1, wt = score) %>% dplyr::ungroup(),
  keep.extra.columns = TRUE
)
as.character(seqnames(uorfs))

# determining corresponding genes
enst_to_ensg <- AnnotationDbi::select(x = gencode_db, keys = as.character(seqnames(uorfs)), 
                       columns = c('GENEID', 'TXNAME'), keytype = 'TXNAME')
ensg_to_enst <- AnnotationDbi::select(x = gencode_db, keys = unique(enst_to_ensg$GENEID),
                       columns = c('GENEID', 'TXNAME'), keytype = 'GENEID')

# filtering out only uORF genes
not_keep <- sapply(ensg_to_enst$GENEID, function(gene_id) {
  tx_names <- ensg_to_enst[ensg_to_enst$GENEID == gene_id, ]$TXNAME
  return(all(tx_names %in% seqnames(uorfs)))
})
ensg_to_enst <- ensg_to_enst[!not_keep, ]

# Now we have all target genes and their transcripts
# Preparing frags for quantification
five_tx <- five_tx[names(five_tx) %in% ensg_to_enst$TXNAME]
cds_tx <- cds_tx[names(cds_tx) %in% ensg_to_enst$TXNAME]
frag_ranges <- mclapply(
  names(five_tx),
  function(tx_name) {
    res <- reduce(c(five_tx[[tx_name]], cds_tx[[tx_name]][1]))
    res <- sort.GenomicRanges(res, decreasing = as.character(strand(res))[1] == '-')
    return(res)
  },
  mc.cores = cores
)
frag_ranges <- GRangesList(frag_ranges)
names(frag_ranges) <- names(five_tx)
export.bed(frag_ranges, 'vkr/raw_ranges.bed')

# filtering
ensg_to_enst <- ensg_to_enst[ensg_to_enst$TXNAME %in% names(frag_ranges), ]
uorfs <- uorfs[seqnames(uorfs) %in% ensg_to_enst$TXNAME]
frag_ranges <- frag_ranges[names(frag_ranges) %in% ensg_to_enst$TXNAME]
export.bed(frag_ranges, 'filtered_ranges.bed')

# seqs
fragmentome <- extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19, 
                                     transcripts = frag_ranges)
extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19,
                      transcripts = five_tx[names(five_tx) %in% names(frag_ranges)])
fragmentome <- fragmentome[!duplicated(fragmentome)]
writeXStringSet(x = fragmentome, filepath = 'vkr/fragmentome.fa')
frag_ranges <- frag_ranges[names(frag_ranges) %in% names(fragmentome)]
export.bed(object = frag_ranges, con = 'vkr/filtered_ranges.bed')

# filtering uorfs
filtered_ranges <- import.bed('vkr/filtered_ranges.bed')
filt_uorfs <- uorfs[seqnames(uorfs) %in% filtered_ranges$name]
export.bed(filt_uorfs, 'vkr/target_uorfs.bed')
