library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(parallel)
cores <- detectCores() - 2

# importing gencode
gencode <- import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- makeTxDbFromGRanges(gencode)

# loading transcript ranges
full_tx <- exonsBy(x = gencode_db, by = 'tx', use.names = TRUE)
five_tx <- fiveUTRsByTranscript(x = gencode_db, use.names = TRUE)
cds_tx <- cdsBy(x = gencode_db, by = 'tx', use.names = TRUE)

# filtering tx to have 5'UTR
full_tx <- full_tx[names(full_tx) %in% names(five_tx)]
cds_tx <- cds_tx[names(cds_tx) %in% names(five_tx)]

# preparing frags for quantification
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
export.bed(frag_ranges, 'vkr/frag_ranges.bed')

# sequences
frag_seqs <- extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19,
                                   transcripts = frag_ranges)
frag_seqs <- frag_seqs[!duplicated(frag_seqs)]
writeXStringSet(x = frag_seqs, filepath = 'vkr/frag_seqs.fa')
frag_ranges <- frag_ranges[names(frag_ranges) %in% names(frag_seqs)]
export.bed(frag_ranges, 'vkr/raw_ranges.bed')
