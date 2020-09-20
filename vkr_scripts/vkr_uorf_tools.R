library(GenomicFeatures)
library(rtracklayer)
library(parallel)
cores <- detectCores()

gencode <- import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- makeTxDbFromGRanges(gencode)
five_utr_by_tx <- fiveUTRsByTranscript(x = gencode_db, use.names = TRUE)
cds_by_tx <- cdsBy(x = gencode_db, use.names = TRUE)
full_tx <- exonsBy(x = gencode_db, by = 'tx', use.names = TRUE)
full_tx <- full_tx[names(full_tx) %in% names(five_utr_by_tx)]


uorf_tools <- readxl::read_excel('Supp/journal.pone.0222459.s003.xlsx')
mean(uorf_tools$ORF_length)
substr(uorf_tools$transcript_id, 1, 15) %in% substr(names(cds_by_tx), 1, 15)
uorf_tools
uorfs_ut <- data.frame(
  chr = uorf_tools$chromosome,
  start = uorf_tools$start,
  stop = uorf_tools$stop,
  name = paste0('ut_', uorf_tools$gene_symbol),
  score = 0,
  strand = uorf_tools$strand,
  len = uorf_tools$ORF_length
)
uorfs_ut <- makeGRangesFromDataFrame(uorfs_ut, keep.extra.columns = TRUE)
ch <- import.chain('Supp/hg38ToHg19.over.chain')
uorfs_ut <- liftOver(uorfs_ut, ch)
uorfs_ut
uorfs_ut <- reduce(uorfs_ut)
names(uorfs_ut) <- paste0('ut_', uorf_tools$gene_symbol, '_', substr(uorf_tools$transcript_id, 1, 15), '_', substr(uorf_tools$gene_id, 1, 15))
uorfs_ut <- uorfs_ut[lengths(uorfs_ut) <= 1]
uorfs_ut <- unlist(uorfs_ut)
uorfs_ut$name <- names(uorfs_ut)
uorfs_ut <- sort.GenomicRanges(uorfs_ut)
any(duplicated(uorfs_ut)) # no duplicates

# preparing BED of uORFs
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
main_codons <- import.bed('BioData/tis/main_codons.bed')
mainorf_kozak_seqs <- stringr::str_split(main_codons$name, ';', simplify = TRUE)[, 2]
kozak_pwm <- PWM(mainorf_kozak_seqs)
uorfs_bed <- uorfs_ut
uorfs_bed$name <- getSeq(BSgenome.Hsapiens.UCSC.hg19, resize(flank(uorfs_ut, 7), 11))
uorfs_bed$score <- sapply(uorfs_bed$name, PWMscoreStartingAt, pwm = kozak_pwm)
uorfs_bed
export.bed(object = uorfs_bed, con = 'vkr/raw_ut_uorfs.bed')

# searching for uoORF
uorfs_ut

res <- mclapply(1:length(uorfs_ut), function(position) {
  enst <- stringr::str_split(uorfs_ut$name[position], '_', simplify = TRUE)[3]
  enst <- substr(enst, 1, 15)
  target_cds <- cds_by_tx[substr(names(cds_by_tx), 1, 15) == enst]
  findOverlaps(uorfs_ut[position], target_cds)
}, mc.cores = cores)
lengths(res)
res[5]
uorfs_ut
# yes, they are present: ZBTB40, LYPLA2, etc

# now to start codons
width(uorfs_bed)
ut_starts <- resize(uorfs_bed, 3)
ut_starts
five_utr_by_tx_with_start <- subsetByOverlaps(x = five_utr_by_tx, ranges = ut_starts, minoverlap = 3L)

# mapping starts to 5'-UTR
uts_in_tx <- mapToTranscripts(x = ut_starts, transcripts = five_utr_by_tx_with_start)
uts_in_tx <- unstrand(uts_in_tx)
ut_starts
mcols(uts_in_tx)$name <- ut_starts$name[uts_in_tx$xHits]
mcols(uts_in_tx)$score <- ut_starts$score[uts_in_tx$xHits]
uts_in_tx

# search space
search_ranges <- mclapply(names(five_utr_by_tx_with_start), function(enst) {
  five_utr <- five_utr_by_tx_with_start[[enst]]
  t_cds <- cds_by_tx[[enst]]
  return(sort.GenomicRanges(reduce(c(five_utr, t_cds)),
                            decreasing = as.character(strand(five_utr))[1] == '-'))
}, mc.cores = cores)
search_ranges <- GRangesList(search_ranges)
names(search_ranges) <- names(five_utr_by_tx_with_start)

search_seqs <- extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19, transcripts = search_ranges)
search_seqs

# splitting
uts_in_tx_by_tx <- split(uts_in_tx, seqnames(uts_in_tx))
uts_in_tx_by_tx

# checking if all start codons are OK
x <- mcmapply(
  function(tx) {
    s <- 0
    for (i in 1:length(tx)) {
      enst <- seqnames(tx)[i]
      search_seq <- search_seqs[[enst]]
      s <- s + (toString(subseq(search_seq, start = start(tx)[i], end = end(tx)[i])) == substr(tx$name[i], 8, 10))
    }
    return(s)
  },
  uts_in_tx_by_tx, mc.cores = cores)
sum(x) # all OK

# searching for uorfs
uorf_uts_in_tx_by_tx <- mclapply(
  uts_in_tx_by_tx,
  function(tx) {
    for (i in 1:length(tx)) {
      enst <- seqnames(tx)[i]
      search_seq <- search_seqs[[enst]]
      s <- start(tx)[i]
      sub_search_seq <- toString(subseq(search_seq, start = s, end = length(search_seq) - 3))
      tx$uORF_end[i] <- stringr::str_locate(
        string = sub_search_seq,
        pattern = '^([ACTG]{3})+?((TAA)|(TAG)|(TGA))'
      )[2] + s - 1
    }
    return(tx)
  },
  mc.cores = cores
)

# filtering
keep <- sapply(uorf_uts_in_tx_by_tx, function(tx) {any(!is.na(tx$uORF_end))})
uorf_uts_in_tx_by_tx <- uorf_uts_in_tx_by_tx[keep]

# enhancing
uorf_uts_in_tx_by_tx <- GRangesList(uorf_uts_in_tx_by_tx)
uorf_uts_in_tx_by_tx <- unlist(uorf_uts_in_tx_by_tx, use.names = FALSE)
uorf_uts_in_tx_by_tx

# filtering NA
uorf_uts_in_tx_by_tx <- uorf_uts_in_tx_by_tx[!is.na(uorf_uts_in_tx_by_tx$uORF_end)]
uorf_uts_in_tx_by_tx <- uorf_uts_in_tx_by_tx[uorf_uts_in_tx_by_tx$uORF_end - start(uorf_uts_in_tx_by_tx) + 1 <= 298]
sum((uorf_uts_in_tx_by_tx$uORF_end - start(uorf_uts_in_tx_by_tx) + 1) %% 3 == 0)

# tx coords
ir <- IRanges(
  start = start(uorf_uts_in_tx_by_tx),
  end = uorf_uts_in_tx_by_tx$uORF_end
)
uorfs <- GRanges(
  seqnames = seqnames(uorf_uts_in_tx_by_tx),
  ranges = ir,
  strand = strand(uorf_uts_in_tx_by_tx),
  name = uorf_uts_in_tx_by_tx$name,
  score = uorf_uts_in_tx_by_tx$score
)

# saving
export.bed(object = uorfs, '~/vkr/ut_uorfs_tx.bed')

# genomic
ut_gene <- mapFromTranscripts(x = uorfs, transcripts = search_ranges)
ut_gene$xHits
ut_gene$name <- uorfs$name[ut_gene$xHits]
ut_gene$score <- uorfs$score[ut_gene$xHits]
ut_gene <- ut_gene[!duplicated(ut_gene)]
export.bed(object = ut_gene, con = 'vkr/ut_uorfs_genomic.bed')
