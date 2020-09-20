library(GenomicFeatures)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(parallel)
cores <- detectCores() - 2

gencode <- import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- makeTxDbFromGRanges(gencode)
five_tx <- fiveUTRsByTranscript(x = gencode_db, use.names = TRUE)
cds_tx <- cdsBy(x = gencode_db, use.names = TRUE)
cds_tx <- cds_tx[names(cds_tx) %in% names(five_tx)]
full_tx <- exonsBy(x = gencode_db, by = 'tx', use.names = TRUE)
full_tx <- full_tx[names(full_tx) %in% names(five_tx)]

# uorf prediction data
uorf_pred <- readxl::read_excel('Supp/Supplemental_Data_Tables_.xlsx', sheet = 4, skip = 2)
str(uorf_pred)
uorf_pred$uORF_start_coordinate < uorf_pred$uORF_end_coordinate
uorf_pred

# making granges
upred <- data.frame(
  chr = uorf_pred$chromosome,
  start = ifelse(uorf_pred$uORF_start_coordinate < uorf_pred$uORF_end_coordinate,
                 uorf_pred$uORF_start_coordinate, uorf_pred$uORF_end_coordinate),
  stop = ifelse(uorf_pred$uORF_start_coordinate < uorf_pred$uORF_end_coordinate,
                uorf_pred$uORF_end_coordinate, uorf_pred$uORF_start_coordinate),
  name = uorf_pred$uORF_ID,
  score = 0,
  strand = uorf_pred$strand
)
upred <- makeGRangesFromDataFrame(upred, keep.extra.columns = TRUE)
export.bed(object = upred, con = 'vkr/upred/raw_upred_genomic.bed')

# taking out start codons
upred_starts <- resize(upred, 3)
upred_starts <- upred_starts[!duplicated(upred_starts)]

# 5utrs with starts
utr_upred <- subsetByOverlaps(x = five_tx, ranges = upred_starts, minoverlap = 3L)

# mapping starts to 5'-UTR
upred_starts_tx <- mapToTranscripts(x = upred_starts,
                                    transcripts = utr_upred)
upred_starts_tx <- upred_starts_tx[width(upred_starts_tx) == 3]
upred_starts_tx <- unstrand(upred_starts_tx)

# transcript seqs
tx_seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19,
                      full_tx[names(full_tx) %in% seqnames(upred_starts_tx)])

# scoring starts
main_codons <- import.bed('BioData/tis/main_codons.bed')
mainorf_kozak_seqs <- stringr::str_split(main_codons$name, ';', simplify = TRUE)[, 2]
kozak_pwm <- PWM(mainorf_kozak_seqs)
uorf_kozak_seqs <- mcmapply(
  function(i) {
    s <- start(upred_starts_tx)[i]
    e <- end(upred_starts_tx)[i]
    hit <- upred_starts_tx$xHits[i]
    if (s < 8) {
      kozak_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19,
                          resize(flank(upred_starts[hit], 7), 11))
    } else {
      tx_name <- as.character(upred_starts_tx[i]@seqnames)
      kozak_seq <- tx_seqs[[tx_name]][(s-7):(e+1)]
    }
    return(as.character(kozak_seq))
  },
  1:length(upred_starts_tx),
  mc.cores = cores
)
upred_starts_tx$name <- uorf_kozak_seqs
upred_starts_tx$score <- sapply(upred_starts_tx$name, PWMscoreStartingAt, pwm = kozak_pwm)

# search space
search_ranges <- mclapply(
  as.character(upred_starts_tx@seqnames), 
  function(tx_name) {
  t_utr <- five_tx[[tx_name]]
  t_cds <- cds_tx[[tx_name]]
  return(
    sort.GenomicRanges(
      reduce(c(t_utr, t_cds)),
      decreasing = as.character(strand(t_utr))[1] == '-'
      )
    )
  },
  mc.cores = cores
)
search_ranges <- GRangesList(search_ranges)
names(search_ranges) <- as.character(upred_starts_tx@seqnames)

# search sequences
search_seqs <- extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19,
                                     transcripts = search_ranges)
search_seqs

# check if start codons are OK
sum(`==`(
  subseq(search_seqs, start = start(upred_starts_tx), end = 1 + end(upred_starts_tx)),
  substr(upred_starts_tx$name, 8, 11)
  )
  )
# OK


# searching for uORFs
uorfs <- stringr::str_locate(
  subseq(search_seqs, start = start(upred_starts_tx), end = width(search_seqs) - 3),
  pattern = '^([ACTG]{3})+?((TAA)|(TAG)|(TGA))'
) + start(upred_starts_tx) - 1
uorfs <- as.data.frame(uorfs)
df <- as.data.frame(upred_starts_tx)[, -(2:4)]
uorfs <- cbind(uorfs, df)


# filtering
uorfs <- uorfs[!is.na(uorfs$end), ]
library(magrittr)
uorfs <- uorfs %>%
  dplyr::filter(end - start + 1 <= 300)
uorfs_tx <- makeGRangesFromDataFrame(uorfs, keep.extra.columns = TRUE)
uorfs_tx
export.bed(uorfs_tx, 'vkr/upred/upred_tx.bed')

# genomic
uorfs_genomic <- mapFromTranscripts(uorfs_tx, full_tx)
uorfs_genomic$name <- uorfs_tx$name[uorfs_genomic$xHits]
uorfs_genomic$score <- uorfs_tx$score[uorfs_genomic$xHits]
uorfs_genomic <- uorfs_genomic[!duplicated(uorfs_genomic)]
export.bed(uorfs_genomic, 'vkr/upred/upred_genomic.bed')
