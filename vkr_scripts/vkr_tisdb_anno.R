library(rtracklayer)
library(GenomicFeatures)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg19)
cores <- detectCores()


tis_start_codons <- import.bed('~/BioData/tis/tis_start_codons.bed')

gencode <- import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- makeTxDbFromGRanges(gencode)
five_utr_by_tx <- fiveUTRsByTranscript(x = gencode_db, use.names = TRUE)
cds_by_tx <- cdsBy(x = gencode_db, use.names = TRUE)
cds_by_tx <- cds_by_tx[names(cds_by_tx) %in% names(five_utr_by_tx)]


mcols(tis_start_codons)$context <- getSeq(
  x = BSgenome.Hsapiens.UCSC.hg19,
  names = resize(flank(tis_start_codons, 7), 11)
)


# overlapping 5'UTRs with start codons
five_utr_by_tx_with_tis <- subsetByOverlaps(x = five_utr_by_tx, ranges = tis_start_codons, minoverlap = 3L)
# are any exons removed?
all(lengths(five_utr_by_tx[names(five_utr_by_tx) %in% names(five_utr_by_tx_with_tis)])
    == lengths(five_utr_by_tx_with_tis))
# no

# now let's map these TIS to 5'UTRs
tis_in_txcoords <- mapToTranscripts(x = tis_start_codons, transcripts = five_utr_by_tx_with_tis)
tis_in_txcoords <- unstrand(tis_in_txcoords)
mcols(tis_in_txcoords)$name <- tis_start_codons$name[tis_in_txcoords$xHits]
mcols(tis_in_txcoords)$context <- tis_start_codons$context[tis_in_txcoords$xHits]
tis_in_txcoords

# search_space
names(five_utr_by_tx_with_tis)
five_utr_by_tx_with_tis$ENST00000449309.1
cds_by_tx$ENST00000449309.1
search_ranges <- mclapply(names(five_utr_by_tx_with_tis), function(enst) {
  five_utr <- five_utr_by_tx_with_tis[[enst]]
  t_cds <- cds_by_tx[[enst]]
  return(sort.GenomicRanges(reduce(c(five_utr, t_cds)),
                            decreasing = as.character(strand(five_utr))[1] == '-'))
}, mc.cores = cores)
search_ranges <- GRangesList(search_ranges)
names(search_ranges) <- names(five_utr_by_tx_with_tis)

search_seqs <- extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19, transcripts = search_ranges)
search_seqs
seqs_five_utr_by_tx_with_tis

search_ranges$ENST00000369439.4
five_utr_by_tx$ENST00000369439.4
cds_by_tx$ENST00000369439.4
tis_in_txcoords_by_tx <- split(tis_in_txcoords, seqnames(tis_in_txcoords))
# checking if all start codons are OK
x <- mcmapply(
  function(tx) {
    s <- 0
    for (i in 1:length(tx)) {
      enst <- seqnames(tx)[i]
      search_seq <- search_seqs[[enst]]
      s <- s + (toString(subseq(search_seq, start = start(tx)[i], end = end(tx)[i])) == tx$name[i])
    }
    return(s > 0)
  },
  tis_in_txcoords_by_tx, mc.cores = cores)
sum(x) # OK

(search_seqs$ENST00000378741.3)
# searching for uorfs
tis_with_uorfs_in_txcoords_by_tx <- mclapply(
  tis_in_txcoords_by_tx,
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
keep <- sapply(tis_with_uorfs_in_txcoords_by_tx, function(tx) {any(!is.na(tx$uORF_end))})
tis_with_uorfs_in_txcoords_by_tx <- tis_with_uorfs_in_txcoords_by_tx[keep]
# enhancing
tis_with_uorfs_in_txcoords_by_tx <- GRangesList(tis_with_uorfs_in_txcoords_by_tx)
tis_with_uorfs_in_txcoords <- unlist(tis_with_uorfs_in_txcoords_by_tx, use.names = FALSE)
# filtering NA
tis_with_uorfs_in_txcoords <- tis_with_uorfs_in_txcoords[!is.na(tis_with_uorfs_in_txcoords$uORF_end)]
# filtering out very long ones
tis_with_uorfs_in_txcoords <- tis_with_uorfs_in_txcoords[tis_with_uorfs_in_txcoords$uORF_end - start(tis_with_uorfs_in_txcoords) + 1 <= 298]
sum((tis_with_uorfs_in_txcoords$uORF_end - start(tis_with_uorfs_in_txcoords) + 1) %% 3 == 0) # 5558 uORFs in transcripts

# uORFs granges in transcript coordinates
ir <- IRanges(
  start = start(tis_with_uorfs_in_txcoords),
  end = tis_with_uorfs_in_txcoords$uORF_end
)
uorfs <- GRanges(
  seqnames = seqnames(tis_with_uorfs_in_txcoords),
  ranges = ir,
  strand = strand(tis_with_uorfs_in_txcoords),
  name = tis_with_uorfs_in_txcoords$context
)
rm(tis_with_uorfs_in_txcoords, tis_with_uorfs_in_txcoords_by_tx,
   tis_start_codons, tis_in_txcoords_by_tx, tis_in_txcoords,
   seqs_five_utr_by_tx_with_tis, x, ir)


# scoring uorfs
main_codons <- import.bed('BioData/tis/main_codons.bed')
mainorf_kozak_seqs <- stringr::str_split(main_codons$name, ';', simplify = TRUE)[, 2]
kozak_pwm <- PWM(mainorf_kozak_seqs)
score(uorfs) <- sapply(uorfs$name, PWMscoreStartingAt, pwm = kozak_pwm)

mean(score(uorfs))
mean(score(uorfs[subseq(uorfs$name, 8, 10) != 'ATG']))

# saving
mean(width(uorfs))
export.bed(object = uorfs, '~/vkr/tisdb_uorfs_tx.bed')
rm(search_seqs, search_ranges, mainorf_kozak_seqs, uorfs)

# saving in genomic coords
tisdb_uorfs <- import.bed('~/vkr/tisdb_uorfs_tx.bed')
tisdb_uorfs_genomic <- mapFromTranscripts(x = tisdb_uorfs, transcripts = search_ranges)
tisdb_uorfs_genomic$xHits
tisdb_uorfs_genomic$name <- tisdb_uorfs$name[tisdb_uorfs_genomic$xHits]
tisdb_uorfs_genomic$score <- tisdb_uorfs$score[tisdb_uorfs_genomic$xHits]
tisdb_uorfs_genomic$xHits <- NULL
tisdb_uorfs_genomic$transcriptsHits <- NULL
tisdb_uorfs_genomic <- tisdb_uorfs_genomic[!duplicated(tisdb_uorfs_genomic)]
export.bed(tisdb_uorfs_genomic, '~/vkr/tisdb_uorfs_genomic.bed')

#### Integrating all the data ####

# fragments preparation
library(pbapply)

frag_ranges <- mclapply(
  names(cds_by_tx),
  function(tx_name) {
    utr_exons <- five_utr_by_tx[[tx_name]]
    cds_exon <- cds_by_tx[[tx_name]][1]
    return(reduce(c(utr_exons, cds_exon)))
  },
  mc.cores = cores
)
names(frag_ranges) <- names(cds_by_tx)
frag_ranges <- GRangesList(frag_ranges)
frag_seqs <- extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19, transcripts = frag_ranges)
not_keep <- !duplicated(frag_seqs)
frag_seqs <- frag_seqs[not_keep]
frag_ranges <- frag_ranges[not_keep]
writeXStringSet(x = frag_seqs, filepath = '~/BioData/transcriptomes/target_fragments.fa')
export.bed(object = frag_ranges, con = '~/BioData/transcriptomes/target_fragments.bed')

# shRNA genes quant preparation 
shRNA_table <- data.frame(readr::read_tsv('BioData/shRNA/shRNA_table.tsv', col_names = FALSE))[-7]
colnames(shRNA_table) <- c('knockdown1', 'knockdown2', 'control1', 'control2', 'cell_line', 'gene_name')
quant_targets <- shRNA_table[shRNA_table$gene_name %in% target_genes$gene_name, ]
readr::write_tsv(x = quant_targets, path = '~/BioData/shRNA/quant_targets.tsv')
writeLines(unlist(quant_targets[1:4]), '~/BioData/shRNA/bam_ids.txt')

# uORF data
uorf_genomic <- mapFromTranscripts(x = uorfs, transcripts = five_utr_by_tx)
mcols(uorf_genomic)$name <- uorfs$name
mcols(uorf_genomic)$score <- uorfs$score
uorf_genomic <- uorf_genomic[!duplicated(uorf_genomic)]
export.bed(object = uorf_genomic, '~/BioData/tis/uorfs_in_genomic_coords.bed')

# target_counts <- target_counts %>% 
#   dplyr::group_by(gene_id) %>% 
#   dplyr::mutate(
#     geneTPM_control = sum(fragTPM_control),
#     geneTPM_xrn1smg6 = sum(fragTPM_xrn1smg6),
#     geneTPM_xrn1upf1 = sum(fragTPM_xrn1upf1)
#   ) %>% dplyr::ungroup()
# y$width <- NULL
# y$strand <- NULL
# colnames(y) <- c('tx_name', 'start', 'end', 'context', 'score')
# target_counts <- dplyr::left_join(target_counts, y, by = 'tx_name')
# # adding main kozak score
# kozak_data <- import.bed('~/BioData/tis/main_codons.bed')
# z <- data.frame(
#   tx_name = substr(kozak_data$name, 1, 17),
#   main_score = kozak_data$score
# )
# target_counts <- dplyr::left_join(target_counts, z, by = 'tx_name')

# calculating delta
target_counts$deltaUPF1 <- `-`(
  log10(target_counts$fragTPM_xrn1upf1 / target_counts$geneTPM_xrn1upf1),
  log10(target_counts$fragTPM_control / target_counts$geneTPM_control) 
)
target_counts$deltaSMG6 <- `-`(
  log10(target_counts$fragTPM_xrn1smg6 / target_counts$geneTPM_xrn1smg6),
  log10(target_counts$fragTPM_control / target_counts$geneTPM_control)
)

# shRNA
shRNA_data <- data.frame(readr::read_tsv('~/BioData/shRNA/total.tsv'))
rownames(shRNA_data) <- shRNA_data$tx_name
shRNA_data$tx_name <- NULL
shRNA_meta <- readr::read_tsv('~/BioData/shRNA/quant_targets.tsv')
eps <- min(apply(shRNA_data, 2, function(x) min(x[x > 0])))
shRNA_data <- shRNA_data + eps
target_counts$fragTPM_gene_control <- rep(NA, nrow(target_counts))
target_counts$fragTPM_gene_knockdown <- rep(NA, nrow(target_counts))
target_counts$geneTPM_gene_control <- rep(NA, nrow(target_counts))
target_counts$geneTPM_gene_knockdown <- rep(NA, nrow(target_counts))

for (gene_name in unique(shRNA_meta$gene_name)) {
  pos <- which(shRNA_meta$gene_name == gene_name)[1]
  ids <- as.character(shRNA_meta[pos, 1:4])
  gene_id <- target_genes$gene_id[which(target_genes$gene_name == gene_name)[1]]
  target_rows <- which(target_counts$gene_id == gene_id)
  target_tx <- target_counts$tx_name[target_rows]
  df <- shRNA_data[target_tx, ids]
  target_counts$fragTPM_gene_knockdown[target_rows] <- rowSums(df[1:2])
  target_counts$fragTPM_gene_control[target_rows] <- rowSums(df[3:4])
  target_counts$geneTPM_gene_control[target_rows] <- sum(rowSums(df[1:2]))
  target_counts$geneTPM_gene_knockdown[target_rows] <- sum(rowSums(df[3:4]))
}

target_counts$deltaGene <- `-`(
  log10(target_counts$fragTPM_gene_knockdown / target_counts$geneTPM_gene_knockdown),
  log10(target_counts$fragTPM_gene_control / target_counts$geneTPM_gene_control)
)
sum(!is.na(target_counts$deltaGene))

# targeting by uorf
target_frags <- frag_ranges[names(frag_ranges) %in% target_counts$tx_name]
uorf_targets <- subset(target_counts, contains_uORF)
non_uorf_targets <- subset(target_counts, !contains_uORF)

uorf_frags <- target_frags[names(target_frags) %in% uorf_targets$tx_name]
uorf_frags <- sort.GRangesList(uorf_frags)
uorf_frags
non_uorf_frags <- target_frags[!(names(target_frags) %in% uorf_targets$tx_name)]
non_uorf_frags <- sort.GRangesList(non_uorf_frags)
non_uorf_frags

mcols(uorf_frags)$thick <- ranges(mapFromTranscripts(x = uorfs, transcripts = uorf_frags))

old_names <- names(uorf_frags)
new_names <- character(length(uorf_frags))
for (g_pos in 1:length(uorf_frags)) {
  tx_name <- old_names[g_pos]
  table_pos <- which(tx_name == uorf_targets$tx_name)
  new_names[g_pos] <- paste0('deltaUPF1=', round(uorf_targets$deltaUPF1[table_pos], 2),
                             '_z=', round(uorf_targets$score[table_pos], 2))
}
names(uorf_frags) <- new_names
export.bed(uorf_frags, '~/BioData/tracks/uORF_deltaUPF1.bed')
for (g_pos in 1:length(uorf_frags)) {
  tx_name <- old_names[g_pos]
  table_pos <- which(tx_name == uorf_targets$tx_name)
  new_names[g_pos] <- paste0('deltaSMG6=', round(uorf_targets$deltaSMG6[table_pos], 2),
                             '_z=', round(uorf_targets$score[table_pos], 2))
}
names(uorf_frags) <- new_names
export.bed(uorf_frags, '~/BioData/tracks/uORF_deltaSMG6.bed')

mcols(non_uorf_frags)$thick <- IRanges(start = min(start(non_uorf_frags)),
                                       end = min(start(non_uorf_frags)) - 1)
old_names <- names(non_uorf_frags)
new_names <- character(length(non_uorf_frags))
for (tx_name in old_names) {
  pos <- which(tx_name == non_uorf_targets$tx_name)
  new_names[pos] <- paste0('deltaUPF1=', round(non_uorf_targets$deltaUPF1[pos], 2))
}
names(non_uorf_frags) <- new_names
export.bed(non_uorf_frags, '~/BioData/tracks/non_uORF_deltaUPF1.bed')
for (tx_name in old_names) {
  pos <- which(tx_name == non_uorf_targets$tx_name)
  new_names[pos] <- paste0('deltaSMG6=', round(non_uorf_targets$deltaSMG6[pos], 2))
}
names(non_uorf_frags) <- new_names
export.bed(non_uorf_frags, '~/BioData/tracks/non_uORF_deltaSMG6.bed')

sub_uorf_targets <- uorf_targets[!is.na(uorf_targets$deltaGene), ]
sub_non_uorf_targets <- non_uorf_targets[!is.na(non_uorf_targets$deltaGene), ]
uorf_frags <- uorf_frags[names(uorf_frags) %in% sub_uorf_targets$tx_name]
non_uorf_frags <- non_uorf_frags[names(non_uorf_frags) %in% sub_non_uorf_targets$tx_name]
names(uorf_frags) == sub_uorf_targets$tx_name
names(non_uorf_frags) == sub_non_uorf_targets$tx_name

old_names <- names(uorf_frags)
new_names <- sapply(old_names, function(tx_name) {
  pos <- which(tx_name == sub_uorf_targets$tx_name)
  paste0('deltaGene=', round(sub_uorf_targets$deltaGene[pos], 2),
         '_z=', round(sub_uorf_targets$score[pos], 2))
}) 
names(uorf_frags) <- new_names
export.bed(uorf_frags, '~/BioData/tracks/uORF_deltaGene.bed')

old_names <- names(non_uorf_frags)
new_names <- sapply(old_names, function(tx_name) {
  pos <- which(tx_name == sub_non_uorf_targets$tx_name)
  paste0('deltaGene=', round(sub_non_uorf_targets$deltaGene[pos], 2))
}) 
names(non_uorf_frags) <- new_names
export.bed(non_uorf_frags, '~/BioData/tracks/non_uORF_deltaGene.bed')

# eCLIP
peaks_merged <- import.bed('BioData/eclip/peaks_merged.bed')
peaks_merged <- split(peaks_merged, peaks_merged$name)
gene_by_name <- gr_gencode[gr_gencode$type == 'gene' & gr_gencode$gene_name %in% names(peaks_merged)]
old_names <- gene_by_name$gene_name
mcols(gene_by_name) <- NULL
gene_by_name$name <- old_names
gene_by_name <- split(gene_by_name, gene_by_name$name)

result <- pbapply::pblapply(names(gene_by_name), function(gene) {
  subsetByOverlaps(peaks_merged[[gene]], gene_by_name[[gene]])
})
names(result) <- names(gene_by_name)
result <- GRangesList(result)
export.bed(unlist(result), '~/BioData/tracks/RBP_peaks.bed')

seqlengths(BSgenome.Hsapiens.UCSC.hg19)

# plots
target_counts
library(tidyverse)
target_table <- map_dfr(unique(target_genes$gene_id), function(gene) {
  sub_df <- target_counts[target_counts$gene_id == gene, ]
  z <- tibble(
    gene_id = gene,
    gene_name = target_genes$gene_name[target_genes$gene_id == gene][1],
    n_uorf = sum(sub_df$contains_uORF),
    n_non_uorf = nrow(sub_df) - n_uorf,
    tpm_uorf_control = sum(sub_df$fragTPM_control[sub_df$contains_uORF]),
    tpm_uorf_smg6 = sum(sub_df$fragTPM_xrn1smg6[sub_df$contains_uORF]),
    tpm_uorf_upf1 = sum(sub_df$fragTPM_xrn1upf1[sub_df$contains_uORF]),
    tpm_gene_control = sub_df$geneTPM_control[1],
    tpm_gene_smg6 = sub_df$geneTPM_xrn1smg6[1],
    tpm_gene_upf1 = sub_df$geneTPM_xrn1upf1[1],
    tpm_uorf_shRNA_control = sum(sub_df$fragTPM_gene_control[sub_df$contains_uORF]),
    tpm_uorf_shRNA_KD = sum(sub_df$fragTPM_gene_knockdown[sub_df$contains_uORF]),
    tpm_gene_shRNA_control = sub_df$geneTPM_gene_control[1],
    tpm_gene_shRNA_KD = sub_df$geneTPM_gene_knockdown[1],
    shRNA = all(!is.na(sub_df$deltaGene)),
    average_zu = mean(sub_df$score, na.rm = TRUE),
  )
  return(z)
})
qplot(
  x = log10(target_table$tpm_gene_smg6 + target_table$tpm_gene_control),
  y = log10(target_table$tpm_uorf_smg6 / target_table$tpm_gene_smg6) -
    log10(target_table$tpm_uorf_control / target_table$tpm_gene_control),
  main = 'SMG6',
  xlab = 'Avg. expression in KD and control',
  ylab = 'Delta'
)
qplot(
  x = log10(target_table$tpm_gene_upf1 + target_table$tpm_gene_control),
  y = log10(target_table$tpm_uorf_upf1 / target_table$tpm_gene_upf1) -
    log10(target_table$tpm_uorf_control / target_table$tpm_gene_control),
  main = 'UPF1',
  xlab = 'Avg. expression in KD and control',
  ylab = 'Delta'
)
qplot(
  x = log10(target_table$tpm_gene_shRNA_control + target_table$tpm_gene_shRNA_KD),
  y = log10(target_table$tpm_uorf_shRNA_KD / target_table$tpm_gene_shRNA_KD) -
    log10(target_table$tpm_uorf_shRNA_control / target_table$tpm_gene_shRNA_control),
  main = 'shRNA',
  xlab = 'Avg. expression in KD and control',
  ylab = 'Delta'
)

qplot(
  x = log10(target_table$tpm_gene_smg6 / target_table$tpm_gene_control),
  y = log10(target_table$tpm_uorf_smg6 / target_table$tpm_gene_smg6) -
    log10(target_table$tpm_uorf_control / target_table$tpm_gene_control),
  main = 'SMG6',
  xlab = 'Expression change in KD and control',
  ylab = 'Delta'
)
qplot(
  x = log10(target_table$tpm_gene_upf1 / target_table$tpm_gene_control),
  y = log10(target_table$tpm_uorf_upf1 / target_table$tpm_gene_upf1) -
    log10(target_table$tpm_uorf_control / target_table$tpm_gene_control),
  main = 'UPF1',
  xlab = 'Expression change in KD and control',
  ylab = 'Delta'
)
qplot(
  x = log10(target_table$tpm_gene_shRNA_KD / target_table$tpm_gene_shRNA_control),
  y = log10(target_table$tpm_uorf_shRNA_KD / target_table$tpm_gene_shRNA_KD) -
    log10(target_table$tpm_uorf_shRNA_control / target_table$tpm_gene_shRNA_control),
  main = 'shRNA',
  xlab = 'Expression change in KD and control',
  ylab = 'Delta'
)

x <- log10(target_table$tpm_uorf_upf1 / target_table$tpm_gene_upf1) -
  log10(target_table$tpm_uorf_control / target_table$tpm_gene_control)
y <- log10(target_table$tpm_uorf_smg6 / target_table$tpm_gene_smg6) -
  log10(target_table$tpm_uorf_control / target_table$tpm_gene_control)

rbp_genes <- names(peaks_merged)

qplot(y = x, x = target_table$average_zu, xlab = 'Z-score', ylab = 'Delta', main = 'UPF1')
qplot(y = y, x = target_table$average_zu, xlab = 'Z-score', ylab = 'Delta', main = 'SMG6')
names(peaks_merged)
p1 <- order(x, decreasing = TRUE)[1:200]
p2 <- order(y, decreasing = TRUE)[1:200]
min(x[p1])
min(y[p2])
interesting <- target_table[intersect(p1, p2), ]$gene_name
cat(
  interesting, sep = '\n'
)
intersect(interesting, rbp_genes)
target_table



mapFromTranscripts(uorfs, five_utr_by_tx_with_tis)


duplicated(unlist(uorf_frags))
test <- GenomicRanges::distanceToNearest(unlist(uorf_frags), unlist(result))
test <- data.frame(test)
k <- unlist(uorf_frags)
mcols(k) <- NULL
test
k <- k[test$queryHits]
k$distance <- test$distance
k
target_genes
