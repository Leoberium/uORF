library(rtracklayer)
library(GenomicFeatures)
library(parallel)
cores <- detectCores()

#### Working with start codons from TISdb ####

tisdb <- readxl::read_excel(path = '~/BioData/tis/human_tisdb_data_1.0.xlsx')
dplyr::glimpse(tisdb)
tis_start_codons <- dplyr::tibble(
  chr = tisdb$Chr,
  start = tisdb$`Coordinate start codon`,
  stop = tisdb$`Coordinate start codon`,
  strand = tisdb$Strand,
  name = tisdb$`Start Codon`
)
tis_start_codons <- makeGRangesFromDataFrame(df = tis_start_codons, keep.extra.columns = TRUE)
tis_start_codons <- resize(tis_start_codons, 3)
tis_start_codons <- tis_start_codons[!duplicated(tis_start_codons)]


# Now we have all start codons from Ribo-Seq experiment without duplicates
tis_start_codons <- sort.GenomicRanges(tis_start_codons)
export.bed(object = tis_start_codons, con = '~/BioData/tis/tis_start_codons.bed')
rm(tisdb)

#### Working with gencode annotation ####

gencodev19 <- import.gff3(con = 'BioData/anno/gz/gencode.v19.annotation.gff3.gz')

# filtering types
table(gencodev19$type)
keep <- gencodev19$type %in% c('gene', 'transcript', 'exon', 'CDS', 'UTR')
sum(keep)
gencodev19 <- gencodev19[keep]

# filtering gene types
table(gencodev19$gene_type)
keep <- gencodev19$gene_type == 'protein_coding'
sum(keep)
gencodev19 <- gencodev19[keep]

# filtering gene status
table(gencodev19$gene_status)
keep <- gencodev19$gene_status == 'KNOWN'
sum(keep)
gencodev19 <- gencodev19[keep]

# filtering out tags mRNA_end_NF, mRNA_start_NF, cds_end_NF, cds_start_NF
keep_out <- sum(gencodev19$tag == 'mRNA_start_NF')
sum(keep_out)
gencodev19 <- gencodev19[!keep_out]
keep_out <- sum(gencodev19$tag == 'cds_start_NF')
sum(keep_out)
keep_out <- sum(gencodev19$tag == 'mRNA_end_NF')
sum(keep_out)
gencodev19 <- gencodev19[!keep_out]
keep_out <- sum(gencodev19$tag == 'cds_end_NF')
sum(keep_out)

# filtering transcript status
table(gencodev19$transcript_status)
keep <- gencodev19$transcript_status == 'KNOWN'
sum(keep)
gencodev19 <- gencodev19[keep]

# filtering transcript type
table(gencodev19$transcript_type)
keep_out <- gencodev19$transcript_type %in% c('non_stop_decay', 'nonsense_mediated_decay', 'processed_transcript')
sum(keep_out)
gencodev19 <- gencodev19[!keep_out]

# now we have filtered annotation
export.gff3(object = gencodev19, con = 'BioData/anno/gencode.v19.filtered.annotation.gff3')

#### uORFs ####

tis_start_codons <- import.bed(con = '~/BioData/tis/tis_start_codons.bed')
gencodev19_db <- makeTxDbFromGFF('~/BioData/anno/gencode.v19.filtered.annotation.gff3')

# context for tis
library(BSgenome.Hsapiens.UCSC.hg19)
mcols(tis_start_codons)$context <- getSeq(
  x = BSgenome.Hsapiens.UCSC.hg19,
  names = reduce(c(flank(tis_start_codons, 7), resize(tis_start_codons, 4)))
)

# transcripts, their 5'UTRs and CDS
five_utr_by_tx <- fiveUTRsByTranscript(x = gencodev19_db, use.names = TRUE)
cds_by_tx <- cdsBy(x = gencodev19_db, use.names = TRUE)

# filtering transcripts which have 5'UTR
cds_by_tx <- cds_by_tx[names(cds_by_tx) %in% names(five_utr_by_tx)]

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

# upstream sequences
seqs_five_utr_by_tx_with_tis <- extractTranscriptSeqs(x = BSgenome.Hsapiens.UCSC.hg19,
                                                      transcripts = five_utr_by_tx_with_tis)
tis_in_txcoords_by_tx <- split(tis_in_txcoords, seqnames(tis_in_txcoords))

# checking if all start codons are OK
x <- mcmapply(
  function(tx) {
  utr_seq <- seqs_five_utr_by_tx_with_tis[[seqnames(tx)[1]]]
  s <- 0
  for (i in 1:length(tx)) {
    s <- s + (toString(subseq(utr_seq, start = start(tx)[i], end = end(tx)[i])) == tx$name[i])
  }
  return(s > 0)
  },
  tis_in_txcoords_by_tx, mc.cores = cores)
sum(x) # OK

# searching for uorfs
tis_with_uorfs_in_txcoords_by_tx <- mclapply(
  tis_in_txcoords_by_tx,
  function(tx) {
    utr_seq <- seqs_five_utr_by_tx_with_tis[[seqnames(tx)[1]]]
    for (i in 1:length(tx)) {
      s <- start(tx)[i]
      sub_utr_seq <- toString(subseq(utr_seq, start = s))
      tx$uORF_end[i] <- stringr::str_locate(
        string = sub_utr_seq,
        pattern = '^([ACTG]{3})+?((TAA)|(TAG)|(TGA))'
      )[2] + s - 1
  }
  return(tx)
  },
  mc.cores = cores
)
tis_with_uorfs_in_txcoords_by_tx <- GRangesList(tis_with_uorfs_in_txcoords_by_tx)
tis_with_uorfs_in_txcoords <- unlist(tis_with_uorfs_in_txcoords_by_tx, use.names = FALSE)
names(tis_with_uorfs_in_txcoords) <- NULL
tis_with_uorfs_in_txcoords <- tis_with_uorfs_in_txcoords[!is.na(tis_with_uorfs_in_txcoords$uORF_end)]

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

# consensus kozak context
context_ranges <- mclapply(
  cds_by_tx,
  function(cds_ranges) {
    first_exon <- cds_ranges[1]
    context_range <- reduce(c(flank(first_exon, 7), resize(first_exon, 4)))
    return(context_range)
  },
  mc.cores = cores
)
context_ranges <- GRangesList(context_ranges)
context_ranges <- unlist(context_ranges)
context_seqs <- getSeq(
  x = BSgenome.Hsapiens.UCSC.hg19,
  names = context_ranges
)
kozak_pwm <- PWM(context_seqs)
mcols(context_ranges)$name <- paste(names(context_ranges), context_seqs, sep=';')
score(context_ranges) <- sapply(context_seqs, PWMscoreStartingAt, pwm = kozak_pwm)
mean(score(context_ranges))
export.bed(object = context_ranges, '~/BioData/tis/main_codons.bed')

# scoring uorfs
score(uorfs) <- sapply(uorfs$name, PWMscoreStartingAt, pwm = kozak_pwm)
mean(score(uorfs))

# saving
mean(width(uorfs))
export.bed(object = uorfs, '~/BioData/tis/uorfs_in_tx_coords.bed')

rm(context_seqs, context_ranges)

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
