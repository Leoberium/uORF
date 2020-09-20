#### Annotations ####

gencode <- rtracklayer::import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- GenomicFeatures::makeTxDbFromGRanges(gencode)
five_utr_by_tx <- GenomicFeatures::fiveUTRsByTranscript(x = gencode_db, use.names = TRUE)
cds_by_tx <- GenomicFeatures::cdsBy(x = gencode_db, use.names = TRUE)

setwd('Supp')
#### uORF track ####
library(tidyverse)

uorf_tools <- readxl::read_excel('journal.pone.0222459.s003.xlsx')
uorfs_ut <- tibble(
  chr = uorf_tools$chromosome,
  start = uorf_tools$start,
  stop = uorf_tools$stop,
  name = 'ut',
  score = 0,
  strand = uorf_tools$strand,
)
uorfs_ut <- GenomicRanges::makeGRangesFromDataFrame(uorfs_ut, keep.extra.columns = TRUE)
ch <- rtracklayer::import.chain('hg38ToHg19.over.chain')
uorfs_ut <- IRanges::reduce(rtracklayer::liftOver(uorfs_ut, ch))
names(uorfs_ut) <- uorf_tools$gene_symbol
uorfs_ut <- uorfs_ut[lengths(uorfs_ut) <= 1]
uorfs_ut <- unlist(uorfs_ut)
uorfs_ut$name <- paste0('ut_', names(uorfs_ut))
length(uorfs_ut)
uorfs_ut$score = 0
uorfs_ut # all uORFs from uORF tools study
any(duplicated(uorfs_ut)) # no duplicates

uorfs_ut

mean(sapply(getSeq(BSgenome.Hsapiens.UCSC.hg19, reduce(c(flank(uorfs_ut, 7), resize(uorfs_ut, 4)))),
       PWMscoreStartingAt, pwm = kozak_pwm))

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

ribotaper <- readxl::read_excel('41592_2016_BFnmeth3688_MOESM242_ESM.xlsx', sheet = 2, na = "NA")
uorf_ribotaper <- ribotaper[ribotaper$category == 'uORF', ]
rm(ribotaper)
str(uorf_ribotaper)
x <- unlist(str_split(uorf_ribotaper$ORF_id_gen, '_'))
uorfs_rt <- tibble(
  chr = x[seq(1, length(x), by = 3)],
  start = x[seq(2, length(x), by = 3)],
  stop = x[seq(3, length(x), by = 3)],
  name = paste0('rt_', uorf_ribotaper$gene_symbol),
  score = 0,
  strand = uorf_ribotaper$strand
)
length(unique(uorf_ribotaper$gene_symbol))
uorfs_rt <- GenomicRanges::makeGRangesFromDataFrame(uorfs_rt, keep.extra.columns = TRUE)
names(uorfs_rt) <- uorf_ribotaper$gene_symbol
uorfs_rt <- uorfs_rt[!duplicated(uorfs_rt)]
uorfs_ut['M6PR']
uorfs_rt['M6PR']

intersect(names(uorfs_ut), names(uorfs_rt))
uorfs_ut['NIPAL3']
uorfs_rt['NIPAL3']

uorfs_ut['DGCR8']
uorf_ribotaper$gene_symbol %in% uorf_tools$gene_symbol

#### test ####

x <- rtracklayer::import.bed('~/BioData/tis/uorfs_in_genomic_coords.bed')

#### other ####

library(tidyverse)

g_table <- read_tsv('uORF/gene_table.tsv')

keep <- g_table$uorf_tx_number < g_table$all_tx_number
sum(keep)

wilcox.test(g_table[keep, ]$nmd_gene_upf1_tpm - g_table[keep, ]$nmd_gene_control_tpm, alternative = "greater")
wilcox.test(g_table[keep, ]$nmd_gene_smg6_tpm - g_table[keep, ]$nmd_gene_control_tpm, alternative = "greater")
wilcox.test(g_table[keep, ]$nmd_gene_upf1_tpm - g_table[keep, ]$nmd_gene_xrn1_tpm, alternative = "greater")
wilcox.test(g_table[keep, ]$nmd_gene_smg6_tpm - g_table[keep, ]$nmd_gene_xrn1_tpm, alternative = "greater")

kg_table <- g_table[keep, ]
kg_table <- arrange(kg_table, desc(kg_table$nmd_gene_upf1_tpm - kg_table$nmd_gene_control_tpm))

#######################################################################################################################

library(Rsubread)
library(edgeR)
x <- featureCounts(files = list.files('~/BioData/nmd/filtered_alignments/', full.names = TRUE),
                   annot.ext = '~/BioData/anno/gencode.v19.filtered.annotation.gff3', 
                   isGTFAnnotationFile = TRUE, 
                   GTF.featureType = 'exon', 
                   GTF.attrType = 'gene_id',
                   GTF.attrType.extra = 'transcript_id',
                   useMetaFeatures = TRUE,
                   countMultiMappingReads = FALSE,
                   isPairedEnd = TRUE,
                   countChimericFragments = FALSE,
                   nthreads = 30)
x <- DGEList(counts = x$counts, genes = x$annotation,
             group = c('CONTROL', 'XRN1', 'XRN1_SMG6', 'XRN1_UPF1'), remove.zeros = TRUE)
keep <- rowSums(cpm(x) > 0) >= 2
x$counts <- cpm(x)
gene_cpm <- x$counts[keep, ]
gene_anno <- x$genes[keep, ]

#######################################################################################################################

library(tidyverse)

uorf_tools <- readxl::read_excel('journal.pone.0222459.s003.xlsx')
id_to_symbol <- rtracklayer::mcols(gencode)[, c('gene_id', 'gene_name')]
id_to_symbol <- id_to_symbol[!duplicated(id_to_symbol), ]
id_to_symbol <- data.frame(id_to_symbol)
gene_anno <- left_join(gene_anno, id_to_symbol, by = c('GeneID' = 'gene_id'))

gene_cpm
gene_anno$GeneID
  uorf_tools <- uorf_tools[substr(uorf_tools$gene_id, 1, 15) %in% substr(gene_anno$GeneID, 1, 15), ]
gene_cpm <- as.data.frame(gene_cpm)
gene_cpm$uORF <- substr(rownames(gene_cpm), 1, 15) %in% substr(uorf_tools$gene_id, 1, 15)

gene_cpm$gene_name <- gene_anno$gene_name
keep <- gene_cpm$control.bam > 0
gene_cpm <- gene_cpm[keep, ]
gene_cpm$fc_smg6 <- gene_cpm$xrn1smg6.bam / gene_cpm$control.bam
gene_cpm$fc_upf1 <- gene_cpm$xrn1upf1.bam / gene_cpm$control.bam

ggplot(gene_cpm, aes(y = log2(fc_smg6), fill = uORF)) +
  geom_boxplot() + ylim(c(-1, 1))
ggplot(gene_cpm, aes(y = log2(fc_upf1), fill = uORF)) +
  geom_boxplot() + ylim(c(-1, 1))

wilcox.test(gene_cpm$fc_smg6[gene_cpm$uORF], gene_cpm$fc_smg6[!gene_cpm$uORF], alternative = 'greater')
wilcox.test(gene_cpm$fc_upf1[gene_cpm$uORF], gene_cpm$fc_upf1[!gene_cpm$uORF], alternative = 'greater')

rbp_ranges <- rtracklayer::import.bed('~/BioData/tracks/input/RBP_peaks.bed')
self_rbp <- unique(str_extract(rbp_ranges$name, '\\w+(?=_)'))



uorf_prediction <- readxl::read_excel('Supplemental_Data_Tables_.xlsx', sheet = 4, skip = 2)


gene_cpm$self_rbp <- gene_cpm$gene_name %in% self_rbp



target_genes <- gene_cpm[(gene_cpm$uORF |
                          gene_cpm$gene_name %in% g_table$gene_name |
                          gene_cpm$gene_name %in% uorf_ribotaper$gene_symbol |
                          gene_cpm$gene_name %in% uorf_prediction[uorf_prediction$UTRonly_vs_CDSpartial == 'UTRonly', ]$gene) &
                          gene_cpm$self_rbp, ]
y1 <- head(arrange(gene_cpm[gene_cpm$uORF, ], desc(fc_upf1)), 50)
y2 <- head(arrange(gene_cpm[gene_cpm$uORF, ], desc(fc_smg6)), 50)
intersect(y1$gene_name, y2$gene_name)


uorf_tools[uorf_tools$gene_symbol == 'U2AF1', ]
uorf_prediction[uorf_prediction$gene == 'U2AF1', ]
uorf_ribotaper[uorf_ribotaper$gene_symbol == 'U2AF1', ]

uorf_tools[uorf_tools$gene_symbol == 'LRP5L', ]
