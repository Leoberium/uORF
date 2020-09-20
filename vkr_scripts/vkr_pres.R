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
keep <- rowSums(cpm(x) > 0.1) >= 4
x$counts <- cpm(x)
gene_cpm <- data.frame(x$counts[keep, ])
gene_anno <- data.frame(x$genes[keep, ])

# saving
gene_cpm <- dplyr::tibble(
  gene_id = rownames(gene_cpm),
  control = gene_cpm$control.bam,
  xrn1 = gene_cpm$xrn1.bam,
  xrn1_smg6 = gene_cpm$xrn1smg6.bam,
  xrn1_upf1 = gene_cpm$xrn1upf1.bam
)
write_tsv(gene_cpm, 'vkr/gene_cpm.tsv')

# loading uorfs
library(rtracklayer)
tu <- import.bed('vkr/upred/upred_tx.bed')
gu <- import.bed('vkr/upred/upred_genomic.bed')
uorf_tx <- as.character(seqnames(tu))

# loading features
library(GenomicFeatures)
gencode <- import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- makeTxDbFromGRanges(gencode)

# genes of interest
enst_to_ensg <- AnnotationDbi::select(
  x = gencode_db, keys = unique(uorf_tx),
  columns = c('GENEID', 'TXNAME'), keytype = 'TXNAME'
)
length(unique(enst_to_ensg$GENEID)) # 198 unique genes with uORFs
ensg_to_enst <- AnnotationDbi::select(
  x = gencode_db, keys = unique(enst_to_ensg$GENEID),
  columns = c('GENEID', 'TXNAME'), keytype = 'GENEID'
)

# adding info on uORF
ensg_to_enst$uORF <- ensg_to_enst$TXNAME %in% uorf_tx
library(tidyverse)
df <- as.data.frame(tu)
df
df <- df %>% 
  group_by(seqnames) %>% 
  summarise(score = max(score))
uorf_genes <- left_join(ensg_to_enst, df, by = c('TXNAME' = 'seqnames'))
main_codons <- import.bed('BioData/tis/main_codons.bed')
x <- stringr::str_split(main_codons$name, ';', simplify = TRUE)[, 1]
df <- tibble(TXNAME = x, main_score = main_codons$score)
uorf_genes <- left_join(uorf_genes, df)

# counting number of uORF txs
uorf_genes <- uorf_genes %>% 
  group_by(GENEID) %>% 
  summarise(u_tx_number = sum(uORF),
            nonu_tx_number = n() - sum(uORF),
            max_kozak = max(score, na.rm = TRUE),
            max_kozak_diff = max(score - main_score, na.rm = TRUE))
uorf_genes$alluORF <- uorf_genes$nonu_tx_number == 0

# making table
allgenes <- gene_cpm$gene_id
gene_data <- tibble(
  gene_id = allgenes,
  control = gene_cpm$control,
  xrn1 = gene_cpm$xrn1,
  xrn1_smg6 = gene_cpm$xrn1_smg6,
  xrn1_upf1 = gene_cpm$xrn1_upf1,
  uORF = allgenes %in% uorf_genes$GENEID,
  alluORF = allgenes %in% uorf_genes$GENEID[uorf_genes$alluORF],
  fc_xrn1_control = gene_cpm$xrn1 / gene_cpm$control,
  fc_smg6_control = gene_cpm$xrn1_smg6 / gene_cpm$control,
  fc_upf1_control = gene_cpm$xrn1_upf1 / gene_cpm$control,
  fc_smg6_xrn1 = gene_cpm$xrn1_smg6 / gene_cpm$xrn1,
  fc_upf1_xrn1 = gene_cpm$xrn1_upf1 / gene_cpm$xrn1
)
gene_data <- left_join(gene_data, uorf_genes[, c(1, 4:5)], by = c('gene_id' = 'GENEID'))

# log-fold
lgc <- pivot_longer(gene_data, cols = 8:12,
                    names_to = 'comparison',
                    values_to = 'fold_change',
                    values_drop_na = TRUE)
lgc$lfc <- log2(lgc$fold_change)
sum(lgc$lfc == Inf | lgc$lfc == -Inf)
lgc <- lgc[!(lgc$lfc == Inf | lgc$lfc == -Inf), ]

# plots
library(ggplot2)
addline_format <- function(x, ...) {gsub('_', '\n', x)}
ggplot(lgc, aes(x = comparison, y = lfc, fill = uORF)) +
  geom_boxplot(alpha = 1, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 4.5, color = 'seagreen', size = 1) +
  xlab('Comparison') + ylab('log2(FC)') +
  ggtitle('uORFs validated by MS') +
  scale_x_discrete(labels = addline_format(c(
    'XRN1 & SMG6_vs_Control',
    'XRN1 & SMG6_vs_XRN1',
    'XRN1 & UPF1_vs_Control',
    'XRN1 & UPF1_vs_XRN1',
    'XRN1_vs_Control'
  ))) +
  scale_fill_brewer(name = 'Contains\nuORF', breaks = c(TRUE, FALSE),
                    labels = c('YES', 'NO'), palette = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y.left = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.size = unit(1.5, 'cm'), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  ylim(c(-1.5, 1.5))

# only UPF1 and SMG6
specific <- lgc[lgc$comparison == 'fc_upf1_control', ]
ggplot(specific[specific$control > 1, ], aes(x = comparison, y = lfc, fill = uORF)) +
  geom_boxplot() +
  ylim(c(-2, 2)) + ylab('log2(FC)') +
  scale_x_discrete(labels = c('UPF1 KD vs Control')) +
  ggtitle('Genes with uORFs are upregulated\nupon NMD inactivation in HEK293\ncells') +
  theme_bw() +
  scale_fill_manual(values = c('#c27ba0ff', '#93c47dff'),
                    name = 'Has uORF',
                    breaks = c(TRUE, FALSE),
                    labels = c('YES', 'NO')) +
  theme(axis.text.x = element_text(family = 'Helvetica', size = 24),
        axis.text.y.left = element_text(family = 'Helvetica', size = 20),
        legend.text = element_text(family = 'Helvetica', size = 20),
        legend.title = element_text(family = 'Helvetica', size = 20),
        legend.key.size = unit(1.5, 'cm'), 
        axis.title.y = element_text(family = 'Helvetica', size = 24),
        axis.title.x = element_blank(),
        plot.title = element_text(family = 'Helvetica', size = 24))

# tests
sapply(
  unique(lgc$comparison), function(comp) {
    x <- wilcox.test(
      x = lgc$lfc[lgc$comparison == comp & lgc$uORF],
      y = lgc$lfc[lgc$comparison == comp & !lgc$uORF],
      alternative = 'greater'
    )
    return(x$p.value)
  }
)
# fail

unique(lgc$comparison)
# let's try kozak
ulgc <- lgc[lgc$uORF & lgc$comparison != 'fc_xrn1_control', ]
ggplot(ulgc, aes(x = max_kozak, y = lfc)) +
  geom_point() + facet_wrap(~ comparison)
lapply(
  unique(ulgc$comparison), function(comp) {
    x <- cor.test(
      x = ulgc$max_kozak[ulgc$comparison == comp],
      y = ulgc$lfc[ulgc$comparison == comp],
      method = 'pearson'
    )
    return(x)
  }
)
unique(ulgc$comparison)
# kozak diff
ggplot(ulgc, aes(x = max_kozak, y = lfc)) +
  geom_point() + facet_wrap(~ comparison)
lapply(
  unique(ulgc$comparison), function(comp) {
    x <- cor.test(
      x = ulgc$max_kozak_diff[ulgc$comparison == comp],
      y = ulgc$lfc[ulgc$comparison == comp],
      method = 'pearson'
    )
    return(x)
  }
)
unique(ulgc$comparison)
comp.labels <- c(
  'fc_smg6_control' = 'XRN1 & SMG6 vs Control',
  'fc_upf1_control' = 'XRN1 & UPF1 vs Control'
)

ulgc2 <- ulgc[ulgc$comparison %in% c("fc_smg6_control", "fc_upf1_control"), ]
ulgc2$comparison <- factor(ulgc2$comparison,
                           levels = c('fc_smg6_control', 'fc_upf1_control'),
                           labels = c('XRN1 & SMG6 vs Control', 'XRN1 & UPF1 vs Control'))
ggplot(ulgc2, aes(x = max_kozak, y = lfc, colour = max_kozak_diff)) +
  geom_point(size = 2.5) +
  facet_wrap(~ comparison) +
  geom_smooth(method = 'lm') +
  theme_bw() + xlab('Max Kozak score') + ylab('log2(FC)') +
  scale_colour_continuous(type = 'viridis', name = 'Kozak\ndifference') +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y.left = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.75, 'cm'))
  
  

ggplot(ulgc2, aes(x = max_kozak > 0.8, y = lfc)) +
  geom_boxplot() + facet_wrap(~ comparison)

sapply(
  unique(ulgc2$comparison), function(comp) {
    x <- wilcox.test(
      x = ulgc2$lfc[ulgc2$max_kozak > 0.8 & ulgc2$comparison == comp],
      y = ulgc2$lfc[ulgc2$max_kozak <= 0.8 & ulgc2$comparison == comp],
      alternative = 'greater'
    )
    return(x$p.value)
  }
)



######################################################

# transcripts with uorfs by upred
upred <- readxl::read_excel('Supp/Supplemental_Data_Tables_.xlsx', sheet = 4, skip = 2)
uorf_tx <- upred$uORF_ID
uorf_tx <- c(
  Filter(function(x) x != '-', unlist(strsplit(upred$additional_transcripts, ','))),
  uorf_tx)
uorf_tx <- unique(str_sub(uorf_tx, 1, 17))
gt <- as.data.frame(gencode)[, c('gene_id', 'transcript_id', 'gene_name')]
gt <- gt[str_sub(gt$transcript_id, 1, 4) == 'ENST', ]
gt <- gt[!duplicated(gt), ]

gt$uORF <- str_sub(gt$transcript_id, 1, 15) %in% str_sub(uorf_tx, 1, 15)
gt <- gt %>% 
  group_by(gene_id) %>% 
  summarise(
    uORF = sum(uORF) > 0,
    alluORF = sum(uORF) == n(),
    gene_name = unique(gene_name)
  )
gd <- left_join(gene_cpm, gt)

# data and plots
gd$fc_xrn1_control <- gd$xrn1 / gd$control
gd$fc_smg6_control <- gd$xrn1_smg6 / gd$control
gd$fc_upf1_control <- gd$xrn1_upf1 / gd$control
gd$fc_smg6_xrn1 <- gd$xrn1_smg6 / gd$xrn1
gd$fc_upf1_xrn1 <- gd$xrn1_upf1 / gd$xrn1

lgc2 <- tidyr::pivot_longer(gd, cols = 9:13,
                           names_to = 'comparison', values_to = 'fold_change',
                           values_drop_na = TRUE)
lgc2$lfc <- log2(lgc2$fold_change)
sum(lgc2$lfc == Inf | lgc2$lfc == -Inf)
lgc2 <- lgc2[!(lgc2$lfc == Inf | lgc2$lfc == -Inf), ]

# plots
addline_format <- function(x, ...) {gsub('_', '\n', x)}
ggplot(lgc2, aes(x = comparison, y = lfc, fill = uORF)) +
  geom_boxplot(alpha = 1, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 4.5, color = 'seagreen', size = 1) +
  xlab('Comparison') + ylab('log2(FC)') +
  scale_x_discrete(labels = addline_format(c(
    'XRN1 & SMG6_vs_Control',
    'XRN1 & SMG6_vs_XRN1',
    'XRN1 & UPF1_vs_Control',
    'XRN1 & UPF1_vs_XRN1',
    'XRN1_vs_Control'
  ))) +
  scale_fill_brewer(name = 'Contains\nuORF', breaks = c(TRUE, FALSE),
                    labels = c('YES', 'NO'), palette = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y.left = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.8, 'cm'), 
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank()) +
  ylim(c(-2, 2))


res <- sapply(
  unique(lgc2$comparison), function(comp) {
    x <- wilcox.test(
      x = lgc2$lfc[lgc2$comparison == comp & lgc2$uORF],
      y = lgc2$lfc[lgc2$comparison == comp & !lgc2$uORF],
      alternative = 'greater'
    )
    return(x$p.value)
  }
)
round(5 * res, 6)


ggplot(lgc2[lgc2$uORF, ], aes(x = comparison, y = lfc, fill = alluORF)) +
  geom_boxplot(alpha = 1, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 4.5, color = 'seagreen', size = 1) +
  xlab('Comparison') + ylab('log2(FC)') +
  ggtitle('uORFs validated by MS') +
  scale_x_discrete(labels = addline_format(c(
    'XRN1 & SMG6_vs_Control',
    'XRN1 & SMG6_vs_XRN1',
    'XRN1 & UPF1_vs_Control',
    'XRN1 & UPF1_vs_XRN1',
    'XRN1_vs_Control'
  ))) +
  scale_fill_brewer(name = 'All\nuORF', breaks = c(TRUE, FALSE),
                    labels = c('YES', 'NO'), palette = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y.left = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.size = unit(1.5, 'cm'), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  ylim(c(-2, 2))
