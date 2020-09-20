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
gene_cpm <- data.frame(x$counts[keep, ])
gene_anno <- data.frame(x$genes[keep, ])

# loading uorfs
library(rtracklayer)
rm(tisdb_uorfs_genomic, tisdb_uorfs, uorfs_ut, uorfs_ut_genomic, uorfs_ut_tx, five_utr_by_tx_with_tis)
tisdb_tx <- import.bed('vkr/tisdb_uorfs_tx.bed')
tisdb_gene <- import.bed('vkr/tisdb_uorfs_genomic.bed')
ut_tx <- import.bed('vkr/ut_uorfs_tx.bed')
ut_gene <- import.bed('vkr/ut_uorfs_genomic.bed')
uorf_prediction <- readxl::read_excel('Supp/Supplemental_Data_Tables_.xlsx', sheet = 4, skip = 2)
length(unique(uorf_prediction$gene))

mean(tisdb_gene$score)
mean(ut_gene$score)
table(substr(ut_gene$name, 8, 10))
table(substr(tisdb_gene$name, 8, 10))

ut_gene
mean(width(ut_tx))
mean(width(tisdb_tx))

# genes for each uORF set
library(GenomicFeatures)
gencode <- import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- makeTxDbFromGRanges(gencode)

tisdb_tx$gene_id <- select(
  x = gencode_db,
  keys = as.character(tisdb_tx@seqnames),
  columns = c('GENEID', 'TXNAME'), 
  keytype = 'TXNAME'
)$GENEID
length(unique(tisdb_tx$gene_id))

ut_tx$gene_id <- select(
  x = gencode_db,
  keys = as.character(ut_tx@seqnames),
  columns = c('GENEID', 'TXNAME'), 
  keytype = 'TXNAME'
)$GENEID
length(unique(ut_tx$gene_id))

# uORFs
gene_cpm$tisdb <- rownames(gene_cpm) %in% tisdb_tx$gene_id
gene_cpm$tisdb_atg <- rownames(gene_cpm) %in% tisdb_tx$gene_id[substr(tisdb_tx$name, 8, 10) == 'ATG']
gene_cpm$ut <- rownames(gene_cpm) %in% ut_tx$gene_id

id_to_symbol <- as.data.frame(mcols(gencode)[, c('gene_id', 'gene_name')])
id_to_symbol <- id_to_symbol[!duplicated(id_to_symbol), ]
upred <- id_to_symbol$gene_name %in% uorf_prediction$gene
id_to_symbol[upred, ]

library(VennDiagram)
library(ggplot2)
venn.diagram(
  x = list(id_to_symbol[upred, ]$gene_id, unique(ut_tx$gene_id), unique(tisdb_tx$gene_id)),
  category.names = c('uORFs validated by MS', 'uORF-Tools annotation', 'TISdb annotation'),
  filename = 'venn.png',
  output = TRUE,
  imagetype = 'png',
  height = 800,
  width = 800,
  resolution = 300,
  compression = 'lzw',
  lwd = 1,
  col=c("#256EFFFF", '#3DDC97FF', '#FF495CFF'),
  fill = c(alpha("#256EFFFF",0.3), alpha('#3DDC97FF',0.3), alpha('#FF495CFF',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#256EFFFF", '#3DDC97FF', '#FF495CFF'),
  rotation = 1
)

sum(duplicated(c(ut_gene, tisdb_gene)))
sum(gene_cpm$tisdb & gene_cpm$ut)


gene_cpm$upred <- rownames(gene_cpm) %in% id_to_symbol$gene_id[upred]
ic <- intersect(tisdb_tx$gene_id, ut_tx$gene_id)
gene_cpm$ic <- rownames(gene_cpm) %in% ic
gene_cpm$gene_id <- rownames(gene_cpm)
gene_cpm <- dplyr::left_join(gene_cpm, id_to_symbol, 'gene_id')
gene_cpm$fc_xrn1_control <- gene_cpm$xrn1.bam / gene_cpm$control.bam
gene_cpm$fc_smg6_control <- gene_cpm$xrn1smg6.bam / gene_cpm$control.bam
gene_cpm$fc_upf1_control <- gene_cpm$xrn1upf1.bam / gene_cpm$control.bam
gene_cpm$fc_smg6_xrn1 <- gene_cpm$xrn1smg6.bam / gene_cpm$xrn1.bam
gene_cpm$fc_upf1_xrn1 <- gene_cpm$xrn1upf1.bam / gene_cpm$xrn1.bam
gene_cpm
lgc <- tidyr::pivot_longer(gene_cpm, cols = 12:16,
                    names_to = 'comparison', values_to = 'fold_change',
                    values_drop_na = TRUE)
lgc$lfc <- log2(lgc$fold_change)
sum(lgc$lfc == Inf | lgc$lfc == -Inf)
lgc <- lgc[!(lgc$lfc == Inf | lgc$lfc == -Inf), ]


# Plots
library(ggplot2)
addline_format <- function(x, ...) {gsub('_', '\n', x)}

ggplot(lgc, aes(x = comparison, y = lfc, fill = upred)) +
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
                    labels = c('YES', 'NO'), palette = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.text.y.left = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.size = unit(1.5, 'cm'), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  ylim(c(-2, 2))


# Tests
unique(lgc$comparison)

comp <- 'fc_upf1_xrn1'

wilcox.test(
  x = lgc$lfc[lgc$comparison == comp & lgc$upred],
  y = lgc$lfc[lgc$comparison == comp & !lgc$upred],
  alternative = 'greater'
)

wilcox.test(
  x = lgc$lfc[lgc$comparison == comp & lgc$ut],
  y = lgc$lfc[lgc$comparison == comp & !lgc$ut],
  alternative = 'greater'
)

wilcox.test(
  x = lgc$lfc[lgc$comparison == comp & lgc$tisdb],
  y = lgc$lfc[lgc$comparison == comp & !lgc$tisdb],
  alternative = 'greater'
)

wilcox.test(
  x = lgc$lfc[lgc$comparison == comp & lgc$tisdb_atg],
  y = lgc$lfc[lgc$comparison == comp & !lgc$tisdb_atg],
  alternative = 'greater'
)

wilcox.test(
  x = lgc$lfc[lgc$comparison == comp & lgc$ic],
  y = lgc$lfc[lgc$comparison == comp & !lgc$ic],
  alternative = 'greater'
)

sum(gene_cpm$upred)
sum(gene_cpm$ut)
sum(gene_cpm$tisdb)
sum(gene_cpm$tisdb_atg)
sum(gene_cpm$ic)

# Enrichment tests
unique(lgc$comparison)
comp <- 'fc_xrn1_control'

x <- lgc[lgc$comparison == comp, ]

fisher.test(
  matrix(c(
    sum(x$upred & x$lfc > 0),
    sum(x$upred & !(x$lfc > 0)),
    sum(!(x$upred) & x$lfc > 0),
    sum(!(x$upred) & !(x$lfc > 0))),
    nrow = 2
  )
)

fisher.test(
  matrix(c(
    sum(x$ut & x$lfc > 0),
    sum(x$ut & !(x$lfc > 0)),
    sum(!(x$ut) & x$lfc > 0),
    sum(!(x$ut) & !(x$lfc > 0))),
    nrow = 2
  )
)

fisher.test(
  matrix(c(
    sum(x$tisdb & x$lfc > 0),
    sum(x$tisdb & !(x$lfc > 0)),
    sum(!(x$tisdb) & x$lfc > 0),
    sum(!(x$tisdb) & !(x$lfc > 0))),
    nrow = 2
  )
)

fisher.test(
  matrix(c(
    sum(x$tisdb_atg & x$lfc > 0),
    sum(x$tisdb_atg & !(x$lfc > 0)),
    sum(!(x$tisdb_atg) & x$lfc > 0),
    sum(!(x$tisdb_atg) & !(x$lfc > 0))),
    nrow = 2
  )
)

fisher.test(
  matrix(c(
    sum(x$ic & x$lfc > 0),
    sum(x$ic & !(x$lfc > 0)),
    sum(!(x$ic) & x$lfc > 0),
    sum(!(x$ic) & !(x$lfc > 0))),
    nrow = 2
  )
)

