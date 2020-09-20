x1 <- readr::read_tsv('vkr/quants/control/quant.sf')
x2 <- readr::read_tsv('vkr/quants/xrn1//quant.sf')
x3 <- readr::read_tsv('vkr/quants/xrn1_smg6/quant.sf')
x4 <- readr::read_tsv('vkr/quants/xrn1_upf1//quant.sf')
tx_counts <- data.frame(
  tx_name = x1$Name,
  control = x1$TPM,
  xrn1 = x2$TPM,
  xrn1_smg6 = x3$TPM,
  xrn1_upf1 = x4$TPM
)
rm(x1, x2, x3, x4)
readr::write_tsv(tx_counts, 'vkr/concatenated.tsv')

# importing
library(tidyverse)
tx_counts <- read_tsv('vkr/concatenated.tsv')
gencode <- rtracklayer::import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- GenomicFeatures::makeTxDbFromGRanges(gencode)
df <- data.frame(gencode)
df <- df[, c('gene_id', 'gene_name')]
df <- df[!duplicated(df), ]

# getting gene ids
ensg_to_enst <- AnnotationDbi::select(
  x = gencode_db, keys = tx_counts$tx_name,
  columns = c('GENEID', 'TXNAME'), keytype = 'TXNAME'
)

# loading uorfs
uorfs <- rtracklayer::import.bed('vkr/ut_tisdb_uorfs.bed')
tx_counts$uORF <- tx_counts$tx_name %in% as.character(uorfs@seqnames)
d_uorfs <- data.frame(uorfs)

# score from uORF
tx_counts <- left_join(tx_counts, d_uorfs[, c(1, 7)], by = c('tx_name' = 'seqnames'))

# score from main ORF
tx_data <- rtracklayer::import.bed('BioData/tis/main_codons.bed')
tx_data <- data.frame(tx_data)
tx_data$tx_name <- substr(tx_data$name, 1, 17)
colnames(tx_data)[7] <- 'main_score'
tx_counts <- left_join(tx_counts, tx_data[, 7:8], by = 'tx_name')

# adding gene id
tx_counts <- left_join(ensg_to_enst, tx_counts, by = c('TXNAME' = 'tx_name'))
tx_counts <- left_join(tx_counts, df, by = c('GENEID' = 'gene_id'))

# self-rbp
self_rbp <- rtracklayer::import('vkr/rbp/self_peaks.bed')


# counting for genes
target_genes <- unique(tx_counts$GENEID[tx_counts$TXNAME %in% as.character(uorfs@seqnames)])
g_table <- map_dfr(
  target_genes,
  function(gene_id) {
    sub_df <- tx_counts[tx_counts$GENEID == gene_id, ]
    tibble(
      gene_id = gene_id,
      gene_name = sub_df$gene_name[1],
      self_rbp = gene_name %in% self_rbp$name,
      avg_score_diff = mean(sub_df$score[sub_df$uORF] - sub_df$main_score[sub_df$uORF]),
      n_uorf = sum(sub_df$uORF),
      n_non_uorf = nrow(sub_df) - n_uorf,
      g_control = sum(sub_df$control),
      g_xrn1 = sum(sub_df$xrn1),
      g_xrn1_smg6 = sum(sub_df$xrn1_smg6),
      g_xrn1_upf1 = sum(sub_df$xrn1_upf1),
      tu_control = sum(sub_df$control[sub_df$uORF]),
      tu_xrn1 = sum(sub_df$xrn1[sub_df$uORF]),
      tu_xrn1_smg6 = sum(sub_df$xrn1_smg6[sub_df$uORF]),
      tu_xrn1_upf1 = sum(sub_df$xrn1_upf1[sub_df$uORF]),
      d_smg6 = (tu_xrn1_smg6 / g_xrn1_smg6) - (tu_control / g_control),
      d_upf1 = (tu_xrn1_upf1 / g_xrn1_upf1) - (tu_control / g_control),
      dx_smg6 = (tu_xrn1_smg6 / g_xrn1_smg6) - (tu_xrn1 / g_xrn1),
      dx_upf1 = (tu_xrn1_upf1 / g_xrn1_upf1) - (tu_xrn1 / g_xrn1),
      dc = (tu_xrn1 / g_xrn1) - (tu_control / g_control)
    )
  }
)

# filtering
sum(g_table$n_non_uorf == 0)
fg_table <- g_table[g_table$n_non_uorf != 0, ]

# plots
lfg_table <- pivot_longer(fg_table, 15:19, names_to = 'comparison',
                          values_to = 'delta', values_drop_na = TRUE)

lfg_table$comparison <- factor(lfg_table$comparison,
                               levels = c('d_smg6', 'dx_smg6', 'd_upf1',
                                          'dx_upf1', 'dc'), ordered = TRUE)
addline_format <- function(x, ...) {gsub('_', '\n', x)}
ggplot(lfg_table[lfg_table$delta != 0, ], aes(x = comparison, y = delta)) +
  geom_boxplot(alpha = 1, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 4.5, color = 'seagreen', size = 1) +
  xlab('Comparison') + ylab(expression(paste(Delta, ' expression'))) +
  ggtitle('') +
  scale_x_discrete(labels = addline_format(c(
    'XRN1 & SMG6_vs_Control',
    'XRN1 & SMG6_vs_XRN1',
    'XRN1 & UPF1_vs_Control',
    'XRN1 & UPF1_vs_XRN1',
    'XRN1_vs_Control'
  ))) +
  scale_fill_brewer(palette = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y.left = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank()) +
  ylim(c(-0.5, 0.5))
lfg_table$comparison
wilcox.test(lfg_table$delta[lfg_table$comparison == 'd_smg6' & lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'dx_smg6' & lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'd_upf1' & lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'dx_upf1'& lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'dc' & lfg_table$delta != 0])

up_genes <- intersect(fg_table$gene_name[fg_table$d_upf1 > 0],
          fg_table$gene_name[fg_table$d_smg6 > 0])
down_genes <- intersect(fg_table$gene_name[fg_table$d_upf1 < 0],
                        fg_table$gene_name[fg_table$d_smg6 < 0])

cat(up_genes, file = 'up_genes.txt', sep = '\n')
cat(down_genes, file = 'down_genes.txt', sep = '\n')

library(ggrepel)
ggplot(fg_table[fg_table$self_rbp, ], aes(x = log10(g_xrn1_upf1 / g_control), y = d_upf1)) +
  geom_point()

fg_table$gene_name[fg_table$d_smg6 > 0 & fg_table$d_upf1 > 0]
cat(
  fg_table$gene_name[fg_table$d_smg6 > 0 & fg_table$d_upf1 > 0],
  file = 'gene_list.txt',
  sep = '\n'
)
addline_format <- function(x, ...) {gsub('\\s', '\n', x)}
enr <- read_tsv('msigdb_results')
colnames(enr) <- c('gene_set_name', 'genes_in_gene_set', 'description',
                   'genes_in_overlap', 'k/K', 'p_value', 'q_value')
ggplot(data = enr, mapping = aes(x = genes_in_overlap,
                                 y = reorder(addline_format(gene_set_name), -q_value),
                                 fill = q_value)) +
  geom_col() + scale_fill_continuous(type = 'gradient', name = 'Q-value') + 
  xlab('Genes in Overlap') + ylab('GO category') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y.left = element_text(size = 13),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.25, 'cm'),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))

fg_table$is_rbp <- fg_table$gene_name %in% merged$name

ggplot(fg_table[fg_table$is_rbp, ],
       aes(x = log10(g_xrn1_upf1 / g_control), y = dx_upf1,
           color = self_rbp, label = gene_name)) +
  geom_point(size = 1) + ggtitle('XRN1 & UPF1 vs Control') +
  geom_label_repel(
    data = fg_table[fg_table$is_rbp, ],
    mapping = aes(x = log10(g_xrn1_upf1 / g_control),
                  y = d_upf1,
                  label = gene_name),
    size = 3
    )

gd <- read_tsv('ME/uORF/gene_table.tsv')
gd
keep_upf1 <- (gd$nmd_gene_upf1_tpm > 0) & (gd$deltaC_upf1 != 0)
keep_smg6 <- (gd$nmd_gene_smg6_tpm > 0) & (gd$deltaC_smg6 != 0)
library(ggrepel)
ggplot(
  data = gd[keep_upf1 & gd$is_rbp, ],
  mapping = aes(
    x = log10(nmd_gene_upf1_tpm / nmd_gene_control_tpm),
    y = deltaC_upf1,
    color = self_rbp,
    label = gene_name
  )
) + 
  geom_point(size = 2) + ggtitle('UPF1 & XRN1 vs Control') +
  geom_label_repel(
    data = gd[keep_upf1 & gd$is_rbp & gd$deltaC_upf1 > 0.1, ],
    mapping = aes(
      x = log10(nmd_gene_upf1_tpm / nmd_gene_control_tpm),
      y = deltaC_upf1,
      label = gene_name
    ),
    size = 7
  ) +
  geom_hline(yintercept = 0.1, size = 0.25) +
  xlab('log2(FC)') + 
  ylab(expression(paste(Delta, ' expression'))) +
  scale_color_hue(name = 'Binds to itself', breaks = c(TRUE, FALSE),
                       labels = c('YES', 'NO'), l = 50, c = 75) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y.left = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(1.5, 'cm'),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))
fg_table$gene_name[fg_table$self_rbp]                        
cat(
  fg_table$gene_name[fg_table$self_rbp], sep = ', '
)
