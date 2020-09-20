# counts
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
tx_counts <- readr::read_tsv('vkr/concatenated.tsv')

# transcripts with uorfs by upred
uorf_tx <- uorf_pred$uORF_ID
uorf_tx <- c(
  Filter(function(x) x != '-', unlist(strsplit(uorf_pred$additional_transcripts, ','))),
  uorf_tx)
uorf_tx <- unique(substr(uorf_tx, 1, 17))
enst_to_ensg <- select(
  x = gencode_db, keys = uorf_tx,
  columns = c('GENEID', 'TXNAME'), keytype = 'TXNAME'
)
ensg_to_enst <- select(
  x = gencode_db, keys = enst_to_ensg$GENEID,
  columns = c('GENEID', 'TXNAME'), keytype = 'GENEID'
)
ensg_to_enst$uORF <- ensg_to_enst$TXNAME %in% uorf_tx

# combining
ensg_to_enst <- dplyr::left_join(ensg_to_enst, tx_counts, by = c('TXNAME' = 'tx_name'))
ensg_to_enst <- ensg_to_enst[!is.na(ensg_to_enst$control), ]
x <- as.data.frame(mcols(gencode)[c('gene_id', 'gene_name')])
x <- x[!duplicated(x), ]
ensg_to_enst <- dplyr::left_join(ensg_to_enst, x, by = c('GENEID' = 'gene_id'))

# RBP
self_peaks <- import.bed('vkr/self_peaks.bed')

# gene table
g_table <- purrr::map_dfr(
  unique(ensg_to_enst$GENEID),
  function(gene_id) {
    sub_df <- ensg_to_enst[ensg_to_enst$GENEID == gene_id, ]
    dplyr::tibble(
      gene_id = gene_id,
      gene_name = sub_df$gene_name[1],
      self_rbp = gene_name %in% self_peaks$name,
      # avg_score_diff = mean(sub_df$score[sub_df$uORF] - sub_df$main_score[sub_df$uORF]),
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
g_table
fg_table <- g_table[g_table$n_non_uorf != 0, ]

# plots
fg_table
lfg_table <- tidyr::pivot_longer(fg_table[fg_table$g_control > 0 & fg_table$g_xrn1 > 0, ],
                                 14:18, names_to = 'comparison', 
                                 values_to = 'delta', values_drop_na = TRUE)
lfg_table$comparison <- factor(lfg_table$comparison,
                               levels = c('d_smg6', 'dx_smg6', 'd_upf1',
                                          'dx_upf1', 'dc'), ordered = TRUE)
addline_format <- function(x, ...) {gsub('_', '\n', x)}
library(ggplot2)
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
  theme(axis.text.x = element_text(size = 12), axis.text.y.left = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.size = unit(1.5, 'cm'),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank()) +
  ylim(c(-0.25, 0.25))
lfg_table$comparison
wilcox.test(lfg_table$delta[lfg_table$comparison == 'd_smg6' & lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'dx_smg6' & lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'd_upf1' & lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'dx_upf1'& lfg_table$delta != 0])
wilcox.test(lfg_table$delta[lfg_table$comparison == 'dc' & lfg_table$delta != 0],
            alternative = 'less')

ggplot(fg_table[fg_table$self_rbp, ], aes(x = log10(g_xrn1_upf1 / g_control), y = d_upf1)) +
  geom_point()
fg_table[fg_table$self_rbp, c('gene_name', 'd_upf1', 'd_smg6')]