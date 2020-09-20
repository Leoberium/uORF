library(tidyverse)

# eCLIP data
hepg2 <- read_tsv('BioData/HepG2.tsv', col_names = FALSE)
k562 <- read_tsv('BioData/K562.tsv', col_names = FALSE)

hepg2 <- hepg2 %>% filter(X7 > 3 & X8 > 3)
k562 <- k562 %>% filter(X7 > 3 & X8 > 3)

genes_hepg2 <- str_extract(hepg2$X4, '^\\w+(?=_HepG2)')
genes_k562 <- str_extract(k562$X4, '^\\w+(?=_K562)')

hepg2$X4 <- str_extract(hepg2$X4, '^\\w+(?=_rep)')
k562$X4 <- str_extract(k562$X4, '^\\w+(?=_rep)')

peaks_merged <- bind_rows(hepg2, k562)
peaks_merged <- peaks_merged %>% select(X1, X2, X3, X4, X5, X6)
colnames(peaks_merged) <- c('chr', 'start', 'end', 'name', 'score', 'strand')
peaks_merged <- GenomicRanges::makeGRangesFromDataFrame(peaks_merged, keep.extra.columns = TRUE)
peaks_merged <- GenomicRanges::sort.GenomicRanges(peaks_merged)
peaks_merged$gene <- str_extract(peaks_merged$name, '^\\w+(?=_(HepG2|K562))')
peaks_merged_list <- GenomicRanges::GRangesList(split(peaks_merged, peaks_merged$gene))
eclip_genes <- unique(peaks_merged$gene)

gencode <- rtracklayer::import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
genes_of_interest <- gencode[gencode$type == 'gene' & gencode$gene_name %in% names(peaks_merged_list)]
x <- genes_of_interest$gene_name
GenomicRanges::mcols(genes_of_interest) <- NULL
genes_of_interest$name <- x

result <- lapply(names(peaks_merged_list), function(gene) {
  IRanges::subsetByOverlaps(peaks_merged_list[[gene]], genes_of_interest[genes_of_interest$name == gene])
})
names(result) <- names(peaks_merged_list)
result <- GenomicRanges::GRangesList(result)
result <- unlist(result)
rtracklayer::export.bed(result, '~/BioData/tracks/input/RBP_peaks.bed')

t_table <- read_tsv('ME/uORF/transcripts.tsv')
g_table <- read_tsv('ME/uORF/gene_table.tsv')
amb_genes <- g_table$gene_name[g_table$uorf_tx_number < g_table$all_tx_number]
ft_table <- t_table[t_table$gene_name %in% amb_genes, ]

filtered_transcripts <- map_dfr(
  ft_table$tx_name,
  function(tx_name) {
    sub_tx <- ft_table[ft_table$tx_name == tx_name, ]
    sub_gene <- g_table[g_table$gene_name == sub_tx$gene_name, ]
    
    res <- tibble(
      tx_name = tx_name,
      gene_id = sub_tx$gene_id,
      gene_name = sub_tx$gene_name,
      contains_uORF = sub_tx$contains_uORF,
      is_rbp_GO = sub_gene$is_rbp,
      is_rbp_eCLIP = sub_tx$gene_name %in% eclip_genes,
      binds_itself = sub_gene$self_rbp,
      is_cancer_gene = sub_gene$is_cancer_gene,
      deltaC_smg6 = sub_tx$delta_smg6,
      deltaC_upf1 = sub_tx$delta_upf1,
      deltaX_smg6 = `-`(sub_tx$tpm_smg6 / sub_gene$nmd_gene_smg6_tpm, sub_tx$tpm_xrn1 / sub_gene$nmd_gene_xrn1_tpm),
      deltaX_upf1 = `-`(sub_tx$tpm_upf1 / sub_gene$nmd_gene_upf1_tpm, sub_tx$tpm_xrn1 / sub_gene$nmd_gene_xrn1_tpm),
      deltaGene = sub_tx$delta_sh
    )
  }
)
x <- write_tsv(filtered_transcripts, 'ME/uORF/filtered_transcripts.tsv')
sum(filtered_transcripts$is_rbp_GO)
sum(filtered_transcripts$is_rbp_eCLIP)
sum(filtered_transcripts$binds_itself)
x <- filtered_transcripts[filtered_transcripts$is_rbp_eCLIP, ]

x <- `-`(g_table$nmd_uorf_tx_upf1_tpm / (g_table$nmd_gene_upf1_tpm - g_table$nmd_uorf_tx_upf1_tpm),
         g_table$nmd_uorf_tx_xrn1_tpm / (g_table$nmd_gene_xrn1_tpm - g_table$nmd_uorf_tx_xrn1_tpm))
wilcox.test(x)

x <- read_tsv('ME/uORF/filtered_transcripts.tsv')
y <- x %>% filter(binds_itself)

pept_uorfs <- readxl::read_xlsx('ME/vkr/uorf_prediction/Supplemental_Data_Tables_.xlsx', sheet = 4, skip = 2)
pept_uorfs <- pept_uorfs %>% filter(UTRonly_vs_CDSpartial == 'UTRonly')
target_genes <- unique(pept_uorfs$gene)
g_table <- read_tsv('ME/uORF/gene_table.tsv')
wilcox.test(g_table[g_table$gene_name %in% target_genes, ]$deltaC_smg6)
wilcox.test(g_table[g_table$gene_name %in% target_genes, ]$deltaC_upf1)
wilcox.test(g_table[g_table$gene_name %in% target_genes, ]$deltaX_smg6)
wilcox.test(g_table[g_table$gene_name %in% target_genes, ]$deltaX_upf1)
