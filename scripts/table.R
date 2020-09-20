library(tidyverse)

# loading data
tx_counts <- read_tsv('BioData/transcriptomes/concatenated.tsv')
uorfs <- rtracklayer::import.bed('BioData/tis/uorfs_in_tx_coords.bed')
gencode <- rtracklayer::import.gff3('BioData/anno/gencode.v19.filtered.annotation.gff3')
gencode_db <- GenomicFeatures::makeTxDbFromGRanges(gencode)

# filtering tx counts
tx_counts <- tx_counts[rowSums(tx_counts[2:5]) > 0, ]

# determining target genes
x <- AnnotationDbi::select(
  x = gencode_db, keys = unique(as.character(uorfs@seqnames)),
  columns = c('GENEID', 'TXNAME'), keytype = 'TXNAME'
)
target_genes <- AnnotationDbi::select(
  x = gencode_db, keys = unique(x$GENEID),
  columns = c('GENEID', 'TXNAME'), keytype = 'GENEID'
)
target_genes <- target_genes[target_genes$TXNAME %in% tx_counts$tx_name, ]
colnames(target_genes) <- c('gene_id', 'tx_name')

gencode_genes <- gencode[gencode$type == 'gene' & gencode$gene_id %in% target_genes$gene_id]

target_genes <- target_genes %>% 
  left_join(
    as_tibble(S4Vectors::mcols(gencode_genes)[, c('gene_id', 'gene_name')]),
    by = 'gene_id'
)

# subsetting counts
tx_counts <- tx_counts[tx_counts$tx_name %in% target_genes$tx_name, ]
tx_counts$contains_uORF <- tx_counts$tx_name %in% as.character(uorfs@seqnames)

# updating tx_counts
target_genes <- target_genes[, c(1, 3, 2)]
tx_counts <- target_genes %>% left_join(tx_counts, by = 'tx_name')
rm(target_genes)
tx_counts <- as_tibble(tx_counts)

# removing some duplicate gene
duplicated_genes <- tx_counts %>% 
  group_by(gene_name) %>% 
  summarise(gene_ids = n_distinct(gene_id)) %>% 
  filter(gene_ids > 1)
tx_counts[tx_counts$gene_name == duplicated_genes$gene_name, ]
# this ENSG00000270011.2 is not of interest cause it has only uORF transcripts
tx_counts <- tx_counts[tx_counts$gene_id != 'ENSG00000270011.2', ]

# scoring
uorfs <- GenomicRanges::makeGRangesFromDataFrame(
  df = data.frame(uorfs) %>% 
    group_by(seqnames) %>% 
    top_n(n = 1, wt = score) %>% ungroup(),
  keep.extra.columns = TRUE
)
x <- data.frame(uorfs)[, c(1, 7)]
colnames(x) <- c('tx_name', 'score')
tx_counts <- tx_counts %>% left_join(x, by = 'tx_name')
rm(x)

# network of cancer genes data
known_cancer_genes <- read_lines('~/BioData/lists/known_cancer_genes.txt')

# shRNA data
quant_targets <- read_tsv('BioData/shRNA/quant_targets.tsv')
quant_shRNA <- read_tsv('BioData/shRNA/total.tsv')

# RPB data
rbp_ranges <- rtracklayer::import.bed('BioData/tracks/RBP_peaks.sorted.bed')
self_rbp <- unique(rbp_ranges$name)

# variable eps
gene_ids <- unique(tx_counts$gene_id)
eps <- c(0, 10 ** (1:6) / 10**6)

target_counts_by_eps <- map_dfr(eps, function(eps) {
  df <- tx_counts
  df[4:7] <- df[4:7] + eps
  
  y <- map_dfr(gene_ids, function(gene_id) {
    sub_df <- df[df$gene_id == gene_id, ]
    x <- tibble(
      gene_id = gene_id,
      gene_name = tx_counts$gene_name[tx_counts$gene_id ==  gene_id][1],
      cancer_gene = gene_name %in% known_cancer_genes,
      self_rbp = gene_name %in% self_rbp,
      shRNA_KD = gene_name %in% quant_targets$gene_name,
      uorftx_number = sum(sub_df$contains_uORF),
      alltx_number = nrow(sub_df),
      eps = eps,
      avg_kozak_score = mean(sub_df$score, na.rm = TRUE),
      max_kozak_score = max(sub_df$score, na.rm = TRUE),
      uorftx_tpm_control = sum(sub_df$tpm_control[sub_df$contains_uORF]),
      uorftx_tpm_smg6 = sum(sub_df$tpm_smg6[sub_df$contains_uORF]),
      uorftx_tpm_upf1 = sum(sub_df$tpm_upf1[sub_df$contains_uORF]),
      gene_tpm_control = sum(sub_df$tpm_control),
      gene_tpm_smg6 = sum(sub_df$tpm_smg6),
      gene_tpm_upf1 = sum(sub_df$tpm_upf1),
      frac_control = uorftx_tpm_control / gene_tpm_control,
      frac_smg6 = uorftx_tpm_smg6 / gene_tpm_smg6,
      frac_upf1 = uorftx_tpm_upf1 / gene_tpm_upf1,
    )
    return(x)
  })
  
  return(y)
})

target_counts_by_eps <- target_counts_by_eps[target_counts_by_eps$gene_tpm_control != 0, ]
write_tsv(x = target_counts_by_eps, path = '~/ME/uORF/genes_vareps.tsv')

# eps = 0

target_counts_by_eps <- read_tsv('~/ME/uORF/genes_vareps.tsv')
target_counts <- filter(target_counts_by_eps, eps == 0)

keep_upf1 <- (target_counts$gene_tpm_upf1 > 0) & (target_counts$frac_upf1 - target_counts$frac_control != 0)
keep_smg6 <- (target_counts$gene_tpm_smg6 > 0) & (target_counts$frac_smg6 - target_counts$frac_control != 0)
sum(target_counts$self_rbp | target_counts$shRNA_KD)

target_counts[keep_upf1 &
                (target_counts$frac_upf1 - target_counts$frac_control < -0.9) &
                log10(target_counts$gene_tpm_upf1 / target_counts$gene_tpm_control) > 2, ]


ggplot(
  data = target_counts[keep_upf1, ],
  mapping = aes(
    x = log10(gene_tpm_upf1 / gene_tpm_control),
    y = frac_upf1 - frac_control,
    color = self_rbp
  )
) +
  geom_point(size = 0.75) + 
  geom_text(
    data = target_counts[keep_upf1 & (target_counts$frac_upf1 - target_counts$frac_control > 0.5), ],
    mapping = aes(
      x = log10(gene_tpm_upf1 / gene_tpm_control),
      y = frac_upf1 - frac_control,
      label = gene_name
    ),
    check_overlap = TRUE,
    size = 4,
    nudge_y = 0.03
  ) +
  xlab('logFC(Gene Expression)') + 
  ylab('Expression Delta') + 
  ggtitle('UPF1') +
  scale_color_discrete(name = 'Binds to itself', label = c('No', 'Yes')) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

peaks_merged <- rtracklayer::import.bed('~/BioData/eclip/peaks_merged.bed')
keep_rbp <- target_counts$shRNA_KD | target_counts$gene_name %in% unique(peaks_merged$name)

library(ggrepel)
ggplot(
  data = target_counts[keep_upf1 & keep_rbp, ],
  mapping = aes(
    x = log10(gene_tpm_upf1 / gene_tpm_control),
    y = frac_upf1 - frac_control,
    color = self_rbp
  )
) +
  geom_point(size = 0.75) + 
  geom_label_repel(
    data = target_counts[keep_upf1 & keep_rbp & (target_counts$frac_upf1 - target_counts$frac_control > 0.05), ],
    mapping = aes(
      x = log10(gene_tpm_upf1 / gene_tpm_control),
      y = frac_upf1 - frac_control,
      label = gene_name
    ),
    arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
    size = 4,
    nudge_y = 0.1
  ) +
  xlab('logFC(Gene Expression)') + 
  ylab('Expression Delta') + 
  ggtitle('UPF1 (RBP)') +
  scale_color_discrete(name = 'Binds to itself', label = c('No', 'Yes')) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

library(clusterProfiler)

qplot(
  x = avg_kozak_score,
  y = frac_upf1 - frac_control,
  data = target_counts[keep_upf1, ]
)
cor.test(
  x = target_counts$avg_kozak_score,
  y = target_counts$frac_upf1 - target_counts$frac_control,
  method = 'spearman'
)

positive_smg6 <- target_counts %>% filter(keep_smg6, target_counts$frac_smg6 - target_counts$frac_control > 0) %>% select(gene_name, gene_id)
positive_upf1 <- target_counts %>% filter(keep_upf1, target_counts$frac_upf1 - target_counts$frac_control > 0) %>% select(gene_name, gene_id)
top_genes <- intersect(positive_smg6$gene_name, positive_upf1$gene_name)

cat(top_genes, sep = '\n')

qplot(
  x = frac_upf1 - frac_control,
  data = filter(target_counts, keep_upf1, target_counts$frac_upf1 - target_counts$frac_control > 0)
)

library(org.Hs.eg.db)
eg = bitr(geneID = top_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

x <- enrichGO(
  gene = eg$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENTREZID'
)
eg[eg$ENTREZID %in% c('55651', '3183', '5073', '55135', '6628'), ]

barplot(x)

library(DOSE)
  y <- enrichDGN(
  gene = eg$ENTREZID
)
head(y)

y <- enrichNCG(
  gene = eg$ENTREZID
)
head(y)

kk <- enrichKEGG(
  gene = eg$ENTREZID,
  organism = 'hsa'
)
head(kk)
