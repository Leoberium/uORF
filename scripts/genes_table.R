library(tidyverse)

# loading data
target_counts_by_eps <- read_tsv('~/ME/uORF/genes_vareps.tsv')
tx_counts <- read_tsv('~/ME/uORF/transcripts.tsv')
target_counts <- filter(target_counts_by_eps, eps == 0)
go_rna_binding <- read_tsv('~/BioData/geneset.grp', comment = '#')

# description
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)
head(keys(org.Hs.eg.db, ''))
genes_data <- select(
  x = org.Hs.eg.db,
  keys = target_counts$gene_name,
  columns = c('SYMBOL', 'GENENAME', 'ENTREZID'),
  keytype = 'SYMBOL'
)

# test
gene_name <- unique(target_counts$gene_name)[2]
sub_gene <- target_counts[target_counts$gene_name == gene_name, ]
sub_tx <- tx_counts[tx_counts$gene_name == gene_name, ]
sub_data = genes_data[genes_data$SYMBOL == gene_name, ]
`-`(sum(sub_tx[sub_tx$contains_uORF, ]$sh_kd) / sub_tx$sh_gene_kd,
    sum(sub_tx[sub_tx$contains_uORF, ]$sh_control) / sub_tx$sh_gene_control)

# preparing final table
tt <- map_dfr(
  unique(target_counts$gene_name),
  
  function(gene_name) {
    sub_gene <- target_counts[target_counts$gene_name == gene_name, ]
    sub_tx <- tx_counts[tx_counts$gene_name == gene_name, ]
    sub_data = genes_data[genes_data$SYMBOL == gene_name, ]

    x <- tibble(
      gene_id = sub_gene$gene_id,
      gene_name = gene_name,
      entrez_id = sub_data$ENTREZID,
      desc = sub_data$GENENAME,
      tx_names = paste(sub_tx$tx_name, collapse = ', '),
      is_rbp = gene_name %in% go_rna_binding$GO_RNA_BINDING,
      self_rbp = sub_gene$self_rbp,
      is_cancer_gene = sub_gene$cancer_gene,
      uorf_tx_number = sub_gene$uorftx_number,
      all_tx_number = sub_gene$alltx_number,
      max_uorf_kozak_score = max(sub_tx$score, na.rm = TRUE),
      max_main_kozak_score = max(sub_tx$main_score),
      diff_score = max_uorf_kozak_score - max_main_kozak_score,
      
      deltaC_smg6 = sub_gene$frac_smg6 - sub_gene$frac_control,
      deltaC_upf1 = sub_gene$frac_upf1 - sub_gene$frac_control,
      
      deltaX_smg6 = `-`( sum(sub_tx$tpm_smg6[sub_tx$contains_uORF]) / sum(sub_tx$tpm_smg6),
                         sum(sub_tx$tpm_xrn1[sub_tx$contains_uORF]) / sum(sub_tx$tpm_xrn1) ),
      deltaX_upf1 = `-`( sum(sub_tx$tpm_upf1[sub_tx$contains_uORF]) / sum(sub_tx$tpm_upf1),
                         sum(sub_tx$tpm_xrn1[sub_tx$contains_uORF]) / sum(sub_tx$tpm_xrn1) ),
      delta_KD = `-`(sum(sub_tx[sub_tx$contains_uORF, ]$sh_kd) / mean(sub_tx$sh_gene_kd),
                     sum(sub_tx[sub_tx$contains_uORF, ]$sh_control) / mean(sub_tx$sh_gene_control)),
      
      nmd_gene_control_tpm = sub_gene$gene_tpm_control,
      nmd_gene_xrn1_tpm = sum(sub_tx$tpm_xrn1),
      nmd_gene_smg6_tpm = sub_gene$gene_tpm_smg6,
      nmd_gene_upf1_tpm = sub_gene$gene_tpm_upf1,
      
      nmd_uorf_tx_control_tpm = sub_gene$uorftx_tpm_control,
      nmd_uorf_tx_xrn1_tpm = sum(sub_tx$tpm_xrn1[sub_tx$contains_uORF]),
      nmd_uorf_tx_smg6_tpm = sub_gene$uorftx_tpm_smg6,
      nmd_uorf_tx_upf1_tpm = sub_gene$uorftx_tpm_upf1,
      
      nmd_frac_control = sub_gene$frac_control,
      nmd_frac_xrn1 = nmd_uorf_tx_xrn1_tpm / nmd_gene_xrn1_tpm,
      nmd_frac_smg6 = sub_gene$frac_smg6,
      nmd_frac_upf1 = sub_gene$frac_upf1,
      
      kd_gene_control_tpm = mean(sub_tx$sh_gene_control),
      kd_gene_knockdown_tpm = mean(sub_tx$sh_gene_kd),
      kd_uorf_tx_control_tpm = sum(sub_tx[sub_tx$contains_uORF, ]$sh_control),
      kd_uorf_tx_knockdown_tpm = sum(sub_tx[sub_tx$contains_uORF, ]$sh_kd),
      kd_frac_control =  kd_uorf_tx_control_tpm / kd_gene_control_tpm,
      kd_frac_knockdown = kd_uorf_tx_knockdown_tpm / kd_gene_knockdown_tpm
    )
  }
)
write_tsv(x = tt, path = '~/ME/uORF/gene_table.tsv')

# tests for delta XRN1
keepC_smg6 <- !is.na(tt$deltaC_smg6) & tt$deltaC_smg6 != 0
keepC_upf1 <- !is.na(tt$deltaC_upf1) & tt$deltaC_upf1 != 0
keepX_smg6 <- !is.na(tt$deltaX_smg6) & tt$deltaX_smg6 != 0
keepX_upf1 <- !is.na(tt$deltaX_upf1) & tt$deltaX_upf1 != 0
# keep_target_genes <- tt$uorf_tx_number == 1 & tt$all_tx_number > 1
keep_target_genes <- tt$gene_name %in% uorfP_genes

boxplot(tt$deltaC_smg6[keepC_smg6 & keep_target_genes])
wilcox.test(tt$deltaC_smg6[keepC_smg6 & keep_target_genes])

boxplot(tt$deltaC_upf1[keepC_upf1 & keep_target_genes])
wilcox.test(tt$deltaC_upf1[keepC_upf1 & keep_target_genes])

boxplot(tt$deltaX_smg6[keepX_smg6 & keep_target_genes])
wilcox.test(tt$deltaX_smg6[keepX_smg6 & keep_target_genes])

boxplot(tt$deltaX_upf1[keepX_upf1 & keep_target_genes])
wilcox.test(tt$deltaX_upf1[keepX_upf1 & keep_target_genes])

# additional uORF data
uorf_paper <- readxl::read_xlsx('ME/vkr/uorf_prediction/Supplemental_Data_Tables_.xlsx',
                                sheet = 6, skip = 2)
uorf_genes <- unique(uorf_paper$gene[uorf_paper$type == 'UTRonly'])

peptide_mapped <- readxl::read_xlsx('ME/vkr/uorf_prediction/Supplemental_Data_Tables_.xlsx',
                                    sheet = 4, skip = 2)
uorfP_genes <- unique(peptide_mapped$gene[peptide_mapped$UTRonly_vs_CDSpartial == 'UTRonly'])

# genes with both isoforms
both <- tt$uorf_tx_number != tt$all_tx_number
both_genes <- tt$gene_name[both]
rbp_genes <- tt$gene_name[tt$self_rbp]

filtered <- tx_counts[tx_counts$gene_name %in% both_genes & tx_counts$gene_name %in% rbp_genes, ]
