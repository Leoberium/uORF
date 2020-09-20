library(tidyverse)
cores <- detectCores()

# loading main codons
main_codons <- rtracklayer::import.bed('BioData/tis/main_codons.bed')
mean_kozak <- mean(main_codons$score)
sd_kozak <- sd(main_codons$score)

# loading uorfs
uorfs <- rtracklayer::import.bed('BioData/tis/uorfs_in_tx_coords.bed')

# loading counts
tx_counts <- read_tsv('BioData/transcriptomes/concatenated.tsv')
quant_targets <- read_tsv('BioData/shRNA/quant_targets.tsv')
quant_shRNA <- read_tsv('BioData/shRNA/total.tsv')

# loading rbp data
rbp_ranges <- rtracklayer::import.bed('~/BioData/tracks/RBP_peaks.sorted.bed')
self_rbp <- unique(rbp_ranges$name)

# loading annotation
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
rm(x)

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

# preparing other data
target_counts_by_eps <- read_tsv('~/ME/uORF/genes_vareps.tsv')
target_counts <- filter(target_counts_by_eps, eps == 0)
rm(target_counts_by_eps)
x <- tibble(
  tx_name = str_sub(main_codons$name, 1, 17), 
  main_score = main_codons$score
)
tx_counts <- tx_counts %>% left_join(x, by = 'tx_name')

# preparing frags
library(GenomicFeatures)
five_utr_by_tx <- fiveUTRsByTranscript(x = gencode_db, use.names = TRUE)
cds_by_tx <- cdsBy(x = gencode_db, use.names = TRUE)
cds_by_tx <- cds_by_tx[names(cds_by_tx) %in% names(five_utr_by_tx)]
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
sum(names(frag_ranges) %in% tx_counts$tx_name)
frag_ranges <- frag_ranges[names(frag_ranges) %in% tx_counts$tx_name, ]
uorfs <- uorfs[seqnames(uorfs) %in% names(frag_ranges), ]


tx_name <- tx_counts$tx_name[51]
gene_name <- tx_counts$gene_name[51]

nmd_df <- target_counts[target_counts$gene_name == gene_name, ]
nmd_df
quant_df <- quant_targets[quant_targets$gene_name == gene_name, ]
quant_df
sh_df <- quant_shRNA[quant_shRNA$tx_name == tx_name, ]
sh_df
sh_df[, quant_df$control1]
sh_df[, quant_df$control2]
sh_control <- as.numeric(sh_df[, quant_df$control1] + sh_df[, quant_df$control2])
sh_control

# preparing final table
row
z <- map2_dfr(
  tx_counts$tx_name, tx_counts$gene_name,
  
  function(tx_name, gene_name) {
    
    nmd_df <- target_counts[target_counts$gene_name == gene_name, ]
    quant_df <- quant_targets[quant_targets$gene_name == gene_name, ]
    if (nrow(quant_df) == 2) {
      quant_df <- quant_df[quant_df$cell_line == 'HepG2', ]
    }
    sh_df <- quant_shRNA[quant_shRNA$tx_name == tx_name, ]
    
    sh_control <- as.numeric(sh_df[, quant_df$control1] + sh_df[, quant_df$control2])
    sh_kd <- as.numeric(sh_df[, quant_df$knockdown1] + sh_df[, quant_df$knockdown2])
    
    x <- tibble(
      tx_name = tx_name,
      gene_name = gene_name,
      gene_tpm_control = nmd_df$gene_tpm_control,
      gene_tpm_smg6 = nmd_df$gene_tpm_smg6,
      gene_tpm_upf1 = nmd_df$gene_tpm_upf1,
      sh_control = ifelse(length(sh_control) > 0, sh_control, NA),
      sh_kd = ifelse(length(sh_kd) > 0, sh_kd, NA)
    )
    
    return(x)
  }
)
any(z$gene_tpm_smg6 == 0)
any(z$gene_tpm_upf1 == 0)
any(z$gene_tpm_control == 0)
zz <- z %>% 
  group_by(gene_name) %>% 
  mutate(
    sh_gene_control = sum(sh_control),
    sh_gene_kd = sum(sh_kd)
  ) %>% 
  ungroup()
tx_counts <- left_join(tx_counts, zz, by = 'tx_name')
tx_counts$gene_name.y <- NULL
colnames(tx_counts)[2] <- 'gene_name'

# subsetting rules
tx_counts <- tx_counts[!is.na(tx_counts$gene_tpm_control), ]
any(tx_counts$gene_tpm_control == 0)
keep_smg6 <- tx_counts$gene_tpm_smg6 != 0
keep_upf1 <- tx_counts$gene_tpm_upf1 != 0
keep_sh <- !is.na(tx_counts$sh_gene_control)
frag_ranges <- frag_ranges[names(frag_ranges) %in% tx_counts$tx_name, ]

# additional calculations
tx_counts$delta_smg6 <- `-`(tx_counts$tpm_smg6 / tx_counts$gene_tpm_smg6,
                            tx_counts$tpm_control / tx_counts$gene_tpm_control)
tx_counts$delta_upf1 <- `-`(tx_counts$tpm_upf1 / tx_counts$gene_tpm_upf1,
                            tx_counts$tpm_control / tx_counts$gene_tpm_control)
tx_counts$delta_sh <- `-`(tx_counts$sh_kd / tx_counts$sh_gene_kd,
                          tx_counts$sh_control / tx_counts$sh_gene_control)
write_tsv(tx_counts, '~/ME/uORF/transcruots.tsv')

# splitting
u_frags <- frag_ranges[names(frag_ranges) %in% tx_counts$tx_name[tx_counts$contains_uORF]]
nu_frags <- frag_ranges[names(frag_ranges) %in% tx_counts$tx_name[!tx_counts$contains_uORF]]

u_frags
mcols(u_frags)$thick <- ranges(mapFromTranscripts(x = uorfs, transcripts = frag_ranges))
mcols(nu_frags)$thick <- IRanges(start = min(start(nu_frags)),
                                       end = min(start(nu_frags)) - 1)

# SMG6 uORF
z <- mcmapply(
  function(name) {
    sub_counts <- tx_counts[keep_smg6 & tx_counts$tx_name == name, ]
    if (nrow(sub_counts) == 0) return(NA)
    else return(paste0(
      'dSMG6=', round(sub_counts$delta_smg6, 2), ';',
      'control=', round(sub_counts$tpm_control / sub_counts$gene_tpm_control, 2), ';',
      'dScore=', round(sub_counts$score - sub_counts$main_score, 2)
    ))
  },
  names(u_frags),
  mc.cores = cores
)
uu_frags <- u_frags
names(uu_frags) <- z
uu_frags <- uu_frags[!is.na(names(uu_frags))]
rtracklayer::export.bed(uu_frags, '~/BioData/tracks/input/uORF_deltaSMG6.bed')

# UPF1 uORF
z <- mcmapply(
  function(name) {
    sub_counts <- tx_counts[keep_upf1 & tx_counts$tx_name == name, ]
    if (nrow(sub_counts) == 0) return(NA)
    else return(paste0(
      'dUPF1=', round(sub_counts$delta_upf1, 2), ';',
      'control=', round(sub_counts$tpm_control / sub_counts$gene_tpm_control, 2), ';',
      'dScore=', round(sub_counts$score - sub_counts$main_score, 2)
    ))
  },
  names(u_frags),
  mc.cores = cores
)
uu_frags <- u_frags
names(uu_frags) <- z
uu_frags <- uu_frags[!is.na(names(uu_frags))]
rtracklayer::export.bed(uu_frags, '~/BioData/tracks/input/uORF_deltaUPF1.bed')

# shRNA uORF
z <- mcmapply(
  function(name) {
    sub_counts <- tx_counts[keep_sh & tx_counts$tx_name == name, ]
    if (nrow(sub_counts) == 0) return(NA)
    else return(paste0(
      'dKD=', round(sub_counts$delta_sh, 2), ';',
      'control=', round(sub_counts$sh_control / sub_counts$sh_gene_control, 2), ';',
      'FC=', round(sub_counts$sh_gene_kd / sub_counts$sh_gene_control, 2), ';',
      'dScore=', round(sub_counts$score - sub_counts$main_score, 2)
    ))
  },
  names(u_frags),
  mc.cores = cores
)
uu_frags <- u_frags
names(uu_frags) <- z
uu_frags <- uu_frags[!is.na(names(uu_frags))]
rtracklayer::export.bed(uu_frags, '~/BioData/tracks/input/uORF_deltaKD.bed')

# SMG6 nuORF
z <- mcmapply(
  function(name) {
    sub_counts <- tx_counts[keep_smg6 & tx_counts$tx_name == name, ]
    if (nrow(sub_counts) == 0) return(NA)
    else return(paste0(
      'dSMG6=', round(sub_counts$delta_smg6, 2), ';',
      'control=', round(sub_counts$tpm_control / sub_counts$gene_tpm_control, 2)
    ))
  },
  names(nu_frags),
  mc.cores = cores
)
nnu_frags <- nu_frags
names(nnu_frags) <- z
nnu_frags <- nnu_frags[!is.na(names(nnu_frags))]
rtracklayer::export.bed(nnu_frags, '~/BioData/tracks/input/non_uORF_deltaSMG6.bed')

# UPF1 nuORF
z <- mcmapply(
  function(name) {
    sub_counts <- tx_counts[keep_upf1 & tx_counts$tx_name == name, ]
    if (nrow(sub_counts) == 0) return(NA)
    else return(paste0(
      'dUPF1=', round(sub_counts$delta_upf1, 2), ';',
      'control=', round(sub_counts$tpm_control / sub_counts$gene_tpm_control, 2)
    ))
  },
  names(nu_frags),
  mc.cores = cores
)
nnu_frags <- nu_frags
names(nnu_frags) <- z
nnu_frags <- nnu_frags[!is.na(names(nnu_frags))]
rtracklayer::export.bed(nnu_frags, '~/BioData/tracks/input/non_uORF_deltaUPF1.bed')

# shRNA nuORF
z <- mcmapply(
  function(name) {
    sub_counts <- tx_counts[keep_sh & tx_counts$tx_name == name, ]
    if (nrow(sub_counts) == 0) return(NA)
    else return(paste0(
      'dKD=', round(sub_counts$delta_sh, 2), ';',
      'control=', round(sub_counts$sh_control / sub_counts$sh_gene_control, 2), ';',
      'FC=', round(sub_counts$sh_gene_kd / sub_counts$sh_gene_control, 2)
    ))
  },
  names(nu_frags),
  mc.cores = cores
)
nnu_frags <- nu_frags
names(nnu_frags) <- z
nnu_frags <- nnu_frags[!is.na(names(nnu_frags))]
rtracklayer::export.bed(nnu_frags, '~/BioData/tracks/input/non_uORF_deltaKD.bed')

# RBP
rtracklayer::export.bed(rbp_ranges, '~/BioData/tracks/input/RBP_peaks.bed')

target_counts
