library(tidyverse)

x1 <- read_tsv(file = 'BioData/transcriptomes/quants/control_quant/quant.sf')
x2 <- read_tsv(file = 'BioData/transcriptomes/quants/xrn1_quant/quant.sf')
x3 <- read_tsv(file = 'BioData/transcriptomes/quants/xrn1_smg6_quant/quant.sf')
x4 <- read_tsv(file = 'BioData/transcriptomes/quants/xrn1_upf1_quant/quant.sf')

fragments_tpm <- tibble(
  tx_name = x1$Name,
  tpm_control = x1$TPM,
  tpm_xrn1 = x2$TPM,
  tpm_smg6 = x3$TPM,
  tpm_upf1 = x4$TPM
)
rm(list = paste0('x', 1:4))

write_tsv(x = fragments_tpm, '~/BioData/transcriptomes/concatenated.tsv')
fragments_tpm
