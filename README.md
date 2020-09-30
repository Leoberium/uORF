# Regulatory Role of Upstream Open Reading Frames in RNA Degradation by Nonsense-Mediated Decay

### Further ideas:
- Prepare new quantification of transcripts with StringTie or something similar
- Prior filtering of genes of with 1 exon, genes with long 3'UTRs, etc. Generally all groups of genes which may be confounding factors

### Folders:
- `diagnostic_plots_files` and `genes_plots_files` are auxiliry folders with images for .Rmd files
- `scripts` - contains scripts with all preliminary data analysis
- `vkr_scripts` - scripts which were used for the preparation of thesis document
- `UCSCHub` - contains track and description files for UCSC Genome Browser Trackhub

### Files:
- `diagnostic_plots` - boxplots and scatterplots of gene expression with different pseudocounts
- `genes_plots` - scatterplots of gene expression
- `filtered_transcripts.tsv` & `transcripts_tsv` & `gene_table.tsv` - tables of expression (including delta) of transcripts and genes with some metadata
- `genes_vareps.tsv` - table of gene expression and delta with different pseudocount
- `msigdb_result.tsv` - enrichment of upregulated genes
- `positive_genes.txt` - list of upregulated genes
