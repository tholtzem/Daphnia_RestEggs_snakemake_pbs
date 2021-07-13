detach("package:stats", unload=TRUE)
library(tidyverse)
#getwd()
#setwd("/mnt/data/snakemake_mach2/Daphnia_RestEggs_snakemake_pbs/")

basedir="/home/uibk/c7701178/scratch/DAPHNIA/Daphnia_RestEggs_snakemake_pbs/"

#args <- commandArgs(trailingOnly=TRUE)

pruned_position <- read_lines(paste0(basedir, "ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id"))
#pruned_position <- read_lines(args[1])
pruned_position <- sub('.*:', '', pruned_position) %>%
  as.integer()

pruned_snp_list <- read_tsv(paste0(basedir, "angsd/angsd_LC_GL2_cutoff_maf0018_nInd55.mafs.gz")) %>%
  dplyr::select(1:4) %>%
  filter(pruned_snp_list$position %in% pruned_position)

write_tsv(pruned_snp_list, "list/LDpruned_snps.list", col_names = F)
