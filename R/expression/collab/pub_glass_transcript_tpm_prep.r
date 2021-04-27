###################################################
# Compare TPM for NUPR1a tpm between primary and recurrent samples
# Collaboration with H.K. Ng group
# Updated: 2020.09.23
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)
library(survival)

#######################################################
rm(list=ls())

#Establish connection
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/kallisto/strand_unadjust/kallisto/final/old/transcript_count_matrix_all_samples.tsv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/pub_glass/glass_transcript_counts.txt"

old_dat <- read.delim(myinf1,row.names=1)
pub_dat <- read.delim(myinf2,row.names=1)

old_dat <- old_dat[,colnames(pub_dat)]

myinf3 <- "/projects/verhaak-lab/GLASS-III/results/kallisto/strand_unadjust/kallisto/final/old/transcript_count_matrix_all_samples.tsv"
pub_tpm <- read.delim(myinf3,row.names=1)

pub_tpm  <-  pub_tpm[,colnames(pub_dat)]
pub_tpm <- pub_tpm %>% rownames_to_column(var = "target_id")

myoutf <- "/projects/verhaak-lab/GLASS-III/data/dataset/pub_glass/glass_transcript_tpm.txt"
write.table(pub_tpm, myoutf, sep="\t", row.names=FALSE, quote=FALSE)