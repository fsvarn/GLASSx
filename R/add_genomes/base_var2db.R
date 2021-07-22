library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "oligo")

dat <- read.delim("/projects/verhaak-lab/yeg/snakemake/GLASS-oligo/results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.tsv", colClasses = c("character","character",NULL,"character",rep(NULL,14)),header=FALSE)
dat[which(dat[,1]=='X'),1] <- 23
dat[,1] <- as.numeric(dat[,1])

colnames(dat) <- c("chrom", "pos", "ref","alt","gene_symbol","variant_classification","secondary_variant_classification","variant_type","genome_change","transcript","transcript_strand","transcript_exon","transcript_position","cdna_change","cds_change","protein_change","gc_content","reference_context")

## script modified because consensus.normalized.sorted.funcotated.tsv will serve as the base variants.anno table
dbWriteTable(con, Id(schema="variants", table="anno"), dat)

EOF
