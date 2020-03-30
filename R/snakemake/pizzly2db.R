#!/usr/bin/env Rscript

#######################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

fusions <- args[1]
#fusions <- "results/kallisto/pizzly/final/fusions_all_samples.tsv"

cat("Reading in inputs...")
tfus <- read.delim(fusions,stringsAsFactors=FALSE)
totalcount <- tfus[,"paircount"] + tfus[,"splitcount"]
tfus <- cbind(tfus,totalcount)
tfus <- tfus[,c(1:7,9,8)]
tfus <- tfus[order(tfus[,"totalcount"],decreasing=TRUE),]
colnames(tfus) <- c("aliquot_barcode","gene_symbol_a","ensembl_gene_id_a","gene_symbol_b","ensembl_gene_id_b","paircount","splitcount","totalcount","ensembl_transcript_list")

dbWriteTable(con, Id(schema="analysis", table="pizzly_fusions"), tfus, overwrite=TRUE)

