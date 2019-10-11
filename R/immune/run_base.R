#!/usr/bin/env Rscript


library(DBI)
library(odbc)
library(reshape2)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT * FROM analysis.gene_tpm"
tpm <- dbGetQuery(con,q)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Please provide an input and a mapping file", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}

sig_matrix <- args[2]
output <- args[3]

#Convert tpm long format into expression matrix
tpm_matrix <- dcast(tpm, aliquot_barcode ~ gene_symbol, value.var=tpm)
tpm_matrix <- log10(tpm_matrix+1)

#Filter out genes that are not expressed in at least 3 samples
tag = apply(tpm_matrix>0,1,sum)
tpm_matrix = tpm_matrix[tag>=3,]

source("/R/immune/base.R")

res <- base(tpm_matrix, sig_matrix, perm=1000, median.norm=T))

write.table(res, myoutf, sep="\t", row.names=T, quote=F)



