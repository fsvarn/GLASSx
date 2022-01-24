###################################################
# Select mesenchymal only cells for CIBERSORTx
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)


#######################################################
rm(list=ls())

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Load previous PRADA fusion run on GLASS
q <- "
SELECT rn.*, cs.idh_codel_subtype
FROM analysis.rna_silver_set rn
JOIN clinical.subtypes cs ON cs.case_barcode = rn.case_barcode
WHERE idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)
initial <- dat[,"tumor_barcode_a"]
recurrent <- dat[,"tumor_barcode_b"]

myinf1 <- "/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"
tpm <- read.delim(myinf1,stringsAsFactor=FALSE)
gene_symbol <- tpm[,"Gene_symbol"]

colnames(tpm) <- gsub("\\.","-",colnames(tpm))
init_tpm <- tpm[,initial]
rec_tpm <- tpm[,recurrent]

init_tpm <- data.frame(gene_symbol, init_tpm)
rec_tpm <- data.frame(gene_symbol, rec_tpm)

myoutf1 <- "/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_init_idhwt_samples.tsv"
myoutf2 <- "/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_rec_idhwt_samples.tsv"
write.table(init_tpm, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(rec_tpm, myoutf2, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
