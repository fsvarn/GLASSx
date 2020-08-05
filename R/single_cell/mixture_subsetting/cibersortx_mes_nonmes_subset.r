###################################################
# Get initial samples that are mesenchymal and then recurrent samples that are mesenchymal
# Updated: 2020.06.17
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

##################################################
rm(list=ls())
# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
SELECT ss.*, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b, cs.idh_codel_subtype
FROM analysis.rna_silver_set ss
LEFT JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
LEFT JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.subtypes cs ON ss.case_barcode = cs.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

myinf1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"
tpm <- read.delim(myinf1, row.names = 1)
colnames(tpm) <- gsub("\\.","-",colnames(tpm))

# Remove genes with 0 expression across cohort
sums <- apply(tpm, 1, function(x) sum(x == 0))
tpm <- tpm[-which(sums == ncol(tpm)),]


# Initial
#-------------------------------

# Pull out initial samples that are mesenchymal
init_mes_samples <- dat[which(dat[,"subtype_a"] == "Mesenchymal"), "tumor_barcode_a"]
init_mes <- tpm[,init_mes_samples]

# Pull out initial samples that are not mesenchymal
init_non_samples <- dat[which(dat[,"subtype_a"] != "Mesenchymal"), "tumor_barcode_a"]
init_non <- tpm[,init_non_samples]

# Write tables
# Save file as a .txt file for input into CIBERSORTx
myoutf1 <- "/projects/verhaak-lab/GLASS-III/data/subset/glass_idhwt_init_mes_tpm.txt"
myoutf2 <- "/projects/verhaak-lab/GLASS-III/data/subset/glass_idhwt_init_non_tpm.txt"

write.table(init_mes, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(init_non, myoutf2, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)


# Recurrent
#-------------------------------

# Pull out recurrent samples that are mesenchymal
rec_mes_samples <- dat[which(dat[,"subtype_b"] == "Mesenchymal"), "tumor_barcode_b"]
rec_mes <- tpm[,rec_mes_samples]

# Pull out recurrent samples that are not mesenchymal
rec_non_samples <- dat[which(dat[,"subtype_b"] != "Mesenchymal"), "tumor_barcode_b"]
rec_non <- tpm[,rec_non_samples]

# Write tables
# Save file as a .txt file for input into CIBERSORTx
myoutf1 <- "/projects/verhaak-lab/GLASS-III/data/subset/glass_idhwt_rec_mes_tpm.txt"
myoutf2 <- "/projects/verhaak-lab/GLASS-III/data/subset/glass_idhwt_rec_non_tpm.txt"

write.table(rec_mes, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(rec_non, myoutf2, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

