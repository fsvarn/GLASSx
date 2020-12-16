library(estimate)
library(DBI)
library(odbc)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())

#Load and prepare gene expression matrix

filterCommonGenes("results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv", "data/sigs/gene_tpm_matrix_all_samples_estimate_filtered.tsv",id="GeneSymbol")
estimateScore("data/sigs/gene_tpm_matrix_all_samples_estimate_filtered.tsv","data/sigs/gene_tpm_matrix_all_samples_estimate_scores.txt","affymetrix")

estimate <- t(read.delim("data/sigs/gene_tpm_matrix_all_samples_estimate_scores.txt"))
colnames(estimate) <- estimate[2,]
estimate <- estimate[-c(1,2),]
aliquot_barcode <- gsub("\\.","-",estimate[,2])
estimate <- estimate[,3:5]
estimate <- apply(estimate,2,as.numeric)
purity <- cos(0.6049872018+0.0001467884 * estimate[,"ESTIMATEScore"])
estimate <- data.frame(aliquot_barcode,estimate,purity)
rownames(estimate) <- NULL
estimate <- estimate[order(estimate[,"aliquot_barcode"]),]
colnames(estimate) <- c("aliquot_barcode","stromal_score","immune_score","estimate_score","purity")

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Write ESTIMATE results to db
dbWriteTable(con, Id(schema="analysis",table="estimate"), estimate, overwrite = TRUE)
