library(estimate)
library(DBI)
library(odbc)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())

#Load and prepare gene expression matrix


filterCommonGenes("/projects/verhaak-lab/GLASS-III/data/estimate/ivygap/ivygap_fpkm_clean.txt", "/projects/verhaak-lab/GLASS-III/data/estimate/ivygap/ivygap_fpkm_clean_estimate_filtered.tsv",id="GeneSymbol")
estimateScore("/projects/verhaak-lab/GLASS-III/data/estimate/ivygap/ivygap_fpkm_clean_estimate_filtered.tsv","/projects/verhaak-lab/GLASS-III/data/estimate/ivygap/ivygap_fpkm_clean_estimate_scores.tsv","affymetrix")

estimate <- t(read.delim("/projects/verhaak-lab/GLASS-III/data/estimate/ivygap/ivygap_fpkm_clean_estimate_scores.tsv"))
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
estimate$aliquot_barcode <- gsub("^X","",estimate$aliquot_barcode)


myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/columns-samples.csv"

samps <- read.csv(myinf2,stringsAsFactor=FALSE)


#Simplify analysis for now: Reference histology
samp_struct <- samps$structure_name
names(samp_struct) <- samps$rna_well_id

tumor_name <- samps$tumor_name
names(tumor_name) <- samps$rna_well_id

estimate[,"structure"] <- samp_struct[as.character(estimate[,"aliquot_barcode"])]
estimate[,"tumor_name"] <- tumor_name[as.character(estimate[,"aliquot_barcode"])]
estimate <- estimate[grep("reference histology", estimate[,"structure"]),]

dat[,"structure"] <- gsub(" sampled by reference histology","",dat[,"structure"])
