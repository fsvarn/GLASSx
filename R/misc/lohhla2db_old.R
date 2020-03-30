library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
myDir1 <- "results/lohhla/final/"
myDir1 <- paste(myDir1,dir(myDir1),sep="")
myinf1 <- unlist(sapply(myDir1,function(x)paste(x,dir(x),sep="/")))

myinf1 <- myinf1[grep("DNA.HLAlossPrediction_CI.txt",myinf1)]

tmp <- list()
for(i in 1:length(myinf1))
{
	tmp[[i]] <- read.delim(myinf1[i],header=TRUE,row.names=NULL)
}

full_lohhla <- do.call(rbind,tmp)
full_lohhla <- full_lohhla[,c("pair_barcode","coverage_filter","HLA_A_type1","HLA_A_type2","HLAtype1Log2MedianCoverage","HLAtype2Log2MedianCoverage",
	"HLAtype1Log2MedianCoverageAtSites","HLAtype2Log2MedianCoverageAtSites","HLA_type1copyNum_withBAFBin","HLA_type1copyNum_withBAFBin_lower","HLA_type1copyNum_withBAFBin_upper",
	"HLA_type2copyNum_withBAFBin","HLA_type2copyNum_withBAFBin_lower","HLA_type2copyNum_withBAFBin_upper","PVal_unique","UnPairedPval_unique",
	"LossAllele","KeptAllele","numMisMatchSitesCov","propSupportiveSites")]
colnames(full_lohhla) <- c("pair_barcode","coverage_filter","hla_type1","hla_type2","hla_type1_log2_median_cov","hla_type2_log2_median_cov",
	"hla_type1_log2_median_cov_at_sites","hla_type2_log2_median_cov_at_sites","hla_type1_copy_number","hla_type1_copy_number_lower","hla_type1_copy_number_upper",
	"hla_type2_copy_number","hla_type2_copy_number_lower","hla_type2_copy_number_upper","pval","unpaired_pval",
	"loss_allele","kept_allele","num_mismatch_sites_cov","prop_supportive_sites")

mycols <- c("hla_type1","hla_type2","loss_allele","kept_allele")
for(i in 1:length(mycols))
{
	full_lohhla[,mycols[i]] <- toupper(full_lohhla[,mycols[i]])
	full_lohhla[,mycols[i]] <- gsub("HLA_A_","HLA-A\\*",full_lohhla[,mycols[i]])
	full_lohhla[,mycols[i]] <- gsub("HLA_B_","HLA-B\\*",full_lohhla[,mycols[i]])
	full_lohhla[,mycols[i]] <- gsub("HLA_C_","HLA-C\\*",full_lohhla[,mycols[i]])
	full_lohhla[,mycols[i]] <- gsub("_",":",full_lohhla[,mycols[i]])

}

dbWriteTable(con,Id(scheme="variants",table="lohhla"), full_lohhla, overwrite=TRUE)