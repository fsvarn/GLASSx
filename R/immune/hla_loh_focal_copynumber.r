library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT hl.*, loss_allele, pval, cs.idh_codel_subtype
FROM analysis.gatk_cnv_by_hla hl
JOIN biospecimen.aliquots bs ON bs.aliquot_barcode = hl.aliquot_barcode
JOIN biospecimen.samples sa On sa.sample_barcode = bs.sample_barcode
JOIN clinical.subtypes cs ON cs.case_barcode = sa.case_barcode
JOIN analysis.pairs pa ON pa.tumor_barcode = hl.aliquot_barcode
JOIN variants.lohhla lo ON lo.pair_barcode = pa.pair_barcode AND SUBSTRING(lo.loss_allele,1,5) = hl.gene_symbol
WHERE coverage_filter = 20
ORDER BY 3, 1
"
dat <- dbGetQuery(con,q)

nrow(dat[which(dat[,"hlvl_call"]==-2 & dat[,"pval"] < 0.05),])/nrow(dat[which(dat[,"hlvl_call"]==-2),])		#1
nrow(dat[which(dat[,"hlvl_call"]==-1 & dat[,"pval"] < 0.05),])/nrow(dat[which(dat[,"hlvl_call"]==-1),])		#0.54
nrow(dat[which(dat[,"hlvl_call"] < 0 & dat[,"pval"] < 0.05),])/nrow(dat[which(dat[,"hlvl_call"] < 0),])		#0.57


nrow(dat[which(dat[,"hlvl_call"] >= 0 & dat[,"pval"] < 0.05),])/nrow(dat[which(dat[,"hlvl_call"] >= 0),])	#0.16
nrow(dat[which(dat[,"hlvl_call"] >= 1 & dat[,"pval"] < 0.05),])/nrow(dat[which(dat[,"hlvl_call"] >= 1),])	#0.38
nrow(dat[which(dat[,"hlvl_call"] >= 2 & dat[,"pval"] < 0.05),])/nrow(dat[which(dat[,"hlvl_call"] >= 2),])	#0.29
