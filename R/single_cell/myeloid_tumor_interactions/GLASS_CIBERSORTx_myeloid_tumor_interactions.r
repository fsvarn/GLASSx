###################################################
# Immune cell interactions
# Updated: 2020.07.28
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(ggdendro)
library(grid)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS/CIBERSORTxHiRes_GLASS_differentiated_tumor_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]

activation_molecules <- c("IFNB1", "IL10", "IL4", "IL13", "IFNG", 
						  "TNF", "LTA", "LTB", "TNFSF4", "CD40LG", "FASLG", "CD70", "TNFSF8", "TNFSF9", "TNFSF10", "TNFSF11", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNFSF18", "EDA")
						  
# Macrophage modules
myinf2 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS/CIBERSORTxHiRes_GLASS_myeloid_Window48.txt"

macs <- read.delim(myinf2, row.names=1)
colnames(macs) <- gsub("\\.","-",colnames(macs))
rem <- apply(macs,1,function(x)sum(is.na(x)))
macs <- macs[-which(rem==ncol(macs)),]

myinf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"
mac_module_df <- read.delim(myinf3,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- paste("module_",1:49,sep="")

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)

res <- t(gsva(data.matrix(macs), mac_modules, method="ssgsea",parallel.sz=1))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT ps.*, ts1.signature_name AS signature_name_a, ts2.signature_name AS signature_name_b
FROM analysis.rna_silver_set ps
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
"
dat <- dbGetQuery(con,q)

act_gep <- t(geps[activation_molecules,c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"])])
act_res <- res[c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"]),]

cor(act_gep, act_res)

tmp <- as.numeric(res[,"module_15"])

sub_res <- res[,]