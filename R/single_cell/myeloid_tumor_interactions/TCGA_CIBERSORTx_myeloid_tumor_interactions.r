###################################################
# Compare macrophage activation profiles in TCGA between transcriptional subtypes
# Updated: 2020.07.22
# Author: Frederick Varn
##################################################
library(tidyverse)
library(ssgsea.GBM.classification)
library(GSVA)
library(RColorBrewer)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.txt"

expr <- read.delim(myinf1,stringsAsFactors=FALSE)

unigene <- table(expr[,"GeneSymbol"])
expr1 <- expr[which(unigene==1),]

# Only one gene symbol, can throw out as its mostly low compared to its duplicate
expr2 <- expr[which(unigene>1),]

rownames(expr1) <- expr1[,1]
expr <- expr1[,-1]
expr <- data.matrix(expr)

activation_molecules <- c("IFNB1", "IL10", "IL4", "IL13", "IFNG", 
						  "TNF", "LTA", "LTB", "TNFSF4", "CD40LG", "FASLG", "CD70", "TNFSF8", "TNFSF9", "TNFSF10", "TNFSF11", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNFSF18", "EDA")

##################################################
# Step 1: Assess macrophage activation activity on CIBERSORTx macrophage profiles
##################################################

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/CIBERSORTxHiRes_TCGA_myeloid_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
#colnames(geps) <- gsub("\\.","-",colnames(geps))
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]

myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"
mac_module_df <- read.delim(myinf2,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- paste("module_",1:49,sep="")

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)

res <- t(gsva(data.matrix(geps), mac_modules, method="ssgsea",parallel.sz=1))

##################################################
# Step 1: Correlate macrophage activity with different activation molecules in BULK TCGA
##################################################

act_bulk <- t(expr[activation_molecules,])
cor(act_bulk, res, method="s")
