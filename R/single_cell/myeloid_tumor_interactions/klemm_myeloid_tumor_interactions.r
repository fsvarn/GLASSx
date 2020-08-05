###################################################
# Compare macrophage activation profiles in Klemm between transcriptional subtypes
# RNAseq count data was variance stabilized using DESeq2 and then batch corrected to get final expression values
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.07.22
# Author: Frederick Varn
##################################################

library(tidyverse)
library(DESeq2)
library(ssgsea.GBM.classification)
library(GSVA)
library(limma)
library(RColorBrewer)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

cts <- read.csv(myinf1,row.names=1)
info <- read.csv(myinf2,row.names=1,stringsAsFactor=FALSE)

# Add a batch column to info
info[,"batch"] <- sapply(strsplit(rownames(info),"_"),function(x)x[1])

sort_name <- colnames(cts)
sample_name <- sapply(strsplit(sort_name, "_"),function(x)paste(x[1:3],collapse="_"))
sub_info <- info[sample_name,]
rownames(sub_info) <- NULL
coldata <- data.frame(sort_name, sample_name, sub_info)

# Create DESeq2 object and perform variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch)                     
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)

# Batch correction
expr <- limma::removeBatchEffect(mat, vsd$batch)

activation_molecules <- c("IFNB1", "IL10", "IL4", "IL13", "IFNG", 
						  "TNF", "LTA", "LTB", "TNFSF4", "CD40LG", "FASLG", "CD70", "TNFSF8", "TNFSF9", "TNFSF10", "TNFSF11", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNFSF18", "EDA")
		
##################################################
# Step 1: Get subtype of each CD45- cells from glioma
##################################################

glioma <- expr[,grep("glioma",colnames(expr))]
glioma <- glioma[,grep("cd45n",colnames(glioma))]

##################################################
# Step 2: Examine whether macrophages or microglia activation modules are enriched in different subtypes (ssGSEA)
##################################################

# Load macrophage modules and convert to signature lists:
module_inf <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"

mac_module_df <- read.delim(module_inf,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- gsub("X","module_",colnames(mac_module_df))

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)

myeloid_vst <- cbind(expr[,grep("glioma_\\d\\d\\d\\d_mdm",colnames(expr))], expr[,grep("glioma_\\d\\d\\d\\d_mg",colnames(expr))])

# Run GSVA with macrophage modules on myeloid cell vst:
res <- t(gsva(myeloid_vst, mac_modules, method="ssgsea",parallel.sz=1))

# MDM
mdm_res <- res[grep("_mdm", rownames(res)),]
act_glioma <- t(glioma[activation_molecules,])

samp1 <- sapply(strsplit(rownames(mdm_res),"_"),function(x)paste(x[1:3],collapse="_"))
samp2 <- sapply(strsplit(rownames(act_glioma),"_"),function(x)paste(x[1:3],collapse="_"))
comxx <- intersect(samp1, samp2)
idhwt <- rownames(info[which(info[,"IDH.Status"] == "wild type"),])
comxx <- intersect(comxx, idhwt)

sub_mdm <- paste(comxx,"_mdm",sep="")
sub_cd45n <- paste(comxx, "_cd45n",sep="")

dat1 <- mdm_res[sub_mdm,]
dat2 <- act_glioma[sub_cd45n,]

cor(dat1, dat2, method = "s")

