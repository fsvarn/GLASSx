###################################################
# Attempt to create synthetic mixtures from Neftel data to validate macrophage associations.... Did not work
# Reference for dataset: An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma: Cell
# Updated: 2020.07.21
# Author: Frederick Varn
##################################################

library(tidyverse)
library(ssgsea.GBM.classification)
library(GSVA)


rm(list=ls())
myinf1 <- "/projects/varnf/Data/Neftel_IDHwt_GBM/IDHwtGBM.processed.SS2.logTPM.txt"
myinf2 <- "/projects/varnf/Data/Neftel_IDHwt_GBM/IDHwt.GBM.Metadata.SS2.txt"

# Loads expression matrix where values = log2((tpm/10) + 1)
tpm <- read.delim(myinf1,row.names=1)
colnames(tpm) <- gsub("\\.","-",colnames(tpm))

# Remove malignant cells, only interested in macrophages
info <- read.delim(myinf2,stringsAsFactor=FALSE)
info <- info[-1,]
info[,8:15] <- apply(info[,8:15],2,as.numeric)
info <- info[-which(info[,"GBMType"]=="Pediatric"),]

malig <- info[which(info[,"CellAssignment"] == "Malignant"),]
myelo <- info[which(info[,"CellAssignment"] == "Macrophage"),]

# Anti-log data to get it out of log2 space
tpm <-  (2^tpm) - 1

malig_tpm <- tpm[,malig[,"NAME"]]
myelo_tpm <- tpm[,myelo[,"NAME"]]

mysamp <- unique(info[,"Sample"])

malig_mix <- myelo_mix <- matrix(0, nrow = nrow(malig_tpm), ncol = length(mysamp))
rownames(malig_mix) <- rownames(myelo_mix) <- rownames(malig_tpm)
colnames(malig_mix) <- colnames(myelo_mix) <- mysamp
for(i in 1:length(mysamp))
{
	cat("\r",i)
	sub_malig <- malig_tpm[,grep(mysamp[i],colnames(malig_tpm))]
	sub_myelo <- data.frame(myelo_tpm[,grep(mysamp[i],colnames(myelo_tpm))])	# Not all samples have myeloid cells so will force data frame format for now and remove later
	
	my_malig_mix <- apply(sub_malig, 1, sum)
	my_myelo_mix <- apply(sub_myelo, 1, sum)
	
	malig_mix[,i] <- my_malig_mix
	myelo_mix[,i] <- my_myelo_mix
}

# Remove samples without myeloid cells
myelo_sums <- apply(myelo_mix,2,sum)
malig_mix <- malig_mix[,-which(myelo_sums == 0)]
myelo_mix <- myelo_mix[,-which(myelo_sums == 0)]


##################################################
# Run ssGSEA on tumor and macrophage data
##################################################

##################################################
# Helper functions
##################################################


# GCT creator function (for transcriptional classifier)
#-------------------------------------------

gct <- function(data, gct_out)
{
	Description <- rep(NA,nrow(data))
	data2 <- cbind(rownames(data),Description,data[,1:ncol(data)])
	colnames(data2)[1] <- "NAME"
	write.table(data2, gct_out, sep="\t", quote=FALSE, row.names=FALSE)

	conIn <- file(gct_out, "r")
	rawfile = readLines(conIn)
	close(conIn)

	mytext <- c("#1.2", paste(nrow(data2),"\t",(ncol(data)-1),sep=""),rawfile)
	conOut = file(gct_out, "w")
	writeLines(mytext, conOut)
	close(conOut)
}	


# Step 1: Bulk transcriptional subtype the malignant cells

# Save glioma-specific TPM file in GCT format for classifier

tumor_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/neftel_2019/neftel_malignant_mixture.gct"
gct(malig_mix, tumor_path)

	
# Run Qianghu's transcriptional classifier
runSsGSEAwithPermutation(tumor_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/neftel_2019/p_result_neftel_malignant_mixture.gct.txt")

aliquot_barcode <- rep(rownames(subtype_ssgsea),3)
signature_name <- c(rep("Proneural",nrow(subtype_ssgsea)),rep("Classical",nrow(subtype_ssgsea)),rep("Mesenchymal",nrow(subtype_ssgsea)))
enrichment_score <- c(subtype_ssgsea[,"Proneural"],subtype_ssgsea[,"Classical"],subtype_ssgsea[,"Mesenchymal"])
p_value <- c(subtype_ssgsea[,"Proneural_pval"],subtype_ssgsea[,"Classical_pval"],subtype_ssgsea[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(aliquot_barcode,signature_name,enrichment_score,p_value,stringsAsFactors=FALSE)
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"aliquot_barcode"]),]


# Step 2: Assess macrophage modules in the macrophages

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

# Run GSVA with macrophage modules on myeloid cell tpm:
myelo_res <- t(gsva(myelo_mix, mac_modules, method="ssgsea",parallel.sz=1))

# Step 3: Identify up-regulated modules in each bulk subtype

mysubtypes <- as.character(unique(transcriptional_subtype[,"signature_name"]))

p_val <- matrix(0, nrow = ncol(myelo_res), ncol = length(mysubtypes))
rownames(p_val) <- rownames(myelo_res)
colnames(p_val) <- mysubtypes
for(i in 1:length(mysubtypes))
{
	g1 <- transcriptional_subtype %>%
			filter(signature_name == mysubtypes[i] & p_value < 0.05) %>%
			dplyr::select(aliquot_barcode) %>%
			pull()
	g2 <- transcriptional_subtype %>%
			filter(signature_name == mysubtypes[i] & p_value > 0.05) %>%
			dplyr::select(aliquot_barcode) %>%
			pull()
	p_val[,i] <- apply(myelo_res, 2, function(x) wilcox.test(x[g1], x[g2])$p.value)
}