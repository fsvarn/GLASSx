###################################################
# Validate CIBERSORTx cell type-specific GEPs
# Updated: 2020.07.14
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/"
mytag <- dir(myDir1)
myinf1 <- paste(myDir1, mytag, sep ="")
myinf1 <- myinf1[grep("_Window48.txt", myinf1)]

myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
#myinf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

facs <- read.delim(myinf2, sep=",",row.names=1)

##################################################
# Helper functions
##################################################

# TPM calculation function
# Source: https://www.biostars.org/p/307603/
# Modified from: https://gist.github.com/slowkow/c6ab0348747f86e2748b
#-------------------------------------------

Counts_to_tpm <- function(counts, featureLength) {

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))

  # Compute effective lengths of features in each library.
  effLen <- featureLength

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

##################################################
# Step 1: Calculate TPM from count data 
##################################################
	
# Get gene lengths
geneLengths <- read.delim("/projects/verhaak-lab/verhaak_ref/star/GRCh38/GRCh38.d1.vd1/gencode.gene.info.v22.tsv",stringsAsFactor=FALSE)

# Clean up gene lengths so that there is one length per gene symbol. In cases where there are multiple isoforms, use the one with the fewest number of exons
featureLengths <- geneLengths %>%
	filter(gene_name %in% rownames(facs)) %>%
	filter(gene_type == "protein_coding") %>%
	group_by(gene_name) %>%
	slice_min(order_by = exon_num) %>%
	slice_min(order_by = exon_length) %>%
	slice_min(order_by = start) %>%
	select(gene_name, exon_length)
	
featureVector <- featureLengths$exon_length
names(featureVector) <- featureLengths$gene_name
featureVector <- featureVector[rownames(facs)]

tpm <- Counts_to_tpm(facs, featureVector)


glioma <- tpm[,grep("glioma_",colnames(tpm))]
norm <- tpm[,grep("nonTumor_",colnames(tpm))]

glioma_cd45n <- apply(glioma[,grep("_cd45n$", colnames(glioma))], 1, mean)
glioma_mdm <- apply(glioma[,grep("_mdm$", colnames(glioma))], 1, mean)
glioma_mg <- apply(glioma[,grep("_mg$", colnames(glioma))], 1, mean)
glioma_neutrophils <- apply(glioma[,grep("_neutrophils$", colnames(glioma))], 1, mean)
glioma_cd4 <- apply(glioma[,grep("_cd4$", colnames(glioma))], 1, mean)
glioma_cd8 <- apply(glioma[,grep("_cd8$", colnames(glioma))], 1, mean)
norm_cd45n <- apply(norm[,grep("_cd45n$", colnames(norm))], 1, mean)
norm_mg <- apply(norm[,grep("_mg$", colnames(norm))], 1, mean)
#norm_mdm <- apply(norm[,grep("_mdm$", colnames(norm))], 1, mean)						#does not exist
#norm_neutrophils <- apply(norm[,grep("_neutrophils$", colnames(norm))], 1, mean)		#does not exist
#norm_cd4 <- apply(norm[,grep("_cd4$", colnames(norm))], 1, mean)						#does not exist
#norm_cd8 <- apply(norm[,grep("_cd8$", colnames(norm))], 1, mean)						#does not exist

groundtruth <- data.frame(glioma_cd45n, glioma_mdm, glioma_mg, glioma_neutrophils, glioma_cd4, glioma_cd8,
						  norm_cd45n, norm_mg)

#sig <- c("prolif_stemcell_tumor","stemcell_tumor","differentiated_tumor","oligodendrocyte", "myeloid", "myeloid", "dendritic_cell", "dendritic_cell", "t_cell", "t_cell", "granulocyte")
sig <- gsub(myDir1, "", myinf1)
sig <- gsub("CIBERSORTxHiRes_GLASS_", "", sig)
sig <- gsub("_Window48.txt", "", sig)

res <- matrix(0, nrow = length(myinf1), ncol = ncol(groundtruth))
rownames(res) <- sig
colnames(res) <- colnames(groundtruth)
for(i in 1:length(myinf1))
{
	cat("\r", i)
	csx <- read.delim(myinf1[i],row.names=1)
	
	#vars <- apply(csx,1,function(x)var(x,na.rm=TRUE))
	#csx <- csx[which(vars > 0),]
	
	comxx <- intersect(rownames(csx), rownames(groundtruth))
	csx <- csx[comxx,]
	test <- groundtruth[comxx,]
	
	cor_matrix <- cor(csx, test, method="p",use="pairwise.complete.obs")
	avg_cor <- apply(cor_matrix, 2, mean)
	res[i,] <- avg_cor
}

res <- res[c("prolif_stemcell_tumor","stemcell_tumor","differentiated_tumor","fibroblast", "pericyte", "endothelial", "oligodendrocyte", "myeloid", "dendritic_cell", "t_cell", "granulocyte", "b_cell"),]
sub_res <- res[c("prolif_stemcell_tumor","stemcell_tumor","differentiated_tumor","oligodendrocyte", "myeloid"),]
