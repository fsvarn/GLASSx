###################################################
# Create ground truth file in CPM format for benchmarking with CIBERSORTx
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.07.29
# Author: Frederick Varn
##################################################

library(tidyverse)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

counts <- read.csv(myinf1,row.names=1)
info <- read.csv(myinf2,row.names=1,stringsAsFactor=FALSE)

# Glioma samples:
glioma <- counts[,grep("_glioma_",colnames(counts))]

#Non-tumor samples (for normal cells):
normal <- counts[,grep("_nonTumor_",colnames(counts))]

# Tumor average profile
tumor_sub <- glioma[,grep("_cd45n", colnames(glioma))]
stemcell_tumor <- differentiated_tumor <- prolif_stemcell_tumor <- apply(tumor_sub,1,mean)

# Myeloid average profile
myeloid_sub <- glioma[,grep("_mg|_mdm", colnames(glioma))]
myeloid <- dendritic_cell <- apply(myeloid_sub,1,mean)

# Granulocyte profile
gran_sub <- glioma[,grep("_neutrophils", colnames(glioma))]
granulocyte <- apply(gran_sub,1,mean)

# Oligodendrocyte, fibroblast, endothelial, pericyte average profile
oligodendrocyte_sub <- normal[,grep("_cd45n", colnames(normal))]
oligodendrocyte <- fibroblast <- endothelial <- pericyte <- apply(oligodendrocyte_sub,1,mean)

# T cell and B cell average profile
t_cell_sub <- glioma[,grep("_cd4|_cd8", colnames(glioma))]
t_cell <- b_cell <- apply(t_cell_sub,1,mean)

ground_truth <- data.frame(stemcell_tumor, myeloid, differentiated_tumor, prolif_stemcell_tumor, oligodendrocyte, t_cell, granulocyte, pericyte, endothelial, dendritic_cell, fibroblast, b_cell)
ground_truth_subset <- data.frame(differentiated_tumor, myeloid, oligodendrocyte, t_cell, granulocyte, oligodendrocyte)
colnames(ground_truth_subset) <- c("tumor", "myeloid", "oligodendrocyte","lymphocyte","granulocyte", "stroma")
ground_truth_tumor_subset <- data.frame(differentiated_tumor, myeloid, oligodendrocyte, t_cell, granulocyte, pericyte, endothelial, dendritic_cell, fibroblast, b_cell)
colnames(ground_truth_tumor_subset) <- c("tumor", "myeloid", "oligodendrocyte","t_cell","granulocyte", "pericyte","endothelial", "dendritic_cell","fibroblast", "b_cell")

##################################################
# Step 1: Calculate CPM
##################################################

cpm <- apply(ground_truth,2, function(x) (x/sum(x))*1000000) 
cpm_subset <- apply(ground_truth_subset,2, function(x) (x/sum(x))*1000000) 
cpm_subset_tumor <- apply(ground_truth_tumor_subset,2, function(x) (x/sum(x))*1000000) 

# Save file as a .txt file for input into CIBERSORTx
myoutf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/klemm_purified_avg_ground_truth_cpm.txt"
myoutf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/klemm_purified_avg_ground_truth_cpm_subset.txt"
myoutf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/klemm_purified_avg_ground_truth_cpm_tumor_subset.txt"

write.table(cpm, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(cpm_subset, myoutf2, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(cpm_subset_tumor, myoutf3, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)