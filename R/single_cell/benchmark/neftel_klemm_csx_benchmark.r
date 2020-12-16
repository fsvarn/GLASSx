###################################################
# Create ground truth file in CPM format for benchmarking with CIBERSORTx with Neftel profiles
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
AClike <- OPClike <- Mesenchymal <- NPClike <- apply(tumor_sub,1,mean)

# Myeloid average profile
myeloid_sub <- glioma[,grep("_mg|_mdm", colnames(glioma))]
Macrophage <- apply(myeloid_sub,1,mean)

# Oligodendrocyte, fibroblast, endothelial, pericyte average profile
oligodendrocyte_sub <- normal[,grep("_cd45n", colnames(normal))]
Oligodendrocyte <- apply(oligodendrocyte_sub,1,mean)

# T cell average profile
t_cell_sub <- glioma[,grep("_cd4|_cd8", colnames(glioma))]
T <- apply(t_cell_sub,1,mean)

ground_truth <- data.frame(AClike, OPClike, Mesenchymal, NPClike, Macrophage, Oligodendrocyte, T)

##################################################
# Step 1: Calculate CPM
##################################################

cpm <- apply(ground_truth,2, function(x) (x/sum(x))*1000000) 

# Save file as a .txt file for input into CIBERSORTx
myoutf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/neftel_2019/klemm_purified_avg_ground_truth_cpm_neftel.txt"

write.table(cpm, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
