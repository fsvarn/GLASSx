###################################################
# Create a reference matrix from the Venteicher IDHmut data that CIBERSORTx can use to make a signature matrix
# Updated: 2020.12.11
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/varnf/data/venteicher_idh-a/IDH_A_processed_data_portal.txt"
myinf2 <- "/projects/verhaak-lab/varnf/data/tirosh_idh-o/OG_processed_data_portal.txt"
myinf3 <- "/projects/verhaak-lab/varnf/data/venteicher_idh-a/IDH_A_cell_type_assignment_portal_v2.txt"
myinf4 <- "/projects/verhaak-lab/varnf/data/tirosh_idh-o/cell_type_assignment_portal.txt"
myinf5 <- "/projects/verhaak-lab/varnf/data/venteicher_idh-a/suva_idhmut_cell_states.txt"

idha_mat <- read.delim(myinf1,row.names=1)
idho_mat <- read.delim(myinf2,row.names=1)
idha_class <- read.delim(myinf3,skip=1,stringsAsFactor=FALSE)
idho_class <- read.delim(myinf4,skip=1,stringsAsFactor=FALSE)
malig_class <- read.delim(myinf5,stringsAsFactor=FALSE)

idha_class[,"TYPE"] <-trimws(idha_class[,"TYPE"]) 
idho_class[,"TYPE"] <-trimws(idho_class[,"TYPE"]) 
malig_class[,"cell_name"] <-trimws(malig_class[,"cell_name"]) 

# Remove duplicate cell state
idha_class <- idha_class[-which(idha_class[,"TYPE"] == "MGH45_P4_F05"),]

# Fix column names
colnames(idho_mat) <- gsub("^X","",colnames(idho_mat))
colnames(idha_mat) <- gsub("^X","",colnames(idha_mat))

# Prepare the IDH-O matrix for combining with IDH-A
tmp_idho <- idho_mat[rownames(idha_mat),]

# Combine the matrices
comb_mat <- cbind(idha_mat, tmp_idho)

# Combine the classifications
tmp_idha_class <- idha_class[,1:3]
comb_class <- rbind(tmp_idha_class, idho_class)		# Several IDH-A cells are missing from the cell type assignment (n = 439 cells)

comb_class <- comb_class %>%
mutate(group = recode(group, "Microglia/Macrophage" = "microglia_macrophage", "microglia/macrophage" = "microglia_macrophage", "Oligodendrocytes (non-malignant)" = "oligodendrocytes"))

# Build final classification table
nonmalig <- comb_class[which(comb_class[,"group"] != "malignant"),]
nonmalig_class <- data.frame(NA, NA, NA, nonmalig$group, nonmalig$TYPE)
colnames(nonmalig_class) <- c("scgs_oligo", "scgs_astro", "scgs_stem","cell_state","cell_name")

final_class <- rbind(malig_class, nonmalig_class)

malig_vect <- malig_class$cell_state
names(malig_vect) <- malig_class$cell_name

comb_class[,"final_class"] <- comb_class$group
comb_class[which(comb_class$final_class == "malignant"),"final_class"] <- malig_vect[comb_class[which(comb_class$final_class == "malignant"),"TYPE"]]

final_class$cell_name[which(!final_class$cell_name %in% colnames(comb_mat))]

final_class <- final_class %>% filter(cell_state!="0")
final_class <- final_class[which(final_class$cell_name %in% colnames(comb_mat)),]

# Build final matrix
final_mat <- comb_mat[,final_class[,"cell_name"]]

#nans <- apply(final_mat,2,function(x)(sum(is.na(x))))
colnames(final_mat) <- final_class[,"cell_state"]

# Step 2: Randomly downsample to 5000 cells using seed 11
set.seed(11)
mysamps <- sample(1:ncol(final_mat),5000)
final_mat <- final_mat[,mysamps]

table(sapply(strsplit(colnames(final_mat),"\\."), "[[", 1))
#           Astro-like microglia_macrophage           Oligo-like 
#                  823                  668                  903 
#     oligodendrocytes              T cells     Undifferentiated 
#                   70                    5                 2531 


# Anti-log data to get it out of log2 space
final_mat_nolog <- (2^final_mat) - 1
final_mat_nolog[is.na(final_mat_nolog)] <- 0
final_mat_nolog <- round(final_mat_nolog,4)

# Add gene symbol as an independent column

gene <- rownames(final_mat_nolog)
final_mat_nolog <- data.frame(gene, final_mat_nolog)
rownames(final_mat_nolog) <- NULL

# Save file as a .txt file for input into CIBERSORTx
myoutf <- "/projects/verhaak-lab/GLASS-III/data/dataset/venteicher_2017/SS2_venteicher_cibersortx_ref_12112020.txt"
write.table(final_mat_nolog, myoutf, sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#md5sum: c8f94d8fac6579bfb27a4222d3d6870d


# Date: 2020-12-11 12:32:10
# Job type: Create Signature Matrix
# Single cell reference matrix file: SS2_venteicher_cibersortx_ref_12112020.txt
# Disable quantile normalization: true
# kappa: 999
# q-value: 0.01
# No. barcode genes: 300 to 500
# Min. Expression: 1
# Replicates: 3
# Sampling: 0
# Filter non-hematopoietic genes from signature matrix during construction: false


