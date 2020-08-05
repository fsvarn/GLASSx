###################################################
# Create purified mixture in CPM format for testing with CIBERSORTx
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.06.12
# Author: Frederick Varn
##################################################

library(tidyverse)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

#Establish connection to db
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
ref <- dbReadTable(con, Id(schema="ref",table="genes"))

counts <- read.csv(myinf1,row.names=1)
info <- read.csv(myinf2,row.names=1,stringsAsFactor=FALSE)

#We want glioma samples only:

counts <- counts[,grep("glioma",colnames(counts))]


##################################################
# Step 1: Calculate CPM
##################################################

cpm <- apply(counts,2, function(x) (x/sum(x))*1000000) 

# Save file as a .txt file for input into CIBERSORTx
myoutf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/klemm_purified_cell_mixture.txt"

write.table(cpm, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

# Job Parameters used for this run:
# 
#     Date: 2020-06-12 12:57:16
#     Job type: Impute Cell Fractions
#     Signature matrix file: CIBERSORTx_Job17_10x_scgp_cibersortx_ref_06092020_inferred_phenoclasses.CIBERSORTx_Job17_10x_scgp_cibersortx_ref_06092020_inferred_refsample.bm.K999.txt
#     Mixture file: klemm_purified_cell_mixture.txt
#     Batch correction: enabled
#     Batch correction mode: S-mode
#     Single cell reference matrix file used for S-mode batch correction: 10x_scgp_cibersortx_ref_06092020.txt
#     Disable quantile normalization: true
#     Run mode (relative or absolute): relative
#     Permutations: 100
# 
