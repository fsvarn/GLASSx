###################################################
# Create synthetic mixtures from the single cell data in the Neftel et al dataset to validate the SCGP signature matrix
# Reference for dataset: An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma: Cell
# Updated: 2020.06.11
# Author: Frederick Varn
##################################################

library(tidyverse)

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
info <- info[-which(info[,"CellAssignment"]=="Malignant" & is.na(info[,"MESlike2"])),]

malig <- info[which(info[,"CellAssignment"] == "Malignant"),]

# Assign each malignant cell a subtype
annot <- malig[,c("NAME","CellAssignment")]
annot[,"subtype"] <- apply(malig[,8:13],1,function(x)which(x==max(x))) %>%
		   			 recode("1" = "Mesenchymal", "2" = "Mesenchymal",
		   	      	 	"3" = "AClike", "4" = "OPClike",
		   	      		"5" = "NPClike", "6" = "NPClike")
		   	      		
# Identify cycling malignant cells and convert them
annot[,"cyc_dist"] <- malig[,"G1S"] > 1 | malig[,"G2M"] > 1
annot[which(annot[,"cyc_dist"] & annot[,"subtype"] %in% c("OPClike","NPClike")),"subtype"] <- "prolif_stemcell_tumor"

# Add the non-malignant cells
nonmalig <- info[which(info[,"CellAssignment"] != "Malignant"),c("NAME","CellAssignment")]
nonmalig[,"subtype"] <- nonmalig[,"CellAssignment"]
nonmalig[,"cyc_dist"] <- NA
annot <- rbind(annot,nonmalig)

# Convert the subtypes to SCGP subtypes
annot[,"subtype"] <- annot[,"subtype"] %>%
					 recode("Mesenchymal" = "differentiated_tumor", "AClike" = "differentiated_tumor", 
					 "OPClike" = "stemcell_tumor", "NPClike" = "stemcell_tumor",
					 "Macrophage" = "myeloid", "Oligodendrocyte" = "oligodendrocyte", "T-cell" = "t_cell")

table(annot[,"subtype"])
#  differentiated_tumor               myeloid       oligodendrocyte 
#                  2838                   536                   210 
# prolif_stemcell_tumor        stemcell_tumor                t_cell 
#                   440                  1595                    80  

sub_tpm <- tpm[,annot[,"NAME"]]
	   	      		

	   	      		
# Assign each cell a subtype
colnames(sub_tpm) <- annot[,"subtype"]

# Anti-log data to get it out of log2 space
tpm_mix <- (2^sub_tpm) - 1

# Build mixture matrix
uni_cell <- unique(annot[,"subtype"])
myprop <- seq(0,1,by=0.1)

synth_mixtures <- matrix(0, nrow=nrow(tpm_mix), ncol=length(uni_cell) * length(myprop))
rownames(synth_mixtures) <- rownames(tpm_mix)
proportions <- matrix(0, nrow=length(uni_cell), ncol=length(uni_cell) * length(myprop))
rownames(proportions) <- uni_cell

ind <- 1
for(i in 1:length(uni_cell))
{
	sc_cell <- tpm_mix[,which(sapply(strsplit(colnames(tpm_mix),"\\."), "[[", 1) == uni_cell[i]),]
	outer_cell <- tpm_mix[,-which(sapply(strsplit(colnames(tpm_mix),"\\."), "[[", 1) == uni_cell[i]),]
	
	for(j in 1:length(myprop))
	{
		cat("\r", i, "-->", j)
	
		ctrl_num <- min(c(round(myprop[j] * ncol(sc_cell)), round(myprop[j] * 250)))
		ctrl_cell <- sc_cell[,sample(1:ncol(sc_cell), ctrl_num)]
		
		noise_num <- 250 - ctrl_num
		noise_cell <- outer_cell[,sample(1:ncol(outer_cell), noise_num)]
		
		mixture_set <- cbind(ctrl_cell, noise_cell)
		fracs <- table(sapply(strsplit(colnames(mixture_set),"\\."), "[[", 1))/ncol(mixture_set)
		fracs <- fracs[uni_cell]
		proportions[,ind] <- fracs
		
		synth_mixtures[,ind] <- apply(mixture_set, 1, sum)
		
		ind <- ind + 1
	}
}

proportions[which(is.na(proportions))] <- 0
mix_names <- paste("mixture_",1:ncol(synth_mixtures),sep="")
colnames(synth_mixtures) <- colnames(proportions) <- mix_names


# Save file as a .txt file for input into CIBERSORTx
myoutf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/SS2_neftel_synthetic_mix_gep_06092020.txt"
myoutf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/SS2_neftel_synthetic_mix_prop_06092020.txt"

write.table(synth_mixtures, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(proportions, myoutf2, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

#mix md5sum: 1088e6469d42f9bdb9074f78cada1f5f
#prop md5sum: d861c646e5c72f7709ea3ca674bf7c0a

# Job Parameters used for this run:
# 
#     Date: 2020-06-11 14:52:56
#     Job type: Impute Cell Fractions
#     Signature matrix file: CIBERSORTx_Job17_10x_scgp_cibersortx_ref_06092020_inferred_phenoclasses.CIBERSORTx_Job17_10x_scgp_cibersortx_ref_06092020_inferred_refsample.bm.K999.txt
#     Mixture file: SS2_neftel_synthetic_mix_gep_06092020.txt
#     Batch correction: enabled
#     Batch correction mode: S-mode
#     Single cell reference matrix file used for S-mode batch correction: 10x_scgp_cibersortx_ref_06092020.txt
#     Disable quantile normalization: true
#     Run mode (relative or absolute): relative
#     Permutations: 100
# 

