###################################################
# Prepare single cell reference matrix for Neftel et al 
# Updated: 2020.08.03
# Author: Frederick Varn
##################################################

library(tidyverse)

rm(list=ls())

# Location of the Neftel et al single cell data
myinf1 <- "/projects/varnf/Data/Neftel_IDHwt_GBM/IDHwtGBM.processed.SS2.logTPM.txt"
myinf2 <- "/projects/varnf/Data/Neftel_IDHwt_GBM/IDHwt.GBM.Metadata.SS2.txt"

# Loads expression matrix where values = log2((tpm/10) + 1)
expr <- read.delim(myinf1,row.names=1)
colnames(expr) <- gsub("\\.","-",colnames(expr))
tpm <- ((2^expr) -1) * 10

# Check to make sure all cells have expression in genes
sums <- apply(tpm, 2, sum)
sum(sums == 0) # 0

# Read in info, remove pediatric tumors
info <- read.delim(myinf2,stringsAsFactor=FALSE)
info <- info[-1,]
info[,8:15] <- apply(info[,8:15],2,as.numeric)
info <- info[-which(info[,"GBMType"]=="Pediatric"),]
info <- info[-which(info[,"CellAssignment"]=="Malignant" & is.na(info[,"MESlike2"])),]
tpm <- tpm[,info[,"NAME"]]


# Subtype malignant cells
malig_info <- info[which(info[,"CellAssignment"] == "Malignant"),]
malig_tpm <- tpm[,malig_info[,"NAME"]]

malig_annot <- malig_info[,c("NAME","CellAssignment")]
malig_annot[,"subtype"] <- unlist(apply(malig_info[,8:13],1,function(x)which(x==max(x,na.rm=TRUE)))) %>%
		   			 recode("1" = "Mesenchymal", "2" = "Mesenchymal",
		   	      	 	"3" = "AClike", "4" = "OPClike",
		   	      		"5" = "NPClike", "6" = "NPClike")
		   	      		

# Add the non-malignant cells
nonmalig_info <- info[which(info[,"CellAssignment"] != "Malignant"),]
nonmalig_tpm <- tpm[,nonmalig_info[,"NAME"]]

nonmalig_annot <- nonmalig_info[,c("NAME","CellAssignment")]
nonmalig_annot[,"subtype"] <- nonmalig_annot[,"CellAssignment"]

# Combine the datasets together
annot <- rbind(malig_annot,nonmalig_annot)
full_tpm <- cbind(malig_tpm, nonmalig_tpm)
colnames(full_tpm) <- annot[,"subtype"]


# Add gene symbol as an independent column
gene <- rownames(full_tpm)

full_tpm <- data.frame(gene, full_tpm)

# Save file as a .txt file for input into CIBERSORTx
myoutf <- "/projects/verhaak-lab/GLASS-III/data/dataset/neftel_2019/neftel_cibersortx_ref_08032020.txt"
write.table(full_tpm, myoutf, sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#md5sum: ec0cafac25213b8653583a176c6e1d83


#     Date: 2020-08-03 10:18:42
#     Job type: Create Signature Matrix
#     Single cell reference matrix file: neftel_cibersortx_ref_08032020.txt
#     Disable quantile normalization: true
#     kappa: 999
#     q-value: 0.01
#     No. barcode genes: 300 to 500
#     Min. Expression: 1 (actually 0.75, CIBERSORTx rounds up for some reason)
#     Replicates: 0
#     Sampling: 0
#     Filter non-hematopoietic genes from signature matrix during construction: false


