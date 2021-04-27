###################################################
# Create synthetic mixtures from the single cell data to validate the single cell signature matrix
# Updated: 2020.06.10
# Author: Frederick Varn
##################################################

library(tidyverse)

# Location of the SCGP single cell data
myinf1 <- "/projects/verhaak-lab/scgp/data/10x/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds"

load(myinf1)

# Remove all cells with no expression
#sums <- apply(log2cpm, 2, sum)	
#range(sums)
# [1]  2453.769 45315.686
# No cells with 0 expression

# Remove QC genes:
qc_genes <- c("ENSGGENES","ENSGUMI","ENSGMITO", "ENSGSEQSAT","ENSGSAMP") 
log2cpm <- log2cpm[-which(rownames(log2cpm) %in% qc_genes),]
featuredata <- featuredata[-which(rownames(featuredata) %in% qc_genes),]

# Annotate clusters using previous definitions
clust_annot = tsne.data %>%
	rownames_to_column('cell') %>%
	mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                              `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                              `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) %>%
	column_to_rownames('cell')

# clust_annot file and log2cpm files are in the same order so no need to match the order between them
cell_type <- clust_annot[,"cell_type"]

# Assign each cell a subtype
log2cpm_annot <- log2cpm
colnames(log2cpm_annot) <- cell_type

# Downsample to 5000 per CIBERSORTx instructions

# Step 1: Remove low quality sample SM008:
sample_id <- sapply(strsplit(rownames(clust_annot), "-"), "[[", 3)
sample_id <- gsub("^6$", "SM008", sample_id)

log2cpm_annot <- log2cpm_annot[-which(sample_id == "SM008")]

# Step 2: Randomly downsample to 5000 cells using seed 11
set.seed(11)
mysamps <- sample(1:ncol(log2cpm_annot),5000)
log2cpm_sig <- log2cpm_annot[,mysamps]

table(sapply(strsplit(colnames(log2cpm_sig),"\\."), "[[", 1))
#                b_cell        dendritic_cell  differentiated_tumor 
#                     4                    23                  1669 
#           endothelial            fibroblast           granulocyte 
#                    14                     7                    87 
#               myeloid       oligodendrocyte              pericyte 
#                  1474                   374                    12 
# prolif_stemcell_tumor        stemcell_tumor                t_cell 
#                   357                   952                    27 

# Now that we have the cells that were involved in creating the signature matrix, we can use the rest to make mixtures

log2cpm_mix <- log2cpm_annot[,-mysamps]

table(sapply(strsplit(colnames(log2cpm_mix),"\\."), "[[", 1))
#                b_cell        dendritic_cell  differentiated_tumor 
#                    47                   302                 14688 
#           endothelial            fibroblast           granulocyte 
#                   124                    68                   684 
#               myeloid       oligodendrocyte              pericyte 
#                 13501                  2941                    74 
# prolif_stemcell_tumor        stemcell_tumor                t_cell 
#                  3286                  9368                   295 

# Anti-log data to get it out of log2 space
cpm_mix <- (2^log2cpm_mix) - 1

# Convert ensmebl id to gene_symbol
featuredata <- featuredata[rownames(cpm_mix),]
rownames(cpm_mix) <- featuredata[,"Associated.Gene.Name"]

# Build mixture matrix
uni_cell <- unique(cell_type)
myprop <- seq(0,1,by=0.1)

synth_mixtures <- matrix(0, nrow=nrow(cpm_mix), ncol=length(uni_cell) * length(myprop))
rownames(synth_mixtures) <- rownames(cpm_mix)
proportions <- matrix(0, nrow=length(uni_cell), ncol=length(uni_cell) * length(myprop))
rownames(proportions) <- uni_cell

ind <- 1
for(i in 1:length(uni_cell))
{
	sc_cell <- cpm_mix[,which(sapply(strsplit(colnames(cpm_mix),"\\."), "[[", 1) == uni_cell[i]),]
	outer_cell <- cpm_mix[,-which(sapply(strsplit(colnames(cpm_mix),"\\."), "[[", 1) == uni_cell[i]),]
	
	for(j in 1:length(myprop))
	{
		cat("\r", i, "-->", j)
	
		ctrl_num <- min(c(round(myprop[j] * ncol(sc_cell)), round(myprop[j] * 5000)))
		ctrl_cell <- sc_cell[,sample(1:ncol(sc_cell), ctrl_num)]
		
		noise_num <- 5000 - ctrl_num
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
myoutf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_synthetic_mix_gep_06092020.txt"
myoutf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_synthetic_mix_prop_06092020.txt"

#write.table(synth_mixtures, myoutf1, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
#write.table(proportions, myoutf2, sep = "\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

#mix md5sum: d509fa7dc09aeaabc54a160d2a420b49
#prop md5sum: b0d0a6d649e2fc41fbec0a5958258d32

# Job Parameters used for this run:
# 
#     Date: 2020-06-10 13:45:33
#     Job type: Impute Cell Fractions
#     Signature matrix file: CIBERSORTx_Job17_10x_scgp_cibersortx_ref_06092020_inferred_phenoclasses.CIBERSORTx_Job17_10x_scgp_cibersortx_ref_06092020_inferred_refsample.bm.K999.txt
#     Mixture file: 10x_scgp_synthetic_mix_gep_06092020.txt
#     Batch correction: disabled
#     Disable quantile normalization: true
#     Run mode (relative or absolute): relative
#     Permutations: 100
