###################################################
# Prepare single cell reference matrix for CIBERSORTx
# Updated: 2020.02.11
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
log2cpm_annot <- log2cpm_annot[,mysamps]

table(sapply(strsplit(colnames(log2cpm_annot),"\\."), "[[", 1))
#                b_cell        dendritic_cell  differentiated_tumor 
#                     4                    23                  1669 
#           endothelial            fibroblast           granulocyte 
#                    14                     7                    87 
#               myeloid       oligodendrocyte              pericyte 
#                  1474                   374                    12 
# prolif_stemcell_tumor        stemcell_tumor                t_cell 
#                   357                   952                    27 

# Anti-log data to get it out of log2 space

# cpm_annot <- (2^log2cpm_annot) - 1 							Original way of doing things (6/10/2020)
cpm_annot <- exp(log2cpm_annot) - 1 							#Per Kevin, this is actually proper (2/11/21)


# Add gene symbol as an independent column
featuredata <- featuredata[rownames(cpm_annot),]

gene <- featuredata[,"Associated.Gene.Name"]
cpm_annot <- data.frame(gene, cpm_annot)

# Save file as a .txt file for input into CIBERSORTx
#myoutf <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_cibersortx_ref_06092020.txt"
myoutf <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_cibersortx_ref_02112021_TEST.txt"
write.table(cpm_annot, myoutf, sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#md5sum: 608182eebc83e2d38991103ff796204c

# Job parameters used to make signature matrix on CIBERSORTx web browser
# Date: 2020-06-10 12:03:07
# Job type: Create Signature Matrix
# Single cell reference matrix file: 10x_scgp_cibersortx_ref_06092020.txt
# Disable quantile normalization: true
# kappa: 999
# q-value: 0.01
# No. barcode genes: 300 to 500
# Min. Expression: 1
# Replicates: 3
# Sampling: 0
# Filter non-hematopoietic genes from signature matrix during construction: false

