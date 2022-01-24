library(tidyverse)
library(data.table)


#######################################################
rm(list=ls())

myinf1 <- "/projects/verhaak-lab/varnf/data/xenahub/toil/TcgaTargetGtex_rsem_gene_tpm.txt"
myinf2 <- "/projects/verhaak-lab/varnf/data/xenahub/toil/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/subset/cibersortx_hires_tcga/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.txt"

toil <- fread(myinf1)
toil <- data.frame(toil)
gtex_anno <- read.delim(myinf2,stringsAsFactors=FALSE)
tcga_anno <- read.delim(myinf3)


# Get brain samples (cortex)
brain <- gtex_anno %>% 
		filter(SMTSD == "Brain - Frontal Cortex (BA9)" | SMTSD == "Brain - Cortex")
		
		
pull_normal <- brain[,1]
pull_normal <- pull_normal[which(pull_normal %in% colnames(toil))]

# Get TCGA samples
colnames(tcga_anno) <- gsub("\\.","-",colnames(tcga_anno))
pull_tumor <- colnames(tcga_anno)
pull_tumor <- sapply(strsplit(pull_tumor, "-"),function(x) paste(x[1:4],collapse="-"))
pull_tumor <- substr(pull_tumor, 1, 15)
pull_tumor <- pull_tumor[which(pull_tumor %in% colnames(toil))]

full_pull <- c("sample", pull_normal, pull_tumor)

# Assemble dataset
glioma_toil <- toil[,..full_pull]

# Save sub-dataset for now
write.table(glioma_toil, "/projects/verhaak-lab/varnf/data/xenahub/toil/TcgaTargetGtex_GBMLGGBRAIN_rsem_gene_tpm.txt", sep = "\t", quote=FALSE,row.names=FALSE)

#---------------------

# Map Ensembl to gene symbol
dat <- read.delim("/projects/verhaak-lab/varnf/data/xenahub/toil/TcgaTargetGtex_GBMLGGBRAIN_rsem_gene_tpm.txt",stringsAsFactor=FALSE)
myinf4 <- "/projects/verhaak-lab/varnf/data/xenahub/toil/probeMap_gencode.v23.annotation.gene.probemap"

gene_anno <- read.delim(myinf4,stringsAsFactor=FALSE)

gene_counts <- table(gene_anno$gene)
uni_gene <- names(gene_counts[which(gene_counts == 1)])
mul_gene <- names(gene_counts[which(gene_counts > 1)])

# One gene one ensembl ID
uni_anno <- gene_anno %>% filter(gene %in% uni_gene)
uni_map <- uni_anno$gene
names(uni_map) <- uni_anno$id

uni_mat <- dat[which(dat$sample %in% uni_anno$id),]
gene_symbol <- uni_map[uni_mat$sample]
uni_mat$sample <- gene_symbol

# Multiple genes multiple ensembl IDs (average?)
mul_anno <- gene_anno %>% filter(gene %in% mul_gene)
mul_map <- mul_anno$gene
names(mul_map) <- mul_anno$id

one_copy <- unique(mul_map)
mul_mat <- list()
for(i in 1:length(one_copy))
{
	mygene <- mul_map[which(mul_map == one_copy[i])]
	pull <- dat[which(dat$sample %in% names(mygene)),]
	avg_pull <- apply(pull[,2:ncol(pull)],2, mean)
	
	mul_mat[[i]] <- avg_pull
}
mul_mat <- do.call(rbind, mul_mat)
mul_mat <- data.frame(one_copy, mul_mat)
colnames(mul_mat)[1] <- "sample"

full_mat <- rbind(uni_mat, mul_mat)
rownames(full_mat) <- NULL
full_mat <- full_mat %>% column_to_rownames(var = "sample")

# Anti-log the values (downloaded them in log2(tpm + 0.001) format)
full_mat <- data.matrix(full_mat)
full_mat <- 2^full_mat - 0.001
full_mat[which(full_mat < 0)] <- 0	# Minimal value is 1e-8

full_mat <- full_mat %>% 
			data.frame() %>%
			rownames_to_column(var = "gene_symbol")
		
# Save this processed table
write.table(full_mat, "/projects/verhaak-lab/varnf/data/xenahub/toil/TcgaTargetGtex_GBMLGGBRAIN_rsem_gene_symbol_tpm_antilog.txt", sep = "\t", quote=FALSE,row.names=FALSE)

#---------------------

# Pull out common genes (saving to the CIBERSORTx running directory):
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/subset/cibersortx_hires_toil/TcgaTargetGtex_GBMLGGBRAIN_rsem_gene_symbol_tpm_antilog.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/subset/cibersortx_hires_toil/10x_scgp_cibersortx_ref_06092020.txt"

dat1 <- read.delim(myinf1,stringsAsFactors=FALSE)
dat2 <- read.delim(myinf2,stringsAsFactors=FALSE)

gene1 <- dat1[,1]
gene2 <- dat2[,1]

common <- intersect(gene1, gene2)
writeLines(common, "/projects/verhaak-lab/GLASS-III/data/subset/cibersortx_hires_toil/common_genes.txt")



