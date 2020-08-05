###################################################
# Calculate TPM for purified macrophage profiles
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.05.29
# Author: Frederick Varn
##################################################

library(tidyverse)
library(DBI)
library(ssgsea.GBM.classification)
library(DESeq2)
library(GSVA)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

#Establish connection to db
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
ref <- dbReadTable(con, Id(schema="ref",table="genes"))

expr <- read.csv(myinf1,row.names=1)
info <- read.csv(myinf2,row.names=1,stringsAsFactor=FALSE)

# Add a batch column to info
info[,"batch"] <- sapply(strsplit(rownames(info),"_"),function(x)x[1])

##################################################
# Helper functions
##################################################

# overrepresentation test of macrophage modules
#-------------------------------------------

or_test <- function(dif, sig, bg)
{
	# Differentially expressed and in signature
	g1 <- length(bg[which(bg %in% dif & bg %in% sig)])
	# Differentially expressed and not in signature
	g2 <- length(bg[which(bg %in% dif & !(bg %in% sig))])
	# Not differentially expressed, in signature
	g3 <- length(bg[which(!(bg %in% dif) & bg %in% sig)])
	# Not differentially expressed, not in signature
	g4 <- length(bg[which(!(bg %in% dif) & !(bg %in% sig))])
	
	p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2))$p.value
	
	return(p.val)
}
	
##################################################
# Step 4: Differential expression analysis between tumor macrophages and normal bmdm, tumor microglia and normal microglia
##################################################

# Macrophages
mac_cts <- expr[,grep("_mdm",colnames(expr))]
#mac_cts <- mac_cts[,-grep("brm",colnames(mac_cts))]

batch <- sapply(strsplit(colnames(mac_cts),"_"),function(x)x[1])
sample_type <- sapply(strsplit(colnames(mac_cts),"_"),function(x)x[2])
tumor <- !grepl("healthyBlood",sample_type)

mac_info <- data.frame(batch, sample_type, tumor, row.names=colnames(mac_cts))

# Use DESeq2 to compare subtype versus no subtype
mac_dds <- DESeqDataSetFromMatrix(countData = mac_cts,
							  colData = mac_info,
							  design= ~ batch + tumor)
dds <- DESeq(mac_dds)
resultsNames(dds)
mac_res <- results(dds, name="tumorTRUE")
mac_up <- rownames(mac_res[which(mac_res[,"padj"]< 0.01 & mac_res[,"log2FoldChange"]> 1),])
length(mac_up)	#3063
mac_genes <- as.integer(mac_res[,"padj"]<.05 & mac_res[,"log2FoldChange"]>1)
names(mac_genes) <- rownames(mac_res)
mac_genes <- mac_genes[-which(is.na(mac_genes))]

# Microglia
mg_cts <- expr[,grep("_mg",colnames(expr))]
#mg_cts <- mg_cts[,-grep("brm",colnames(mg_cts))]

batch <- sapply(strsplit(colnames(mg_cts),"_"),function(x)x[1])
sample_type <- sapply(strsplit(colnames(mg_cts),"_"),function(x)x[2])
tumor <- !grepl("nonTumor",sample_type)

mg_info <- data.frame(batch, sample_type, tumor, row.names=colnames(mg_cts))

# Use DESeq2 to compare subtype versus no subtype
mg_dds <- DESeqDataSetFromMatrix(countData = mg_cts,
							  colData = mg_info,
							  design= ~ batch + tumor)
dds <- DESeq(mg_dds)
resultsNames(dds)
mg_res <- results(dds, name="tumorTRUE")
mg_up <- rownames(mg_res[which(mg_res[,"padj"]< 0.01 & mg_res[,"log2FoldChange"]> 1),])
length(mg_up)	#1186


##################################################
# Step 5: Overrepresentation analysis of different macrophage activation states
##################################################

pwf=nullp(mac_genes,"hg19","geneSymbol")

GO.wall=goseq(pwf,"hg19","geneSymbol")

module_inf <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"

mac_module_df <- read.delim(module_inf,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- gsub("X","module",colnames(mac_module_df))

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)


measured_genes <- rownames(expr)
mac_pval <- mg_pval <- rep(0, length(mac_modules))
names(mac_pval) <- names(mg_pval) <- names(mac_modules)
for(i in 1:length(mac_modules))
{
	mymodule <- mac_modules[[i]]
	mymodule <- mymodule[which(mymodule %in% measured_genes)]
	
	mac_pval[i] <- or_test(mac_up, mymodule, measured_genes)
	mg_pval[i] <- or_test(mg_up, mymodule, measured_genes)
}

mac_qval <- p.adjust(mac_pval,"BH")
mg_qval <- p.adjust(mac_pval,"BH")

which(mac_qval[,3] < 0.1)
# module41 
#       41 
which(mg_qval[,3] < 0.1)
#  module3  module6 module15 
#        3        6       15 