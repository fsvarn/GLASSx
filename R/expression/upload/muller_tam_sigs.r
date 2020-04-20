library(tidyverse)
library(GSVA)
library(DBI)
library(odbc)
library(reshape)

rm(list=ls())

#Load and prepare gene expression matrix
myinf1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"

expr <- read.delim(myinf1,row.names=1)
colnames(expr) <- gsub("\\.","-",colnames(expr))
expr <- as.matrix(expr)

#only use genes expressed in at least half of the cohort
sums <- apply(expr,1,function(x)sum(x>0))
#means <- apply(expr,1,mean)
expr <- expr[which(sums>(0.5*max(sums))),]

#Load Davoli gene sets
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT * FROM ref.immune_signatures WHERE signature_set = 'Muller'"

imm_sigs <- dbGetQuery(con, q)

#Convert table to signature list
immune_cell <- unique(imm_sigs[,"signature_name"])
gene_list <- list()
for(i in 1:length(immune_cell))
{
	sub_imm_sig <- imm_sigs[which(imm_sigs[,"signature_name"] == immune_cell[i]),]
	sig_genes <- sub_imm_sig[,"gene_symbol"]
	gene_list[[i]] <- sig_genes
}
names(gene_list) <- immune_cell

res <- gsva(expr, gene_list, method="ssgsea")

cor(t(res),method="s")

#             Microglia Macrophages
# Microglia    1.000000    0.415856
# Macrophages  0.415856    1.000000

dbWriteTable(con, Id(schema="analysis", table="muller_tam_score"), res, overwrite=TRUE)

