###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_rec_061920/CIBERSORTxGEP_Job28_GEPs_Filtered.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_rec_061920/CIBERSORTxGEP_Job28_GEPs_StdErrs.txt"
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_init_061920/CIBERSORTxGEP_Job27_GEPs_Filtered.txt"
myinf4 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_init_061920/CIBERSORTxGEP_Job27_GEPs_StdErrs.txt"

#Recurrent gene expression profiles
geps1 <- read.delim(myinf1, row.names=1)
stderr1 <- read.delim(myinf2, row.names=1)
stderr1 <- stderr1[rownames(geps1),]

#Initial gene expression profiles
geps2 <- read.delim(myinf3, row.names=1)
stderr2 <- read.delim(myinf4, row.names=1)
stderr2 <- stderr2[rownames(geps2),]

gene_qvals <- gene_diffs <- list()
for (i in 1:ncol(geps1))
{
	cat("\r",i)
	vBetaZ <- sapply(1:nrow(geps1), function(j) (geps1[j,i]-geps2[j,i])/sqrt(stderr1[j,i]^2+stderr2[j,i]^2))
	fc <- sapply(1:nrow(geps1), function(j) log2(geps1[j,i] +1 ) - log2(geps2[j,i] + 1))
	names(vBetaZ) <- names(fc) <- rownames(geps1)
	ZPs <- 2*pnorm(-abs(vBetaZ))
	Zqvals <- p.adjust(ZPs, method="BH")
	fc[which(Zqvals < 0.05 & !is.na(Zqvals))]
	
	gene_qvals[[i]] <- Zqvals[which(Zqvals < 0.05 & !is.na(Zqvals))]
	gene_diffs[[i]] <- fc[which(Zqvals < 0.05 & !is.na(Zqvals))]
	
}
names(gene_qvals) <- names(gene_diffs) <- colnames(geps1)

sapply(gene_diffs, length)

gene_diffs$myeloid[which(gene_qvals$myeloid < 0.05)]
gene_diffs$t_cell[which(gene_qvals$t_cell < 0.05)]