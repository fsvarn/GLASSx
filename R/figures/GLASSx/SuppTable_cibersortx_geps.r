library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)
library(survival)
library(topGO)

#######################################################
rm(list=ls())

#  Read in post-treatment stem cell signature
myDir1 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/analysis/"
mysigs <- dir(myDir1)
mysigs <- paste(myDir1, mysigs, sep = "")[1:6]

idh_status <- c(rep("IDHmut",3), rep("IDHwt",3))
cell_state <- rep(c("Diff.-like","Prolif. stem-like","Stem-like"),2)
res <- list()
for(i in 1:length(mysigs))
{
	dat <- read.delim(mysigs[i])
	
	dat <- dat %>% filter(sig) %>% rownames_to_column(var="gene_symbol")
	dat <- dat[,1:4]
	dat$cell_state <- cell_state[i]
	dat$idh_status <- idh_status[i]
	dat$eff <- 2^dat$eff
	
	res[[i]] <- dat
}

long_res <- do.call(rbind, res)

ord <- rep(0, nrow(long_res))
ord[which(long_res$idh_status=="IDHwt" & long_res$cell_state == "Diff.-like")] <- 1
ord[which(long_res$idh_status=="IDHwt" & long_res$cell_state == "Stem-like")] <- 2
ord[which(long_res$idh_status=="IDHwt" & long_res$cell_state == "Prolif. stem-like")] <- 3
ord[which(long_res$idh_status=="IDHmut" & long_res$cell_state == "Diff.-like")] <- 4
ord[which(long_res$idh_status=="IDHmut" & long_res$cell_state == "Stem-like")] <- 5
ord[which(long_res$idh_status=="IDHmut" & long_res$cell_state == "Prolif. stem-like")] <- 6

long_res <- long_res[order(ord),]

write.table(long_res, "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/analysis/SuppTable_posttreatment_sigs.txt", sep = "\t", quote=FALSE, row.names=FALSE)

#---------------------------------------------


#  Read in post-treatment stem cell signature
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v3.txt"
mysig <- read.delim(myinf1)

mysig <- mysig %>% filter(q.val < 0.1, eff > log2(1.1)) %>% rownames_to_column(var="gene_symbol")
mysig <- mysig[,1:4]
mysig$eff <- 2^mysig$eff

write.table(mysig, "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/analysis/SuppTable_mes_myeloid_sig.txt", sep = "\t", quote=FALSE, row.names=FALSE)

