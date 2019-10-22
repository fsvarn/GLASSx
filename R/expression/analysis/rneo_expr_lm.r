library(DBI)
library(odbc)
library(tidyverse)

rm(list=ls())

#Load and prepare gene expression matrix
myinf1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"

genes <- read.delim(myinf1,row.names=1)
colnames(genes) <- gsub("\\.","-",colnames(genes))

#Load clinical data
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT al.aliquot_barcode, sa.sample_type, su.idh_codel_subtype, aliquot_batch, ROUND(nd.rneo,2) AS rneo
FROM biospecimen.aliquots al 
LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode
LEFT JOIN biospecimen.samples sa ON al.sample_barcode = sa.sample_barcode
JOIN analysis.analyte_sets an ON an.rna_barcode = al.aliquot_barcode
LEFT JOIN analysis.neoantigen_depletion nd ON nd.aliquot_barcode = an.dna_barcode
INNER JOIN analysis.gold_set gs ON an.dna_barcode = gs.tumor_barcode_a OR an.dna_barcode = gs.tumor_barcode_b
WHERE aliquot_analyte_type = 'R'"

clin <- dbGetQuery(con,q)

sub_genes <- genes[,clin[,"aliquot_barcode"]]
gene_sums <- apply(sub_genes,1,sum)
sub_genes <- sub_genes[which(gene_sums>0),]
sub_genes <- log10(sub_genes + 1)

#Create model comparing neoantigen depletion to expression when adjusting for subtype and batch
pval <- eff <- rep(0,nrow(sub_genes))
for(i in 1:nrow(sub_genes))
{
	cat("\r", i)
	mygene = as.numeric(sub_genes[i,])
	xx = cbind(mygene, clin)
	formula_vec = paste("rneo ~ mygene +", paste(colnames(clin)[3:4],collapse=" + ",sep=" "))
	mylm <- lm(formula_vec, data = xx)	
	mylm = summary(mylm)
	pval[i] =  mylm[["coefficients"]][2, 4]
	eff[i] =  mylm[["coefficients"]][2, 1]
}

#Multiple testing correction
fdr <- p.adjust(pval,"BH")

#Create results
res <- data.frame(rownames(sub_genes),pval,fdr,eff)
colnames(res)[1] <- "gene_symbol"
res <- res[order(res[,"pval"]),]

#Get the genes negatively associated with the ratio
res_neg <- res[which(res[,"eff"] < 0 & res[,"pval"] < 0.05),]