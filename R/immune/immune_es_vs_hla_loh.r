library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
WITH collapse AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, string_agg(hla_loh_change,'; ') AS hla_change
	FROM variants.hla_loh_change 
	GROUP BY case_barcode, tumor_barcode_a, tumor_barcode_b
	ORDER BY 1
)
SELECT co.*, im1.signature_name, im1.enrichment_score AS es_a, im2.enrichment_score AS es_b
FROM collapse co
JOIN analysis.rna_dna_pairs rd ON dna_barcode_a = tumor_barcode_a AND dna_barcode_b = tumor_barcode_b
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
"

dat <- dbGetQuery(con,q)

mycells <- unique(dat[,"signature_name"])

gain_pval <- gain_eff <- loss_pval <- loss_eff <- change_pval <- change_eff <- rep(0,length(mycells))
for(i in 1:length(mycells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==mycells[i]),]
	
	loh_dat <- sub_dat[grep("gain",sub_dat[,"hla_change"]),]
	gain_pval[i] <- wilcox.test(loh_dat[,"es_a"], loh_dat[,"es_b"],paired=TRUE)$p.value
	gain_eff[i] <- median(loh_dat[,"es_b"]) - median(loh_dat[,"es_a"])

	noloh_dat <- sub_dat[grep("loss",sub_dat[,"hla_change"]),]
	loss_pval[i] <- wilcox.test(noloh_dat[,"es_a"], noloh_dat[,"es_b"],paired=TRUE)$p.value
	loss_eff[i] <- median(noloh_dat[,"es_b"]) - median(noloh_dat[,"es_a"])
	
	
	sub_dat[,"infiltration_change"] <- sub_dat[,"es_b"] - sub_dat[,"es_a"]
	g1 <- sub_dat[grep("gain",sub_dat[,"hla_change"]),"infiltration_change"]
	g2 <- sub_dat[grep("gain",sub_dat[,"hla_change"],invert=TRUE),"infiltration_change"]

	change_pval[i] <- wilcox.test(g1,g2)$p.value		#P = 0.02
	change_eff[i] <- median(g2) - median(g1)
}

res <- data.frame(gain_pval, gain_eff, loss_pval, loss_eff, change_pval, change_eff)
rownames(res) <- mycells