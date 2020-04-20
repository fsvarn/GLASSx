###################################################
# Test % of genome with CNV associates with immune levels
# Updated: 2020.04.16
# Author: Frederick Varn
##################################################

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
SELECT ps.case_barcode, 
an1.prop_aneuploidy AS prop_aneuploidy_a, 
an2.prop_aneuploidy AS prop_aneuploidy_b,
im1.signature_name,
im1.enrichment_score AS enrichment_score_a,
im2.enrichment_score AS enrichment_score_b,
pu1.purity AS purity_a,
pu2.purity AS purity_b,
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a 
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b 
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im1.signature_name = im2.signature_name
JOIN analysis.estimate pu1 ON pu1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate pu2 ON pu2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"
dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

#Immune score

mycor1 <- mycor2 <- mycor3 <- mycor4 <- matrix(0,nrow=length(cells),ncol=length(subtypes))
rownames(mycor1) <- rownames(mycor2) <- rownames(mycor3) <- rownames(mycor4) <- cells
colnames(mycor1) <- colnames(mycor2) <- colnames(mycor3) <- colnames(mycor4) <- subtypes
for(i in 1:length(cells))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"signature_name"] == cells[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		mycor1[i,j] <- cor(sub_dat[,"enrichment_score_a"],sub_dat[,"prop_aneuploidy_a"],method="s")
		mycor2[i,j] <- cor(sub_dat[,"enrichment_score_b"],sub_dat[,"prop_aneuploidy_b"],method="s")
		mycor3[i,j] <- cor(sub_dat[,"enrichment_score_a"],sub_dat[,"prop_aneuploidy_b"],method="s")
		mycor4[i,j] <- cor(sub_dat[,"enrichment_score_b"],sub_dat[,"prop_aneuploidy_a"],method="s")
	}
}

#ESTIMATE purity

pur_cor1 <- pur_cor2 <- rep(0,length(subtypes))
names(pur_cor1) <- names(pur_cor2) <- subtypes

for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"signature_name"] == cells[i] & dat[,"idh_codel_subtype"]==subtypes[i]),]
	pur_cor1[i] <- cor(sub_dat[,"purity_a"],sub_dat[,"prop_aneuploidy_a"],method="s")
	pur_cor2[i] <- cor(sub_dat[,"purity_b"],sub_dat[,"prop_aneuploidy_b"],method="s")
}


