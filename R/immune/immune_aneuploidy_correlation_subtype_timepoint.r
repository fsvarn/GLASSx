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
SELECT gs.case_barcode, 
an1.prop_aneuploidy AS prop_aneuploidy_a, 
an2.prop_aneuploidy AS prop_aneuploidy_b,
im1.signature_name,
im1.enrichment_score AS enrichment_score_a,
im2.enrichment_score AS enrichment_score_b,
idh_codel_subtype
FROM analysis.gold_set gs
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = gs.tumor_barcode_a 
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = gs.tumor_barcode_b 
JOIN analysis.rna_dna_pairs rd ON an1.aliquot_barcode = rd.dna_barcode_a AND an2.aliquot_barcode = rd.dna_barcode_b
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im1.signature_name = im2.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = gs.case_barcode
"
dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

mycor1 <- mycor2 <- matrix(0,nrow=length(cells),ncol=length(subtypes))
rownames(mycor1) <- rownames(mycor2) <- cells
colnames(mycor1) <- colnames(mycor2) <- subtypes
for(i in 1:length(cells))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"signature_name"] == cells[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		mycor1[i,j] <- cor(sub_dat[,"enrichment_score_a"],sub_dat[,"prop_aneuploidy_a"],method="s")
		mycor2[i,j] <- cor(sub_dat[,"enrichment_score_b"],sub_dat[,"prop_aneuploidy_b"],method="s")
	}
}

