library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggbeeswarm)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "SELECT ps.case_barcode, 
ps.rna_barcode_a, 
ps.rna_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b,
CASE WHEN mf1.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf1.coverage_adj_mut_freq < 10 THEN 0 END AS hm_a,
CASE WHEN mf2.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf2.coverage_adj_mut_freq < 10 THEN 0 END AS hm_b
FROM analysis.platinum_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = ps.dna_barcode_b
ORDER BY 1, 2, 6"

hm <- dbGetQuery(con, q)

#See if there are significant changes across samples over time

#Where are those missing TCGA samples???? (Microarray??)

#hm_only <- hm[which(hm[,"hm_b"] ==1 & hm[,"subtype_a"]!="IDHwt"),]

cells <- unique(hm[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_hm <- hm_only[which(hm_only[,"signature_name"]==cells[i]),]
	eff[i] <- median(sub_hm[,"es_b"]) - median(sub_hm[,"es_a"])
	p.value[i] <- wilcox.test(sub_hm[,"es_a"], sub_hm[,"es_b"],paired=TRUE)$p.value
}
time_res <- data.frame(cells,eff,p.value)