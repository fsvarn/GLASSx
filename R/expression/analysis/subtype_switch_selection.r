library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(ggbeeswarm)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
WITH subtype_label AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	ts1.signature_name, 
	ts1.p_value AS pval_a,
	ts2.p_value AS pval_b,
	CASE WHEN ts1.p_value < 0.05 THEN ts1.signature_name ELSE NULL END AS subtype_a,
	CASE WHEN ts2.p_value < 0.05 THEN ts2.signature_name ELSE NULL END AS subtype_b
	FROM analysis.rna_silver_set ss
	JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
	WHERE ts1.p_value < 0.05 OR ts2.p_value < 0.05
),
subtype_switch AS
(
	SELECT sl.tumor_pair_barcode,
	sl.case_barcode,
	string_agg(subtype_a,',') AS subtype_a,
	string_agg(subtype_b,',') AS subtype_b
	FROM subtype_label sl
	GROUP BY tumor_pair_barcode, sl.case_barcode
	ORDER BY 1
)
SELECT ss.*, 
si1.simplicity_score AS ss_a,
si2.simplicity_score AS ss_b,
CASE WHEN sw.subtype_a IS NULL THEN 'Mixed' ELSE sw.subtype_a END AS subtype_a,
CASE WHEN sw.subtype_b IS NULL THEN 'Mixed' ELSE sw.subtype_b END AS subtype_b,
CASE WHEN sw.subtype_a = sw.subtype_b THEN 'None' ELSE 'Switch' END AS switch,
su.treatment_tmz,
su.treatment_radiotherapy,
CASE WHEN su.treatment_chemotherapy_other LIKE '%Nivolumab%' OR su.treatment_chemotherapy_other LIKE '%Pembrolizumab%' THEN true ELSE false END AS treatment_pd1 ,
su.idh_codel_subtype,
se1.most_probable_classification AS select_a,
se2.most_probable_classification AS select_b
FROM analysis.rna_silver_set ss
JOIN analysis.rna_dna_pairs rd ON ss.tumor_barcode_a = rd.rna_barcode_a AND ss.tumor_barcode_b = rd.rna_barcode_b
JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.surgeries su ON su.sample_barcode = al1.sample_barcode
JOIN analysis.simplicity_score si1 ON ss.tumor_barcode_a = si1.aliquot_barcode
JOIN analysis.simplicity_score si2 ON ss.tumor_barcode_b = si2.aliquot_barcode
JOIN analysis.subclonalselection se1 ON rd.dna_barcode_a = se1.aliquot_barcode
JOIN analysis.subclonalselection se2 ON rd.dna_barcode_b = se2.aliquot_barcode
JOIN subtype_switch sw ON ss.tumor_pair_barcode = sw.tumor_pair_barcode
WHERE idh_codel_subtype IS NOT NULL -- ts1.p_value < 0.05 OR ts2.p_value < 0.05
ORDER BY 1, 2 
"

dat <- dbGetQuery(con, q)
dat <- dat[-which(duplicated(dat)),]

#dat <- dat[-which(dat[,"subtype_b"]=="Mixed"),]
g1 <- nrow(dat[which(dat[,"select_b"]=="S" & dat[,"switch"]=="Switch"),])
g2 <- nrow(dat[which(dat[,"select_b"]=="N" & dat[,"switch"]=="Switch"),])
g3 <- nrow(dat[which(dat[,"select_b"]=="S" & dat[,"switch"]!="Switch"),])
g4 <- nrow(dat[which(dat[,"select_b"]=="N" & dat[,"switch"]!="Switch"),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2)
fisher.test(ct)

#dat <- dat[-which(dat[,"subtype_b"]=="Mixed"),]
g1 <- nrow(dat[which(dat[,"select_b"]=="S" & dat[,"subtype_a"]=="Mesenchymal"),])
g2 <- nrow(dat[which(dat[,"select_b"]=="N" & dat[,"subtype_a"]=="Mesenchymal"),])
g3 <- nrow(dat[which(dat[,"select_b"]=="S" & dat[,"subtype_a"]!="Mesenchymal"),])
g4 <- nrow(dat[which(dat[,"select_b"]=="N" & dat[,"subtype_a"]!="Mesenchymal"),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2)
fisher.test(ct)

g1 <- nrow(dat[which(dat[,"select_b"]=="S" & dat[,"switch"]=="Switch" & dat[,"subtype_b"]=="Mesenchymal"),])
g2 <- nrow(dat[which(dat[,"select_b"]=="S" & (dat[,"switch"]!="Switch" | dat[,"subtype_b"]!="Mesenchymal")),])
g3 <- nrow(dat[which(dat[,"select_b"]=="N" & dat[,"switch"]=="Switch" & dat[,"subtype_b"]=="Mesenchymal"),])
g4 <- nrow(dat[which(dat[,"select_b"]=="N" & (dat[,"switch"]!="Switch" | dat[,"subtype_b"]!="Mesenchymal")),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2)
fisher.test(ct)

g1 <- nrow(dat[which(dat[,"select_b"]=="S" & dat[,"switch"]=="Switch" & dat[,"subtype_b"]=="Proneural"),])
g2 <- nrow(dat[which(dat[,"select_b"]=="S" & (dat[,"switch"]!="Switch" | dat[,"subtype_b"]!="Proneural")),])
g3 <- nrow(dat[which(dat[,"select_b"]=="N" & dat[,"switch"]=="Switch" & dat[,"subtype_b"]=="Proneural"),])
g4 <- nrow(dat[which(dat[,"select_b"]=="N" & (dat[,"switch"]!="Switch" | dat[,"subtype_b"]!="Proneural")),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2)
fisher.test(ct)

g1 <- nrow(dat[which(dat[,"select_b"]=="S" & dat[,"switch"]=="Switch" & dat[,"subtype_b"]=="Classical"),])
g2 <- nrow(dat[which(dat[,"select_b"]=="S" & (dat[,"switch"]!="Switch" | dat[,"subtype_b"]!="Classical")),])
g3 <- nrow(dat[which(dat[,"select_b"]=="N" & dat[,"switch"]=="Switch" & dat[,"subtype_b"]=="Classical"),])
g4 <- nrow(dat[which(dat[,"select_b"]=="N" & (dat[,"switch"]!="Switch" | dat[,"subtype_b"]!="Classical")),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2)
fisher.test(ct)


