library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


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
	ss.tumor_barcode_a,
	ss.tumor_barcode_b,
	ts1.signature_name, 
	ts1.p_value AS pval_a,
	ts2.p_value AS pval_b,
	CASE WHEN ts1.p_value < 0.05 THEN ts1.signature_name ELSE NULL END AS subtype_a,
	CASE WHEN ts2.p_value < 0.05 THEN ts2.signature_name ELSE NULL END AS subtype_b,
	su.idh_codel_subtype
	FROM analysis.rna_silver_set ss
	JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
	JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
	JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
	JOIN clinical.surgeries su ON su.sample_barcode = al1.sample_barcode
	WHERE idh_codel_subtype IS NOT NULL AND
	ts1.p_value < 0.05 OR ts2.p_value < 0.05
),
subtype_collapse AS
(
	SELECT sl.tumor_pair_barcode,
	sl.case_barcode,
	sl.tumor_barcode_a,
	sl.tumor_barcode_b,
	string_agg(subtype_a,',') AS subtype_a,
	string_agg(subtype_b,',') AS subtype_b,
	idh_codel_subtype
	FROM subtype_label sl
	GROUP BY tumor_pair_barcode, sl.case_barcode, sl.tumor_barcode_a, sl.tumor_barcode_b, idh_codel_subtype
	ORDER BY 1
)
SELECT sc.*,
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
im2.enrichment_score - im1.enrichment_score AS es_dif
FROM subtype_collapse sc
JOIN analysis.xue_macrophage_score im1 ON im1.aliquot_barcode = sc.tumor_barcode_a
JOIN analysis.xue_macrophage_score im2 ON im2.aliquot_barcode = sc.tumor_barcode_b AND im2.signature_name = im1.signature_name
WHERE idh_codel_subtype = 'IDHwt'
"
dat <- dbGetQuery(con, q)

modules <- unique(dat[,"signature_name"])

p.val <- eff <- rep(0,length(modules))
for(i in 1:length(modules))
{
	sub_dat <- dat %>% 
		filter(signature_name == modules[i])
	
	g1 <- sub_dat %>%
		filter(subtype_a == "Mesenchymal") %>%
		pull(es_a)
	g2 <- sub_dat %>%
		filter(subtype_a != "Mesenchymal") %>%
		pull(es_a)
	
	p.val[i] <- wilcox.test(g1,g2)$p.value
	eff[i] <- median(g1) - median(g2)
}
fdr <- p.adjust(p.val,"BH")

res <- data.frame(modules, eff, p.val, fdr)

p.val <- eff <- rep(0,length(modules))
for(i in 1:length(modules))
{
	sub_dat <- dat %>% 
		filter(signature_name == modules[i])
	
	g1 <- sub_dat %>%
		filter(subtype_a == "Mesenchymal" & subtype_b == "Mesenchymal") %>%
		pull(es_a)
	g2 <- sub_dat %>%
		filter(subtype_a == "Mesenchymal" & subtype_b == "Mesenchymal") %>%
		pull(es_b)
	
	p.val[i] <- wilcox.test(g1,g2,paired=TRUE)$p.value
	eff[i] <- median(g1) - median(g2)
}
fdr <- p.adjust(p.val,"BH")

res <- data.frame(modules, eff, p.val, fdr)