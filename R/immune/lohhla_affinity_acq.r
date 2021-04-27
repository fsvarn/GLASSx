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
WITH lohhla_pairs AS
(
	SELECT ps.*,
	lh1.hla_type1, lh1.hla_type2, 
	lh1.hla_type1_copy_number AS hla_type1_copy_number_a, lh2.hla_type1_copy_number AS hla_type1_copy_number_b,
	lh1.hla_type2_copy_number AS hla_type2_copy_number_a, lh2.hla_type2_copy_number AS hla_type2_copy_number_b,
	lh1.pval AS pval_a, lh1.loss_allele AS loss_allele_a, lh1.kept_allele AS kept_allele_a,  
	lh2.pval AS pval_b, lh2.loss_allele AS loss_allele_b, lh2.kept_allele AS kept_allele_b,
	CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
	FROM analysis.silver_set ps
	JOIN analysis.pairs pa1 ON pa1.tumor_barcode = ps.tumor_barcode_a
	JOIN analysis.pairs pa2 ON pa2.tumor_barcode = ps.tumor_barcode_b
	JOIN variants.lohhla lh1 ON lh1.pair_barcode = pa1.pair_barcode
	JOIN variants.lohhla lh2 ON lh2.pair_barcode = pa2.pair_barcode AND lh2.hla_type1 = lh1.hla_type1
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	WHERE lh1.coverage_filter = 20 AND lh2.coverage_filter = 20
),
acq_loss AS
(
	SELECT tumor_pair_barcode, loss_allele_b, kept_allele_b, idh_status
	FROM lohhla_pairs lp
	WHERE pval_a > 0.05 AND pval_b < 0.05 AND (hla_type1_copy_number_b < 0.5 OR hla_type2_copy_number_b < 0.5) 
),
neoag AS
(
	SELECT tumor_pair_barcode, variant_id, transcript, pvacseq_variant_type, gene_name, hla_allele, netmhcpan_mt_score
	FROM analysis.neoantigens_by_pair
	WHERE fraction = 'R'
),
loss_neoag AS
(
	SELECT ac.*, AVG(nl.netmhcpan_mt_score) AS loss_affinity, COUNT(*) AS loss_count
	FROM acq_loss ac
	JOIN neoag nl ON nl.tumor_pair_barcode = ac.tumor_pair_barcode AND nl.hla_allele = substring(ac.loss_allele_b, 1, 11)
	GROUP BY ac.tumor_pair_barcode, loss_allele_b, kept_allele_b, idh_status
),
kept_neoag AS
(
	SELECT ac.*, AVG(nk.netmhcpan_mt_score) AS kept_affinity, COUNT(*) AS kept_count
	FROM acq_loss ac
	JOIN neoag nk ON nk.tumor_pair_barcode = ac.tumor_pair_barcode AND nk.hla_allele = substring(ac.kept_allele_b, 1, 11)
	GROUP BY ac.tumor_pair_barcode, loss_allele_b, kept_allele_b, idh_status
)
SELECT lo.*, ke.kept_affinity, ke.kept_count
FROM loss_neoag lo
JOIN kept_neoag ke ON ke.tumor_pair_barcode = lo.tumor_pair_barcode AND ke.loss_allele_b = lo.loss_allele_b

"

dat <- dbGetQuery(con,q)

dat %>% 
group_by(idh_status) %>%
summarise(pval = wilcox.test(loss_affinity, kept_affinity, paired=TRUE)$p.value,
		  eff = median(loss_affinity - kept_affinity))
