library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(survival)

# Need more IDH mutant samples for this analysis to be good

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "WITH neoag_by_ali AS
(
	SELECT aliquot_barcode, variant_id, gene_name, mutation, pvacseq_protein_position, peptide_length, sub_peptide_position, mt_epitope_seq 
	FROM analysis.neoantigens_by_aliquot neo
	WHERE ssm2_pass_call = TRUE 
	GROUP BY aliquot_barcode, variant_id, gene_name, mutation, pvacseq_protein_position, peptide_length, sub_peptide_position, mt_epitope_seq
),
idh_counts AS
(
	SELECT aliquot_barcode, 
	CASE 
		WHEN gene_name = 'IDH1' AND mutation = 'R/H' AND pvacseq_protein_position = '132' THEN 1
		ELSE 0
	END AS idh1_neoag
	FROM neoag_by_ali
),
idh_neoag AS
(
	SELECT aliquot_barcode, SUM(idh1_neoag) 
	FROM idh_counts
	GROUP BY aliquot_barcode
)
SELECT ps.case_barcode, 
im1.shannon AS shannon_a, 
im2.shannon AS shannon_b, 
im1.evenness AS evenness_a, 
im2.evenness AS evenness_b,
im1.richness AS richness_a, 
im2.richness AS richness_b, 
im1.total_tcr AS total_tcr_a, 
im2.total_tcr AS total_tcr_b, 
na1.sum AS idh_neoag_a, 
na2.sum AS idh_neoag_b,
bs1.sample_type AS sample_type_a
FROM analysis.platinum_set ps
JOIN analysis.tcr_stats im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.tcr_stats im2 ON im2.aliquot_barcode = ps.rna_barcode_b
JOIN idh_neoag na1 ON na1.aliquot_barcode = ps.dna_barcode_a
JOIN idh_neoag na2 ON na2.aliquot_barcode = ps.dna_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
WHERE cs1.idh_status = 'IDHmut" AND bs1.sample_type = 'TP'
ORDER BY 2"

dat <- dbGetQuery(con,q)

wilcox.test(dat[which(dat[,"idh_neoag_a"]==1),"evenness_a"],dat[which(dat[,"idh_neoag_a"]==0),"evenness_a"])