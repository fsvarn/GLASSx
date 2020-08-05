###################################################
# Compare proliferating stem cell fraction in IDHmut receiving radiotherapy
# Updated: 2020.07.06
# Author: Frederick Varn
##################################################
library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
WITH
roman_to_int(grade, grade_int) AS 
(
	VALUES ('I'::text,1), ('II'::text,2), ('III'::text,3), ('IV'::text,4)
)
SELECT ss.*, 
CASE
	WHEN r2.grade_int > r1.grade_int THEN 'Grade up'::text
	WHEN r2.grade_int = r1.grade_int THEN 'Grade stable'::text
	WHEN r2.grade_int < r1.grade_int THEN 'Grade down'::text
	ELSE NULL::text
END AS grade_change,
su1.treatment_tmz,
su1.treatment_radiotherapy,
CASE WHEN su1.treatment_chemotherapy_other LIKE '%Nivolumab%' OR su1.treatment_chemotherapy_other LIKE '%Pembrolizumab%' THEN true ELSE false END AS treatment_pd1 ,
ci1.cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
su1.idh_codel_subtype
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.surgeries su1 ON su1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries su2 ON su2.sample_barcode = al2.sample_barcode
JOIN roman_to_int r1 ON r1.grade = su1.grade::text
JOIN roman_to_int r2 ON r2.grade = su2.grade::text
WHERE su1.idh_codel_subtype LIKE 'IDHmut%' AND ci1.cell_state = 'prolif_stemcell_tumor'
ORDER BY 1, 2 
"
dat <- dbGetQuery(con, q)

# Compare proliferating stem cell fraction in IDHmut

# Radiation receivers
rt_dat <- dat %>%
		  filter(treatment_radiotherapy==1)

p.val1 <- wilcox.test(rt_dat[,"fraction_a"], rt_dat[,"fraction_b"], paired=TRUE)	#0.01
eff1 <- mean(rt_dat[,"fraction_b"] - rt_dat[,"fraction_a"])							#0.06

# Non-radiation receivers
nrt_dat <- dat %>%
		  filter(treatment_radiotherapy==0)

p.val2 <- wilcox.test(nrt_dat[,"fraction_a"], nrt_dat[,"fraction_b"], paired=TRUE)	#0.14
eff2 <- mean(nrt_dat[,"fraction_b"] - nrt_dat[,"fraction_a"])						#0.03


