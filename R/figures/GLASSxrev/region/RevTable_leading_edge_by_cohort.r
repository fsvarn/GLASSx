###################################################
# Compare CIBERSORTx profiles across Ivy GAP histological features
# Author: Frederick Varn
# Date: 2021.10.29
# Revision comment
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
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ci1.cell_state AS region,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_ivygap ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
WHERE ci1.cell_state = 'LE' AND idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

sum(dat$fraction_b > dat$fraction_a)/nrow(dat)				#61%
mean(dat$fraction_b - dat$fraction_a)						#10.5%
t.test(dat$fraction_b, dat$fraction_a, paired=TRUE)			#

#Read in data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ci1.cell_state AS region,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
cc.case_source
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_ivygap ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = ss.case_barcode
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
WHERE ci1.cell_state = 'LE' AND idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)
dat[which(grepl("TCGA-", dat$case_barcode)),"case_source"] <- "TCGA"
dat[which(grepl("CU-P", dat$case_barcode)),"case_source"] <- "PD1"


dat %>%
group_by(case_source) %>%
summarise(inc = sum(fraction_b > fraction_a), n = n(), pct = sum(fraction_b > fraction_a)/n())