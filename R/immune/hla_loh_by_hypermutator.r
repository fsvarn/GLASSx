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
SELECT co.*, mf1.coverage_adj_mut_freq AS mut_freq_a, mf2.coverage_adj_mut_freq AS mut_freq_b,
CASE WHEN mf2.coverage_adj_mut_freq > 10 THEN 1 ELSE 0 END AS hypermutator,
CASE WHEN co.hla_change LIKE '%gain%' OR co.hla_change LIKE '%stable%' THEN 1 ELSE 0 END AS rec_loh,
CASE WHEN co.hla_change LIKE '%loss%' OR co.hla_change LIKE '%stable%' THEN 1 ELSE 0 END AS init_loh,
cs.idh_codel_subtype
FROM collapse co
JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = co.tumor_barcode_a
JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = co.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = co.case_barcode
--WHERE co.hla_change NOT LIKE '%stable%'
ORDER BY 6 DESC
"

dat <- dbGetQuery(con,q)
#dat <- dat[which(dat[,"idh_codel_subtype"]!="IDHwt"),]

g1 <- sum(dat[which(dat[,"hypermutator"]==1),"rec_loh"])
g2 <- sum(dat[which(dat[,"hypermutator"]==0),"rec_loh"])
g3 <- sum(dat[which(dat[,"hypermutator"]==1),"rec_loh"]==0)
g4 <- sum(dat[which(dat[,"hypermutator"]==0),"rec_loh"]==0)
ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)		
#P = 0.08 for having LOH at recurrence (52% in hypermutators versus 33% in non-hypermutators)
#P = 0.08 for having LOH at recurrence in IDHmut (62.5% in hypermutators versus 36% in non-hypermutators)

