library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(qusage)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/res/verhaak_prelim/msigdb_c2_enrichment_idhwt.txt"

es <- read.delim(myinf1)
hippo_es <- es["REACTOME_SIGNALING_BY_HIPPO",]
colnames(hippo_es) <- gsub("\\.","-",colnames(hippo_es))

hippo_es <- t(hippo_es)
hippo_es <- data.frame(rownames(hippo_es), hippo_es)
rownames(hippo_es) <- NULL
colnames(hippo_es) <- c("aliquot_barcode","es")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- 
"SELECT ps.rna_barcode_a AS aliquot_barcode, an1.prop_aneuploidy
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON ps.dna_barcode_a = an1.aliquot_barcode
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE idh_codel_subtype = 'IDHwt' AND an1.aliquot_barcode LIKE '%-TP-%'"

cnv <- dbGetQuery(con,q)
cnv <- cnv %>%
	   inner_join(y = hippo_es, by = "aliquot_barcode")

# Define responders as patients whose tumors do not come back within a year
thr <- median(cnv[,"prop_aneuploidy"])
an_hi <- cnv %>% filter(prop_aneuploidy >= thr) %>% .$es
an_lo <- cnv %>% filter(prop_aneuploidy < thr) %>% .$es


