library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggdendro)
library(grid)

#######################################################
rm(list=ls())
#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")


#Read in data
q <- "
SELECT gs.*, mf1.coverage_adj_mut_freq AS mut_freq_a, mf2.coverage_adj_mut_freq AS mut_freq_b, 
ss1.simplicity_score AS simplicity_score_a, ss2.simplicity_score AS simplicity_score_b, 
CASE WHEN cs.idh_codel_subtype LIKE 'IDHmut%' THEN 'IDHmut' ELSE 'IDHwt' END AS idh_codel_subtype
FROM analysis.gold_set gs
JOIN analysis.mut_freq mf1 ON gs.tumor_barcode_a = mf1.aliquot_barcode
JOIN analysis.mut_freq mf2 ON gs.tumor_barcode_b = mf2.aliquot_barcode
JOIN analysis.analyte_sets an1 ON mf1.aliquot_barcode = an1.dna_barcode
JOIN analysis.analyte_sets an2 ON mf2.aliquot_barcode = an2.dna_barcode
JOIN analysis.simplicity_score ss1 ON ss1.aliquot_barcode = an1.rna_barcode
JOIN analysis.simplicity_score ss2 ON ss2.aliquot_barcode = an2.rna_barcode
JOIN clinical.subtypes cs ON cs.case_barcode = gs.case_barcode
ORDER BY 4 DESC
"

dat <- dbGetQuery(con,q)

#Initial tumors
#--------------------------

cor.test(dat[,"mut_freq_a"],dat[,"simplicity_score_a"],method="s")		#-0.21, P = 0.04

g1 <- dat[which(dat[,"simplicity_score_a"] < median(dat[,"simplicity_score_a"])),"mut_freq_a"]
g2 <- dat[which(dat[,"simplicity_score_a"] >= median(dat[,"simplicity_score_a"])),"mut_freq_a"]
wilcox.test(g1,g2)		#P = 0.05


#Recurrent tumors
#--------------------------

cor.test(dat[,"mut_freq_b"],dat[,"simplicity_score_b"],method="s")		#-0.03, P = 0.74

g1 <- dat[which(dat[,"simplicity_score_b"] < median(dat[,"simplicity_score_b"])),"mut_freq_b"]
g2 <- dat[which(dat[,"simplicity_score_b"] >= median(dat[,"simplicity_score_b"])),"mut_freq_b"]
wilcox.test(g1,g2)		#P = 0.83


#PRIMARY tumors
#--------------------------
sub_dat <- dat[grep("-TP-",dat[,"tumor_barcode_a"]),]
cor.test(sub_dat[,"mut_freq_a"],sub_dat[,"simplicity_score_a"],method="s")		#-0.20, P = 0.05

g1 <- dat[which(sub_dat[,"simplicity_score_a"] < median(sub_dat[,"simplicity_score_a"])),"mut_freq_a"]
g2 <- dat[which(sub_dat[,"simplicity_score_a"] >= median(sub_dat[,"simplicity_score_a"])),"mut_freq_a"]
wilcox.test(g1,g2)		#P = 0.42


#Recurrent tumors ignoring hypermutators
#--------------------------
sub_dat <- dat[-which(dat[,"mut_freq_b"]>=10),]
cor.test(sub_dat[,"mut_freq_b"],sub_dat[,"simplicity_score_b"],method="s")		#-0.11, P = 0.30

g1 <- sub_dat[which(sub_dat[,"simplicity_score_b"] < median(sub_dat[,"simplicity_score_b"])),"mut_freq_b"]
g2 <- sub_dat[which(sub_dat[,"simplicity_score_b"] >= median(sub_dat[,"simplicity_score_b"])),"mut_freq_b"]
wilcox.test(g1,g2)		#P = 0.38

g1 <- sub_dat[which(sub_dat[,"simplicity_score_b"] < quantile(0.25,sub_dat[,"simplicity_score_b"])),"mut_freq_b"]
g2 <- sub_dat[which(sub_dat[,"simplicity_score_b"] >= quantile(0.75,sub_dat[,"simplicity_score_b"])),"mut_freq_b"]
wilcox.test(g1,g2)		#P = 0.37

#Multivariate linear regression:
#----------------------
summary(lm(simplicity_score_a ~ mut_freq_a + idh_codel_subtype, data=dat))
summary(lm(simplicity_score_b ~ mut_freq_b + idh_codel_subtype, data=dat))

