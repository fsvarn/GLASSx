###################################################
# Test how macrophage and microglia infiltration levels differ between mesenchymal and non mesenchymal samples
# Updated: 2020.04.10
# Author: Frederick Varn
##################################################
library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(reshape)


##################################################
rm(list=ls())

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
WITH subtype_rank AS
(
	SELECT *,
	RANK() OVER (PARTITION BY aliquot_barcode ORDER BY p_value ASC) AS p_rank
	FROM analysis.transcriptional_subtype
),
top_rank AS
(
	SELECT *
	FROM subtype_rank
	WHERE p_rank = 1
),
agg AS
(
	SELECT aliquot_barcode, 
	string_agg(signature_name,',') AS subtype 
	FROM top_rank
	GROUP BY aliquot_barcode
)
SELECT ps.*, 
tm1.signature_name, tm1.enrichment_score AS es_a, tm2.enrichment_score AS es_b, 
ts1.subtype AS subtype_a, ts2.subtype AS subtype_b,
cs.idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.muller_tam_score tm1 ON tm1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.muller_tam_score tm2 ON tm2.aliquot_barcode = ps.rna_barcode_b AND tm1.signature_name = tm2.signature_name
JOIN agg ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
JOIN agg ts2 ON ts2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE idh_codel_subtype LIKE 'IDHwt'
"

dat <- dbGetQuery(con,q)


dat[which(dat[,"subtype_a"]=="Mesenchymal,Classical"),"subtype_a"]  <- "Classical,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Classical"),"subtype_b"]  <- "Classical,Mesenchymal"

dat[which(dat[,"subtype_a"]=="Mesenchymal,Proneural"),"subtype_a"]  <- "Proneural,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Proneural"),"subtype_b"]  <- "Proneural,Mesenchymal"

#Convert for the purposes of plotting
dat[which(dat[,"subtype_a"]=="Proneural,Classical"),"subtype_a"]  <- "Proneural"
dat[which(dat[,"subtype_b"]=="Proneural,Classical"),"subtype_b"]  <- "Proneural"

dat[which(dat[,"subtype_a"]=="Proneural,Mesenchymal"),"subtype_a"]  <- "Mesenchymal"
dat[which(dat[,"subtype_b"]=="Proneural,Mesenchymal"),"subtype_b"]  <- "Mesenchymal"

dat[which(dat[,"subtype_a"]=="Classical,Mesenchymal"),"subtype_a"]  <- "Classical"
dat[which(dat[,"subtype_b"]=="Classical,Mesenchymal"),"subtype_b"]  <- "Classical"

##################################################

# Test how initial mesenchymal samples differ in macrophage and microglia levels

g1 <- dat[which(dat[,"subtype_a"] == "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_a"]
g2 <- dat[which(dat[,"subtype_a"] != "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_a"]
wilcox.test(g1,g2) #P = 2.7e-6

g1 <- dat[which(dat[,"subtype_a"] == "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_a"]
g2 <- dat[which(dat[,"subtype_a"] != "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_a"]
wilcox.test(g1,g2) #P = 0.41


# Test how recurrent mesenchymal samples differ in macrophage and microglia levels

g1 <- dat[which(dat[,"subtype_b"] == "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_b"]
g2 <- dat[which(dat[,"subtype_b"] != "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_b"]
wilcox.test(g1,g2) #P = 2.8e-8

g1 <- dat[which(dat[,"subtype_b"] == "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_b"]
g2 <- dat[which(dat[,"subtype_b"] != "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_b"]
wilcox.test(g1,g2) #P = 0.41

##################################################

# Test how macrophage levels change in samples switching TO mesenchymal

g1 <- dat[which(dat[,"subtype_a"] != "Mesenchymal" & dat[,"subtype_b"] == "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_a"]
g2 <- dat[which(dat[,"subtype_a"] != "Mesenchymal" & dat[,"subtype_b"] == "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_b"]
wilcox.test(g1, g2, paired=TRUE) #P = 0.01

# Test how macrophage levels change in samples switching AWAY from mesenchymal

g1 <- dat[which(dat[,"subtype_a"] == "Mesenchymal" & dat[,"subtype_b"] != "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_a"]
g2 <- dat[which(dat[,"subtype_a"] == "Mesenchymal" & dat[,"subtype_b"] != "Mesenchymal" & dat[,"signature_name"] == "Macrophages"),"es_b"]
wilcox.test(g1, g2, paired=TRUE) #P = 5e-4


# Test how microglia levels change in samples switching TO mesenchymal

g1 <- dat[which(dat[,"subtype_a"] != "Mesenchymal" & dat[,"subtype_b"] == "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_a"]
g2 <- dat[which(dat[,"subtype_a"] != "Mesenchymal" & dat[,"subtype_b"] == "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_b"]
wilcox.test(g1, g2, paired=TRUE) #P = 0.53

# Test how microglia levels change in samples switching AWAY from mesenchymal

g1 <- dat[which(dat[,"subtype_a"] == "Mesenchymal" & dat[,"subtype_b"] != "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_a"]
g2 <- dat[which(dat[,"subtype_a"] == "Mesenchymal" & dat[,"subtype_b"] != "Mesenchymal" & dat[,"signature_name"] == "Microglia"),"es_b"]
wilcox.test(g1, g2, paired=TRUE) #P = 0.1

