library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q1 <- "SELECT ic.case_barcode, tumor_barcode_a, tumor_barcode_b, change, 
s1.surgery_extent_of_resection AS extent1, s2.surgery_extent_of_resection AS extent2,
CASE
WHEN s1.surgery_location::text = s2.surgery_location::text AND (s1.surgery_laterality::text = s2.surgery_laterality::text OR s1.surgery_laterality IS NULL AND s2.surgery_laterality IS NULL) THEN 'Local'::text
WHEN s1.surgery_location::text <> s2.surgery_location::text OR s1.surgery_laterality::text <> s2.surgery_laterality::text THEN 'Distal'::text
ELSE NULL::text
END AS recurrence_location,
s1.idh_codel_subtype
FROM analysis.immune_cluster_change ic
JOIN analysis.rna_silver_set ss ON ic.case_barcode = ss.case_barcode
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.surgeries s1 ON s1.sample_barcode = al2.sample_barcode
JOIN clinical.surgeries s2 ON s2.sample_barcode = al2.sample_barcode
ORDER BY 1"

dat <- dbGetQuery(con,q1)

sub_dat <- dat[which(dat[,"recurrence_location"] == "Local"),]

g1 <- nrow(dat[which(dat[,"extent1"] == "Subtotal" & dat[,"change"] == "None"),])
g2 <- nrow(dat[which(dat[,"extent1"] == "Subtotal" & dat[,"change"] != "None"),])
g3 <- nrow(dat[which(dat[,"extent1"] != "Subtotal" & dat[,"change"] == "None"),])
g4 <- nrow(dat[which(dat[,"extent1"] != "Subtotal" & dat[,"change"] != "None"),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2)

fisher.test(ct)

g1 <- nrow(dat[which(dat[,"extent1"] == "Subtotal" & dat[,"change"] == "Increase"),])
g2 <- nrow(dat[which(dat[,"extent1"] == "Subtotal" & dat[,"change"] != "Increase"),])
g3 <- nrow(dat[which(dat[,"extent1"] != "Subtotal" & dat[,"change"] == "Increase"),])
g4 <- nrow(dat[which(dat[,"extent1"] != "Subtotal" & dat[,"change"] != "Increase"),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2)

fisher.test(ct)

g1 <- nrow(dat[which(dat[,"extent1"] == "Subtotal" & dat[,"change"] == "Decrease"),])
g2 <- nrow(dat[which(dat[,"extent1"] == "Subtotal" & dat[,"change"] != "Decrease"),])
g3 <- nrow(dat[which(dat[,"extent1"] != "Subtotal" & dat[,"change"] == "Decrease"),])
g4 <- nrow(dat[which(dat[,"extent1"] != "Subtotal" & dat[,"change"] != "Decrease"),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2)

fisher.test(ct)