library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggbeeswarm)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Radiotherapy

#Read in data
q <- "
SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b,
cs1.treatment_radiotherapy
FROM analysis.rna_silver_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.tumor_barcode_b AND im2.signature_name = im1.signature_name
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
WHERE cs1.treatment_radiotherapy IS true
ORDER BY 1, 2, 6
"

dat <- dbGetQuery(con,q)

#See if there are significant changes across samples over time

cells <- unique(dat[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i] & dat[,"subtype_a"]=="IDHwt"),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}
time_res <- data.frame(cells,eff,p.value)

#                    cells           eff    p.value
# 1                B.cells -0.0008587021 0.38086771
# 2             CD4.mature -0.0222623790 0.02194575	***
# 3           CD8.effector  0.0127815552 0.42946626
# 4  CD8.effector.NK.cells  0.0100555999 0.86878354
# 5              Dendritic -0.0248952236 0.30435600
# 6            Macrophages  0.0069865275 0.76289862
# 7         Macrophages.M1 -0.0096128605 0.88010328
# 8         Macrophages.M2  0.0010396352 0.90281413
# 9               NK.cells  0.0133690382 0.74109194
# 10                 T.reg -0.0025225089 0.34306512

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i] & dat[,"subtype_a"]!="IDHwt"),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}
time_res <- data.frame(cells,eff,p.value)

#                    cells          eff    p.value
# 1                B.cells  0.027164136 0.74665833
# 2             CD4.mature  0.028906497 0.37782288
# 3           CD8.effector  0.007395984 0.43067932
# 4  CD8.effector.NK.cells  0.016129885 0.05688477
# 5              Dendritic  0.001184473 0.57905579
# 6            Macrophages -0.019272284 0.15937805
# 7         Macrophages.M1 -0.015036040 0.43067932
# 8         Macrophages.M2 -0.018527086 0.10888672
# 9               NK.cells  0.034586018 0.01499939
# 10                 T.reg  0.007727198 0.40376282

#----------------------------------------------------------------

#Chemotherapy

#Read in data
q <- "
SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b,
cs1.treatment_tmz
FROM analysis.rna_silver_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.tumor_barcode_b AND im2.signature_name = im1.signature_name
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
WHERE cs1.treatment_tmz IS true
ORDER BY 1, 2, 6
"

dat <- dbGetQuery(con,q)

#See if there are significant changes across samples over time

cells <- unique(dat[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i] & dat[,"subtype_a"]=="IDHwt"),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}
time_res <- data.frame(cells,eff,p.value)

#                    cells          eff    p.value
# 1                B.cells  0.002806734 0.08975723
# 2             CD4.mature -0.015695569 0.02974777
# 3           CD8.effector  0.011719480 0.45314340
# 4  CD8.effector.NK.cells  0.004036001 0.93295079
# 5              Dendritic -0.023511266 0.51438195
# 6            Macrophages -0.001700453 0.48762186
# 7         Macrophages.M1 -0.014405783 0.88848502
# 8         Macrophages.M2 -0.002358022 0.72592323
# 9               NK.cells  0.008047323 0.58447238
# 10                 T.reg -0.011044058 0.18747629

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i] & dat[,"subtype_a"]!="IDHwt"),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}
time_res <- data.frame(cells,eff,p.value)

#                    cells         eff    p.value
# 1                B.cells  0.02101366 0.55664062
# 2             CD4.mature -0.01407242 0.43164062
# 3           CD8.effector  0.04520202 0.19335938
# 4  CD8.effector.NK.cells  0.01347239 0.06445312
# 5              Dendritic  0.03898694 0.76953125
# 6            Macrophages -0.01167744 0.32226562
# 7         Macrophages.M1  0.01362770 0.84570312
# 8         Macrophages.M2  0.01239615 0.49218750
# 9               NK.cells  0.03884037 0.13085938
# 10                 T.reg  0.04687770 0.03710938