library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Query for grade increase
q <- "SELECT ps.*, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs.idh_codel_subtype
FROM analysis.tumor_clinical_comparison ps
JOIN analysis.rna_dna_pairs rd ON rd.dna_pair_barcode = ps.tumor_pair_barcode
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE grade_change = 'Grade up' AND idh_codel_subtype LIKE 'IDHmut%'
ORDER BY 2, 1, 16"

dat <- dbGetQuery(con, q)


cells <- unique(dat[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}

sig <- rep("",length(p.value))
sig[which(p.value < 0.05)] <- "*"
sig[which(p.value < 0.01)] <- "**"
time_res <- data.frame(cells,eff,p.value, sig)

#                    cells          eff     p.value sig
# 1                B.cells  0.018932614 0.241210938    
# 2             CD4.mature  0.028611303 0.426269531    
# 3           CD8.effector  0.043551361 0.357543945    
# 4  CD8.effector.NK.cells  0.040157104 0.067626953    
# 5              Dendritic -0.021196160 0.241210938    
# 6            Macrophages -0.006956297 0.357543945    
# 7         Macrophages.M1 -0.007395265 0.714843750    
# 8         Macrophages.M2 -0.015133314 0.295776367    
# 9               NK.cells  0.043426059 0.001708984  **
# 10                 T.reg -0.000926465 1.000000000   


#Query for grade stable
q <- "SELECT ps.*, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs.idh_codel_subtype
FROM analysis.tumor_clinical_comparison ps
JOIN analysis.rna_dna_pairs rd ON rd.dna_pair_barcode = ps.tumor_pair_barcode
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE grade_change = 'Grade stable'
ORDER BY 2, 1, 16"

dat <- dbGetQuery(con, q)


cells <- unique(dat[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}

sig <- rep("",length(p.value))
sig[which(p.value < 0.05)] <- "*"
sig[which(p.value < 0.01)] <- "**"
time_res <- data.frame(cells,eff,p.value,sig)

#                    cells          eff    p.value sig
# 1                B.cells  0.009590539 0.11489563    
# 2             CD4.mature -0.028543629 0.03105225   *
# 3           CD8.effector  0.014181102 0.23850949    
# 4  CD8.effector.NK.cells  0.013230272 0.47052665    
# 5              Dendritic -0.035416718 0.57493960    
# 6            Macrophages -0.004560506 0.23850949    
# 7         Macrophages.M1 -0.020332797 0.69856104    
# 8         Macrophages.M2 -0.006347169 0.52967000    
# 9               NK.cells  0.010901454 0.24019957    
# 10                 T.reg -0.008478395 0.58072445


#######################################################

q <- "SELECT ps.*, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs.idh_codel_subtype
FROM analysis.tumor_clinical_comparison ps
JOIN analysis.rna_dna_pairs rd ON rd.dna_pair_barcode = ps.tumor_pair_barcode
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE received_alk IS true AND idh_codel_subtype ='IDHwt'
ORDER BY 2, 1, 16"


dat <- dbGetQuery(con, q)

cells <- unique(dat[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}

sig <- rep("",length(p.value))
sig[which(p.value < 0.05)] <- "*"
sig[which(p.value < 0.01)] <- "**"
time_res <- data.frame(cells,eff,p.value,sig)

#                    cells           eff   p.value
# 1                B.cells  0.0049032093 0.2729101
# 2             CD4.mature -0.0284350572 0.0442438
# 3           CD8.effector  0.0092045707 0.4729129
# 4  CD8.effector.NK.cells -0.0002294454 0.4806453
# 5              Dendritic -0.0237492488 0.1963760
# 6            Macrophages -0.0094369017 0.9188035
# 7         Macrophages.M1 -0.0054385640 0.8076878
# 8         Macrophages.M2 -0.0141705342 0.5957590
# 9               NK.cells  0.0163263808 0.1241865
# 10                 T.reg -0.0220467205 0.5176235

#######################################################

q <- "
SELECT ps.*, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs.idh_codel_subtype
FROM analysis.tumor_clinical_comparison ps
JOIN analysis.rna_dna_pairs rd ON rd.dna_pair_barcode = ps.tumor_pair_barcode
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE received_rt IS true AND idh_codel_subtype ='IDHwt'
ORDER BY 2, 1, 16
"


dat <- dbGetQuery(con, q)

cells <- unique(dat[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}

sig <- rep("",length(p.value))
sig[which(p.value < 0.05)] <- "*"
sig[which(p.value < 0.01)] <- "**"
time_res <- data.frame(cells,eff,p.value,sig)

