library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(survival)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
im1.aaseqcdr3 AS tcr_amino_acid,
im1.clonecount AS clonecount_a, 
im2.clonecount AS clonecount_b, 
im1.clonefraction AS clonefraction_a, 
im2.clonefraction AS clonefraction_b, 
sc1.fraction AS t_cell_a,
sc2.fraction AS t_cell_b,
cs.idh_codel_subtype,
ts1.signature_name AS subtype_a,
ts2.signature_name AS subtype_b
FROM analysis.rna_silver_set ps
JOIN analysis.mixcr_tcr im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.mixcr_tcr im2 ON im2.aliquot_barcode = ps.tumor_barcode_b AND im1.aaseqcdr3 = im2.aaseqcdr3
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ps.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE sc1.cell_state = 't_cell'"

dat <- dbGetQuery(con, q)

#Need to filter for 

q <- "SELECT cc.*, cs.idh_codel_subtype
FROM clinical.cases cc
JOIN clinical.subtypes cs ON cc.case_barcode = cs.case_barcode
JOIN analysis.rna_silver_set rs ON cc.case_barcode = rs.case_barcode
WHERE idh_codel_subtype = 'IDHwt' 
ORDER BY 2"

clin <- dbGetQuery(con, q)

shared <- as.numeric(clin[,"case_barcode"] %in% dat[,"case_barcode"])
clin[,"shared"] <- shared
clin[,"vital_status"] <- ifelse(clin[,"case_vital_status"] == "dead",1,0)

shared_cox <- coxph(Surv(case_overall_survival_mo, vital_status) ~ shared, data=clin)
summary(shared_cox)
# Call:
# coxph(formula = Surv(case_overall_survival_mo, vital_status) ~ 
#     shared, data = clin)
# 
#   n= 98, number of events= 83 
#    (1 observation deleted due to missingness)
# 
#           coef exp(coef) se(coef)      z Pr(>|z|)
# shared -0.5136    0.5983   0.3737 -1.374    0.169
# 
#        exp(coef) exp(-coef) lower .95 upper .95
# shared    0.5983      1.671    0.2877     1.245
# 
# Concordance= 0.517  (se = 0.025 )
# Likelihood ratio test= 2.16  on 1 df,   p=0.1
# Wald test            = 1.89  on 1 df,   p=0.2
# Score (logrank) test = 1.93  on 1 df,   p=0.2


#Do shared patients have a longer surgical interval

q <- "SELECT tc.*, cs.idh_codel_subtype
FROM analysis.tumor_tcr_comparison tc
JOIN clinical.subtypes cs ON cs.case_barcode = tc.case_barcode
WHERE union_ab > 0 AND idh_codel_subtype = 'IDHwt'"

clin <- dbGetQuery(con, q)

clin[,"shared"] <- as.numeric(clin[,"intersection_ab"] > 0)
clin[,"recur_status"] <- rep(1,nrow(clin))
shared_cox <- coxph(Surv(surgical_interval_mo, recur_status) ~ shared, data=clin)
summary(shared_cox)

# Call:
# coxph(formula = Surv(surgical_interval_mo, recur_status) ~ shared, 
#     data = clin)
# 
#   n= 101, number of events= 101 
# 
#           coef exp(coef) se(coef)     z Pr(>|z|)
# shared 0.05869   1.06045  0.28992 0.202     0.84
# 
#        exp(coef) exp(-coef) lower .95 upper .95
# shared      1.06      0.943    0.6008     1.872
# 
# Concordance= 0.501  (se = 0.021 )
# Rsquare= 0   (max possible= 0.999 )
# Likelihood ratio test= 0.04  on 1 df,   p=0.8407
# Wald test            = 0.04  on 1 df,   p=0.8396
# Score (logrank) test = 0.04  on 1 df,   p=0.8395