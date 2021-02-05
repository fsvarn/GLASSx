library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(survival)
library(topGO)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT * 
FROM analysis.cibersortx_scgp
WHERE cell_state LIKE '%tumor'
"

dat <- dbGetQuery(con, q)

stats = dat %>%
group_by(aliquot_barcode) %>%
summarise(cell_state, fraction = fraction/sum(fraction)) %>%
filter(fraction > 0) %>%
summarise(shannon = -1 * sum(fraction * log(fraction)), 
		  evenness = (-1 * sum(fraction * log(fraction)))/ log(n()))

q <- "SELECT tc.*, cc.case_age_diagnosis_years
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
WHERE idh_codel_subtype = 'IDHwt'
"

surv <- dbGetQuery(con,q)

comb = surv %>%
	   inner_join(stats, c("tumor_barcode_a" = "aliquot_barcode")) %>%
	   inner_join(stats, c("tumor_barcode_b" = "aliquot_barcode"), suffix=c("_a","_b")) 
	   

#Mark surgery (all happened)
comb$event <- 1
n = 4
shannon_res <- comb %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ shannon_a + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ shannon_a + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ shannon_b + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ shannon_b + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(1*n) + 1])

summary(coxph(Surv(surgical_interval, event) ~ shannon_a + case_age_diagnosis_years + received_rt + received_tmz, data = comb))
summary(coxph(Surv(surgical_interval, event) ~ shannon_b + case_age_diagnosis_years + received_rt + received_tmz, data = comb))


summary(coxph(Surv(surgical_interval, event) ~ evenness_a + case_age_diagnosis_years + received_rt + received_tmz, data = comb))
summary(coxph(Surv(surgical_interval, event) ~ evenness_b + case_age_diagnosis_years + received_rt + received_tmz, data = comb))


