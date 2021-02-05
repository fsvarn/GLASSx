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

q <- "SELECT tc.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status,
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
"

dat <- dbGetQuery(con,q)

dat$mes_init <- as.factor(dat$subtype_a == "Mesenchymal")
dat$mes_rec <- dat$subtype_b == "Mesenchymal"
dat$mes_trans <- dat$subtype_a != "Mesenchymal" & dat$subtype_b == "Mesenchymal"
dat$nomes_trans <- dat$subtype_a == "Mesenchymal" & dat$subtype_b != "Mesenchymal"

dat$received_treatment = as.factor(dat$received_treatment)

#Mark surgery (all happened)
dat$event <- 1
dat$case_vital_status <- recode(dat$case_vital_status, 'alive' = 0, 'dead' = 1)

# Examine mesenchymal subtype associations (none significant)
summary(coxph(Surv(surgical_interval, event) ~ mes_init + idh_status + case_age_diagnosis_years + received_treatment, data = dat))
summary(coxph(Surv(surgical_interval, event) ~ mes_rec +  idh_status + case_age_diagnosis_years + received_treatment, data = dat))
summary(coxph(Surv(surgical_interval, event) ~ mes_trans + idh_status + case_age_diagnosis_years + received_treatment, data = dat))
summary(coxph(Surv(surgical_interval, event) ~ nomes_trans + idh_status + case_age_diagnosis_years + received_treatment, data = dat))
    
summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ mes_init + idh_status + case_age_diagnosis_years + received_treatment, data = dat))
summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ mes_rec +  idh_status + case_age_diagnosis_years + received_treatment, data = dat))
summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ mes_trans + idh_status + case_age_diagnosis_years + received_treatment, data = dat))
summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ nomes_trans + idh_status + case_age_diagnosis_years + received_treatment, data = dat))
    

#######################################################
# Examine whether anatomic features associate with survival

q <- "SELECT tc.*, cc.case_age_diagnosis_years,
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_ivygap sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
"

dat <- dbGetQuery(con,q)

#Mark surgery (all happened)
dat$event <- 1
n = 4
ivygap_res <- dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

#######################################################
# Examine whether SCGP cell states associate with survival

q <- "SELECT tc.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status,
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
"

dat <- dbGetQuery(con,q)

#Mark surgery (all happened)
dat$event <- 1
dat$case_vital_status <- recode(dat$case_vital_status, 'alive' = 0, 'dead' = 1)

n = 5
scgp_res <- dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

scgp_os <- dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + idh_codel_subtype + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

#-------------------
# Tumor only

tumor_dat <- dat %>%
	filter(base::grepl("_tumor", cell_state)) %>%
	group_by(tumor_pair_barcode) %>%
	summarise(cell_state,
	fraction_a = fraction_a/sum(fraction_a), 
	fraction_b = fraction_b/sum(fraction_b),
	diff = (fraction_b/sum(fraction_b)) - (fraction_a/sum(fraction_a)),
	surgical_interval, event, case_overall_survival_mo, case_vital_status, case_age_diagnosis_years, idh_status, received_treatment) 

tumor_int_res <- tumor_dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

tumor_os_res <- tumor_dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

#-------------------
# Non-tumor only

non_dat <- dat %>%
	filter(!base::grepl("_tumor", cell_state)) %>%
	group_by(tumor_pair_barcode) %>%
	summarise(cell_state,
	fraction_a = fraction_a/sum(fraction_a), 
	fraction_b = fraction_b/sum(fraction_b),
	diff = (fraction_b/sum(fraction_b)) - (fraction_a/sum(fraction_a)),
	surgical_interval, event, case_overall_survival_mo, case_vital_status, case_age_diagnosis_years, idh_status, received_treatment) 

non_int_res <- non_dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment, data = non_dat %>% filter(cell_state == "myeloid")))

non_os_res <- non_dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

#######################################################
# Adjust for anatomic feature proportions 

q <- "SELECT tc.*, cc.case_age_diagnosis_years,
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS diff,
iv1.cell_state AS region, iv1.fraction AS region_fraction_a, iv2.fraction AS region_fraction_b, iv2.fraction - iv1.fraction AS region_diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN analysis.cibersortx_ivygap iv1 ON iv1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap iv2 ON iv2.aliquot_barcode = ss.tumor_barcode_b AND iv2.cell_state = iv1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
WHERE sc1.cell_state = 'myeloid'
"

dat <- dbGetQuery(con,q)

#Mark surgery (all happened)
dat$event <- 1
n = 5
mvp_adj <- dat %>%
filter(region == 'CTmvp') %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment + region_diff))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment +region_diff))$coefficient[(1*n) + 1])

pan_adj <- dat %>%
filter(region == 'CTpan') %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment + region_diff))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment +region_diff))$coefficient[(1*n) + 1])


ct_adj <- dat %>%
filter(region == 'CT') %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment + region_fraction_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment + region_diff))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + idh_status + case_age_diagnosis_years + received_treatment +region_diff))$coefficient[(1*n) + 1])
