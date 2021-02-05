library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(survival)
library(topGO)
library(gridExtra)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT tc.*, gr1.grade AS grade_a, gr2.grade AS grade_b, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status,
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
ts1.signature_name, ts1.enrichment_score/1000 AS score_a, ts2.enrichment_score/1000 AS score_b,
(ts2.enrichment_score/1000) - (ts1.enrichment_score/1000) AS score_diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
JOIN clinical.surgeries gr1 ON gr1.sample_barcode = substring(ss.tumor_barcode_a, 1,15)
JOIN clinical.surgeries gr2 ON gr2.sample_barcode = substring(ss.tumor_barcode_b, 1,15)
"

dat <- dbGetQuery(con,q)

#Mark surgery (all happened)
dat$event <- 1
dat$case_vital_status <- recode(dat$case_vital_status, 'alive' = 0, 'dead' = 1)

# IDHwt (all were treated so no need to adjust)
n = 3
idhwt_subtype_si <- dat %>%
filter(idh_status == "IDHwt") %>%
group_by(signature_name) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ score_a + case_age_diagnosis_years + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ score_a + case_age_diagnosis_years + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ score_b + case_age_diagnosis_years + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ score_b + case_age_diagnosis_years + grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ score_diff + case_age_diagnosis_years + grade_change))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ score_diff + case_age_diagnosis_years + grade_change))$coefficient[(1*n) + 1]) %>%
pivot_longer(-signature_name) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="subtype", analysis = "si", subtype = "IDHwt")

# IDHmut
n = 5
idhmut_subtype_si <- dat %>%
filter(idh_status == "IDHmut") %>%
group_by(signature_name) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ score_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ score_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ score_b + case_age_diagnosis_years + received_treatment + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ score_b + case_age_diagnosis_years + received_treatment + grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ score_diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(4*(n-1)) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ score_diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(1*(n-1)) + 1]) %>%
pivot_longer(-signature_name) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="subtype", analysis = "si", subtype = "IDHmut")

# IDHwt (all were treated so no need to adjust)
n = 3
idhwt_subtype_os <- dat %>%
filter(idh_status == "IDHwt") %>%
group_by(signature_name) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_a + case_age_diagnosis_years + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_a + case_age_diagnosis_years + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_b + case_age_diagnosis_years + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_b + case_age_diagnosis_years + grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_diff + case_age_diagnosis_years + grade_change))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_diff + case_age_diagnosis_years + grade_change))$coefficient[(1*n) + 1]) %>%
pivot_longer(-signature_name) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="subtype", analysis = "os", subtype = "IDHwt")


# IDHmut
n = 5
idhmut_subtype_os <- dat %>%
filter(idh_status == "IDHmut") %>%
group_by(signature_name) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_b + case_age_diagnosis_years + received_treatment + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_b + case_age_diagnosis_years + received_treatment+ grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(4*(n-1)) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ score_diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(1*(n-1)) + 1]) %>%
pivot_longer(-signature_name) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="subtype", analysis = "os", subtype = "IDHmut")


plot_res <- rbind(idhwt_subtype_si, idhmut_subtype_si, idhwt_subtype_os, idhmut_subtype_os)
colnames(plot_res)[1] <- "cell_state"
#######################################################
# Examine whether anatomic features associate with survival

q <- "SELECT tc.*, gr1.grade AS grade_a, gr2.grade AS grade_b, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status,
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction * 100 AS fraction_a, sc2.fraction *100 AS fraction_b, (sc2.fraction - sc1.fraction) * 100 AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_ivygap sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
JOIN clinical.surgeries gr1 ON gr1.sample_barcode = substring(ss.tumor_barcode_a, 1,15)
JOIN clinical.surgeries gr2 ON gr2.sample_barcode = substring(ss.tumor_barcode_b, 1,15)
"

dat <- dbGetQuery(con,q)

#Mark surgery (all happened)
dat$event <- 1
dat$case_vital_status <- recode(dat$case_vital_status, 'alive' = 0, 'dead' = 1)

# IDHwt (all were treated so no need to adjust)
n = 3
idhwt_region_si <- dat %>%
filter(idh_status == "IDHwt") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(1*n) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="region", analysis = "si", subtype = "IDHwt")

# IDHmut
n = 5
idhmut_region_si <- dat %>%
filter(idh_status == "IDHmut") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + received_treatment + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + received_treatment+ grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(4*(n-1)) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(1*(n-1)) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="region", analysis = "si", subtype = "IDHmut")

# IDHwt (all were treated so no need to adjust)
n = 3
idhwt_region_os <- dat %>%
filter(idh_status == "IDHwt") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(1*n) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="region", analysis = "os", subtype = "IDHwt")


# IDHmut
n = 5
idhmut_region_os <- dat %>%
filter(idh_status == "IDHmut") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + received_treatment + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + received_treatment+ grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(4*(n-1)) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(1*(n-1)) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="region", analysis = "os", subtype = "IDHmut")


plot_res <- rbind(plot_res, idhwt_region_si, idhmut_region_si, idhwt_region_os, idhmut_region_os)

#######################################################
# Examine whether SCGP cell states associate with survival

q <- "SELECT tc.*,  gr1.grade AS grade_a, gr2.grade AS grade_b, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status,
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction * 100 AS fraction_a, sc2.fraction * 100 AS fraction_b, (sc2.fraction - sc1.fraction)*100 AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
JOIN clinical.surgeries gr1 ON gr1.sample_barcode = substring(ss.tumor_barcode_a, 1,15)
JOIN clinical.surgeries gr2 ON gr2.sample_barcode = substring(ss.tumor_barcode_b, 1,15)
"

dat <- dbGetQuery(con,q)

#Mark surgery (all happened)
dat$event <- 1
dat$case_vital_status <- recode(dat$case_vital_status, 'alive' = 0, 'dead' = 1)

# IDHwt (all were treated so no need to adjust)
n = 3
idhwt_cell_si <- dat %>%
filter(idh_status == "IDHwt") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(1*n) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="cell_state", analysis = "si", subtype = "IDHwt")

# IDHmut
n = 5
idhmut_cell_si <- dat %>%
filter(idh_status == "IDHmut") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + received_treatment + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + received_treatment+ grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(4*(n-1)) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(1*(n-1)) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="cell_state", analysis = "si", subtype = "IDHmut")

# IDHwt (all were treated so no need to adjust)
n = 3
idhwt_cell_os <- dat %>%
filter(idh_status == "IDHwt") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + grade_change))$coefficient[(1*n) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="cell_state", analysis = "os", subtype = "IDHwt")


# IDHmut
n = 5
idhmut_cell_os <- dat %>%
filter(idh_status == "IDHmut") %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_a + case_age_diagnosis_years + received_treatment + grade_a))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + received_treatment + grade_b))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ fraction_b + case_age_diagnosis_years + received_treatment+ grade_b))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(4*(n-1)) + 1],
		  diff_hr = summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ diff + case_age_diagnosis_years + received_treatment + grade_change))$coefficient[(1*(n-1)) + 1]) %>%
pivot_longer(-cell_state) %>% 
separate(name,c("timepoint","stat")) %>% 
pivot_wider(names_from = stat, values_from = value) %>%
mutate(type="cell_state", analysis = "os", subtype = "IDHmut")


plot_res <- rbind(plot_res, idhwt_cell_si, idhmut_cell_si, idhwt_cell_os, idhmut_cell_os)
eff <- rep("none",nrow(plot_res))
eff[which(plot_res$pval < 0.05 & plot_res$hr < 1)] <- "prolonged"
eff[which(plot_res$pval < 0.05 & plot_res$hr > 1)] <- "short"
plot_res$eff <- eff

plot_res <- plot_res %>%
mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		mutate(cell_state = fct_relevel(cell_state, "Prolif. stem-like", "Stem-like","Diff.-like",
										"Fibroblast", "Pericyte","Endothelial", "Oligodendrocyte",
										"Myeloid", "Dendritic cell", "T cell", "Granulocyte","B cell")) %>%
		mutate(type = recode(type, "cell_state" = "Cell state", "region" = "Region", "subtype" = "Subtype")) %>%
		mutate(type = fct_relevel(type, "Subtype", "Region", "Cell state")) %>%
		mutate(timepoint = recode(timepoint, "init" = "Initial", "rec" = "Recurrent")) %>%
		mutate(subtype = fct_relevel(subtype, "IDHwt", "IDHmut")) 


os_plot <- ggplot(plot_res %>% filter(timepoint!="diff",analysis=="os"), aes(x=timepoint, y = cell_state, fill=eff)) +
geom_tile() + 
scale_fill_manual(values=c("white","#CD4F39")) +
theme_classic() +
facet_grid(type~subtype, scales = "free_y", space = "free_y") +
theme(axis.text.x = element_text(size=7, angle = 45, hjust=1),
	axis.text.y = element_text(size=7),
	axis.title = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_blank(),
	legend.position = "none")
	
si_plot <- ggplot(plot_res %>% filter(timepoint!="diff",analysis=="si"), aes(x=timepoint, y = cell_state, fill=eff)) +
geom_tile() + 
scale_fill_manual(values=c("white","#27408B","#CD4F39")) +
theme_classic() +
facet_grid(type~subtype, scales = "free_y", space = "free_y") +
theme(axis.text.x = element_text(size=7, angle = 45, hjust=1),
	axis.text.y = element_text(size=7),
	axis.title = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_blank(),
	legend.position = "none")

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_surv_heatmap.pdf",width=4,height=4.65)
grid.arrange(os_plot, si_plot, ncol=2)
dev.off()