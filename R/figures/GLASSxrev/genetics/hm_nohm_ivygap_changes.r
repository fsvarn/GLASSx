###################################################
# Compare microvascular proliferation fraction in samples with hypermutation/cell cycle alterations
# Author: Frederick Varn
# Date: 2022.01.21
# Figure S3B
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(survival)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
SELECT cc.*, ci1.cell_state,
CASE WHEN cc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
mf1.coverage_adj_mut_freq AS mut_freq_a,
mf2.coverage_adj_mut_freq AS mut_freq_b,
mf2.coverage_adj_mut_freq > 10 AS hm,
CASE WHEN cn.cnv_driver_change_a LIKE '%+CDKN2A del%' OR cn.cnv_driver_change_a LIKE '%+CCND2 amp%' THEN TRUE 
WHEN (cn.cnv_driver_change_a NOT LIKE '%+CDKN2A del%'  AND cn.cnv_driver_change_a NOT LIKE '%+CCND2 amp%') OR (cn.cnv_driver_change_a IS NULL AND mf2.coverage_adj_mut_freq IS NOT NULL) THEN FALSE 
ELSE NULL END AS new_cell_cycle,
ca.case_overall_survival_mo, ca.case_vital_status, ca.case_age_diagnosis_years
FROM analysis.rna_silver_set ss
JOIN analysis.tumor_rna_clinical_comparison cc ON cc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_ivygap ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
LEFT JOIN analysis.platinum_set rd ON ss.tumor_barcode_a = rd.rna_barcode_a AND ss.tumor_barcode_b = rd.rna_barcode_b
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = rd.dna_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = rd.dna_barcode_b
LEFT JOIN analysis.driver_status_cnv cn ON cn.tumor_barcode_a = rd.dna_barcode_a AND cn.tumor_barcode_b = rd.dna_barcode_b
JOIN clinical.cases ca ON ca.case_barcode = ss.case_barcode
WHERE cc.idh_codel_subtype IS NOT NULL 
"

dat <- dbGetQuery(con, q)

# Ignore acquired cell cycle alterations from IDH-wild-type tumors as we only care about them in IDH-mutant
dat[which(dat$idh_status == "IDHwt"),"new_cell_cycle"] = 0

dat$alt_status = as.numeric(dat$new_cell_cycle) | as.numeric(dat$hm)

cells <- unique(dat$cell_state)

hyp <- dat %>% 
		filter(!(is.na(alt_status))) %>%
		group_by(alt_status, cell_state, idh_status) %>% 
		summarise(pval = t.test(fraction_b,fraction_a, paired=TRUE)$p.value, eff=mean(fraction_b) - mean(fraction_a), n = n()) %>%
		as.data.frame()

hm_only <- dat %>% filter(alt_status == 1, cell_state == "CTmvp")
case_barcode <- rep(hm_only$case_barcode, 2)
fraction <- c(hm_only$fraction_a, hm_only$fraction_b)
idh_codel_subtype <- rep(hm_only$idh_codel_subtype, 2)
idh_status <- rep(hm_only$idh_status, 2)
timepoint <- rep(c("Initial", "Recurrent"),each = nrow(hm_only))

plot_dat <- data.frame(case_barcode, fraction, idh_codel_subtype, idh_status, timepoint)
plot_dat <- plot_dat %>% 
mutate(idh_status = fct_relevel(idh_status, "IDHwt", "IDHmut"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/hm_mvp.pdf",width=2,height=1.8)
ggplot(data = plot_dat, aes(x = timepoint, y = fraction * 100)) +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size = 1, aes(colour = idh_codel_subtype)) +
scale_colour_manual(values = c("IDHmut-codel" = "#F8766D", "IDHmut-noncodel" = "#00BA38", "IDHwt" = "#619CFF")) +
labs(y = "MVP fraction (%)") +
facet_grid(.~idh_status, scale = "free_y") +
theme_classic() +
theme(
	axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim = c(0,60))
dev.off()

		