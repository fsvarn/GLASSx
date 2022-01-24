##################################################
# Compare CIBERSORTx profiles in samples with acquired hypermutation, test survival associations with proliferating neoplastic cells
# Author: Frederick Varn
# Date: 2022.01.06
# Figure 3C-F
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(survival)
library(survminer)

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
CASE WHEN cn.cnv_driver_change_b LIKE '%-EGFR amp%' THEN TRUE 
WHEN (cn.cnv_driver_change_b NOT LIKE '%-EGFR amp%') OR (cn.cnv_driver_change_b IS NULL AND mf2.coverage_adj_mut_freq IS NOT NULL) THEN FALSE 
ELSE NULL END AS lost_egfr,
ca.case_overall_survival_mo, ca.case_vital_status, ca.case_age_diagnosis_years
FROM analysis.rna_silver_set ss
JOIN analysis.tumor_rna_clinical_comparison cc ON cc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
LEFT JOIN analysis.platinum_set rd ON ss.tumor_barcode_a = rd.rna_barcode_a AND ss.tumor_barcode_b = rd.rna_barcode_b
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = rd.dna_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = rd.dna_barcode_b
LEFT JOIN analysis.driver_status_cnv cn ON cn.tumor_barcode_a = rd.dna_barcode_a AND cn.tumor_barcode_b = rd.dna_barcode_b
JOIN clinical.cases ca ON ca.case_barcode = ss.case_barcode
WHERE cc.idh_codel_subtype IS NOT NULL 
"
dat <- dbGetQuery(con, q)

cells <- unique(dat$cell_state)

# P-values

hyp <- dat %>% 
		filter(!(is.na(hm))) %>%
		group_by(hm, cell_state, idh_status) %>% 
		summarise(pval = t.test(fraction_b,fraction_a, paired=TRUE)$p.value, eff=mean(fraction_b) - mean(fraction_a), n = n()) %>%
		as.data.frame()


plot_dat <- dat %>% 
			pivot_longer(cols = starts_with("fraction"),
						 names_to = "timepoint",
						 names_prefix = "fraction_",
						 values_to = "fraction") %>%
			mutate(timepoint = recode(timepoint, "a" = "Init.", "b" = "Rec.")) %>%
			mutate(fraction = fraction * 100) %>%
			mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
					"differentiated_tumor" = "Differentiated tumor", "endothelial" = "Endothelial",
					"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
					"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
					"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem cell tumor",
					"stemcell_tumor" = "Stem cell tumor","t_cell" = "T cell")) %>%
			mutate(cell_state = as_factor(cell_state)) %>%
			mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
											"Oligodendrocyte", 
											"Endothelial", "Pericyte",
											"Fibroblast", 
											"Differentiated tumor", "Stem cell tumor", "Prolif. stem cell tumor")) %>%
			mutate(hm = recode(hm, "1" = "Yes", "0" = "No")) %>%
			mutate(hm = as_factor(hm)) %>%
			mutate(hm = fct_relevel(hm, "Yes", "No"))

# Plot the effect of hypermutators on IDHwt and IDHmut tumors
# IDHwt
#----------------------------------
# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Prolif. stem cell tumor", idh_status == "IDHwt", !is.na(hm)),
	   aes(x = timepoint, y = fraction, colour= idh_codel_subtype)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~hm) +
scale_colour_manual(values="#619CFF") +
labs(y = "Proportion (%)") +
theme_classic() +
theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim=c(0,100))

plot_bars <- plot_dat %>%
			 group_by(cell_state, timepoint, idh_status, hm) %>%
			 summarise(mean = mean(fraction))

bars <- ggplot(plot_bars %>% filter(idh_status == "IDHwt", !is.na(hm)), aes(x=timepoint, y = mean, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
facet_grid(.~hm) +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Prolif. stem cell tumor" = "#a50f15")) +
labs(y = "Proportion (%)") +
theme(axis.text = element_text(size=7),
	axis.title.x = element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_blank(),
	strip.background = element_blank(),
	legend.position = "none") +
    coord_cartesian(ylim=c(0,100))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_hypermutator_prolif_stem_cell_v2.pdf",width=1.45,height=2)
grid.arrange(boxes, bars, nrow = 2, ncol = 1, heights = c(1.25, 1))
dev.off()


# IDHmut
#----------------------------------
# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Prolif. stem cell tumor", idh_status == "IDHmut", !is.na(hm)),
	   aes(x = timepoint, y = fraction, colour= idh_codel_subtype)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~hm) +
scale_colour_manual(values=c("#F8766D","#00BA38")) +
labs(y = "Proportion (%)") +
theme_classic() +
theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim=c(0,100))

plot_bars <- plot_dat %>%
			 group_by(cell_state, timepoint, idh_status, hm) %>%
			 summarise(mean = mean(fraction))

bars <- ggplot(plot_bars %>% filter(idh_status == "IDHmut", !is.na(hm)), aes(x=timepoint, y = mean, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
facet_grid(.~hm) +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Prolif. stem cell tumor" = "#a50f15")) +
labs(y = "Proportion (%)") +
theme(axis.text = element_text(size=7),
	axis.title.x = element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_blank(),
	strip.background = element_blank(),
	legend.position = "none") +
    coord_cartesian(ylim=c(0,100))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_hypermutator_prolif_stem_cell_v2.pdf",width=1.45,height=2)
grid.arrange(boxes, bars, nrow = 2, ncol = 1, heights = c(1.25, 1))
dev.off()


# Waterfall plot
sub_dat <- dat %>% filter(cell_state == "prolif_stemcell_tumor", idh_status == "IDHmut")
sub_dat$diff <- (sub_dat[,"fraction_b"] - sub_dat[,"fraction_a"]) * 100

sub_dat <- sub_dat[order(sub_dat$diff, decreasing=TRUE),]
sub_dat <- sub_dat %>% 
			mutate(tumor_pair_barcode = as_factor(tumor_pair_barcode)) %>%
			mutate(tumor_pair_barcode = fct_relevel(tumor_pair_barcode, sub_dat$tumor_pair_barcode))
sub_dat$fill <- sub_dat$diff > 0

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_prolif_stem_cell_waterfall_v2.pdf",width=2.4,height=1.1)
ggplot(sub_dat, aes(x=tumor_pair_barcode, y = diff, fill = factor(fill))) +
geom_bar(stat="identity") +
theme_classic() +
scale_fill_manual(values=c("#27408B","#CD4F39")) +
labs(title = "Prolif. stem cell", y = "Change (%)") +
theme(plot.title = element_text(size=7, hjust = 0.5),
	axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),
	axis.ticks.x = element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_blank(),
	strip.background = element_blank(),
	legend.position = "none")
dev.off()

plot_dat <- sub_dat %>% 
			pivot_longer(cols = c("idh_codel_subtype", "received_tmz","received_rt","grade_change","hm","new_cell_cycle"),
						 names_to = "var",
						 values_to = "value") %>%
			mutate(value = recode(value, "Grade up" = "1", "Grade stable" = "0")) %>%
			mutate(value = recode(value, "IDHmut-noncodel" = "0", "IDHmut-codel" = "1")) %>%
			mutate(var = recode(var, "idh_codel_subtype" = "1p/19q co-deletion", "received_tmz" = "TMZ",
									 "received_rt" = "Radiotherapy", "grade_change" = "Grade increase",
									 "hm" = "Hypermutator", "new_cell_cycle" = "Acquired cell cycle alt.")) %>%
			mutate(var = as_factor(var)) %>%
			mutate(var = fct_relevel(var,"Acquired cell cycle alt.", "Hypermutator","Grade increase","Radiotherapy","TMZ","1p/19q co-deletion")) %>%
			mutate(value = as_factor(value))
type <- rep("Clinical", nrow(plot_dat))
type[which(plot_dat$var %in% c("Acquired cell cycle alt.", "Hypermutator"))] <- "Mutation"
plot_dat$type <- type
plot_dat <- plot_dat %>% mutate(type = fct_relevel(type,"Mutation", "Clinical"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_prolif_stem_cell_tile_v2.pdf",width=2.2,height=1.1)
ggplot(plot_dat, aes(x=tumor_pair_barcode, y = var, fill = factor(value))) +
geom_tile() +
theme_classic() +
scale_fill_manual(values=c("white","black"),na.value="#E5E5E5") +
facet_grid(type~.,scales="free_y",space="free") +
theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.ticks.x = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_blank(),
	strip.background = element_blank(),
	legend.position = "none")
dev.off()

# Survival analysis

sub_dat <- sub_dat %>%
mutate(case_vital_status = recode(case_vital_status, 'alive' = 0, 'dead' = 1))

mytag <- rep(0, nrow(sub_dat))
mytag[which(sub_dat$diff <= 0)] <- 1
mytag[which(sub_dat$diff > 0)] <- 2
sub_dat$mytag <- mytag

diff = survdiff(Surv(case_overall_survival_mo, case_vital_status)~mytag, data=sub_dat)
p.val <- 1-pchisq(diff$chisq, length(diff$n) - 1)

source("/projects/verhaak-lab/USERS/varnf/SofWar/R/createSurvivalFrame.r")
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/qplot_survival.r")

fit <- survfit(Surv(case_overall_survival_mo, case_vital_status)~mytag, data=sub_dat)
frame <- createSurvivalFrame(fit)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/prolif_stemcell_km_plots.pdf",width=2.4,height=1.9)
qplot_survival(frame, FALSE) +
scale_colour_manual(values = c("#CD4F39", "#27408B"), 
	labels = c(paste("Increase (n=",sum(mytag==2),")",sep=""),paste("Non-increase (n=",sum(mytag==1),")",sep=""))) +
annotate("text", x=0.83*max(sub_dat[,"case_overall_survival_mo"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(p.val )))), parse=TRUE, size=2.5) +
labs(x = "Overall survival (months)", y = "Probability of survival") +
theme_classic() +
theme(axis.text= element_text(size=7),axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.justification=c(1,1), legend.position=c(x=1.03,y=0.96), 
	legend.text=element_text(size=7), legend.title=element_blank(),
	legend.key.size=unit(0.5, "line"),legend.background = element_rect(fill = "transparent",colour=NA) ) +
coord_cartesian(ylim=c(0,1))
dev.off()


# Multivariable forest plot
sub_dat$mytag <- as.factor(mytag)
sub_dat$grade_a <- as.factor(sub_dat$grade_a)
mulvar <- coxph(Surv(case_overall_survival_mo, case_vital_status)~mytag + case_age_diagnosis_years + grade_a + idh_codel_subtype, data=sub_dat)
mulvar_summary <- summary(mulvar)
hr <- mulvar_summary$coefficients[6:10]
low <- mulvar_summary$conf.int[11:15]
high <- mulvar_summary$conf.int[16:20]
vars <- c("Prolif. stem-like", "Age", "Init. grade (3)", "Init. grade (4)","1p/19q co-deletion")

plot_mult <- data.frame(vars, hr, low, high)
plot_mult <- plot_mult %>%
			 mutate(vars = fct_relevel(vars, rev(c("Prolif. stem-like", "Age", "Init. grade (3)", "Init. grade (4)","1p/19q co-deletion") )))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/prolif_stemcell_forest_ggplot2_v2.pdf",width=1.8,height=1.75)
ggplot(data=plot_mult,aes(x=hr,y=vars))+
  geom_point(colour="black")+
  geom_errorbarh(aes(xmin=low,xmax=high),height=0.3)+
  geom_vline(xintercept=1,linetype="dashed") +
  labs(title = "IDH-mutant") + 
  theme_classic() +
	theme(axis.text= element_text(size=7),
	axis.title = element_blank(),
	plot.title = element_text(size=7),	
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position = NULL) +
	scale_x_continuous(breaks=seq(0,12,3))
dev.off()


