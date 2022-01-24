###################################################
# Examine association between survival and mesenchymal subtype in GLASS
# Author: Frederick Varn
# Updated: 2022.01.06
# Figures S5A, S5B
##################################################

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
q <- "
WITH full_data AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	subtype_a, 
	subtype_b,
	grade_a,
	received_rt,
	received_tmz,
	received_treatment,
	surgical_interval,
	CASE WHEN ss.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
	FROM analysis.tumor_rna_clinical_comparison ss
	WHERE idh_codel_subtype IS NOT NULL
)
SELECT fd.case_barcode, subtype_a, subtype_b, grade_a, idh_status, received_rt, received_tmz, received_treatment, surgical_interval, cc.case_overall_survival_mo, cc.case_vital_status, cc.case_age_diagnosis_years
FROM full_data fd
JOIN clinical.cases cc ON cc.case_barcode = fd.case_barcode
"

dat <- dbGetQuery(con,q)

dat <- dat %>%
mutate(case_vital_status = recode(case_vital_status, 'alive' = 0, 'dead' = 1)) %>%
mutate(recurrence = 1) %>% 
filter(idh_status == "IDHwt")

dat$grade_a = factor(dat$grade_a)

mytag <- rep(0, nrow(dat))
mytag[which(dat$subtype_b != "Mesenchymal")] <- 1
mytag[which(dat$subtype_b == "Mesenchymal")] <- 2
dat$mytag <- mytag

#survdiff(Surv(case_overall_survival_mo, case_vital_status)~mytag, data=dat)
diff = survdiff(Surv(surgical_interval, recurrence)~mytag, data=dat)
p.val <- 1-pchisq(diff$chisq, length(diff$n) - 1)

source("/projects/verhaak-lab/USERS/varnf/SofWar/R/createSurvivalFrame.r")
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/qplot_survival.r")

fit <- survfit(Surv(surgical_interval, recurrence)~mytag, data=dat)
frame <- createSurvivalFrame(fit)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/mesenchymal_km_plots.pdf",width=1.7,height=1.35)
qplot_survival(frame, FALSE) +
scale_colour_manual(values = c("#CD4F39", "#27408B"), 
	labels = c(paste("Mes. (n=",sum(mytag==2),")",sep=""),paste("Non-mes. (n=",sum(mytag==1),")",sep=""))) +
annotate("text", x=0.83*max(dat[,"case_overall_survival_mo"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(p.val )))), parse=TRUE, size=2.5) +
labs(x = "Surgical interval (months)", y = "Probability of surgery") +
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
dat$mytag <- as.factor(mytag)
summary(coxph(Surv(surgical_interval, recurrence)~mytag + case_age_diagnosis_years , data=dat))
summary(coxph(Surv(surgical_interval, recurrence)~mytag + case_age_diagnosis_years + grade_a , data=dat))

mulvar <- coxph(Surv(surgical_interval, recurrence)~mytag + case_age_diagnosis_years + grade_a , data=dat)
mulvar_summary <- summary(mulvar)
hr <- mulvar_summary$coefficients[4:6]
low <- mulvar_summary$conf.int[7:9]
high <- mulvar_summary$conf.int[10:12]
vars <- c("Mes. recurrence", "Age", "Init. grade")

plot_mult <- data.frame(vars, hr, low, high)
plot_mult <- plot_mult %>%
			 mutate(vars = fct_relevel(vars, rev(c("Mes. recurrence", "Age", "Init. grade"))))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/mesenchymal_forest_ggplot2.pdf",width=1.8,height=1.35)
ggplot(data=plot_mult,aes(x=hr,y=vars))+
  geom_point(colour="black")+
  geom_errorbarh(aes(xmin=low,xmax=high),height=0.3)+
  geom_vline(xintercept=1,linetype="dashed") +
  labs(title = "IDH-wild-type") + 
  theme_classic() +
	theme(axis.text= element_text(size=7),
	axis.title.y = element_blank(),
	axis.title.x = element_text(size=7),
	plot.title = element_blank(),	
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position = NULL) 
dev.off()
