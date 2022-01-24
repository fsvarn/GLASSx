###################################################
# Examine association between HLA LOH, SCNA, and immune checkpoint blockade response
# Author: Frederick Varn
# Date: 2021.10.28
# Revision comment
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

#Read in data
q <- "
 WITH roman_to_int(grade, grade_int) AS (
         VALUES ('I'::text,1), ('II'::text,2), ('III'::text,3), ('IV'::text,4)
        ),
tumor_pd1_comparison AS (
 SELECT tp.tumor_pair_barcode,
    tp.case_barcode,
    tp.tumor_barcode_a,
    tp.tumor_barcode_b,
        CASE
            WHEN s1.surgery_location::text = s2.surgery_location::text AND (s1.surgery_laterality::text = s2.surgery_laterality::text OR s1.surgery_laterality IS NULL AND s2.surgery_laterality IS NULL) THEN 'Local'::text
            WHEN s1.surgery_location::text <> s2.surgery_location::text OR s1.surgery_laterality::text <> s2.surgery_laterality::text THEN 'Distal'::text
            ELSE NULL::text
        END AS recurrence_location,
    ( SELECT bool_or(ss.treatment_tmz OR ss.treatment_concurrent_tmz) AS bool_or
           FROM clinical.surgeries ss
          WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_tmz,
    (( SELECT sum(ss.treatment_tmz_cycles) AS sum
           FROM clinical.surgeries ss
          WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number))::integer AS received_tmz_sum_cycles,
    ( SELECT bool_or(ss.treatment_radiotherapy) AS bool_or
           FROM clinical.surgeries ss
          WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt,
    ( SELECT sum(ss.treatment_radiation_dose_gy) AS sum
           FROM clinical.surgeries ss
          WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_rt_sum_gy,
    ( SELECT bool_or(ss.treatment_alkylating_agent) AS bool_or
           FROM clinical.surgeries ss
          WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_alk,
    ( SELECT bool_or(ss.treatment_alkylating_agent OR ss.treatment_radiotherapy) AS bool_or
           FROM clinical.surgeries ss
          WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_treatment,
         ( SELECT bool_or(
                CASE
                    WHEN ss.treatment_chemotherapy_other ~~ '%Nivolumab%'::text OR ss.treatment_chemotherapy_other ~~ '%Pembrolizumab%'::text THEN true
                    ELSE false
                END) AS bool_or
           FROM clinical.surgeries ss
          WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number) AS received_pd1,
		  CASE
            WHEN r2.grade_int > r1.grade_int THEN 'Grade up'::text
            WHEN r2.grade_int = r1.grade_int THEN 'Grade stable'::text
            WHEN r2.grade_int < r1.grade_int THEN 'Grade down'::text
            ELSE NULL::text
        END AS grade_change,
    s2.surgical_interval_mo - s1.surgical_interval_mo AS surgical_interval,
        CASE
            WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10::numeric THEN true
            WHEN mf2.coverage_adj_mut_freq < 10::numeric THEN false
            ELSE NULL::boolean
        END AS hypermutator_status,
        CASE
            WHEN mf2.coverage_adj_mut_freq > mf1.coverage_adj_mut_freq AND mf2.coverage_adj_mut_freq > 10::numeric AND (( SELECT bool_or(ss.treatment_alkylating_agent) AS bool_or
               FROM clinical.surgeries ss
              WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number)) IS TRUE THEN true
            WHEN mf2.coverage_adj_mut_freq < 10::numeric AND (( SELECT bool_or(ss.treatment_alkylating_agent) AS bool_or
               FROM clinical.surgeries ss
              WHERE ss.case_barcode = tp.case_barcode AND ss.surgery_number >= s1.surgery_number AND ss.surgery_number < s2.surgery_number)) IS TRUE THEN false
            ELSE NULL::boolean
        END AS alk_assoc_hypermutator_status
   FROM analysis.tumor_pairs tp
     LEFT JOIN biospecimen.aliquots a1 ON a1.aliquot_barcode = tp.tumor_barcode_a
     LEFT JOIN biospecimen.aliquots a2 ON a2.aliquot_barcode = tp.tumor_barcode_b
     LEFT JOIN clinical.surgeries s1 ON s1.sample_barcode = a1.sample_barcode
     LEFT JOIN clinical.surgeries s2 ON s2.sample_barcode = a2.sample_barcode
     LEFT JOIN roman_to_int r1 ON r1.grade = s1.grade::text
     LEFT JOIN roman_to_int r2 ON r2.grade = s2.grade::text
     LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = tp.tumor_barcode_a
     LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = tp.tumor_barcode_b
)
SELECT rc.*, cc.case_vital_status, cc.case_overall_survival_mo, an1.prop_aneuploidy AS prop_aneuploidy_a, an2.prop_aneuploidy AS prop_aneuploidy_b,
idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.gold_set ps
JOIN tumor_pd1_comparison rc ON rc.tumor_barcode_a = ps.tumor_barcode_a AND rc.tumor_barcode_b = ps.tumor_barcode_b
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE received_pd1 --AND ps.case_barcode LIKE 'GLSS-CU-P%'

"

dat <- dbGetQuery(con,q)

# Get responder data from paper (Zhao et al 2019 Nat Med):
resp <- c("GLSS-CU-P101" = "R", "GLSS-CU-P071" = "R", "GLSS-CU-P020" = "NR", "GLSS-CU-P055" = "R", "GLSS-MD-0042" = NA, "GLSS-NS-0001" = NA, "GLSS-MD-0035" = NA)

plot_dat <- dat[,c("case_barcode", "prop_aneuploidy_a", "prop_aneuploidy_b", "idh_status", "case_vital_status", "case_overall_survival_mo")]
plot_dat[,"aneuploidy_change"] <- plot_dat$prop_aneuploidy_b - plot_dat$prop_aneuploidy_a
resp <- resp[plot_dat[,"case_barcode"]]

plot_dat <- plot_dat %>%
			mutate(case_vital_status = recode(case_vital_status, "alive" = 0, "dead" = 1))

plot_dat[,"resp"] <- resp

p1 <- ggplot(plot_dat %>% filter(grepl("GLSS-CU-P", case_barcode)), aes(fill=idh_status, y=prop_aneuploidy_a, x=case_barcode)) + 
geom_bar(position="dodge", stat="identity") +
scale_fill_manual(values=c("#F8766D", "#619CFF")) +
labs(y = "SCNA burden (initial)") +
facet_grid(.~resp, scales = "free_x", space = "free_x") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none",
legend.title = element_blank(),
legend.text = element_text(size=7))

p2 <- ggplot(plot_dat %>% filter(grepl("GLSS-CU-P", case_barcode)), aes(fill=idh_status, y=prop_aneuploidy_b, x=case_barcode)) + 
geom_bar(position="dodge", stat="identity") +
scale_fill_manual(values=c("#F8766D", "#619CFF")) +
labs(y = "SCNA burden (recurrent)") +
facet_grid(.~resp, scales = "free_x", space = "free_x") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none",
legend.title = element_blank(),
legend.text = element_text(size=7))

p3 <- ggplot(plot_dat %>% filter(grepl("GLSS-CU-P", case_barcode)), aes(fill=idh_status, y=aneuploidy_change, x=case_barcode)) + 
geom_bar(position="dodge", stat="identity") +
scale_fill_manual(values=c("#F8766D", "#619CFF")) +
labs(y = "SCNA burden (longitudinal change)") +
facet_grid(.~resp, scales = "free_x", space = "free_x") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none",
legend.title = element_blank(),
legend.text = element_text(size=7))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scna_burden_pd1_legend.pdf", width=5, height = 2) #,width=2.7,height=3)
ggplot(plot_dat %>% filter(grepl("GLSS-CU-P", case_barcode)), aes(fill=idh_status, y=aneuploidy_change, x=case_barcode)) + 
geom_bar(position="dodge", stat="identity") +
scale_fill_manual(values=c("#F8766D", "#619CFF")) +
labs(y = "SCNA burden (longitudinal change)") +
facet_grid(.~resp, scales = "free_x", space = "free_x") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scna_burden_pd1.pdf", width=5, height = 2) #,width=2.7,height=3)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()


# Check survival
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/createSurvivalFrame.r")
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/qplot_survival.r")

mytag1 <- as.numeric(plot_dat[,"prop_aneuploidy_a"] > median(plot_dat[,"prop_aneuploidy_a"]))
mytag2 <- as.numeric(plot_dat[,"prop_aneuploidy_b"] > median(plot_dat[,"prop_aneuploidy_b"]))
mytag3 <- as.numeric(plot_dat[,"aneuploidy_change"] > 0)

diff1 = survdiff(Surv(case_overall_survival_mo, case_vital_status)~mytag1, data=plot_dat)
p.val1= 1-pchisq(diff1$chisq, length(diff1$n)-1)

diff2 = survdiff(Surv(case_overall_survival_mo, case_vital_status)~mytag2, data=plot_dat)
p.val2= 1-pchisq(diff2$chisq, length(diff2$n)-1)

diff3 = survdiff(Surv(case_overall_survival_mo, case_vital_status)~mytag3, data=plot_dat)
p.val3= 1-pchisq(diff3$chisq, length(diff3$n)-1)

fit1 <- survfit(Surv(case_overall_survival_mo, case_vital_status)~mytag1, data=plot_dat)
fit2 <- survfit(Surv(case_overall_survival_mo, case_vital_status)~mytag2, data=plot_dat)
fit3 <- survfit(Surv(case_overall_survival_mo, case_vital_status)~mytag3, data=plot_dat)
frame1 <- createSurvivalFrame(fit1)
frame2 <- createSurvivalFrame(fit2)
frame3 <- createSurvivalFrame(fit3)
p.val1 <- signif(p.val1, digits = 2)
p.val2 <- signif(p.val2, digits = 2)
p.val3 <- signif(p.val3, digits = 2)

se1 <- qplot_survival(frame1, FALSE) +
scale_colour_manual(values = c("tomato3","royalblue4"), 
	labels = c(paste("High (n=",sum(mytag1==1),")",sep=""),paste("Low (n=",sum(mytag1==0),")",sep=""))) +
annotate("text", x=0.83*max(plot_dat[,"case_overall_survival_mo"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(p.val1)))), parse=TRUE, size=2.5) +
labs(x = "Time (months)", y = "Overall survival probability",title="Initial SCNA burden") +
theme_bw() +
theme(axis.text= element_text(size=7),axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.justification=c(1,1), legend.position=c(x=1.03,y=0.96), 
	legend.text=element_text(size=7), legend.title=element_blank(),
	legend.key.size=unit(0.5, "line"),legend.background = element_rect(fill = "transparent",colour=NA) ) +
coord_cartesian(ylim=c(0,1))

se2 <- qplot_survival(frame2, FALSE) +
scale_colour_manual(values = c("tomato3","royalblue4"), 
	labels = c(paste("High (n=",sum(mytag2==1),")",sep=""),paste("Low (n=",sum(mytag2==0),")",sep=""))) +
annotate("text", x=0.83*max(plot_dat[,"case_overall_survival_mo"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(p.val2)))), parse=TRUE, size=2.5) +
labs(x = "Time (months)", y = "Overall survival probability",title="Recurrent SCNA burden") +
theme_bw() +
theme(axis.text= element_text(size=7),axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.justification=c(1,1), legend.position=c(x=1.03,y=0.96), 
	legend.text=element_text(size=7), legend.title=element_blank(),
	legend.key.size=unit(0.5, "line"),legend.background = element_rect(fill = "transparent",colour=NA) ) +
coord_cartesian(ylim=c(0,1))

se3 <- qplot_survival(frame3, FALSE) +
scale_colour_manual(values = c("tomato3","royalblue4"), 
	labels = c(paste("Inc. (n=",sum(mytag3==1),")",sep=""),paste("Dec. (n=",sum(mytag3==0),")",sep=""))) +
annotate("text", x=0.83*max(plot_dat[,"case_overall_survival_mo"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(p.val3)))), parse=TRUE, size=2.5) +
labs(x = "Time (months)", y = "Overall survival probability",title="SCNA burden (change)") +
theme_bw() +
theme(axis.text= element_text(size=7),axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.justification=c(1,1), legend.position=c(x=1.03,y=0.96), 
	legend.text=element_text(size=7), legend.title=element_blank(),
	legend.key.size=unit(0.5, "line"),legend.background = element_rect(fill = "transparent",colour=NA) ) +
coord_cartesian(ylim=c(0,1))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scna_pd1_km_plots.pdf",width=6.1,height=1.85)
grid.arrange(se1, se2, se3, ncol = 3)
dev.off()