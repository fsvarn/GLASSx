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
q <- "SELECT ps.*, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs.idh_codel_subtype,
cc.case_age_diagnosis_years AS age
FROM analysis.tumor_clinical_comparison ps
JOIN analysis.rna_dna_pairs rd ON rd.dna_pair_barcode = ps.tumor_pair_barcode
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
ORDER BY 2, 1, 16"

dat <- dbGetQuery(con, q)

dat[,"recur_status"] <- rep(1,nrow(dat))

cells <- unique(dat[,"signature_name"])

hr <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	mycox <- coxph(Surv(surgical_interval, recur_status) ~ es_a + idh_codel_subtype + age, data=sub_dat)
	mycox <- summary(mycox)
	hr[i] <- mycox$coefficients[5]
	p.value[i] <- mycox$coefficients[17]
}

initial_res <- data.frame(cells,hr,p.value)

hr <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	mycox <- coxph(Surv(surgical_interval, recur_status) ~ es_b + idh_codel_subtype+ age, data=sub_dat)
	mycox <- summary(mycox)
	hr[i] <- mycox$coefficients[5]
	p.value[i] <- mycox$coefficients[17]
}

recurrent_res <- data.frame(cells,hr,p.value)

#When split up by subtypes, only significant for IDHwt

#Log-rank tests split by median (IDHwt only)

#Subset on IDHwt only
dat <- dat[which(dat[,"idh_codel_subtype"]=="IDHwt"),]

lr.p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]

	mytag <- as.numeric(sub_dat[,"es_a"] > median(sub_dat[,"es_a"]))
	diff = survdiff(Surv(surgical_interval, recur_status)~mytag, data=sub_dat)
	lr.p.value[i] = 1-pchisq(diff$chisq, length(diff$n)-1)
}
initial_km_res <- data.frame(cells,lr.p.value)

hr <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]

	mytag <- as.numeric(sub_dat[,"es_b"] > median(sub_dat[,"es_b"]))
	diff = survdiff(Surv(surgical_interval, recur_status)~mytag, data=sub_dat)
	lr.p.value[i] = 1-pchisq(diff$chisq, length(diff$n)-1)
}

recurrent_km_res <- data.frame(cells,lr.p.value)

#######################################################

#Kaplan-Meier plots
source("/projects/varnf/SofWar/R/createSurvivalFrame.r")
source("/projects/varnf/SofWar/R/qplot_survival.r")

#Macrophages
mac_dat <- dat[which(dat[,"signature_name"]=="Macrophages"),]
mytag1 <- as.numeric(mac_dat[,"es_a"] > median(mac_dat[,"es_a"]))
mytag2 <- as.numeric(mac_dat[,"es_b"] > median(mac_dat[,"es_b"]))
fit1 <- survfit(Surv(surgical_interval, recur_status)~mytag1, data=mac_dat)
fit2 <- survfit(Surv(surgical_interval, recur_status)~mytag2, data=mac_dat)
frame1 <- createSurvivalFrame(fit1)
frame2 <- createSurvivalFrame(fit2)
pval1 <- initial_km_res[which(initial_km_res[,"cells"]=="Macrophages"),"lr.p.value"]
pval1 <- signif(pval1, digits = 1)
pval2 <- recurrent_km_res[which(recurrent_km_res[,"cells"]=="Macrophages"),"lr.p.value"]
pval2 <- signif(pval2, digits = 1)
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/macrophage_km_plots.pdf",width=1.8,height=3)
se1 <- qplot_survival(frame1, FALSE) +
scale_colour_manual(values = c("tomato3","royalblue4"), 
	labels = c(paste("High (n=",sum(mytag1==1),")",sep=""),paste("Low (n=",sum(mytag1==0),")",sep=""))) +
annotate("text", x=0.83*max(mac_dat[,"surgical_interval"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval1)))), parse=TRUE, size=2.5) +
labs(x = "Surgical interval (months)", y = "Probability of surgery",title="Initial macrophage content") +
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
annotate("text", x=0.80*max(mac_dat[,"surgical_interval"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval2)))), parse=TRUE, size=2.5) +
labs(x = "Surgical interval (months)", y = "Probability of surgery",title="Recurrent macrophage content") +
theme_bw() +
theme(axis.text= element_text(size=7),axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.justification=c(1,1), legend.position=c(x=1.03,y=0.96), 
	legend.text=element_text(size=7), legend.title=element_blank(),
	legend.key.size=unit(0.5, "line"),legend.background = element_rect(fill = "transparent",colour=NA) ) +
coord_cartesian(ylim=c(0,1))

grid.arrange(se1,se2,nrow=2)
dev.off()

#---------------

#CD4+ T cells
cd4_dat <- dat[which(dat[,"signature_name"]=="CD4.mature"),]
mytag1 <- as.numeric(cd4_dat[,"es_a"] > median(cd4_dat[,"es_a"]))
mytag2 <- as.numeric(cd4_dat[,"es_b"] > median(cd4_dat[,"es_b"]))
fit1 <- survfit(Surv(surgical_interval, recur_status)~mytag1, data=cd4_dat)
fit2 <- survfit(Surv(surgical_interval, recur_status)~mytag2, data=cd4_dat)
frame1 <- createSurvivalFrame(fit1)
frame2 <- createSurvivalFrame(fit2)
pval1 <- initial_km_res[which(initial_km_res[,"cells"]=="CD4.mature"),"lr.p.value"]
pval1 <- signif(pval1, digits = 1)
pval2 <- recurrent_km_res[which(recurrent_km_res[,"cells"]=="CD4.mature"),"lr.p.value"]
pval2 <- signif(pval2, digits = 1)


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cd4_km_plots.pdf",width=1.8,height=3)
se1 <- qplot_survival(frame1, FALSE) +
scale_colour_manual(values = c("tomato3","royalblue4"), 
	labels = c(paste("High (n=",sum(mytag1==1),")",sep=""),paste("Low (n=",sum(mytag1==0),")",sep=""))) +
annotate("text", x=0.83*max(cd4_dat[,"surgical_interval"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval1)))), parse=TRUE, size=2.5) +
labs(x = "Surgical interval (months)", y = "Probability of surgery",title="Initial CD4+ T cell content") +
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
annotate("text", x=0.80*max(cd4_dat[,"surgical_interval"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval2)))), parse=TRUE, size=2.5) +
labs(x = "Surgical interval (months)", y = "Probability of surgery",title="Recurrent CD4+ T cell content") +
theme_bw() +
theme(axis.text= element_text(size=7),axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.justification=c(1,1), legend.position=c(x=1.03,y=0.96), 
	legend.text=element_text(size=7), legend.title=element_blank(),
	legend.key.size=unit(0.5, "line"),legend.background = element_rect(fill = "transparent",colour=NA) ) +
coord_cartesian(ylim=c(0,1))

grid.arrange(se1,se2,nrow=2)
dev.off()


#######################################################