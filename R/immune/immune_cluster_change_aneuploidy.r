library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggdendro)
library(grid)

#######################################################
rm(list=ls())
#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
SELECT ps.*, an1.prop_aneuploidy AS prop_an1, an2.prop_aneuploidy AS prop_an2, an1.aneuploidy_score AS an_score1, an2.aneuploidy_score AS an_score2, ic.change
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b
JOIN analysis.immune_cluster_change ic ON ic.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con,q)

dat_inc <- dat[which(dat[,"change"]=="Increase"),]
wilcox.test(dat_inc[,"an_score1"],dat_inc[,"an_score2"],paired=TRUE)
median(dat_inc[,"an_score2"]) - median(dat_inc[,"an_score1"])
wilcox.test(dat_inc[,"prop_an1"],dat_inc[,"prop_an2"],paired=TRUE)
median(dat_inc[,"prop_an2"]) - median(dat_inc[,"prop_an1"])

dat_dec <- dat[which(dat[,"change"]=="Decrease"),]
wilcox.test(dat_dec[,"an_score1"],dat_dec[,"an_score2"],paired=TRUE)
median(dat_dec[,"an_score2"]) - median(dat_dec[,"an_score1"])
wilcox.test(dat_dec[,"prop_an1"],dat_dec[,"prop_an2"],paired=TRUE)
median(dat_dec[,"prop_an2"]) - median(dat_dec[,"prop_an1"])

dat_non <- dat[which(dat[,"change"]=="None"),]
wilcox.test(dat_non[,"an_score1"],dat_non[,"an_score2"],paired=TRUE)
median(dat_non[,"an_score2"]) - median(dat_non[,"an_score1"])
wilcox.test(dat_non[,"prop_an1"],dat_non[,"prop_an2"],paired=TRUE)
median(dat_non[,"prop_an2"]) - median(dat_non[,"prop_an1"])

case_barcode <- as.factor(rep(dat[,"case_barcode"],2))
prop_an <- as.numeric(c(dat[,"prop_an1"],dat[,"prop_an2"]))
an_score <- as.numeric(c(dat[,"an_score1"],dat[,"an_score2"]))
change <- as.factor(rep(dat[,"change"],2))
status <- as.factor(c(rep("Initial",nrow(dat)),rep("Recurrent",nrow(dat))))

plot_dat <- data.frame(case_barcode,prop_an,an_score,change,status)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_change_aneuploidy_score.pdf",width=3,height=2)
ggplot(data = plot_dat, aes(x = status, y = an_score, colour= status)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~change) +
scale_colour_manual(values=c("royalblue4","tomato3")) +
labs(y = "Aneuploidy score") +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")
dev.off()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_change_cnv_prop.pdf",width=3,height=2)
ggplot(data = plot_dat, aes(x = status, y = prop_an, colour= status)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~change) +
scale_colour_manual(values=c("royalblue4","tomato3")) +
labs(y = "Proportion of genome altered") +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")
dev.off()

an_dif <- as.numeric(dat[,"an_score2"] - dat[,"an_score1"])
prop_dif <- as.numeric(dat[,"prop_an2"] - dat[,"prop_an1"])
change <- as.factor(dat[,"change"])
case_barcode <- as.factor(dat[,"case_barcode"])
plot_dat2  <- data.frame(case_barcode,an_dif,prop_dif,change)

wilcox.test(plot_dat2[which(plot_dat2[,"change"]=="Increase"),"an_dif"], plot_dat2[which(plot_dat2[,"change"]!="Increase"),"an_dif"])
wilcox.test(plot_dat2[which(plot_dat2[,"change"]=="Decrease"),"an_dif"], plot_dat2[which(plot_dat2[,"change"]=="Decrease"),"an_dif"])

wilcox.test(plot_dat2[which(plot_dat2[,"change"]=="Increase"),"prop_dif"], plot_dat2[which(plot_dat2[,"change"]!="Increase"),"prop_dif"])
wilcox.test(plot_dat2[which(plot_dat2[,"change"]=="Decrease"),"prop_dif"], plot_dat2[which(plot_dat2[,"change"]=="Decrease"),"prop_dif"])


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cluster_score_dif.pdf",width=3,height=2)
ggplot(data = plot_dat2, aes(x = change, y = an_dif, fill= change)) +
geom_hline(yintercept=0) +
geom_boxplot() +
scale_fill_manual(values=c("#26ABE2", "#BD1E2D","None")) +
labs(y = "Aneuploidy score difference") +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position = "none")
dev.off()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cluster_prop_dif.pdf",width=3,height=2)
ggplot(data = plot_dat2, aes(x = change, y = prop_dif, fill= change)) +
geom_hline(yintercept=0) +
geom_boxplot() +
scale_fill_manual(values=c("#26ABE2", "#BD1E2D","None")) +
labs(y = "Chromosomal instability proportion") +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position = "none")
dev.off()