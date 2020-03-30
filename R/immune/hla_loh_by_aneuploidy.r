library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
WITH collapse AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, string_agg(hla_loh_change,'; ') AS hla_change
	FROM variants.hla_loh_change 
	GROUP BY case_barcode, tumor_barcode_a, tumor_barcode_b
	ORDER BY 1
)
SELECT co.*, an1.prop_aneuploidy AS prop_aneuploidy_a, an2.prop_aneuploidy AS prop_aneuploidy_b, an2.prop_aneuploidy-an1.prop_aneuploidy AS prop_aneuploidy_change,idh_codel_subtype
FROM collapse co
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = co.tumor_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = co.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = co.case_barcode
--WHERE idh_codel_subtype = 'IDHmut-noncodel'
"

dat <- dbGetQuery(con,q)


loh_dat <- dat[grep("gain",dat[,"hla_change"]),]
wilcox.test(loh_dat[,"prop_aneuploidy_a"], loh_dat[,"prop_aneuploidy_b"],paired=TRUE)	#P = 7e-3

noloh_dat <- dat[grep("loss",dat[,"hla_change"]),]
wilcox.test(noloh_dat[,"prop_aneuploidy_a"], noloh_dat[,"prop_aneuploidy_b"],paired=TRUE)

loh_dat <- dat[grep("stable",dat[,"hla_change"]),]
wilcox.test(loh_dat[,"prop_aneuploidy_a"], loh_dat[,"prop_aneuploidy_b"],paired=TRUE)


g1 <- dat[grep("gain",dat[,"hla_change"]),"prop_aneuploidy_change"]
g2 <- dat[grep("gain",dat[,"hla_change"],invert=TRUE),"prop_aneuploidy_change"]

gain <- rep("None", nrow(dat))
gain[grep("gain",dat[,"hla_change"])] <- "Gain"
dat <- cbind(dat,gain)
dat[,"gain"] <- factor(dat[,"gain"],levels=c("None","Gain"))

wilcox.test(g1,g2)		#P = 0.02

wilcox.test(dat[,"prop_aneuploidy_a"],dat[,"prop_aneuploidy_b"])

plot_dat <- dat
hla_change <- rep("None",nrow(dat))
hla_change[grep("loss",plot_dat[,"hla_change"])] <- "Loss"
hla_change[grep("gain",plot_dat[,"hla_change"])] <- "Gain"
hla_change <- rep(hla_change, 2)

case_barcode <- rep(dat[,"case_barcode"],2)
timepoint <- c(rep("Initial",nrow(dat)),rep("Recurrent",nrow(dat)))
cnv_prop <- c(dat[,"prop_aneuploidy_a"],dat[,"prop_aneuploidy_b"])
plot_dat <- data.frame(case_barcode,timepoint,cnv_change,hla_change)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cnv_prop_vs_lohhla.pdf",width=3,height=2)
ggplot(data = plot_dat, aes(x = timepoint, y = cnv_prop, colour= timepoint)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~hla_change) +
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


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cnv_prop_vs_lohhla_change.pdf",width=1.5,height=2)
ggplot(data = dat, aes(x = gain, y = prop_aneuploidy_change, fill= gain)) +
geom_boxplot(outlier.size=0,colour="black") +
#geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
scale_fill_manual(values=c("royalblue4","tomato3")) +
labs(y = "Delta genome altered proportion") +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")
dev.off()