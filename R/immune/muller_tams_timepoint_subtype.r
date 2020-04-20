###################################################
# How do macrophage and microglia levels differ across timepoint and subtype?
# Updated: 2020.04.10
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(reshape)


##################################################
rm(list=ls())

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
SELECT ps.case_barcode, m1.aliquot_barcode AS aliquot_barcode_a, m2.aliquot_barcode AS aliquot_barcode_b,
m1.signature_name, m1.enrichment_score AS es_a, m2.enrichment_score AS es_b, cs.idh_codel_subtype, cc.grade_change
FROM analysis.platinum_set ps
JOIN analysis.muller_tam_score m1 ON ps.rna_barcode_a = m1.aliquot_barcode
JOIN analysis.muller_tam_score m2 ON ps.rna_barcode_b = m2.aliquot_barcode AND m1.signature_name = m2.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.tumor_clinical_comparison cc ON cc.tumor_barcode_a = ps.dna_barcode_a AND cc.tumor_barcode_b = ps.dna_barcode_b
ORDER BY signature_name, m1.aliquot_barcode
"

dat <- dbGetQuery(con,q)

###################################################
# Step 1: Examine how macrophages and microglia differ across subtypes
##################################################

#Compare macrophage levels across subtypes in initial tumors:
g1 <- dat[which(dat[,"idh_codel_subtype"]=="IDHwt" & dat[,"signature_name"]=="Macrophages"),"es_a"]
g2 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-noncodel" & dat[,"signature_name"]=="Macrophages"),"es_a"]
g3 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-codel" & dat[,"signature_name"]=="Macrophages"),"es_a"]

wilcox.test(g1,g2)			# IDHwt vs IDHmut-noncodel: P = 9e-6				***
wilcox.test(g1,g3)			# IDHwt vs IDHmut-codel: P = 6e-3					***
wilcox.test(g2,g3)			# IDHmut-noncodel vs IDHmut-noncodel: P = 0.38
wilcox.test(g1,c(g2,g3))	# IDHwt vs IDHmut: P = 8e-7							***
init_mac_pval <- formatC(wilcox.test(g1,c(g2,g3))$p.value, format = "e", digits = 0)

#Compare microglia levels across subtypes in initial tumors:
g1 <- dat[which(dat[,"idh_codel_subtype"]=="IDHwt" & dat[,"signature_name"]=="Microglia"),"es_a"]
g2 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-noncodel" & dat[,"signature_name"]=="Microglia"),"es_a"]
g3 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-codel" & dat[,"signature_name"]=="Microglia"),"es_a"]

wilcox.test(g1,g2)			# IDHwt vs IDHmut-noncodel: P = 0.92
wilcox.test(g1,g3)			# IDHwt vs IDHmut-codel: P = 0.13
wilcox.test(g2,g3)			# IDHmut-noncodel vs IDHmut-noncodel: P = 0.26
wilcox.test(g1,c(g2,g3))	# IDHwt vs IDHmut: P = 0.41
init_mg_pval <- round(wilcox.test(g1,c(g2,g3))$p.value, 2)

#Compare macrophage levels across subtypes in recurrent tumors:
g1 <- dat[which(dat[,"idh_codel_subtype"]=="IDHwt" & dat[,"signature_name"]=="Macrophages"),"es_b"]
g2 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-noncodel" & dat[,"signature_name"]=="Macrophages"),"es_b"]
g3 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-codel" & dat[,"signature_name"]=="Macrophages"),"es_b"]

wilcox.test(g1,g2)			# IDHwt vs IDHmut-noncodel: P = 7e-6				***
wilcox.test(g1,g3)			# IDHwt vs IDHmut-codel: P = 2e-3					***
wilcox.test(g2,g3)			# IDHmut-noncodel vs IDHmut-noncodel: P = 0.95
wilcox.test(g1,c(g2,g3))	# IDHwt vs IDHmut: P = 2e-7							***
rec_mac_pval <- formatC(wilcox.test(g1,c(g2,g3))$p.value, format = "e", digits = 0)

#Compare microglia levels across subtypes in recurrent tumors:
g1 <- dat[which(dat[,"idh_codel_subtype"]=="IDHwt" & dat[,"signature_name"]=="Microglia"),"es_b"]
g2 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-noncodel" & dat[,"signature_name"]=="Microglia"),"es_b"]
g3 <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-codel" & dat[,"signature_name"]=="Microglia"),"es_b"]

wilcox.test(g1,g2)			# IDHwt vs IDHmut-noncodel: P = 0.31
wilcox.test(g1,g3)			# IDHwt vs IDHmut-codel: P = 6e-3					***
wilcox.test(g2,g3)			# IDHmut-noncodel vs IDHmut-noncodel: P = 0.21
wilcox.test(g1,c(g2,g3))	# IDHwt vs IDHmut: P = 0.03							**
rec_mg_pval <- round(wilcox.test(g1,c(g2,g3))$p.value, 2)

###################################################
# Step 2: Examine how macrophages and microglia change over time
##################################################

subtypes <- unique(dat[,"idh_codel_subtype"])
eff <- p_val <- rep(0, length(subtypes))
for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"idh_codel_subtype"] == subtypes[i] & dat[,"signature_name"] == "Macrophages"),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p_val[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"], paired=TRUE)$p.value
}
mac_time_res <- data.frame(subtypes, eff, p_val)	#no significannce

subtypes <- unique(dat[,"idh_codel_subtype"])
eff <- p_val <- rep(0, length(subtypes))
for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"idh_codel_subtype"] == subtypes[i] & dat[,"signature_name"] == "Microglia"),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p_val[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"], paired=TRUE)$p.value
}
mg_time_res <- data.frame(subtypes, eff, p_val)		#no significance

# Check whether grade change associates with macrophage or microglia increase

sub_dat <- dat[which(dat[,"grade_change"] == "Grade up" & dat[,"signature_name"] == "Macrophages" & dat[,"idh_codel_subtype"] != "IDHwt"),]
mac_eff <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
mac_p_val <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"], paired=TRUE)$p.value

sub_dat <- dat[which(dat[,"grade_change"] == "Grade up" & dat[,"signature_name"] == "Microglia" & dat[,"idh_codel_subtype"] != "IDHwt"),]
mg_eff <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
mg_p_val <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"], paired=TRUE)$p.value

mac_p_val <- round(mac_p_val, 2)
mg_p_val <- round(mg_p_val, 2)

###################################################
# Step 3: Plot these findings
##################################################

case_barcode <- rep(dat[,"case_barcode"],2)
aliquot_barcode <- c(dat[,"aliquot_barcode_a"], dat[,"aliquot_barcode_b"])
signature_name <- rep(dat[,"signature_name"],2)
enrichment_score <- c(dat[,"es_a"], dat[,"es_b"])
idh_codel_subtype <- rep(dat[,"idh_codel_subtype"],2)
timepoint <- c(rep("Initial",nrow(dat)), rep("Recurrent",nrow(dat)))
grade_change <- rep(dat[,"grade_change"],2)

plot_data <- data.frame(case_barcode, aliquot_barcode, signature_name, enrichment_score, idh_codel_subtype, timepoint, grade_change)

# Plot showing the enrichment scores across all tumors:

p_val <- c(deparse((bquote(italic("P") ~" = " ~ .(init_mac_pval)))),
		   deparse((bquote(italic("P") ~" = " ~ .(init_mg_pval)))),
		   deparse((bquote(italic("P") ~" = " ~ .(rec_mac_pval)))),
		   deparse((bquote(italic("P") ~" = " ~ .(rec_mg_pval)))))
annotation_text <- data.frame(signature_name = factor(rep(c("Macrophages","Microglia"),2)),
							  timepoint = factor(c("Initial","Initial","Recurrent","Recurrent")),
							  idh_codel_subtype = 1.4,
							  enrichment_score = min(plot_data[,"enrichment_score"])-(0.15 * (max(plot_data[,"enrichment_score"]) - min(plot_data[,"enrichment_score"]))) + 
							  0.2 * abs(min(plot_data[,"enrichment_score"])-(0.15 * (max(plot_data[,"enrichment_score"]) - min(plot_data[,"enrichment_score"])))),
							  p_val)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/muller_tams_initial_recurrent.pdf",width=2.5,height=2.5)
ggplot(plot_data, aes(x = idh_codel_subtype, y = enrichment_score)) +
geom_boxplot(outlier.size=0, aes(fill=idh_codel_subtype)) +
geom_point(size=1,colour="black") +
facet_grid(signature_name ~ timepoint) +
labs(y = "Infiltration score") +
geom_text(data=annotation_text,label=p_val, size=2.5, parse=TRUE) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(min(plot_data[,"enrichment_score"])-(0.15 * (max(plot_data[,"enrichment_score"]) - min(plot_data[,"enrichment_score"]))),
max(plot_data[,"enrichment_score"])))
dev.off()



# Plot comparing grade increases in IDHmuts:

sub_plot_data <- plot_data[which(plot_data[,"grade_change"]=="Grade up" & plot_data[,"idh_codel_subtype"] != "IDHwt"),]

p_val <- c(deparse((bquote(italic("P") ~" = " ~ .(mac_p_val)))),
			 deparse((bquote(italic("P") ~" = " ~ .(mg_p_val)))))
annotation_text <- data.frame(signature_name = factor(c("Macrophages","Microglia")),
							  timepoint = 1.1,
							  enrichment_score = min(sub_plot_data[,"enrichment_score"])-(0.15 * (max(sub_plot_data[,"enrichment_score"]) - min(sub_plot_data[,"enrichment_score"]))) + 
							  0.04 * abs(min(sub_plot_data[,"enrichment_score"])-(0.15 * (max(sub_plot_data[,"enrichment_score"]) - min(sub_plot_data[,"enrichment_score"])))),
							  p_val)
							  
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/muller_tams_grade_up.pdf",width=2,height=2)
ggplot(sub_plot_data, aes(x = timepoint, y = enrichment_score)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode,colour=idh_codel_subtype)) +
geom_point(size=1,colour="black") +
scale_colour_manual(values=c("#F8766D","#00BA38")) +
facet_grid(.~signature_name) +
labs(y = "Infiltration score") +
geom_text(data=annotation_text,label=p_val, size=2.5, parse=TRUE) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(min(sub_plot_data[,"enrichment_score"])-(0.15 * (max(sub_plot_data[,"enrichment_score"]) - min(sub_plot_data[,"enrichment_score"]))),
max(sub_plot_data[,"enrichment_score"])))
dev.off()

###################################################
# Step 4: Compare the microglia grade increase results to other immune cells (Davoli)
##################################################

# Read in data
q <- 
"
SELECT ps.case_barcode, m1.aliquot_barcode AS aliquot_barcode_a, m2.aliquot_barcode AS aliquot_barcode_b,
m1.signature_name, m1.enrichment_score AS es_a, m2.enrichment_score AS es_b, 
an1.prop_aneuploidy AS prop_aneuploidy_a, an2.prop_aneuploidy AS prop_aneuploidy_b,
cs.idh_codel_subtype, cc.grade_change
FROM analysis.platinum_set ps
JOIN analysis.davoli_immune_score m1 ON ps.rna_barcode_a = m1.aliquot_barcode
JOIN analysis.davoli_immune_score m2 ON ps.rna_barcode_b = m2.aliquot_barcode AND m1.signature_name = m2.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.tumor_clinical_comparison cc ON cc.tumor_barcode_a = ps.dna_barcode_a AND cc.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b
WHERE grade_change = 'Grade up' AND idh_codel_subtype LIKE 'IDHmut%'
ORDER BY signature_name, m1.aliquot_barcode
"

dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])

p_val <- eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	p_val[i] <- wilcox.test(sub_dat[,"es_a"],sub_dat[,"es_b"],paired=TRUE)$p.value
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
}

data.frame(cells, p_val, eff)

# Davoli immune cells: no difference, though dendritic cell results are close
#                   cells      p_val          eff
# 1                B.cells 0.58789063 -0.024854240
# 2             CD4.mature 0.33959961  0.027151535
# 3           CD8.effector 0.68481445  0.047303768
# 4  CD8.effector.NK.cells 0.24389648  0.033719959
# 5              Dendritic 0.05737305 -0.031101909
# 6            Macrophages 0.49731445 -0.015585682
# 7         Macrophages.M1 0.83935547 -0.005017327
# 8         Macrophages.M2 0.27343750 -0.015578422
# 9               NK.cells 0.06811523  0.043499327
# 10                 T.reg 0.73535156 -0.011032812
