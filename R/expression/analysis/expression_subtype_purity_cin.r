###################################################
# Plot chromosomal instability and purity values in initial and recurrent tumors
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
WITH subtype_rank AS
(
	SELECT *,
	RANK() OVER (PARTITION BY aliquot_barcode ORDER BY p_value ASC) AS p_rank
	FROM analysis.transcriptional_subtype
),
top_rank AS
(
	SELECT *
	FROM subtype_rank
	WHERE p_rank = 1
),
agg AS
(
	SELECT aliquot_barcode, 
	string_agg(signature_name,',') AS subtype 
	FROM top_rank
	GROUP BY aliquot_barcode
)
SELECT ps.*, 
an1.prop_aneuploidy AS prop_aneuploidy_a, an2.prop_aneuploidy AS prop_aneuploidy_b, 
ts1.subtype AS subtype_a, ts2.subtype AS subtype_b,
es1.purity AS purity_a, es2.purity AS purity_b,
cs.idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b
JOIN agg ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
JOIN agg ts2 ON ts2.aliquot_barcode = ps.rna_barcode_b
JOIN analysis.estimate es1 ON es1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate es2 ON es2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con,q)


dat[which(dat[,"subtype_a"]=="Mesenchymal,Classical"),"subtype_a"]  <- "Classical,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Classical"),"subtype_b"]  <- "Classical,Mesenchymal"

dat[which(dat[,"subtype_a"]=="Mesenchymal,Proneural"),"subtype_a"]  <- "Proneural,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Proneural"),"subtype_b"]  <- "Proneural,Mesenchymal"

#Convert for the purposes of plotting
dat[which(dat[,"subtype_a"]=="Proneural,Classical"),"subtype_a"]  <- "Proneural"
dat[which(dat[,"subtype_b"]=="Proneural,Classical"),"subtype_b"]  <- "Proneural"

dat[which(dat[,"subtype_a"]=="Proneural,Mesenchymal"),"subtype_a"]  <- "Mesenchymal"
dat[which(dat[,"subtype_b"]=="Proneural,Mesenchymal"),"subtype_b"]  <- "Mesenchymal"

dat[which(dat[,"subtype_a"]=="Classical,Mesenchymal"),"subtype_a"]  <- "Classical"
dat[which(dat[,"subtype_b"]=="Classical,Mesenchymal"),"subtype_b"]  <- "Classical"

##################################################
# Step 1: Test/plot purity and chromosomal instability values across subtypes
##################################################

# Initial purity
g1 <- dat[which(dat[,"subtype_a"]=="Mesenchymal"),"purity_a"]
g2 <- dat[which(dat[,"subtype_a"]!="Mesenchymal"),"purity_a"]
init_pur_pval <- wilcox.test(g1,g2)$p.value			#9e-9

# Recurrent purity
g1 <- dat[which(dat[,"subtype_b"]=="Mesenchymal"),"purity_b"]
g2 <- dat[which(dat[,"subtype_b"]!="Mesenchymal"),"purity_b"]
rec_pur_pval <- wilcox.test(g1,g2)$p.value			#2e-8

init_pur_pval <- formatC(init_pur_pval, format = "e", digits = 0)
rec_pur_pval <- formatC(rec_pur_pval, format = "e", digits = 0)

# Initial chromosomal instability
g1 <- dat[which(dat[,"subtype_a"]=="Mesenchymal"),"prop_aneuploidy_a"]
g2 <- dat[which(dat[,"subtype_a"]!="Mesenchymal"),"prop_aneuploidy_a"]
wilcox.test(g1,g2)			#9e-9

# Recurrent chromosomal instability
g1 <- dat[which(dat[,"subtype_b"]=="Mesenchymal"),"prop_aneuploidy_b"]
g2 <- dat[which(dat[,"subtype_b"]!="Mesenchymal"),"prop_aneuploidy_b"]
wilcox.test(g1,g2)			#2e-8

# Initial chromosomal instability ANOVA
init_an_pval <- round(summary(aov(dat[,"prop_aneuploidy_a"]~dat[,"subtype_a"]))[[1]][["Pr(>F)"]][1],2)
rec_an_pval <- round(summary(aov(dat[,"prop_aneuploidy_b"]~dat[,"subtype_b"]))[[1]][["Pr(>F)"]][1],2)

##################################################
# Step 2: Test/plot purity and chromosomal instability values across subtypes
##################################################

case_barcode <- rep(dat[,"case_barcode"],2)
aliquot_barcode <- c(dat[,"rna_barcode_a"], dat[,"rna_barcode_b"])
transcriptional_subtype <- c(dat[,"subtype_a"], dat[,"subtype_b"])
purity <- c(dat[,"purity_a"], dat[,"purity_b"])
chromosomal_instability <- c(dat[,"prop_aneuploidy_a"], dat[,"prop_aneuploidy_b"])
idh_codel_subtype <- rep(dat[,"idh_codel_subtype"],2)
timepoint <- c(rep("Initial",nrow(dat)), rep("Recurrent",nrow(dat)))

plot_data <- data.frame(case_barcode, aliquot_barcode, transcriptional_subtype, purity, chromosomal_instability, idh_codel_subtype, timepoint)

plot_data[,"transcriptional_subtype"] <- factor(plot_data[,"transcriptional_subtype"], levels = c("Proneural","Classical","Mesenchymal"))

p_val <- c(deparse((bquote(italic("P") ~" = " ~ .(init_pur_pval)))),
		   deparse((bquote(italic("P") ~" = " ~ .(rec_pur_pval)))))
		   
annotation_text <- data.frame(timepoint = factor(c("Initial","Recurrent")),
							  transcriptional_subtype = 1.4,
							  purity = 0.0,
							  p_val)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/transcriptional_subtype_purity.pdf",width=2.5,height=2.5)
ggplot(plot_data, aes(x = transcriptional_subtype, y = purity)) +
geom_boxplot(outlier.shape=NA,colour="black") +
geom_jitter(size=1,aes(colour = transcriptional_subtype)) +
scale_colour_manual(values=c("#00458a","#008a22","#8a0000")) +
facet_grid(. ~ timepoint) +
labs(y = "Purity") +
geom_text(data=annotation_text,label=p_val, size=2.5, parse=TRUE) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(0,1))
dev.off()



p_val <- c(deparse((bquote(italic("P") ~" = " ~ .(init_an_pval)))),
		   deparse((bquote(italic("P") ~" = " ~ .(rec_an_pval)))))
		   
annotation_text <- data.frame(timepoint = factor(c("Initial","Recurrent")),
							  transcriptional_subtype = 1.4,
							  chromosomal_instability = min(plot_data[,"chromosomal_instability"])-(0.15 * (max(plot_data[,"chromosomal_instability"]) - min(plot_data[,"chromosomal_instability"]))) + 
							  0.04 * abs(min(plot_data[,"chromosomal_instability"])-(0.15 * (max(plot_data[,"chromosomal_instability"]) - min(plot_data[,"chromosomal_instability"])))),
							  p_val)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/transcriptional_subtype_chromosomal_instability.pdf",width=2.5,height=2.5)
ggplot(plot_data, aes(x = transcriptional_subtype, y = chromosomal_instability)) +
geom_boxplot(outlier.shape=NA,colour="black") +
geom_jitter(size=1,aes(colour = transcriptional_subtype)) +
scale_colour_manual(values=c("#00458a","#008a22","#8a0000")) +
facet_grid(. ~ timepoint) +
labs(y = "Purity") +
geom_text(data=annotation_text,label=p_val, size=2.5, parse=TRUE) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(min(plot_data[,"chromosomal_instability"])-(0.15 * (max(plot_data[,"chromosomal_instability"]) - min(plot_data[,"chromosomal_instability"]))),
max(plot_data[,"chromosomal_instability"])))
dev.off()