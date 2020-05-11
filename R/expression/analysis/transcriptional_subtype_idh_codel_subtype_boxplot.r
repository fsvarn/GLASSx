###################################################
# Examine how transcriptional subtype and molecular subtype associate with ESTIMATE immune score
# Updated: 2020.05.07
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(ggalluvial)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
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
SELECT ss.*, 
ag1.subtype AS subtype_a, 
ag2.subtype AS subtype_b,
es1.immune_score AS immune_a,
es2.immune_score AS immune_b,
cs.idh_codel_subtype
FROM analysis.rna_silver_set ss
JOIN agg ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN agg ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN analysis.estimate es1 ON es1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.estimate es2 ON es2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
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

#dat[,"subtype_a"] <- recode(dat[,"subtype_a"], "Proneural" = "Pro.","Classical" = "Cla.","Mesenchymal" = "Mes.")
#dat[,"subtype_b"] <- recode(dat[,"subtype_b"], "Proneural" = "Pro.","Classical" = "Cla.","Mesenchymal" = "Mes.")

#######################################################

aliquot_barcode <- c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"])
timepoint <- rep(c("Initial","Recurrent"),each = nrow(dat))
immune <- c(dat[,"immune_a"], dat[,"immune_b"])
expr_subtype <- c(dat[,"subtype_a"], dat[,"subtype_b"])
idh_codel_subtype <- rep(dat[,"idh_codel_subtype"],2)
idh_codel_subtype[grep("IDHmut",idh_codel_subtype)] <- "IDHmut"
plot_dat <- data.frame(aliquot_barcode, timepoint, immune, expr_subtype, idh_codel_subtype)

g1 <- plot_dat[which(plot_dat[,"timepoint"]=="Initial" & plot_dat[,"expr_subtype"] == "Proneural" & plot_dat[,"idh_codel_subtype"]=="IDHwt"),"immune"]
g2 <- plot_dat[which(plot_dat[,"timepoint"]=="Initial" & plot_dat[,"expr_subtype"] == "Proneural" & plot_dat[,"idh_codel_subtype"]=="IDHmut"),"immune"]
wilcox.test(g1,g2) #3e-3

g1 <- plot_dat[which(plot_dat[,"timepoint"]=="Initial" & plot_dat[,"expr_subtype"] == "Mesenchymal" & plot_dat[,"idh_codel_subtype"]=="IDHwt"),"immune"]
g2 <- plot_dat[which(plot_dat[,"timepoint"]=="Initial" & plot_dat[,"expr_subtype"] == "Mesenchymal" & plot_dat[,"idh_codel_subtype"]=="IDHmut"),"immune"]
wilcox.test(g1,g2) #0.22


g1 <- plot_dat[which(plot_dat[,"timepoint"]=="Recurrent" & plot_dat[,"expr_subtype"] == "Proneural" & plot_dat[,"idh_codel_subtype"]=="IDHwt"),"immune"]
g2 <- plot_dat[which(plot_dat[,"timepoint"]=="Recurrent" & plot_dat[,"expr_subtype"] == "Proneural" & plot_dat[,"idh_codel_subtype"]=="IDHmut"),"immune"]
wilcox.test(g1,g2) #0.18

g1 <- plot_dat[which(plot_dat[,"timepoint"]=="Recurrent" & plot_dat[,"expr_subtype"] == "Mesenchymal" & plot_dat[,"idh_codel_subtype"]=="IDHwt"),"immune"]
g2 <- plot_dat[which(plot_dat[,"timepoint"]=="Recurrent" & plot_dat[,"expr_subtype"] == "Mesenchymal" & plot_dat[,"idh_codel_subtype"]=="IDHmut"),"immune"]
wilcox.test(g1,g2) #0.07

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/idh_codel_transcript_subtype_immune.pdf",width=4,height=3)
ggplot(plot_dat, aes(x = expr_subtype, y = immune)) +
geom_boxplot(outlier.size=0, aes(fill=idh_codel_subtype)) +
#geom_point(size=1,colour="black") +
facet_grid(. ~ timepoint) +
labs(y = "Immune score") +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.text = element_text(size=7),
legend.title = element_blank(),
legend.position="right")
#coord_cartesian(ylim=c(min(plot_data[,"enrichment_score"])-(0.15 * (max(plot_data[,"enrichment_score"]) - min(plot_data[,"enrichment_score"]))),
#max(plot_data[,"enrichment_score"])))
dev.off()
