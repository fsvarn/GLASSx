library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
WITH lohhla_pairs AS
(
	SELECT gs.case_barcode, gs.tumor_barcode_a, gs.tumor_barcode_b, 
	lh1.hla_type1, lh1.hla_type2, 
	lh1.hla_type1_copy_number AS hla_type1_copy_number_a, lh2.hla_type1_copy_number AS hla_type1_copy_number_b,
	lh1.hla_type2_copy_number AS hla_type2_copy_number_a, lh2.hla_type2_copy_number AS hla_type2_copy_number_b,
	lh1.pval AS pval_a, lh1.loss_allele AS loss_allele_a, 
	lh2.pval AS pval_b, lh2.loss_allele AS loss_allele_b
	FROM analysis.lohhla_set gs
	JOIN analysis.pairs pa1 ON pa1.tumor_barcode = gs.tumor_barcode_a
	JOIN analysis.pairs pa2 ON pa2.tumor_barcode = gs.tumor_barcode_b
	JOIN variants.lohhla lh1 ON lh1.pair_barcode = pa1.pair_barcode
	JOIN variants.lohhla lh2 ON lh2.pair_barcode = pa2.pair_barcode AND lh2.hla_type1 = lh1.hla_type1
	WHERE lh1.coverage_filter = 20 AND lh2.coverage_filter = 20
),
any_loss AS
(
	SELECT lp.case_barcode, tumor_barcode_a, tumor_barcode_b, 
	SUM(CASE WHEN pval_a < 0.1 AND (hla_type1_copy_number_a < 0.5 OR hla_type2_copy_number_a < 0.5) THEN 1 ELSE 0 END) AS loss_a, 
	SUM(CASE WHEN pval_b < 0.1 AND (hla_type1_copy_number_b < 0.5 OR hla_type2_copy_number_b < 0.5) THEN 1 ELSE 0 END) AS loss_b,
	idh_codel_subtype
	FROM lohhla_pairs lp
	JOIN clinical.subtypes cs ON cs.case_barcode = lp.case_barcode
	GROUP BY lp.case_barcode, tumor_barcode_a, tumor_barcode_b, idh_codel_subtype
)
SELECT idh_codel_subtype,
CASE WHEN loss_a > 0 THEN 1 ELSE 0 END AS loss_a,
CASE WHEN loss_b > 0 THEN 1 ELSE 0 END AS loss_b
FROM any_loss
"

dat <- dbGetQuery(con,q)
dat[,"loss_a"] <- as.numeric(dat[,"loss_a"])
dat[,"loss_b"] <- as.numeric(dat[,"loss_b"])

subtypes <- unique(dat[,"idh_codel_subtype"])

subtype <- rep(subtypes,each=4)
loss <- rep(c(rep("loss",2),rep("no_loss",2)),3)
timepoint <- rep(c("Initial","Recurrent"),6)
count <- rep(0,12)
plot_res  <- data.frame(subtype,loss,timepoint,count)
p.val <- rep(0,length(subtypes))
names(p.val) <- subtypes

for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"idh_codel_subtype"]==subtypes[i]),]
	
	g1 <- sum(sub_dat[,"loss_a"])
	g2 <- sum(sub_dat[,"loss_b"])
	g3 <- nrow(sub_dat) - sum(sub_dat[,"loss_a"])
	g4 <- nrow(sub_dat) - sum(sub_dat[,"loss_b"])

	plot_res[which(plot_res[,"subtype"]==subtypes[i] & plot_res[,"timepoint"]=="Initial" & plot_res[,"loss"]=="loss"),"count"] <- g1
	plot_res[which(plot_res[,"subtype"]==subtypes[i] & plot_res[,"timepoint"]=="Recurrent" & plot_res[,"loss"]=="loss"),"count"] <- g2
	plot_res[which(plot_res[,"subtype"]==subtypes[i] & plot_res[,"timepoint"]=="Initial" & plot_res[,"loss"]=="no_loss"),"count"] <- g3
	plot_res[which(plot_res[,"subtype"]==subtypes[i] & plot_res[,"timepoint"]=="Recurrent" & plot_res[,"loss"]=="no_loss"),"count"] <- g4
	
	ct <- matrix(c(g1,g2,g3,g4), nrow=2)
	p.val[i] <- fisher.test(ct)$p.value
}

#Check IDHmut specifically:
sub_dat <- dat[grep("IDHmut",dat[,"idh_codel_subtype"]),]

g1 <- sum(sub_dat[,"loss_a"])
g2 <- sum(sub_dat[,"loss_b"])
g3 <- nrow(sub_dat) - sum(sub_dat[,"loss_a"])
g4 <- nrow(sub_dat) - sum(sub_dat[,"loss_b"])

ct <- matrix(c(g1,g2,g3,g4), nrow=2)
p.val_mut <- fisher.test(ct)$p.value

#Check all :

g1 <- sum(dat[,"loss_a"])
g2 <- sum(dat[,"loss_b"])
g3 <- nrow(dat) - sum(dat[,"loss_a"])
g4 <- nrow(dat) - sum(dat[,"loss_b"])

ct <- matrix(c(g1,g2,g3,g4), nrow=2)
p.val_all <- fisher.test(ct)$p.value

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/lohhla_prevalance_timepoint_subtype.pdf",width=3,height=2)
p1 <- ggplot(data = plot_res, aes(x = timepoint, y = count)) +
geom_bar(aes(fill=loss),stat="identity") +
scale_fill_manual(values=c("#BD1E2D","gray50")) +
facet_grid(.~subtype) +
labs(y = "Number of samples") +
theme_bw() +
theme(axis.text.x= element_text(size=7),axis.text.y= element_text(size=7),
	axis.title.x = element_blank(),axis.title.y = element_text(size=7),
	plot.title = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position="none")
p1
dev.off()


#Quick SQL query to get table:

q<-
"
WITH lohhla_pairs AS
(
	SELECT gs.case_barcode, gs.tumor_barcode_a, gs.tumor_barcode_b, 
	lh1.hla_type1, lh1.hla_type2, 
	lh1.hla_type1_copy_number AS hla_type1_copy_number_a, lh2.hla_type1_copy_number AS hla_type1_copy_number_b,
	lh1.hla_type2_copy_number AS hla_type2_copy_number_a, lh2.hla_type2_copy_number AS hla_type2_copy_number_b,
	lh1.pval AS pval_a, lh1.loss_allele AS loss_allele_a, 
	lh2.pval AS pval_b, lh2.loss_allele AS loss_allele_b
	FROM analysis.lohhla_set gs
	JOIN analysis.pairs pa1 ON pa1.tumor_barcode = gs.tumor_barcode_a
	JOIN analysis.pairs pa2 ON pa2.tumor_barcode = gs.tumor_barcode_b
	JOIN variants.lohhla lh1 ON lh1.pair_barcode = pa1.pair_barcode
	JOIN variants.lohhla lh2 ON lh2.pair_barcode = pa2.pair_barcode AND lh2.hla_type1 = lh1.hla_type1
	WHERE lh1.coverage_filter = 20 AND lh2.coverage_filter = 20
),
any_loss AS
(
	SELECT lp.case_barcode, tumor_barcode_a, tumor_barcode_b, 
	SUM(CASE WHEN pval_a < 0.1 AND (hla_type1_copy_number_a < 0.5 OR hla_type2_copy_number_a < 0.5) THEN 1 ELSE 0 END) AS loss_count_a, 
	SUM(CASE WHEN pval_b < 0.1 AND (hla_type1_copy_number_b < 0.5 OR hla_type2_copy_number_b < 0.5) THEN 1 ELSE 0 END) AS loss_count_b,
	idh_codel_subtype
	FROM lohhla_pairs lp
	JOIN clinical.subtypes cs ON cs.case_barcode = lp.case_barcode
	GROUP BY lp.case_barcode, tumor_barcode_a, tumor_barcode_b, idh_codel_subtype
),
loss_bin AS
(
	SELECT idh_codel_subtype,
	CASE WHEN loss_count_a > 0 THEN 1 ELSE 0 END AS loss_a,
	CASE WHEN loss_count_b > 0 THEN 1 ELSE 0 END AS loss_b
	FROM any_loss
)
SELECT idh_codel_subtype,
ROUND(CAST(SUM(loss_a) AS DECIMAL)/COUNT(*),2) AS fraction_a,
ROUND(CAST(SUM(loss_b) AS DECIMAL)/COUNT(*),2) AS fraction_b
FROM loss_bin
GROUP BY idh_codel_subtype
"

dat <- dbGetQuery(con,q)

#   idh_codel_subtype fraction_a fraction_b
# 1   IDHmut-noncodel       0.11       0.22
# 2             IDHwt       0.19       0.21
# 3      IDHmut-codel       0.05       0.10
