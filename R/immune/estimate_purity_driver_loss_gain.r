###################################################
# How does the acquisition/loss of copy number variants associate with immune changes
# Updated: 2020.04.20
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(reshape)
library(gridExtra)
library(ggbeeswarm)


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
dc.idh_codel_subtype, dc.cnv_driver_shared, dc.cnv_driver_change_a, dc.cnv_driver_change_b, 
ds.snv_driver_shared, ds.snv_driver_change_a, ds.snv_driver_change_b, 
ts1.subtype AS subtype_a, ts2.subtype AS subtype_b,
eg1.egfrviii_status AS egfrviii_a, eg2.egfrviii_status AS egfrviii_b,
es1.immune_score AS immune_a, es2.immune_score AS immune_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv dc ON dc.tumor_barcode_a = ps.dna_barcode_a AND dc.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN agg ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
JOIN agg ts2 ON ts2.aliquot_barcode = ps.rna_barcode_b
JOIN analysis.estimate es1 ON es1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate es2 ON es2.aliquot_barcode = ps.rna_barcode_b
JOIN variants.prada_egfrviii eg1 ON eg1.aliquot_barcode = ps.rna_barcode_a
JOIN variants.prada_egfrviii eg2 ON eg2.aliquot_barcode = ps.rna_barcode_b
WHERE dc.idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

# Recode the subtypes of samples with joint subtypes
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

###################################################
# Step 1: Identify the driver alterations that associate with purity changes
##################################################

# Start with copy number:
#---------------------------

# Get list of unique driver alterations
cnv_drivers <- unique(unlist(strsplit(c(dat[,"cnv_driver_shared"],dat[,"cnv_driver_change_a"],dat[,"cnv_driver_change_b"]),", ")))
cnv_drivers <- gsub("\\+","",cnv_drivers)
cnv_drivers <- gsub("-","",cnv_drivers)
cnv_drivers <- unique(cnv_drivers)
cnv_drivers <- cnv_drivers[-which(is.na(cnv_drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
cnv_genes <- cnv_drivers

gain_p.val <- gain_eff <- loss_p.val <- loss_eff <- rep(NA, length(cnv_genes))
for(i in 1:length(cnv_genes))
{
	if(nrow(dat[grep(cnv_genes[i],dat[,"cnv_driver_change_a"]),]) > 0)
	{
		g1 <- dat[grep(cnv_genes[i],dat[,"cnv_driver_change_a"]),"immune_a"]
		g2 <- dat[grep(cnv_genes[i],dat[,"cnv_driver_change_a"]),"immune_b"]
		gain_p.val[i] <- wilcox.test(g1,g2,paired=TRUE)$p.value
		gain_eff[i] <- median(g2 - g1)
	}

	if(nrow(dat[grep(cnv_genes[i],dat[,"cnv_driver_change_b"]),]) > 0)
	{
		g1 <- dat[grep(cnv_genes[i],dat[,"cnv_driver_change_b"]),"immune_a"]
		g2 <- dat[grep(cnv_genes[i],dat[,"cnv_driver_change_b"]),"immune_b"]
		loss_p.val[i] <- wilcox.test(g1,g2,paired=TRUE)$p.value
		loss_eff[i] <- median(g2 - g1)
	}
}

cnv_results <- data.frame(cnv_genes, gain_p.val, gain_eff, loss_p.val, loss_eff)

# Check out single nucleotide variants:
#---------------------------

# Get list of unique driver alterations
snv_drivers <- unique(unlist(strsplit(c(dat[,"snv_driver_shared"],dat[,"snv_driver_change_a"],dat[,"snv_driver_change_b"]),", ")))
snv_drivers <- unique(sapply(strsplit(snv_drivers," "),function(x)x[1]))
snv_drivers <- gsub("\\+","",snv_drivers)
snv_drivers <- gsub("-","",snv_drivers)
snv_drivers <- unique(snv_drivers)
snv_drivers <- snv_drivers[-which(is.na(snv_drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
snv_genes <- snv_drivers

gain_p.val <- gain_eff <- loss_p.val <- loss_eff <- rep(NA, length(snv_genes))
for(i in 1:length(snv_genes))
{
	if(nrow(dat[grep(snv_genes[i],dat[,"snv_driver_change_a"]),]) > 0)
	{
		g1 <- dat[grep(snv_genes[i],dat[,"snv_driver_change_a"]),"immune_a"]
		g2 <- dat[grep(snv_genes[i],dat[,"snv_driver_change_a"]),"immune_b"]
		gain_p.val[i] <- wilcox.test(g1,g2,paired=TRUE)$p.value
		gain_eff[i] <- median(g2 - g1)
	}

	if(nrow(dat[grep(snv_genes[i],dat[,"snv_driver_change_b"]),]) > 0)
	{
		g1 <- dat[grep(snv_genes[i],dat[,"snv_driver_change_b"]),"immune_a"]
		g2 <- dat[grep(snv_genes[i],dat[,"snv_driver_change_b"]),"immune_b"] 
		loss_p.val[i] <- wilcox.test(g1,g2,paired=TRUE)$p.value
		loss_eff[i] <- median(g2 - g1)
	}
}

snv_results <- data.frame(snv_genes, gain_p.val, gain_eff, loss_p.val, loss_eff)

# Only two significant results: Gains of CDKN2A associate with increased immune, losses of EGFR associate with decreased immune
# Plot them:

sub_dat <- dat[grep("EGFR",dat[,"cnv_driver_change_b"]),]
case_barcode <- rep(sub_dat[,"case_barcode"],2)
immune <- c(sub_dat[,"immune_a"], sub_dat[,"immune_b"])
subtype <- c(sub_dat[,"subtype_a"],sub_dat[,"subtype_b"])
timepoint <- c(rep("Initial",nrow(sub_dat)),rep("Recurrent",nrow(sub_dat)))
egfr_dat <- data.frame(case_barcode, immune, subtype, timepoint)
pval1 <- round(cnv_results[2,4],2)

sub_dat <- dat[grep("CDKN2A",dat[,"cnv_driver_change_a"]),]
case_barcode <- rep(sub_dat[,"case_barcode"],2)
immune <- c(sub_dat[,"immune_a"], sub_dat[,"immune_b"])
subtype <- c(sub_dat[,"subtype_a"],sub_dat[,"subtype_b"])
timepoint <- c(rep("Initial",nrow(sub_dat)),rep("Recurrent",nrow(sub_dat)))
cdkn2a_dat <- data.frame(case_barcode, immune, subtype, timepoint)
pval2 <- round(cnv_results[1,2],2)

egfr_dat[,"subtype"] <- factor(egfr_dat[,"subtype"],levels=c("Proneural","Classical","Mesenchymal"))
cdkn2a_dat[,"subtype"] <- factor(cdkn2a_dat[,"subtype"],levels=c("Proneural","Classical","Mesenchymal"))

###################################################
# Step 2: Plot the significant results
##################################################

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/ckn2a_egfr_gain_loss_ladder.pdf",width=3,height=2)
egfr <- ggplot(egfr_dat ,aes(x = timepoint, y = immune)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode)) +
geom_point(size=1,aes(colour=subtype)) +
scale_colour_manual(values=c("#00458a","#008a22","#8a0000")) +
labs(title="EGFR amplification loss", y = "Immune score") +
annotate("text", x=1.1, y = min(egfr_dat[,"immune"]) + 0.05 * (max(egfr_dat[,"immune"] - min(egfr_dat[,"immune"]))), 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval1)))), parse=TRUE, size=2.5) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")

cdkn2a <- ggplot(cdkn2a_dat ,aes(x = timepoint, y = immune)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode)) +
geom_point(size=1,aes(colour=subtype)) +
scale_colour_manual(values=c("#00458a","#008a22","#8a0000")) +
labs(title="CDKN2A deletion gain", y = "Immune score") +
annotate("text", x=1.1,  y = min(cdkn2a_dat[,"immune"]) + 0.05 * (max(cdkn2a_dat[,"immune"] - min(cdkn2a_dat[,"immune"]))), 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval2)))), parse=TRUE, size=2.5) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 

grid.arrange(egfr,cdkn2a,nrow=1)
dev.off()

###################################################
# Step 3: Split samples by EGFR alterationn to see which drives effect
##################################################

case_barcode <- rep(dat[,"case_barcode"],2)
aliquot_barcode <- c(dat[,"rna_barcode_a"], dat[,"rna_barcode_b"])
idh_codel_subtype <- rep(dat[,"idh_codel_subtype"],2)
transcriptional_subtype <- c(dat[,"subtype_a"], dat[,"subtype_b"])
immune <- c(dat[,"immune_a"], dat[,"immune_b"])
timepoint <- c(rep("Initial",nrow(dat)), rep("Recurrent",nrow(dat)))
egfr_status <- rep("", length(case_barcode))
egfrviii <- as.logical(as.numeric(c(dat[,"egfrviii_a"], dat[,"egfrviii_b"])))


# Four groups: extracellular domain mutation, EGFRvIII, EGFR amp (other), EGFR wt

# Identify all samples that are EGFR wildtype at each timepoint:
#--------------------------------------------------

# Initial wt
egfr_wt_ini <- !(grepl("EGFR",dat[,"cnv_driver_shared"]) | grepl("EGFR",dat[,"cnv_driver_change_b"]) |
grepl("EGFR",dat[,"snv_driver_shared"]) | grepl("EGFR",dat[,"snv_driver_change_b"]))
 
# Recurrence wt
egfr_wt_rec <- !(grepl("EGFR",dat[,"cnv_driver_shared"]) | grepl("EGFR",dat[,"cnv_driver_change_a"]) |
grepl("EGFR",dat[,"snv_driver_shared"]) | grepl("EGFR",dat[,"snv_driver_change_a"]))

egfr_wt <- c(egfr_wt_ini, egfr_wt_rec)
egfr_status[which(egfr_wt)] <- "WT"

# Identify all samples with an EGFR amplification:
#--------------------------------------------------

# Initial
egfr_amp_ini <- grepl("EGFR amp", dat[,"cnv_driver_shared"]) | grepl("EGFR amp", dat[,"cnv_driver_change_b"])

# Recurrent
egfr_amp_rec <- grepl("EGFR amp", dat[,"cnv_driver_shared"]) | grepl("EGFR amp", dat[,"cnv_driver_change_a"])

egfr_amp <- c(egfr_amp_ini, egfr_amp_rec)
egfr_status[which(egfr_amp)] <- "Amplified"

# Identify all samples with an EGFR extracellular domain mutations (some may overlap with amp and that's okay):
#--------------------------------------------------
# Shared
egfr_mut_shared <- sapply(strsplit(dat[,"snv_driver_shared"],", "), 
	function(x){
		muts <- x[grep("EGFR",x)]
		aa_pos <- as.numeric(unlist(sapply(strsplit(muts, " "), function(x) regmatches(x[2], gregexpr("[[:digit:]]+", x[2])))))
		return(aa_pos)})
egfr_ec_shared <- sapply(egfr_mut_shared, function(x) x[1] < 630)		
egfr_ec_shared[which(is.na(egfr_ec_shared))] <- FALSE

# Initial only
egfr_mut_ini <- sapply(strsplit(dat[,"snv_driver_change_b"],", "), 
	function(x){
		muts <- x[grep("EGFR",x)]
		aa_pos <- as.numeric(unlist(sapply(strsplit(muts, " "), function(x) regmatches(x[2], gregexpr("[[:digit:]]+", x[2])))))
		return(aa_pos)})
egfr_ec_ini <- sapply(egfr_mut_ini, function(x) x[1] < 630)		
egfr_ec_ini[which(is.na(egfr_ec_ini))] <- FALSE
egfr_ec_ini <- egfr_ec_shared | egfr_ec_ini

# Recurrence only
egfr_mut_rec <- sapply(strsplit(dat[,"snv_driver_shared"],", "), 
	function(x){
		muts <- x[grep("EGFR",x)]
		aa_pos <- as.numeric(unlist(sapply(strsplit(muts, " "), function(x) regmatches(x[2], gregexpr("[[:digit:]]+", x[2])))))
		return(aa_pos)})
egfr_ec_rec <- sapply(egfr_mut_rec, function(x) x[1] < 630)	
egfr_ec_rec[which(is.na(egfr_ec_rec))] <- FALSE
egfr_ec_rec <- egfr_ec_shared | egfr_ec_rec

egfr_ec <- c(egfr_ec_ini, egfr_ec_rec)
egfr_status[which(egfr_ec)] <- "Extracellular"

# Identify all samples with an EGFR kinase domain mutations (some may overlap with amp and that's okay):
#--------------------------------------------------

# Shared
egfr_kd_shared <- sapply(egfr_mut_shared, function(x) x[1] > 700 & x[1] < 1000)		
egfr_kd_shared[which(is.na(egfr_kd_shared))] <- FALSE

# Initial
egfr_kd_ini <- sapply(egfr_mut_ini, function(x) x[1] > 700 & x[1] < 1000)	
egfr_kd_ini[which(is.na(egfr_kd_ini))] <- FALSE
egfr_kd_ini <- egfr_kd_shared | egfr_kd_ini

# Recurrent
egfr_kd_rec <- sapply(egfr_mut_rec, function(x) x[1] > 700 & x[1] < 1000)	
egfr_kd_rec[which(is.na(egfr_kd_rec))] <- FALSE
egfr_kd_rec <- egfr_kd_shared | egfr_kd_rec

egfr_kd <- c(egfr_kd_ini, egfr_kd_rec)
egfr_status[which(egfr_kd)] <- "Kinase domain"

# Identify all samples with EGFRvIII fusions:
#--------------------------------------------------
egfr_status[which(egfrviii)] <- "EGFRvIII"



plot_egfr_class <- data.frame(case_barcode, aliquot_barcode, idh_codel_subtype, transcriptional_subtype, immune, timepoint, egfr_status)

# Get p-values comparing each alteration to wildtype:
#--------------------------------------------------

uni_stat <- unique(egfr_status)
egfr_ini_pval <- egfr_rec_pval <- rep(0, length(uni_stat))
for(i in 1:length(uni_stat))
{
	g1 <- plot_egfr_class[which(plot_egfr_class[,"timepoint"] == "Initial" & plot_egfr_class[,"egfr_status"] == uni_stat[i]), "immune"]
	g2 <- plot_egfr_class[which(plot_egfr_class[,"timepoint"] == "Initial" & plot_egfr_class[,"egfr_status"] == "WT"), "immune"]
	egfr_ini_pval[i] <- wilcox.test(g1, g2)$p.value

	g1 <- plot_egfr_class[which(plot_egfr_class[,"timepoint"] == "Recurrent" & plot_egfr_class[,"egfr_status"] == uni_stat[i]), "immune"]
	g2 <- plot_egfr_class[which(plot_egfr_class[,"timepoint"] == "Recurrent" & plot_egfr_class[,"egfr_status"] == "WT"), "immune"]
	egfr_rec_pval[i] <- wilcox.test(g1, g2)$p.value
}
egfr_res <- data.frame(uni_stat, egfr_ini_pval, egfr_rec_pval)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/egfr_status_ini_rec_box.pdf",width=4,height=2)
ggplot(plot_egfr_class ,aes(x = egfr_status, y = immune)) +
geom_boxplot(outlier.size=0,colour="black") +
#geom_line(size=0.8,alpha=0.4,aes(group=case_barcode)) +
geom_quasirandom(size=1,aes(colour=transcriptional_subtype)) +
scale_colour_manual(values=c("#008a22","#8a0000","#00458a")) +
labs(title="EGFR status", y = "Immune score") +
# annotate("text", x=1.1, y = 0, 
# 	label=deparse((bquote(italic("P") ~" = " ~ .(pval1)))), parse=TRUE, size=2.5) +
facet_grid(.~timepoint) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_blank(),
legend.position="none")
dev.off()

###################################################
# Step 4: Examine how each alteration changes immune in a given patient
##################################################

ini_egfr_class <- plot_egfr_class[which(plot_egfr_class[,"timepoint"]=="Initial"),]
rec_egfr_class <- plot_egfr_class[which(plot_egfr_class[,"timepoint"]=="Recurrent"),]

change_egfr_class <- data.frame(ini_egfr_class[,"case_barcode"],
								ini_egfr_class[,"aliquot_barcode"],
								rec_egfr_class[,"aliquot_barcode"],
								ini_egfr_class[,"transcriptional_subtype"],
								rec_egfr_class[,"transcriptional_subtype"],
								ini_egfr_class[,"immune"],
								rec_egfr_class[,"immune"],
								ini_egfr_class[,"egfr_status"],
								rec_egfr_class[,"egfr_status"])
colnames(change_egfr_class) <- c("case_barcode","tumor_barcode_a","tumor_barcode_b",
								 "transcriptional_subtype_a","transcriptional_subtype_b",
								 "immune_a","immune_b","egfr_status_a","egfr_status_b")

# Testing changes in EGFRvIII status

g1 <- change_egfr_class[which(change_egfr_class[,"egfr_status_a"] == "EGFRvIII" & change_egfr_class[,"egfr_status_b"] != "EGFRvIII"),"immune_a"]
g2 <- change_egfr_class[which(change_egfr_class[,"egfr_status_a"] == "EGFRvIII" & change_egfr_class[,"egfr_status_b"] != "EGFRvIII"),"immune_b"]
wilcox.test(g1,g2,paired=TRUE)		#0.875


g1 <- change_egfr_class[which(change_egfr_class[,"egfr_status_a"] != "EGFRvIII" & change_egfr_class[,"egfr_status_b"] == "EGFRvIII"),"immune_a"]
g2 <- change_egfr_class[which(change_egfr_class[,"egfr_status_a"] != "EGFRvIII" & change_egfr_class[,"egfr_status_b"] == "EGFRvIII"),"immune_b"]
wilcox.test(g1,g2,paired=TRUE)		#1.00