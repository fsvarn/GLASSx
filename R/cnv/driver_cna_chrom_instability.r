##################################################
# Exploratory script looking at how different alterations associate with chromosomal instability
# Additionally examines how these features associate with immune infiltrate
# Focus on EGFR and NF1 which are the main alterations associated with immune cells and subtypes
# Updated: 2020.03.30
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(lme4)
library(car)
library(reshape)
library(grid)

##################################################
rm(list=ls())
# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
SELECT da.*, an.prop_aneuploidy
FROM analysis.drivers_by_aliquot da
JOIN analysis.gatk_aneuploidy an ON an.aliquot_barcode = da.aliquot_barcode
ORDER BY 1
"

dat <- dbGetQuery(con,q)

###################################################
# Step 1: Run mixed effect modeling for chromosomal instability vs genomic alteration (most interested in CNV)
##################################################

# Start with copy number alterations because that is what we are more interested in:

# Get list of unique driver alterations
cnv_drivers <- unique(unlist(strsplit(dat[,"cnv_driver"],", ")))
cnv_drivers <- unique(cnv_drivers)
cnv_drivers <- cnv_drivers[-which(is.na(cnv_drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
cnv_genes <- cnv_drivers

# Create results matrix
cnv_results <- rep(NA, length(cnv_genes))
names(cnv_results) <- cnv_genes

cnv_num <- rep(0,length(cnv_genes))
names(cnv_num) <- cnv_genes

for(i in 1:length(cnv_genes))
{
	#Create a driver_status vector 
	driver_status <- rep(0, nrow(dat))
	
	# Identify all samples that have a mutation in this gene
	indices <- grep(cnv_genes[i], dat[,"cnv_driver"])	
	
	# Change driver_status vector values to indicate driver is present at initial tumor
	driver_status[indices] <- 1
	if(sum(driver_status)==0){
		next}
	test_dat <- cbind(dat, driver_status)
	
	lmm <- lmer(prop_aneuploidy ~ driver_status + idh_codel_subtype + (1|case_barcode), data = test_dat)
	eff <- summary(lmm)[["coefficients"]][2]
	p.val <- Anova(lmm)[["Pr(>Chisq)"]][1]
	
	p.val <- ifelse(eff < 0, p.val*-1, p.val)
	cnv_results[i] <- p.val
	cnv_num[i] <- sum(driver_status)
	
}

#     PTEN del   CDKN2A del     EGFR amp     CDK6 amp      MET amp   PDGFRA amp 
# 4.279384e-02 2.013926e-05 2.515659e-03 3.284160e-04 2.799733e-06 5.187614e-05 
#     MDM4 amp     CDK4 amp     MDM2 amp     MYCN amp    CCND2 amp     ATRX del 
# 6.877258e-04 7.444719e-09 1.225771e-01 3.319288e-02 1.221680e-01 2.768344e-01 

#---------------------------------

# Now examine SNVs for completeness:

# Get list of unique driver alterations
snv_drivers <- unique(unlist(strsplit(dat[,"snv_driver"],", ")))
snv_drivers <- sapply(strsplit(snv_drivers," "),function(x)x[1])
snv_drivers <- unique(snv_drivers)
snv_drivers <- snv_drivers[-which(is.na(snv_drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
snv_genes <- snv_drivers

# Create results matrix
snv_results <- rep(NA, length(snv_genes))
names(snv_results) <- snv_genes

snv_num <- rep(0,length(snv_genes))
names(snv_num) <- snv_genes

for(i in 1:length(snv_genes))
{
	#Create a driver_status vector 
	driver_status <- rep(0, nrow(dat))
	
	# Identify all samples that have a mutation in this gene
	indices <- grep(snv_genes[i], dat[,"snv_driver"])	
	
	# Change driver_status vector values to indicate driver is present at initial tumor
	driver_status[indices] <- 1
	if(sum(driver_status)==0){
		next}
	test_dat <- cbind(dat, driver_status)
	
	lmm <- lmer(prop_aneuploidy ~ driver_status + idh_codel_subtype + (1|case_barcode), data = test_dat)
	eff <- summary(lmm)[["coefficients"]][2]
	p.val <- Anova(lmm)[["Pr(>Chisq)"]][1]
	
	p.val <- ifelse(eff < 0, p.val*-1, p.val)
	snv_results[i] <- p.val
	snv_num[i] <- sum(driver_status)
	
}

#         NF1        TP53      PIK3CA        PTEN        EGFR         RB1 
# -0.43884492  0.02161066  0.86297301  0.28984289  0.22743127  0.02598989 
#        ATRX       FUBP1      PIK3R1         CIC 
# -0.85538050  0.40948478  0.38321395 -0.39718103


###################################################
# Step 2: Examine how acquired EGFR amplifications associate with chromosomal instability changes
# Also check how this associates with NF1 mutations
##################################################

# Read in data

q <- "
SELECT dc.tumor_pair_barcode, dc.case_barcode, dc.tumor_barcode_a, dc.tumor_barcode_b, dc.idh_codel_subtype, 
dc.cnv_driver_shared, dc.cnv_driver_change_a, dc.cnv_driver_change_b,
sv.snv_driver_shared, sv.snv_driver_change_a, sv.snv_driver_change_b,
an1.prop_aneuploidy AS prop_a, an2.prop_aneuploidy AS prop_b
FROM analysis.driver_status_cnv dc
JOIN analysis.driver_status_snv sv ON dc.tumor_barcode_a = sv.tumor_barcode_a AND dc.tumor_barcode_b = sv.tumor_barcode_b
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = dc.tumor_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = dc.tumor_barcode_b
JOIN analysis.gold_set gs ON gs.tumor_pair_barcode = dc.tumor_pair_barcode
"
dat <- dbGetQuery(con,q)

# Again starting with copy number alterations:

# Create gain results matrix
cnv_gain_results <- rep(NA, length(cnv_genes))
names(cnv_gain_results) <- cnv_genes

cnv_gain_num <- rep(0,length(cnv_genes))
names(cnv_gain_num) <- cnv_genes

# Create loss results matrix
cnv_loss_results <- rep(NA, length(cnv_genes))
names(cnv_loss_results) <- cnv_genes

cnv_loss_num <- rep(0,length(cnv_genes))
names(cnv_loss_num) <- cnv_genes

# nested for loop, one iteration for each cell, one for each mutation

for(i in 1:length(cnv_genes))
{
	gain_dat <- dat[grep(cnv_genes[i], dat[,"cnv_driver_change_a"]),]
	loss_dat <- dat[grep(cnv_genes[i], dat[,"cnv_driver_change_b"]),]
	
	if(nrow(gain_dat) > 0)
	{
		p.val  <- wilcox.test(gain_dat[,"prop_a"], gain_dat[,"prop_b"], paired=TRUE)$p.value
		eff <- median(gain_dat[,"prop_b"]) - median(gain_dat[,"prop_a"])
	
		p.val <- ifelse(eff < 0, p.val*-1, p.val)
		cnv_gain_results[i] <- p.val
		cnv_gain_num[i] <- nrow(gain_dat)
	}

	if(nrow(loss_dat) > 0)
	{	
		p.val  <- wilcox.test(loss_dat[,"prop_a"], loss_dat[,"prop_b"], paired=TRUE)$p.value
		eff <- median(loss_dat[,"prop_b"]) - median(loss_dat[,"prop_a"])
	
		p.val <- ifelse(eff < 0, p.val*-1, p.val)
		cnv_loss_results[i] <- p.val
		cnv_loss_num[i] <- nrow(loss_dat)
	}
}

# gain results:
#     PTEN del   CDKN2A del     EGFR amp     CDK6 amp      MET amp   PDGFRA amp 
# -0.750000000  0.000625364  0.054687500  0.500000000  0.109375000  0.083251953 
#     MDM4 amp     CDK4 amp     MDM2 amp     MYCN amp    CCND2 amp     ATRX del 
#           NA  0.062500000  1.000000000  0.250000000  0.218750000           NA 

# loss_results:
#   PTEN del CDKN2A del   EGFR amp   CDK6 amp    MET amp PDGFRA amp   MDM4 amp 
# -0.6875000 -0.4422989 -0.1013974 -0.1230469  0.3125000 -0.0078125 -0.1250000 
#   CDK4 amp   MDM2 amp   MYCN amp  CCND2 amp   ATRX del 
# -0.0312500 -0.5000000  1.0000000  0.5000000  1.0000000 

#---------------------------------

# Repeating with SNVs:

# Create gain results matrix
snv_gain_results <- rep(NA, length(snv_genes))
names(snv_gain_results) <- snv_genes

snv_gain_num <- rep(0,length(snv_genes))
names(snv_gain_num) <- snv_genes

# Create loss results matrix
snv_loss_results <- rep(NA, length(snv_genes))
names(snv_loss_results) <- snv_genes

snv_loss_num <- rep(0,length(snv_genes))
names(snv_loss_num) <- snv_genes

# nested for loop, one iteration for each cell, one for each mutation

for(i in 1:length(snv_genes))
{
	gain_dat <- dat[grep(snv_genes[i], dat[,"snv_driver_change_a"]),]
	loss_dat <- dat[grep(snv_genes[i], dat[,"snv_driver_change_b"]),]
	
	if(nrow(gain_dat) > 0)
	{
		p.val  <- wilcox.test(gain_dat[,"prop_a"], gain_dat[,"prop_b"], paired=TRUE)$p.value
		eff <- median(gain_dat[,"prop_b"]) - median(gain_dat[,"prop_a"])
	
		p.val <- ifelse(eff < 0, p.val*-1, p.val)
		snv_gain_results[i] <- p.val
		snv_gain_num[i] <- nrow(gain_dat)
	}

	if(nrow(loss_dat) > 0)
	{	
		p.val  <- wilcox.test(loss_dat[,"prop_a"], loss_dat[,"prop_b"], paired=TRUE)$p.value
		eff <- median(loss_dat[,"prop_b"]) - median(loss_dat[,"prop_a"])
	
		p.val <- ifelse(eff < 0, p.val*-1, p.val)
		snv_loss_results[i] <- p.val
		snv_loss_num[i] <- nrow(loss_dat)
	}
}

# gain results:
#         NF1        TP53      PIK3CA        PTEN        EGFR         RB1 
# -0.65782738  0.51708984  0.04882812  0.24349976 -0.74665833  0.29687500 
#        ATRX       FUBP1      PIK3R1         CIC 
# -0.83105469  1.00000000  0.43750000  0.39273834 

# loss_results:
#        NF1       TP53     PIK3CA       PTEN       EGFR        RB1       ATRX 
#  1.0000000  0.2333984  0.6250000 -0.1562500  0.1474609  0.5625000  0.0625000 
#      FUBP1     PIK3R1        CIC 
#  0.5000000  1.0000000  0.3125000 


###################################################
# Step 3: How often do EGFR amplifications and NF1 mutations co-occur?
##################################################

#Initial tumors only since they have not received treatment
q <- "
SELECT *,
CASE WHEN snv_driver LIKE '%NF1%' THEN 1 ELSE 0 END AS nf1,
CASE WHEN cnv_driver LIKE '%EGFR%' THEN 1 ELSE 0 END AS egfr,
CASE WHEN snv_driver LIKE '%NF1%' AND cnv_driver LIKE '%EGFR%' THEN 1 ELSE 0 END AS co_occur
FROM analysis.drivers_by_aliquot da
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = da.aliquot_barcode
"

dat <- dbGetQuery(con, q)

nf1_prob <- sum(dat[,"nf1"])/nrow(dat)			#0.14
egfr_prob <- sum(dat[,"egfr"])/nrow(dat)		#0.32

expected <- nf1_prob * egfr_prob 				#0.04
observed <- sum(dat[,"co_occur"])/nrow(dat)		#0.01

g1 <- sum(dat[,"nf1"])
g2 <- sum(dat[,"egfr"])
g3 <- sum(dat[,"co_occur"])
g4 <- nrow(dat) - sum(c(g1,g2,g3))

p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value #2.3e-9

sum(dat[,"co_occur"])/sum(dat[,"nf1"])	#0.07
sum(dat[,"co_occur"])/sum(dat[,"egfr"])	#0.03

###################################################
# Step 4: What is the proportion of samples with EGFR amplifications that are classical? 
# What happens to the transcriptional subtype  when the mutations co-occur?
##################################################

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
SELECT da.*,
CASE WHEN snv_driver LIKE '%NF1%' THEN 1 ELSE 0 END AS nf1,
CASE WHEN cnv_driver LIKE '%EGFR%' THEN 1 ELSE 0 END AS egfr,
CASE WHEN snv_driver LIKE '%NF1%' AND cnv_driver LIKE '%EGFR%' THEN 1 ELSE 0 END AS co_occur,
ag.subtype
FROM analysis.drivers_by_aliquot da
JOIN analysis.platinum_set ps ON ps.dna_barcode_a = da.aliquot_barcode
JOIN agg ag ON ps.rna_barcode_a = ag.aliquot_barcode
"

dat <- dbGetQuery(con, q)

dat[which(dat[,"co_occur"]==1),"subtype"] #Classical (n = 1)

egfr_only <- dat[which(dat[,"egfr"]==1),]
prop_class <- length(grep("Classical",egfr_only[,"subtype"]))/nrow(egfr_only)	#0.74

