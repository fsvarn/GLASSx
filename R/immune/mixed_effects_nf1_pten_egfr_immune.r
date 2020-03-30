##################################################
# Mixed effect models that link immune infiltrate and purity to different driver alterations
# Adjusts for different mutations
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
SELECT ds.*,im.signature_name, im.enrichment_score
FROM analysis.drivers_by_aliquot ds
JOIN analysis.analyte_sets an ON an.dna_barcode = ds.aliquot_barcode
JOIN analysis.davoli_immune_score im ON im.aliquot_barcode = an.rna_barcode
ORDER BY ds.aliquot_barcode, signature_name
"

dat <- dbGetQuery(con,q)

###################################################
# Step 1: Model of immune infiltration with each significant driver
##################################################

# Create driver mutation variables:
dat[,"nf1_status"] <- as.numeric(grepl("NF1", dat[,"snv_driver"]))
dat[,"pten_status"] <- as.numeric(grepl("PTEN", dat[,"snv_driver"]))
dat[,"egfr_status"] <- as.numeric(grepl("EGFR", dat[,"cnv_driver"]))

# Set idh_codel_subtype variable to factor so that IDHwt is the reference in the model
dat[,"idh_codel_subtype"] <- factor(dat[,"idh_codel_subtype"], levels = c("IDHwt","IDHmut-noncodel","IDHmut-codel"))

# Testing how specific cell infiltration changes in response to the presence/absence of different genes
cells <- unique(dat[,"signature_name"])

# Results matrix, one row for each variable (IDH/codel subtype, NF1, PTEN, and EGFR)
results <- matrix(0, nrow = 4, ncol = length(cells))
rownames(results) <- c("idh_codel_subtype", "NF1", "PTEN","EGFR_amp")
colnames(results) <- cells

for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	lmm <- lmer(enrichment_score ~  + idh_codel_subtype + nf1_status + pten_status + egfr_status + (1|case_barcode), data = sub_dat)
	eff <- summary(lmm)[["coefficients"]][c(2,4:6),3]
	p.val <- Anova(lmm)[["Pr(>Chisq)"]]
	
	p.val <- ifelse(eff < 0, p.val*-1, p.val)
	results[,i] <- p.val
}


