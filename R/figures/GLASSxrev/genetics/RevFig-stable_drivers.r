##################################################
# Test cell state changes in longitudinal pairs with stable driver genes
# Author: Frederick Varn
# Date: 2021.10.26
# Revision comment
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)

#######################################################
rm(list=ls())
#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#######################################################
# Step 1: Begin analysis in IDHwt tumors:
#######################################################

# 1A: Single nucleotide variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, snv_driver_shared, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
JOIN analysis.mut_freq mf ON mf.aliquot_barcode = ps.dna_barcode_b
WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND coverage_adj_mut_freq <= 10
"

dat <- dbGetQuery(con,q)
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

snv_driv <- driv %>% filter(has_mut==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gene_snv <- list()
for(i in 1:nrow(snv_driv))
{
	mygene <- snv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gene_dat <- dat %>%
		filter(grepl(mygene, snv_driver_shared))
	
	if(length(unique(gene_dat[,"case_barcode"])) < 3){
		next}
	
	gene_pval <- rep(0, length(cell_states))
	gene_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_dat <- gene_dat %>% 
					filter(cell_state == cell_states[j])
		gene_pval[j] <- t.test(sub_dat[,"fraction_a"], sub_dat[,"fraction_b"], paired=TRUE)$p.value
		gene_eff[j] <- mean(sub_dat[,"fraction_b"] - sub_dat[,"fraction_a"])
	}
	gene_snv[[i]]<- data.frame(cell_states, gene_pval, gene_eff, nrow(sub_dat), mygene)
	colnames(gene_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene")
}

snv_res <- do.call(rbind, gene_snv)
snv_res[which(snv_res[,"p.value"] < 0.05),]
 
# 1B: Copy number variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_shared, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt'
"

dat <- dbGetQuery(con,q)
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

cnv_driv <- driv %>% filter(has_cnv==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])


gene_cnv <- list()
for(i in 1:nrow(cnv_driv))
{
	mygene <- cnv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gene_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_shared))
	
	if(length(unique(gene_dat[,"case_barcode"])) < 3){
		next}
	
	gene_pval <- rep(0, length(cell_states))
	gene_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_dat <- gene_dat %>% 
					filter(cell_state == cell_states[j])
		gene_pval[j] <- t.test(sub_dat[,"fraction_a"], sub_dat[,"fraction_b"], paired=TRUE)$p.value
		gene_eff[j] <- mean(sub_dat[,"fraction_b"] - sub_dat[,"fraction_a"])
	}
	gene_cnv[[i]]<- data.frame(cell_states, gene_pval, gene_eff, nrow(sub_dat), mygene)
	colnames(gene_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene")
}

cnv_res <- do.call(rbind, gene_cnv)

snv_res[,"variant"] <- "snv"
cnv_res[,"variant"] <- "cnv"
idhwt_res <- rbind(snv_res, cnv_res)


#######################################################
# Step 2: Examine in IDHmut tumors:
#######################################################

# 2A: Single nucleotide variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, snv_driver_shared, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
JOIN analysis.mut_freq mf ON mf.aliquot_barcode = ps.dna_barcode_b
WHERE ds.idh_codel_subtype LIKE 'IDHmut%' AND coverage_adj_mut_freq <= 10
"

dat <- dbGetQuery(con,q)
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

snv_driv <- driv %>% filter(has_mut==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gene_snv <- list()
for(i in 1:nrow(snv_driv))
{
	mygene <- snv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gene_dat <- dat %>%
		filter(grepl(mygene, snv_driver_shared))
	
	if(length(unique(gene_dat[,"case_barcode"])) < 3){
		next}
	
	gene_pval <- rep(0, length(cell_states))
	gene_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_dat <- gene_dat %>% 
					filter(cell_state == cell_states[j])
		gene_pval[j] <- t.test(sub_dat[,"fraction_a"], sub_dat[,"fraction_b"], paired=TRUE)$p.value
		gene_eff[j] <- mean(sub_dat[,"fraction_b"] - sub_dat[,"fraction_a"])
	}
	gene_snv[[i]]<- data.frame(cell_states, gene_pval, gene_eff, nrow(sub_dat), mygene)
	colnames(gene_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene")
}

snv_res <- do.call(rbind, gene_snv)
snv_res[which(snv_res[,"p.value"] < 0.05),]
 
# 2B: Copy number variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_shared, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%'
"

dat <- dbGetQuery(con,q)
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

cnv_driv <- driv %>% filter(has_cnv==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])


gene_cnv <- list()
for(i in 1:nrow(cnv_driv))
{
	mygene <- cnv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gene_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_shared))
	
	if(length(unique(gene_dat[,"case_barcode"])) < 3){
		next}
	
	gene_pval <- rep(0, length(cell_states))
	gene_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_dat <- gene_dat %>% 
					filter(cell_state == cell_states[j])
		gene_pval[j] <- t.test(sub_dat[,"fraction_a"], sub_dat[,"fraction_b"], paired=TRUE)$p.value
		gene_eff[j] <- mean(sub_dat[,"fraction_b"] - sub_dat[,"fraction_a"])
	}
	gene_cnv[[i]]<- data.frame(cell_states, gene_pval, gene_eff, nrow(sub_dat), mygene)
	colnames(gene_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene")
}

cnv_res <- do.call(rbind, gene_cnv)

snv_res[,"variant"] <- "snv"
cnv_res[,"variant"] <- "cnv"
idhmut_res <- rbind(snv_res, cnv_res)


idhwt_res[,"idh_codel_subtype"] <- "IDHwt"
idhmut_res[,"idh_codel_subtype"] <- "IDHmut"
full_res <- rbind(idhwt_res, idhmut_res)
full_res[which(full_res[,"p.value"] < 0.05),]

#              cell_state     p.value        effect  n   gene variant
# 34        dendritic_cell 0.009714140  0.0035860094 21   PTEN     snv
# 63  differentiated_tumor 0.001464776 -0.1598181249 14    NF1     snv
# 87  differentiated_tumor 0.037601949 -0.0637262402 21 CDKN2A     cnv
# 111 differentiated_tumor 0.047295516  0.0794028367  5   CDK4     cnv
# 125      oligodendrocyte 0.024797132  0.0663622091 21   EGFR     cnv
# 179           fibroblast 0.005489971  0.0008393045  3 CDKN2A     cnv
#     idh_codel_subtype
# 34              IDHwt
# 63              IDHwt
# 87              IDHwt
# 111             IDHwt
# 125             IDHwt
# 179            IDHmut
# 
