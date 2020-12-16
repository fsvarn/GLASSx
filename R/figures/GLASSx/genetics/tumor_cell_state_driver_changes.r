# Code examining whether acquisition or loss of specific drivers results in changes in tumor cell steady states
# Examines both SCGP and Neftel cell states, does not find anything of note
# Date: 11/25/2020

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
SELECT ps.*, ds.idh_codel_subtype, snv_driver_change_a, snv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(grepl("tumor", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), snv_driver_change_a, snv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()

driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

snv_driv <- driv %>% filter(has_mut==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gain_snv <- loss_snv <- list()
for(i in 1:nrow(snv_driv))
{
	mygene <- snv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_a) & !grepl(mygene, snv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_snv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an SNV at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_b) & !grepl(mygene, snv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_snv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_snv <- do.call(rbind, gain_snv)
loss_snv <- do.call(rbind, loss_snv)
snv_res <- rbind(gain_snv, loss_snv)
snv_res[which(snv_res[,"p.value"] < 0.05),]


# 1B: Copy number variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(grepl("tumor", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), cnv_driver_change_a, cnv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()

driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

cnv_driv <- driv %>% filter(has_cnv==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])


gain_cnv <- loss_cnv <- list()
for(i in 1:nrow(cnv_driv))
{
	mygene <- cnv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an cnv at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_a) & !grepl(mygene, cnv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_cnv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an cnv at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_b) & !grepl(mygene, cnv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_cnv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_cnv <- do.call(rbind, gain_cnv)
loss_cnv <- do.call(rbind, loss_cnv)
cnv_res <- rbind(gain_cnv, loss_cnv)
cnv_res[which(cnv_res[,"p.value"] < 0.05),]

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
SELECT ps.*, ds.idh_codel_subtype, snv_driver_change_a, snv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(grepl("tumor", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), snv_driver_change_a, snv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()
		
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

snv_driv <- driv %>% filter(has_mut==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gain_snv <- loss_snv <- list()
for(i in 1:nrow(snv_driv))
{
	mygene <- snv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_a) & !grepl(mygene, snv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_snv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an SNV at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_b) & !grepl(mygene, snv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_snv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_snv <- do.call(rbind, gain_snv)
loss_snv <- do.call(rbind, loss_snv)
snv_res <- rbind(gain_snv, loss_snv)
snv_res[which(snv_res[,"p.value"] < 0.05),]
 
# 2B: Copy number variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(grepl("tumor", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), cnv_driver_change_a, cnv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()
		
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

cnv_driv <- driv %>% filter(has_cnv==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gain_cnv <- loss_cnv <- list()
for(i in 1:nrow(cnv_driv))
{
	mygene <- cnv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an cnv at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_a) & !grepl(mygene, cnv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_cnv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an cnv at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_b) & !grepl(mygene, cnv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_cnv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_cnv <- do.call(rbind, gain_cnv)
loss_cnv <- do.call(rbind, loss_cnv)
cnv_res <- rbind(gain_cnv, loss_cnv)
cnv_res[which(cnv_res[,"p.value"] < 0.05),]

snv_res[,"variant"] <- "snv"
cnv_res[,"variant"] <- "cnv"
idhmut_res <- rbind(snv_res, cnv_res)

idhwt_res[,"idh_codel_subtype"] <- "IDHwt"
idhmut_res[,"idh_codel_subtype"] <- "IDHmut"
full_res <- rbind(idhwt_res, idhmut_res)
full_res <- full_res[which(full_res[,"p.value"] < 0.05),]




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
SELECT ps.*, ds.idh_codel_subtype, snv_driver_change_a, snv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_neftel cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_neftel cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(cell_state %in% c("AClike","OPClike","Mesenchymal","NPClike")) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), snv_driver_change_a, snv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()

driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

snv_driv <- driv %>% filter(has_mut==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gain_snv <- loss_snv <- list()
for(i in 1:nrow(snv_driv))
{
	mygene <- snv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_a) & !grepl(mygene, snv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_snv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an SNV at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_b) & !grepl(mygene, snv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_snv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_snv <- do.call(rbind, gain_snv)
loss_snv <- do.call(rbind, loss_snv)
snv_res <- rbind(gain_snv, loss_snv)
snv_res[which(snv_res[,"p.value"] < 0.05),]


# 1B: Copy number variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_neftel cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_neftel cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(cell_state %in% c("AClike","OPClike","Mesenchymal","NPClike")) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), cnv_driver_change_a, cnv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()

driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

cnv_driv <- driv %>% filter(has_cnv==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])


gain_cnv <- loss_cnv <- list()
for(i in 1:nrow(cnv_driv))
{
	mygene <- cnv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an cnv at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_a) & !grepl(mygene, cnv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_cnv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an cnv at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_b) & !grepl(mygene, cnv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_cnv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_cnv <- do.call(rbind, gain_cnv)
loss_cnv <- do.call(rbind, loss_cnv)
cnv_res <- rbind(gain_cnv, loss_cnv)
cnv_res[which(cnv_res[,"p.value"] < 0.05),]

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
SELECT ps.*, ds.idh_codel_subtype, snv_driver_change_a, snv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_neftel cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_neftel cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(cell_state %in% c("AClike","OPClike","Mesenchymal","NPClike")) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), snv_driver_change_a, snv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()
		
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

snv_driv <- driv %>% filter(has_mut==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gain_snv <- loss_snv <- list()
for(i in 1:nrow(snv_driv))
{
	mygene <- snv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an SNV at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_a) & !grepl(mygene, snv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_snv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an SNV at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, snv_driver_change_b) & !grepl(mygene, snv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_snv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_snv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_snv <- do.call(rbind, gain_snv)
loss_snv <- do.call(rbind, loss_snv)
snv_res <- rbind(gain_snv, loss_snv)
snv_res[which(snv_res[,"p.value"] < 0.05),]
 
# 2B: Copy number variants:
#-----------------------------------

#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_neftel cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_neftel cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%'
"

dat <- dbGetQuery(con,q)

dat <- dat %>% 
		filter(cell_state %in% c("AClike","OPClike","Mesenchymal","NPClike")) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), cnv_driver_change_a, cnv_driver_change_b, idh_codel_subtype) %>%
		as.data.frame()
		
driv <- dbReadTable(con, Id(schema="ref", table="driver_genes"))

cnv_driv <- driv %>% filter(has_cnv==1 & gene_symbol != 'IDH1')
cell_states <- unique(dat[,"cell_state"])

gain_cnv <- loss_cnv <- list()
for(i in 1:nrow(cnv_driv))
{
	mygene <- cnv_driv[i,"gene_symbol"]
	
	# Compare samples gaining an cnv at recurrence
	gain_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_a) & !grepl(mygene, cnv_driver_change_b))
	
	if(length(unique(gain_dat[,"case_barcode"])) < 3){
		next}
	
	gain_pval <- rep(0, length(cell_states))
	gain_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_gain <- gain_dat %>% 
					filter(cell_state == cell_states[j])
		gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
		gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
	}
	gain_cnv[[i]]<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), mygene, "gain")
	colnames(gain_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
	
	# Compare samples losing an cnv at recurrence
	loss_dat <- dat %>%
		filter(grepl(mygene, cnv_driver_change_b) & !grepl(mygene, cnv_driver_change_a))
	
	if(length(unique(loss_dat[,"case_barcode"])) < 3){
		next}
	
	loss_pval <- rep(0, length(cell_states))
	loss_eff <- rep(0, length(cell_states))
	for(j in 1:length(cell_states))
	{
		sub_loss <- loss_dat %>% 
					filter(cell_state == cell_states[j])
		loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
		loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
	}
	loss_cnv[[i]]<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), mygene, "loss")
	colnames(loss_cnv[[i]]) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
}

gain_cnv <- do.call(rbind, gain_cnv)
loss_cnv <- do.call(rbind, loss_cnv)
cnv_res <- rbind(gain_cnv, loss_cnv)
cnv_res[which(cnv_res[,"p.value"] < 0.05),]

snv_res[,"variant"] <- "snv"
cnv_res[,"variant"] <- "cnv"
idhmut_res <- rbind(snv_res, cnv_res)

idhwt_res[,"idh_codel_subtype"] <- "IDHwt"
idhmut_res[,"idh_codel_subtype"] <- "IDHmut"
full_res <- rbind(idhwt_res, idhmut_res)
full_res <- full_res[which(full_res[,"p.value"] < 0.05),]

