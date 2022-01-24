##################################################
# Identify cell state associations for acquired/lost driver genes
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
SELECT ps.*, ds.idh_codel_subtype, snv_driver_change_a, snv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
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
snv_res[which(snv_res[,"p.value"] < 0.1),]
 
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
JOIN analysis.mut_freq mf ON mf.aliquot_barcode = ps.dna_barcode_b
WHERE ds.idh_codel_subtype LIKE 'IDHmut%' AND coverage_adj_mut_freq <= 10
"

dat <- dbGetQuery(con,q)
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
full_res[which(full_res[,"p.value"] < 0.05),]

#                cell_state     p.value        effect  n   gene change variant
# 43            granulocyte 0.027238001  0.0146687407  3    NF1   gain     snv
# 48                 b_cell 0.008688277 -0.0004779909  3    NF1   gain     snv
# 80               pericyte 0.030424886 -0.0044655232  3    NF1   loss     snv
# 122               myeloid 0.037545685 -0.0733913550  4 PDGFRA   gain     cnv
# 126                t_cell 0.024848150 -0.0115399294  4 PDGFRA   gain     cnv
# 135  differentiated_tumor 0.012801403 -0.2218988492  7 CDKN2A   loss     cnv
# 145        stemcell_tumor 0.024128620 -0.0610513645 11   EGFR   loss     cnv
# 156                b_cell 0.013503486 -0.0002648105 11   EGFR   loss     cnv
# 157        stemcell_tumor 0.031054309 -0.1303500652  4 PDGFRA   loss     cnv
# 168                b_cell 0.035852384 -0.0003527813  4 PDGFRA   loss     cnv
# 172 prolif_stemcell_tumor 0.039199162  0.1391666736  4 CDKN2A   gain     cnv
# 184 prolif_stemcell_tumor 0.024972439  0.1534343835  3  CCND2   gain     cnv
# 196 prolif_stemcell_tumor 0.047721529  0.1671272285  3 PDGFRA   gain     cnv
# 201           endothelial 0.041306407  0.0054980399  3 PDGFRA   gain     cnv
# 
#     idh_codel_subtype
# 43              IDHwt
# 48              IDHwt
# 80              IDHwt
# 122             IDHwt
# 126             IDHwt
# 135             IDHwt
# 145             IDHwt
# 156             IDHwt
# 157             IDHwt
# 168             IDHwt
# 172            IDHmut
# 184            IDHmut
# 196            IDHmut
# 201            IDHmut
# 

# 5 plots: NF1 IDHwt, CDKN2A IDHwt, EGFR IDHwt, PIK3CA IDHmut, CDKN2A IDHmut
