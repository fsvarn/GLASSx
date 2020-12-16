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
WHERE ds.idh_codel_subtype LIKE 'IDHmut%'
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

# Note: For IDHmut need to group by function (cell cycle) due to low numbers
# mygene <- cnv_driv[which(cnv_driv[,"pathway"]=="Cell cycle"),"gene_symbol"]
# string <- paste(mygene,collapse="|")
# 
# # Compare samples gaining an cnv at recurrence
# gain_dat <- dat %>%
# 	filter(grepl(string, cnv_driver_change_a) & !grepl(string, cnv_driver_change_b))
# 
# if(length(unique(gain_dat[,"case_barcode"])) < 3){
# 	next}
# 
# gain_pval <- rep(0, length(cell_states))
# gain_eff <- rep(0, length(cell_states))
# for(j in 1:length(cell_states))
# {
# 	sub_gain <- gain_dat %>% 
# 				filter(cell_state == cell_states[j])
# 	gain_pval[j] <- t.test(sub_gain[,"fraction_a"], sub_gain[,"fraction_b"], paired=TRUE)$p.value
# 	gain_eff[j] <- mean(sub_gain[,"fraction_b"] - sub_gain[,"fraction_a"])
# }
# gain_cnv<- data.frame(cell_states, gain_pval, gain_eff, nrow(sub_gain), "cell cycle", "gain")
# colnames(gain_cnv) <- c("cell_state", "p.value", "effect", "n", "gene", "change")
# 
# # Compare samples losing an cnv at recurrence
# loss_dat <- dat %>%
# 	filter(grepl(string, cnv_driver_change_b) & !grepl(string, cnv_driver_change_a))
# 
# if(length(unique(loss_dat[,"case_barcode"])) < 3){
# 	next}
# 
# loss_pval <- rep(0, length(cell_states))
# loss_eff <- rep(0, length(cell_states))
# for(j in 1:length(cell_states))
# {
# 	sub_loss <- loss_dat %>% 
# 				filter(cell_state == cell_states[j])
# 	loss_pval[j] <- t.test(sub_loss[,"fraction_a"], sub_loss[,"fraction_b"], paired=TRUE)$p.value
# 	loss_eff[j] <- mean(sub_loss[,"fraction_b"] - sub_loss[,"fraction_a"])
# }
# loss_cnv<- data.frame(cell_states, loss_pval, loss_eff, nrow(sub_loss), string, "loss")
# colnames(loss_cnv) <- c("cell_state", "p.value", "effect", "n", "gene", "change")

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

#                cell_state     p.value        effect  n   gene change variant	idh_codel_subtype
# 31            granulocyte 0.004364083  0.0342400238  5    NF1   gain     snv				IDHwt
# 35             fibroblast 0.038543585  0.0066686847  5    NF1   gain     snv				IDHwt
# 65        oligodendrocyte 0.045634417  0.0361527787 15 CDKN2A   gain     cnv				IDHwt
# 82         dendritic_cell 0.014716762 -0.0189109810  3   PTEN   gain     cnv				IDHwt
# 111  differentiated_tumor 0.022947092 -0.2280723036  7 CDKN2A   loss     cnv				IDHwt
# 122               myeloid 0.004776607  0.1230122873 10   EGFR   loss     cnv				IDHwt
# 123  differentiated_tumor 0.026351548 -0.1523071242 10   EGFR   loss     cnv				IDHwt
# 147  differentiated_tumor 0.016812167 -0.1123352288  4 PIK3CA   gain     snv			   IDHmut
# 152              pericyte 0.048638383 -0.0007791239  4 PIK3CA   gain     snv			   IDHmut
# 172 prolif_stemcell_tumor 0.010897768  0.1151759197  4 CDKN2A   gain     cnv			   IDHmut

# 5 plots: NF1 IDHwt, CDKN2A IDHwt, EGFR IDHwt, PIK3CA IDHmut, CDKN2A IDHmut

#######################################################
# Step 3: Plot the results: Box plots and then average barplots
#######################################################

# NF1 plot
#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, snv_driver_change_a, snv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND snv_driver_change_a LIKE '%NF1%' AND (snv_driver_change_b NOT LIKE '%NF1%' OR snv_driver_change_b IS NULL)
"

dat <- dbGetQuery(con,q)

plot_dat <- dat %>% 
			pivot_longer(cols = starts_with("fraction"),
						 names_to = "timepoint",
						 names_prefix = "fraction_",
						 values_to = "fraction") %>%
			mutate(timepoint = recode(timepoint, "a" = "Init.", "b" = "Rec.")) %>%
			mutate(fraction = fraction * 100) %>%
			mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
					"differentiated_tumor" = "Differentiated tumor", "endothelial" = "Endothelial",
					"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
					"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
					"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Proliferating stem cell tumor",
					"stemcell_tumor" = "Stem cell tumor","t_cell" = "T cell")) %>%
			mutate(cell_state = as_factor(cell_state)) %>%
			mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
											"Oligodendrocyte", 
											"Endothelial", "Pericyte",
											"Fibroblast", 
											"Differentiated tumor", "Stem cell tumor", "Proliferating stem cell tumor"))
# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Granulocyte" | cell_state == "Fibroblast"),
	   aes(x = timepoint, y = fraction, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~cell_state) +
scale_colour_manual(values=c("gray","tomato3")) +
labs(y = "Proportion (%)") +
theme_classic() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") + 
coord_cartesian(ylim=c(0,10))

plot_bars <- plot_dat %>%
			 group_by(cell_state, timepoint) %>%
			 summarise(mean = mean(fraction))

bars <- ggplot(plot_bars, aes(x=timepoint, y = mean, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
theme(axis.text = element_text(size=7),
	axis.title = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/NF1_cell_state_changes.pdf",width=2.5,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(2, 1))
dev.off()


# EGFR plot
#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND cnv_driver_change_b LIKE '%EGFR%' AND (cnv_driver_change_a NOT LIKE '%EGFR%' OR cnv_driver_change_a IS NULL)
"

dat <- dbGetQuery(con,q)

plot_dat <- dat %>% 
			pivot_longer(cols = starts_with("fraction"),
						 names_to = "timepoint",
						 names_prefix = "fraction_",
						 values_to = "fraction") %>%
			mutate(timepoint = recode(timepoint, "a" = "Init.", "b" = "Rec.")) %>%
			mutate(fraction = fraction * 100) %>%
			mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
					"differentiated_tumor" = "Differentiated tumor", "endothelial" = "Endothelial",
					"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
					"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
					"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Proliferating stem cell tumor",
					"stemcell_tumor" = "Stem cell tumor","t_cell" = "T cell")) %>%
			mutate(cell_state = as_factor(cell_state)) %>%
			mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
											"Oligodendrocyte", 
											"Endothelial", "Pericyte",
											"Fibroblast", 
											"Differentiated tumor", "Stem cell tumor", "Proliferating stem cell tumor"))
# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Differentiated tumor" | cell_state == "Myeloid"),
	   aes(x = timepoint, y = fraction, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~cell_state) +
scale_colour_manual(values=c("gray","royalblue4")) +
labs(y = "Proportion (%)") +
theme_classic() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") + 
coord_cartesian(ylim=c(0,100))

plot_bars <- plot_dat %>%
			 group_by(cell_state, timepoint) %>%
			 summarise(mean = mean(fraction))

bars <- ggplot(plot_bars, aes(x=timepoint, y = mean, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
theme(axis.text = element_text(size=7),
	axis.title = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/egfr_cell_state_changes.pdf",width=2.5,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(2, 1))
dev.off()


# PIK3CA plot
#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, snv_driver_change_a, snv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_snv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%' AND snv_driver_change_a LIKE '%PIK3CA%' AND (snv_driver_change_b NOT LIKE '%PIK3CA%' OR snv_driver_change_b IS NULL)
"

dat <- dbGetQuery(con,q)

plot_dat <- dat %>% 
			pivot_longer(cols = starts_with("fraction"),
						 names_to = "timepoint",
						 names_prefix = "fraction_",
						 values_to = "fraction") %>%
			mutate(timepoint = recode(timepoint, "a" = "Init.", "b" = "Rec.")) %>%
			mutate(fraction = fraction * 100) %>%
			mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
					"differentiated_tumor" = "Differentiated tumor", "endothelial" = "Endothelial",
					"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
					"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
					"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Proliferating stem cell tumor",
					"stemcell_tumor" = "Stem cell tumor","t_cell" = "T cell")) %>%
			mutate(cell_state = as_factor(cell_state)) %>%
			mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
											"Oligodendrocyte", 
											"Endothelial", "Pericyte",
											"Fibroblast", 
											"Differentiated tumor", "Stem cell tumor", "Proliferating stem cell tumor"))
# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Differentiated tumor"),
	   aes(x = timepoint, y = fraction, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~cell_state) +
scale_colour_manual(values=c("gray","tomato3")) +
labs(y = "Proportion (%)") +
theme_classic() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") + 
coord_cartesian(ylim=c(0,100))

plot_bars <- plot_dat %>%
			 group_by(cell_state, timepoint) %>%
			 summarise(mean = mean(fraction))

bars <- ggplot(plot_bars, aes(x=timepoint, y = mean, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
theme(axis.text = element_text(size=7),
	axis.title = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/PIK3CA_cell_state_changes_idhmut.pdf",width=1.9,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()

# CDKN2A plot
#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%' AND cnv_driver_change_a LIKE '%CDKN2A%' AND (cnv_driver_change_b NOT LIKE '%CDKN2A%' OR cnv_driver_change_b IS NULL)
"

dat <- dbGetQuery(con,q)

plot_dat <- dat %>% 
			pivot_longer(cols = starts_with("fraction"),
						 names_to = "timepoint",
						 names_prefix = "fraction_",
						 values_to = "fraction") %>%
			mutate(timepoint = recode(timepoint, "a" = "Init.", "b" = "Rec.")) %>%
			mutate(fraction = fraction * 100) %>%
			mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
					"differentiated_tumor" = "Differentiated tumor", "endothelial" = "Endothelial",
					"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
					"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
					"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Proliferating stem cell tumor",
					"stemcell_tumor" = "Stem cell tumor","t_cell" = "T cell")) %>%
			mutate(cell_state = as_factor(cell_state)) %>%
			mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
											"Oligodendrocyte", 
											"Endothelial", "Pericyte",
											"Fibroblast", 
											"Differentiated tumor", "Stem cell tumor", "Proliferating stem cell tumor"))
# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Proliferating stem cell tumor"),
	   aes(x = timepoint, y = fraction, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~cell_state) +
scale_colour_manual(values=c("gray","tomato3")) +
labs(y = "Proportion (%)") +
theme_classic() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") + 
coord_cartesian(ylim=c(0,25))

plot_bars <- plot_dat %>%
			 group_by(cell_state, timepoint) %>%
			 summarise(mean = mean(fraction))

bars <- ggplot(plot_bars, aes(x=timepoint, y = mean, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
theme(axis.text = element_text(size=7),
	axis.title = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CDKN2A_cell_state_changes_idhmut.pdf",width=1.9,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()