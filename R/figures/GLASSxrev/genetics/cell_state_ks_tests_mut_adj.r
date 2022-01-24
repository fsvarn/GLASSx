###################################################
# Test cell state change distributions for deviation from the normal distribution
# Author: Frederick Varn
# Date: 2021.10.26
# Figure 3B in the initial submission, not in revision; revision comment
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(survival)
library(DescTools)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
SELECT ps.*, CASE WHEN ds.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
snv_driver_change_a, snv_driver_change_b, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state,
cs1.fraction * 100 AS fraction_a,
cs2.fraction * 100 AS fraction_b,
(cs2.fraction - cs1.fraction) * 100 AS fraction_diff,
mf1.coverage_adj_mut_freq AS mut_freq_a,
mf2.coverage_adj_mut_freq AS mut_freq_b,
mf2.coverage_adj_mut_freq > 10 AS hm,
ts1.signature_name AS subtype_a,
ts2.signature_name AS subtype_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = ps.dna_barcode_b
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.rna_barcode_b
"

dat <- dbGetQuery(con, q)

dat <- dat %>% mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligo.",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		mutate(cell_state = fct_relevel(cell_state, "Prolif. stem-like", "Stem-like","Diff.-like",
										"Fibroblast", "Pericyte","Endothelial", "Oligo.",
										"Myeloid", "Dendritic cell", "T cell", "Granulocyte","B cell")) %>%
		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) 


# Analysis for the revision of what happens when we leave out these mutated samples
# sub_dat <- dat %>%
# 		   filter(!grepl("+CDKN2A del", cnv_driver_change_a)) %>%
# 		   filter(!grepl("+CCND2 amp", cnv_driver_change_a)) %>%
# 		   filter(hm == 0)
# 		   
# tumor_dat <- sub_dat %>% 
# 		filter(grepl("-like", cell_state)) %>%
# 		group_by(case_barcode) %>%
# 		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, fraction_diff = (fraction_b - fraction_a), idh_status) %>%
# 	   	mutate(cell_state = fct_relevel(cell_state, c("Prolif. stem-like","Stem-like","Diff.-like"))) %>%
# 		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>% 
# 		as.data.frame()	
# 
# tumor_res <- tumor_dat %>% group_by(cell_state, idh_status) %>% 
# 		summarise(pval = t.test(fraction_a, fraction_b, paired=TRUE)$p.value, n = n()) %>%
# 		as.data.frame()


# Test changes when the noted mutations are gone
sub_dat <- dat %>%
		   filter(!grepl("+NF1", snv_driver_change_a)) %>%
		   filter(!grepl("+EGFR amp", cnv_driver_change_a)) %>%
		   filter(!grepl("+PDGFRA amp", cnv_driver_change_a)) %>%
		   filter(!grepl("+CDKN2A del", cnv_driver_change_a)) %>%
		   filter(!grepl("+CCND2 amp", cnv_driver_change_a)) %>%
		   filter(hm == 0)

full_tumor_dat <- dat %>% 
		filter(grepl("-like", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, fraction_diff = (fraction_b - fraction_a), idh_status) %>%
	   	mutate(cell_state = fct_relevel(cell_state, c("Prolif. stem-like","Stem-like","Diff.-like"))) %>%
		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>% 
		as.data.frame()	


full_non_tumor_dat <- dat %>% 
		filter(!grepl("-like", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, fraction_diff = (fraction_b - fraction_a), idh_status) %>%
	   	mutate(cell_state = fct_relevel(cell_state, c("Prolif. stem-like","Stem-like","Diff.-like"))) %>%
		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>% 
		as.data.frame()		
		
tumor_dat <- sub_dat %>% 
		filter(grepl("-like", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, fraction_diff = (fraction_b - fraction_a), idh_status) %>%
	   	mutate(cell_state = fct_relevel(cell_state, c("Prolif. stem-like","Stem-like","Diff.-like"))) %>%
		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>% 
		as.data.frame()	

non_tumor_dat <- sub_dat %>% 
		filter(!grepl("-like", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, fraction_diff = (fraction_b - fraction_a), idh_status) %>%
	   	mutate(cell_state = fct_relevel(cell_state, c("Prolif. stem-like","Stem-like","Diff.-like"))) %>%
		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>% 
		as.data.frame()		
		
full_res <- dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n(), sd = sd(fraction_diff)) %>%
		as.data.frame()
	
res <- sub_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n(), sd = sd(fraction_diff)) %>%
		as.data.frame()


full_tumor_res <- full_tumor_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n(), sd = sd(fraction_diff)) %>%
		as.data.frame()
		
tumor_res <- tumor_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n(), sd = sd(fraction_diff)) %>%
		as.data.frame()

full_non_tumor_res <- full_non_tumor_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n(), sd = sd(fraction_diff)) %>%
		as.data.frame()
		
non_tumor_res <- non_tumor_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n(), sd = sd(fraction_diff)) %>%
		as.data.frame()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# All samples (main figure)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

p1 <- ggplot(dat, aes(x = fraction_diff, fill=cell_state)) +
geom_density() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligo." = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
facet_grid(cell_state ~ idh_status,scales="free_y",switch="y") +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_blank(),
strip.text.y.left = element_text(size=7,angle = 0,hjust=1),
strip.background = element_blank(),
strip.placement = "outside",
legend.position="none") +
coord_cartesian(xlim=c(-35, 35))

p2 <- ggplot(full_tumor_dat, aes(x = fraction_diff, fill=cell_state)) +
geom_density() +
scale_fill_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
facet_grid(cell_state ~ idh_status,scales="free_y",switch="y") +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_blank(),
strip.text.y.left = element_text(size=7,angle = 0, hjust = 1),
strip.background = element_blank(),
strip.placement = "outside",
legend.position="none")  +
coord_cartesian(xlim=c(-35, 35))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_diff_density_mut_uncorrected.pdf", width=2.75, height = 3.75) #,width=2.7,height=3)
grid.arrange(p1, p2, heights = c(3.1,1))
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Mutations excluded samples (revision figure)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
p1 <- ggplot(sub_dat, aes(x = fraction_diff, fill=cell_state)) +
geom_density() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligo." = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
facet_grid(cell_state ~ idh_status,scales="free_y",switch="y") +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_blank(),
strip.text.y.left = element_text(size=7,angle = 0,hjust=1),
strip.background = element_blank(),
strip.placement = "outside",
legend.position="none") +
coord_cartesian(xlim=c(-35, 35))

p2 <- ggplot(tumor_dat, aes(x = fraction_diff, fill=cell_state)) +
geom_density() +
scale_fill_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
facet_grid(cell_state ~ idh_status,scales="free_y",switch="y") +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_blank(),
strip.text.y.left = element_text(size=7,angle = 0, hjust = 1),
strip.background = element_blank(),
strip.placement = "outside",
legend.position="none")  +
coord_cartesian(xlim=c(-35, 35))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_diff_density_mut_corrected.pdf", width=2.75, height = 3.75) #,width=2.7,height=3)
grid.arrange(p1, p2, heights = c(3.1,1))
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Variance graph (revision figure)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
var_res <- c(res$sd, full_res$sd, tumor_res$sd, full_tumor_res$sd)
#var_res <- var_res ^2
idh_status <- c(as.character(res$idh_status), as.character(full_res$idh_status), as.character(tumor_res$idh_status), as.character(full_tumor_res$idh_status))
cell_state <- c(as.character(res$cell_state), as.character(full_res$cell_state), as.character(tumor_res$cell_state), as.character(full_tumor_res$cell_state))
var_class <- c(rep("Driver excluded", nrow(res)), rep("All", nrow(res)), rep("Driver excluded", nrow(tumor_res)), rep("All", nrow(full_tumor_res)))
fraction_class <- c(rep("All", nrow(res)*2), rep("Tumor only", nrow(tumor_res)*2))

plot_var <- data.frame(var_res, idh_status, cell_state, var_class, fraction_class)
plot_var <- plot_var %>%
			mutate(cell_state = fct_relevel(cell_state, "Prolif. stem-like", "Stem-like","Diff.-like",
			"Fibroblast", "Pericyte","Endothelial", "Oligo.",
			"Myeloid", "Dendritic cell", "T cell", "Granulocyte","B cell")) %>%
			mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) 


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_diff_var_comp.pdf", width=5, height = 3) #,width=2.7,height=3)
ggplot(plot_var, aes(fill=var_class, y=var_res, x=cell_state)) + 
geom_bar(position="dodge", stat="identity") +
scale_fill_manual(values=c("#27408B", "#CD4F39")) +
labs(y = "Standard deviation (cell state change)") +
facet_grid(idh_status ~ fraction_class, scales = "free_x", space = "free_x") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_diff_var_comp2.pdf", width=5, height = 3) #,width=2.7,height=3)
ggplot(plot_var, aes(fill=var_class, y=var_res^2, x=cell_state)) + 
geom_bar(position="dodge", stat="identity") +
scale_fill_manual(values=c("#27408B", "#CD4F39")) +
labs(y = "Standard deviation (cell state change)") +
facet_grid(idh_status ~ fraction_class, scales = "free_x", space = "free_x") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7))
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Proliferating stem-like density change (revision figure)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

d1 <- tumor_dat %>% filter(cell_state == "Prolif. stem-like")
d2 <- full_tumor_dat %>% filter(cell_state == "Prolif. stem-like")
plot_pro <- rbind(d1,d2)
dat_class <- c(rep("Driver excluded", nrow(d1)), rep("All", nrow(d2)))
plot_pro <- data.frame(plot_pro, dat_class)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/prolif_stem_highlight_dens.pdf", width=2, height = 2) #,width=2.7,height=3)
ggplot(plot_pro, aes(x = fraction_diff, fill=cell_state)) +
geom_density() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligo." = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
facet_grid(idh_status~dat_class,scales="free") +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_text(size=7),
strip.text.y = element_text(size=7,hjust=0.5),
strip.background = element_blank(),
strip.placement = "outside",
legend.position="none") +
coord_cartesian(xlim=c(-35, 35))
dev.off()
