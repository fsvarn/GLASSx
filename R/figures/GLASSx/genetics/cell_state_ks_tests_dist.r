###################################################
# Compare cell state fractions between initial and recurrent tumors (used in Figure 2)
# Updated: 2020.11.19 
# Author: Frederick Varn
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
SELECT ss.*, ci1.cell_state,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
ci1.fraction * 100 AS fraction_a,
ci2.fraction * 100 AS fraction_b,
(ci2.fraction - ci1.fraction) * 100 AS fraction_diff
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
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

		
tumor_dat <- dat %>% 
		filter(grepl("-like", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, fraction_diff = (fraction_b - fraction_a), idh_status) %>%
	   	mutate(cell_state = fct_relevel(cell_state, c("Prolif. stem-like","Stem-like","Diff.-like"))) %>%
		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>% 
		as.data.frame()
	
#pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tumor_cell_state_diff_density.pdf", width=2.75, height = .75) #,width=2.7,height=3)
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
#dev.off()

		
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_diff_density.pdf", width=2.75, height = 3.75) #,width=2.7,height=3)
grid.arrange(p1, p2, heights = c(3.1,1))
dev.off()

res <- dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=residuals, y="pnorm", mean=0, sd = sd(residuals))$p.value, n = n()) %>%
		as.data.frame()

tumor_res <- tumor_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=residuals, y="pnorm", mean=0, sd = sd(residuals))$p.value, n = n()) %>%
		as.data.frame()		
		
		
# Get correlation coefficients

# res <- dat %>% group_by(cell_state, idh_status) %>% 
# 		summarise(cor = cor(fraction_a,fraction_b,method="p"), p.value=cor.test(fraction_a, fraction_b)$p.value, n = n()) %>%
# 		as.data.frame()

res <- dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n()) %>%
		as.data.frame()

tumor_res <- tumor_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(norm_pval = ks.test(x=fraction_diff, y="pnorm", mean=0, sd = sd(fraction_diff))$p.value, n = n()) %>%
		as.data.frame()