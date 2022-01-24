##################################################
# Compare CIBERSORTx profiles in samples with lost EGFR/PDGFRA amplifications
# Author: Frederick Varn
# Date: 2022.01.05
# Figure S3E
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
											"Differentiated tumor", "Stem cell tumor", "Proliferating stem cell tumor")) %>%
			mutate(cell_class = recode(cell_state, "B cell" = "Normal", "Dendritic cell" = "Normal",
					"Differentiated tumor" = "Malignant", "Endothelial" = "Normal",
					"Fibroblast" = "Normal", "Granulocyte" = "Normal",
					"Myeloid" = "Normal", "Oligodendrocyte" = "Normal",
					"Pericyte" = "Normal", "Proliferating stem cell tumor" = "Malignant",
					"Stem cell tumor" = "Malignant","T cell" = "Normal"))
		
norm_sum <- plot_dat %>%
			group_by(case_barcode, cell_class, timepoint) %>%
			summarise(fraction = sum(fraction/100))
			
g1 <- norm_sum %>% filter(cell_class == "Normal" & timepoint == "Init.") %>% .$fraction
g2 <- norm_sum %>% filter(cell_class == "Normal" & timepoint == "Rec.") %>% .$fraction
t.test(g1,g2,paired=TRUE)
wilcox.test(g1,g2, paired=TRUE)

# Boxplot
boxes <- ggplot(data = norm_sum %>% 
	   		  filter(cell_class == "Normal"),
	   aes(x = timepoint, y = fraction*100, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
scale_colour_manual(values=c("gray","#27408B")) +
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

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/egfr_cell_state_changes.pdf",width=1.9,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()

# PDGFRA plot
#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND cnv_driver_change_b LIKE '%PDGFRA%' AND (cnv_driver_change_a NOT LIKE '%PDGFRA%' OR cnv_driver_change_a IS NULL)
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
											"Differentiated tumor", "Stem cell tumor", "Proliferating stem cell tumor")) %>%
			mutate(cell_class = recode(cell_state, "B cell" = "Normal", "Dendritic cell" = "Normal",
					"Differentiated tumor" = "Malignant", "Endothelial" = "Normal",
					"Fibroblast" = "Normal", "Granulocyte" = "Normal",
					"Myeloid" = "Normal", "Oligodendrocyte" = "Normal",
					"Pericyte" = "Normal", "Proliferating stem cell tumor" = "Malignant",
					"Stem cell tumor" = "Malignant","T cell" = "Normal"))
		
norm_sum <- plot_dat %>%
			group_by(case_barcode, cell_class, timepoint) %>%
			summarise(fraction = sum(fraction/100))
			
g1 <- norm_sum %>% filter(cell_class == "Normal" & timepoint == "Init.") %>% .$fraction
g2 <- norm_sum %>% filter(cell_class == "Normal" & timepoint == "Rec.") %>% .$fraction
t.test(g1,g2,paired=TRUE)
wilcox.test(g1,g2, paired=TRUE)

# Boxplot
boxes <- ggplot(data = norm_sum %>% 
	   		  filter(cell_class == "Normal"),
	   aes(x = timepoint, y = fraction*100, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
scale_colour_manual(values=c("gray","#27408B")) +
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

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/pdgfra_cell_state_changes.pdf",width=1.9,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()
