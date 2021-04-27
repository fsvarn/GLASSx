#######################################################
# Step 3: Plot the results: Box plots and then average barplots
#######################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

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

g1 <- plot_dat %>% filter(cell_state == "Fibroblast" & timepoint == "Init.") %>% .$fraction
g2 <- plot_dat %>% filter(cell_state == "Fibroblast" & timepoint == "Rec.") %>% .$fraction
t.test(g1,g2,paired=TRUE)

g1 <- plot_dat %>% filter(cell_state == "Granulocyte" & timepoint == "Init.") %>% .$fraction
g2 <- plot_dat %>% filter(cell_state == "Granulocyte" & timepoint == "Rec.") %>% .$fraction
t.test(g1,g2,paired=TRUE)
wilcox.test(g1,g2,paired=TRUE)

# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Granulocyte" | cell_state == "Fibroblast"),
	   aes(x = timepoint, y = fraction, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
facet_grid(.~cell_state) +
scale_colour_manual(values=c("gray","#CD4F39")) +
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

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/NF1_cell_state_changes.pdf",width=2.5,height=1.5)
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

pdf("/projects/verhaak-lab/GLASS-III/egfr_cell_state_changes.pdf",width=1.9,height=1.5)
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

pdf("/projects/verhaak-lab/GLASS-III/pdgfra_cell_state_changes.pdf",width=1.9,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()


# CDKN2A and CCND2 plot- IDHwt
#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND (
	cnv_driver_change_a LIKE '%CDKN2A%' AND (cnv_driver_change_b NOT LIKE '%CDKN2A%' OR cnv_driver_change_b IS NULL) OR
	cnv_driver_change_a LIKE '%CCND2%' AND (cnv_driver_change_b NOT LIKE '%CCND2%' OR cnv_driver_change_b IS NULL))
--WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND cnv_driver_change_a LIKE '%CDKN2A%' AND (cnv_driver_change_b NOT LIKE '%CDKN2A%' OR cnv_driver_change_b IS NULL)
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

dat %>% group_by(cell_state) %>% summarise(p.val = t.test(fraction_a, fraction_b, paired=TRUE)$p.value, eff = median(fraction_b - fraction_a))
# Only oligodendrocyte (P = 0.04)

g1 <- plot_dat %>% filter(cell_state == "Oligodendrocyte" & timepoint == "Init.") %>% .$fraction
g2 <- plot_dat %>% filter(cell_state == "Oligodendrocyte" & timepoint == "Rec.") %>% .$fraction
t.test(g1,g2,paired=TRUE)
wilcox.test(g1,g2, paired=TRUE)

# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Proliferating stem cell tumor"),
	   aes(x = timepoint, y = fraction, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
scale_colour_manual(values=c("gray","#CD4F39")) +
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

#pdf("/projects/verhaak-lab/GLASS-III/cdkn2a_idhwt_cell_state_changes.pdf",width=1.9,height=1.5) Old file for CDKN2A
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_cycle_idhwt_cell_state_changes.pdf",width=1.9,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()


# CDKN2A plot- IDHmut
#Read in data
q <- "
SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b
FROM analysis.platinum_set ps
JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp cs1 ON cs1.aliquot_barcode =  ps.rna_barcode_a 
JOIN analysis.cibersortx_scgp cs2 ON cs2.aliquot_barcode =  ps.rna_barcode_b AND cs1.cell_state = cs2.cell_state
WHERE ds.idh_codel_subtype LIKE 'IDHmut%' AND (
	cnv_driver_change_a LIKE '%CDKN2A%' AND (cnv_driver_change_b NOT LIKE '%CDKN2A%' OR cnv_driver_change_b IS NULL) OR
	cnv_driver_change_a LIKE '%CCND2%' AND (cnv_driver_change_b NOT LIKE '%CCND2%' OR cnv_driver_change_b IS NULL))
--WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND cnv_driver_change_a LIKE '%CDKN2A%' AND (cnv_driver_change_b NOT LIKE '%CDKN2A%' OR cnv_driver_change_b IS NULL)
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
		
g1 <- plot_dat %>% filter(cell_state == "Proliferating stem cell tumor" & timepoint == "Init.") %>% .$fraction
g2 <- plot_dat %>% filter(cell_state == "Proliferating stem cell tumor" & timepoint == "Rec.") %>% .$fraction
t.test(g1,g2,paired=TRUE)
wilcox.test(g1,g2, paired=TRUE)

# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "Proliferating stem cell tumor"),
	   aes(x = timepoint, y = fraction, colour= timepoint)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1) +
scale_colour_manual(values=c("gray","#CD4F39")) +
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

#pdf("/projects/verhaak-lab/GLASS-III/cdkn2a_idhmut_cell_state_changes.pdf",width=1.9,height=1.5) Old file for CDKN2A
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_cycle_idhmut_cell_state_changes.pdf",width=1.9,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()
