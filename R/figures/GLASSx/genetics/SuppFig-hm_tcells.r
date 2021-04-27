###################################################
# Compare cell state fractions in IDHmut receiving therapy (used in Figure 2)
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

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
SELECT cc.*, ci1.cell_state,
CASE WHEN cc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
mf1.coverage_adj_mut_freq AS mut_freq_a,
mf2.coverage_adj_mut_freq AS mut_freq_b,
mf2.coverage_adj_mut_freq > 10 AS hm
FROM analysis.platinum_set ps
JOIN analysis.tumor_rna_clinical_comparison cc ON cc.tumor_barcode_a = ps.rna_barcode_a AND cc.tumor_barcode_b = ps.rna_barcode_b
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ps.rna_barcode_b AND ci2.cell_state = ci1.cell_state
LEFT JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = ps.dna_barcode_a
LEFT JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = ps.dna_barcode_b
LEFT JOIN analysis.driver_status_cnv cn ON cn.tumor_barcode_a = ps.dna_barcode_a AND cn.tumor_barcode_b = ps.dna_barcode_b
"
dat <- dbGetQuery(con, q)

cells <- unique(dat$cell_state)

hyp <- dat %>% 
		filter(!(is.na(hm))) %>%
		group_by(hm, cell_state, idh_status) %>% 
		summarise(pval = t.test(fraction_b,fraction_a, paired=TRUE)$p.value, eff=mean(fraction_b) - mean(fraction_a), n = n()) %>%
		as.data.frame()

hm_only <- dat %>% filter(hm == 1, cell_state == "t_cell")
case_barcode <- rep(hm_only$case_barcode, 2)
fraction <- c(hm_only$fraction_a, hm_only$fraction_b)
idh_codel_subtype <- rep(hm_only$idh_codel_subtype, 2)
idh_status <- rep(hm_only$idh_status, 2)
timepoint <- rep(c("Initial", "Recurrent"),each = nrow(hm_only))

plot_dat <- data.frame(case_barcode, fraction, idh_codel_subtype, idh_status, timepoint)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/hm_t_cells.pdf",width=2,height=1.8)
ggplot(data = plot_dat, aes(x = timepoint, y = fraction * 100)) +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size = 1, aes(colour = idh_codel_subtype)) +
scale_colour_manual(values = c("IDHmut-codel" = "#F8766D", "IDHmut-noncodel" = "#00BA38", "IDHwt" = "#619CFF")) +
labs(y = "T cell fraction (%)") +
facet_grid(.~idh_status, scale = "free_y") +
theme_classic() +
theme(
	axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim = c(0,7))
dev.off()

		