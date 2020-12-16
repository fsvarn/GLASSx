###################################################
# Compare CIBERSORTx profiles in initial vs recurrent tumors across glioma molecular subtypes
# Updated: 2020.06.17
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ag1.signature_name AS subtype_a, 
ag2.signature_name AS subtype_b,
ci1.cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
"

dat <- dbGetQuery(con,q)

cells <- unique(dat[,"cell_state"])
subtypes <- unique(dat[,"idh_status"])

p.val <- matrix(0,nrow = length(cells), ncol = length(subtypes))
colnames(p.val) <- subtypes
rownames(p.val) <- cells
for(i in 1:length(cells))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat %>%
				filter(cell_state == cells[i] & idh_status == subtypes[j])
		p.val[i,j] <- wilcox.test(sub_dat[,"fraction_a"], sub_dat[,"fraction_b"])$p.value
	}
}

plot_fract <- dat %>% 
			  group_by(idh_status, cell_state) %>%
			  summarise(fraction_a = mean(fraction_a), fraction_b = mean(fraction_b)) %>%
			  pivot_longer(-c("idh_status","cell_state"), values_to = "fraction") %>%
			  mutate(name = recode(name, "fraction_a" = "Initial", "fraction_b" = "Recurrent"),
			  cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
			  				"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
			  				"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
			  				"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
			  				"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Proliferating stem-like",
			  				"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
			  mutate(cell_state = as_factor(cell_state)) %>%
			  mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
			  											  "Oligodendrocyte", 
			  											  "Endothelial", "Pericyte",
			  											  "Fibroblast", 
			  											  "Diff.-like", "Stem-like", "Proliferating stem-like"))

ts_fract <- ts_fract %>% mutate(fraction = fraction*100)
			  				
			  							  
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cibersortx_stacked_barplot_subtype_timepoint.pdf",width=1.286,height=1.9127)  #,width=1.287,height=1.9127) 
ggplot(plot_fract, aes(fill=cell_state, y=fraction, x=name)) + 
geom_bar(position="stack", stat="identity") +
facet_grid(.~idh_status) +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "#ndothelial" = "#ffffd4", "Pericyte" = "#fee391",
						 "Fibroblast" = "#feb24c",
						 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Proliferating stem-like" = "#a50f15")) +
labs(y = "Proportion") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "none")
dev.off()

