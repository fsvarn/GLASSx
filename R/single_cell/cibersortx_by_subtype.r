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
WITH subtype_rank AS
(
	SELECT *,
	RANK() OVER (PARTITION BY aliquot_barcode ORDER BY p_value ASC) AS p_rank
	FROM analysis.transcriptional_subtype
),
top_rank AS
(
	SELECT *
	FROM subtype_rank
	WHERE p_rank = 1
),
agg AS
(
	SELECT aliquot_barcode, 
	string_agg(signature_name,',') AS subtype
	FROM top_rank
	GROUP BY aliquot_barcode
)
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ag1.subtype AS subtype_a, 
ag2.subtype AS subtype_b,
ci1.cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN agg ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN agg ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
"

dat <- dbGetQuery(con,q)

dat[which(dat[,"subtype_a"]=="Mesenchymal,Classical"),"subtype_a"]  <- "Classical,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Classical"),"subtype_b"]  <- "Classical,Mesenchymal"

dat[which(dat[,"subtype_a"]=="Mesenchymal,Proneural"),"subtype_a"]  <- "Proneural,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Proneural"),"subtype_b"]  <- "Proneural,Mesenchymal"


cells <- unique(dat[,"cell_state"])
subtypes <- unique(dat[,"idh_codel_subtype"])

p.val <- matrix(0,nrow = length(cells), ncol = length(subtypes))
colnames(p.val) <- subtypes
rownames(p.val) <- cells
for(i in 1:length(cells))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat %>%
				filter(cell_state == cells[i] & idh_codel_subtype == subtypes[j])
		p.val[i,j] <- wilcox.test(sub_dat[,"fraction_a"], sub_dat[,"fraction_b"])$p.value
	}
}

plot_fract <- dat %>% 
			  group_by(idh_codel_subtype, cell_state) %>%
			  summarise(fraction_a = mean(fraction_a), fraction_b = mean(fraction_b)) %>%
			  pivot_longer(-c("idh_codel_subtype","cell_state"), values_to = "fraction") %>%
			  mutate(name = recode(name, "fraction_a" = "Initial", "fraction_b" = "Recurrent"),
			  cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
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
			  							  
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_stacked_barplot_subtype_timepoint.pdf",width=6,height=5)  
ggplot(plot_fract, aes(fill=cell_state, y=fraction, x=name)) + 
geom_bar(position="stack", stat="identity") +
facet_grid(.~idh_codel_subtype) +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "#ndothelial" = "#ffffd4", "Pericyte" = "#fee391",
						 "Fibroblast" = "#feb24c",
						 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
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
legend.position = "right")
dev.off()

