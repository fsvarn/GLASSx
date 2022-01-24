###################################################
# Compare CIBERSORTx Neftel et al profiles in initial vs recurrent tumors across glioma molecular subtypes
# Author: Frederick Varn
# Date: 2021.12.17
# Figure S1G- initial/recurrent plots
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
JOIN analysis.cibersortx_neftel ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_neftel ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
--WHERE idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

cells <- unique(dat[,"cell_state"])
subtypes <- unique(dat[,"idh_status"])

dat %>%
group_by(cell_state, idh_status) %>%
summarise(pval = wilcox.test(fraction_a, fraction_b, paired=TRUE)$p.value,
		  eff = median(fraction_b - fraction_a)) %>%
data.frame()

dat %>%
group_by(cell_state, idh_status) %>%
summarise(pval = t.test(fraction_a, fraction_b, paired=TRUE)$p.value,
		  eff = mean(fraction_b - fraction_a)) %>%
data.frame()

plot_fract <- dat %>% 
			  group_by(cell_state, idh_status) %>%
			  summarise(fraction_a = mean(fraction_a), fraction_b = mean(fraction_b)) %>%
			  pivot_longer(-c("idh_status","cell_state"), values_to = "fraction") %>%
			  mutate(name = recode(name, "fraction_a" = "Initial", "fraction_b" = "Recurrent"),
			  cell_state = recode(cell_state, "AClike" = "AC-like", "Macrophage" = "Macrophage",
			  				"Mesenchymal" = "MES-like", "NPClike" = "NPC-like",
			  				"Oligodendrocyte" = "Oligodendrocyte", "OPClike" = "OPC-like",
			  				"T" = "T cell")) %>%
			  mutate(cell_state = as_factor(cell_state)) %>%
			  mutate(cell_state = fct_relevel(cell_state, "T cell", "Macrophage", 
			  											  "Oligodendrocyte", 
			  											  "MES-like", "AC-like", "OPC-like", "NPC-like"))
			  							  							  
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cibersortx_stacked_barplot_subtype_timepoint_neftel.pdf",width=1.286,height=1.9127)  #,width=1.287,height=1.9127) 
ggplot(plot_fract, aes(fill=cell_state, y=fraction, x=name)) + 
geom_bar(position="stack", stat="identity") +
facet_grid(.~idh_status) +
scale_fill_manual(values=c("T cell" = "#6baed6", "Macrophage" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "AC-like" = "#fcbba1", "OPC-like" = "#fb6a4a", "MES-like" = "#fee391", "NPC-like" = "#a50f15")) +
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

