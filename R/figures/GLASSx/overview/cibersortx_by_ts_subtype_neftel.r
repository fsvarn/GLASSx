###################################################
# Compare CIBERSORTx Neftel et al profiles across glioma transcriptional subtypes
# Updated: 2020.04.13
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
JOIN analysis.cibersortx_neftel ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_neftel ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
WHERE idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

grp_b <- dat %>% 
		 select("tumor_barcode_b","subtype_b", "cell_state", "fraction_b") %>%
		 pivot_longer(-c("tumor_barcode_b","subtype_b", "cell_state"), values_to = "fraction") %>%
		 rename(aliquot_barcode = tumor_barcode_b, subtype = subtype_b)
ts_fract <- dat %>% 
		 select("tumor_barcode_a","subtype_a", "cell_state", "fraction_a") %>%
		 pivot_longer(-c("tumor_barcode_a","subtype_a", "cell_state"), values_to = "fraction") %>%
		 rename(aliquot_barcode = tumor_barcode_a, subtype = subtype_a) %>%
		 bind_rows(grp_b) %>%
		 filter(!grepl(",", subtype)) %>%	 
		 group_by(subtype, cell_state) %>%
		 summarise(fraction = mean(fraction)) %>%
		 mutate(cell_state = recode(cell_state, "AClike" = "AC-like", "Macrophage" = "Macrophage",
			"Mesenchymal" = "MES-like", "NPClike" = "NPC-like",
			"Oligodendrocyte" = "Oligodendrocyte", "OPClike" = "OPC-like",
			"T" = "T cell")) %>%
		 mutate(cell_state = as_factor(cell_state)) %>%
		 mutate(cell_state = fct_relevel(cell_state, "T cell", "Macrophage", 
		    "Oligodendrocyte", 
		    "MES-like", "AC-like", "OPC-like", "NPC-like"))

ts_fract <- ts_fract %>% mutate(fraction = fraction*100)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cibersortx_stacked_barplot_ts_subtype.pdf",width=1,height=1.82)
ggplot(ts_fract, aes(fill=cell_state, y=fraction, x=subtype)) + 
geom_bar(position="stack", stat="identity") +
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