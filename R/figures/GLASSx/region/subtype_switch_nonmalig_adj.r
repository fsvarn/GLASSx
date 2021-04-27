# Malignant cell adjustment only


library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

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
ci1.cell_state AS cell_state,
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

dat <- dat %>%
mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell"))


ss_res <- dat %>%
		  filter(cell_state == "Stem-like" | cell_state =="Prolif. stem-like" | cell_state == "Diff.-like") %>%
		  group_by(tumor_pair_barcode, idh_status) %>%
		  summarise(cell_state, subtype_a, subtype_b, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b), n = n()) %>%
		  ungroup() %>%
		  group_by(cell_state, idh_status, subtype_a, subtype_b) %>%
		  summarise(pval = wilcox.test(fraction_a, fraction_b, paired=TRUE)$p.value, eff = median(fraction_b - fraction_a)) %>% data.frame()
ss_res %>% filter(pval < 0.05, subtype_a != subtype_b)		  
		  
		  mutate(logp = log10(pval), status = case_when(
			  pval < 0.05 & eff > 0 ~ "up",
			  pval < 0.05 & eff < 0 ~ "dn",
			  pval > 0.05 ~ "none")) %>%
		  mutate(logp = case_when(
		  	  pval < 0.05 & eff < 0 ~ logp * 1,
		  	  pval < 0.05 & eff >= 0 ~ logp * -1,
		  	  pval > 0.05 ~ logp * 0)) %>%
		  mutate(region = recode(region, 
		  "CT" = "Cellular tumor", "CTmvp" = "Microvascular proliferation", "CTpan" = "Pseudopalisading cells around necrosis", "LE" = "Leading edge")) %>%
		  mutate(region = as_factor(region)) %>%
		  mutate(region = fct_relevel(region, "Leading edge", "Infiltrating tumor", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation")) %>%
		  as.data.frame()
