
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
ss.grade_a,
ss.grade_b,
ag1.signature_name AS subtype_a, 
ag2.signature_name AS subtype_b,
ci1.cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
ss.idh_codel_subtype,
CASE WHEN ss.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.tumor_rna_clinical_comparison ss
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
WHERE (ss.idh_codel_subtype = 'IDHwt' AND grade_a = 4) OR (grade_a < 4 AND grade_b = 4)
"

dat <- dbGetQuery(con,q)

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
			  											  "Diff.-like", "Stem-like", "Proliferating stem-like")) %>%
			  mutate(fraction = fraction*100)

tumor_fract <- dat %>% 
			  filter(base::grepl("_tumor",cell_state)) %>%
			  group_by(case_barcode) %>%
			  summarise(idh_status, cell_state, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b)) %>%
			  ungroup() %>%
			  group_by(idh_status, cell_state) %>%
			  summarise(fraction_a = mean(fraction_a), fraction_b = mean(fraction_b)) %>%
			  pivot_longer(-c("idh_status","cell_state"), values_to = "fraction") %>%
			  mutate(name = recode(name, "fraction_a" = "Initial", "fraction_b" = "Recurrent"),
			  cell_state = recode(cell_state, "differentiated_tumor" = "Diff.-like", "prolif_stemcell_tumor" = "Proliferating stem-like",
			  				"stemcell_tumor" = "Stem-like")) %>%
			  mutate(cell_state = as_factor(cell_state)) %>%
			  mutate(cell_state = fct_relevel(cell_state, "Diff.-like", "Stem-like", "Proliferating stem-like")) %>%
			  mutate(fraction = fraction * 100)
			  

g1 <- plot_fract %>% filter(idh_status == "IDHwt", name == "Initial") %>% .$fraction
g2 <- plot_fract %>% filter(idh_status == "IDHwt", name == "Recurrent")	%>% .$fraction	
g3 <- plot_fract %>% filter(idh_status == "IDHmut", name == "Initial") %>% .$fraction
g4 <- plot_fract %>% filter(idh_status == "IDHmut", name == "Recurrent") %>% .$fraction		

avg_matrix <- data.frame(g1,g2,g3,g4)
cor_matrix <- cor(avg_matrix,method="p")
			  				
			
g1 <- tumor_fract %>% filter(idh_status == "IDHwt", name == "Initial") %>% .$fraction
g2 <- tumor_fract %>% filter(idh_status == "IDHwt", name == "Recurrent")	%>% .$fraction	
g3 <- tumor_fract %>% filter(idh_status == "IDHmut", name == "Initial") %>% .$fraction
g4 <- tumor_fract %>% filter(idh_status == "IDHmut", name == "Recurrent") %>% .$fraction		

avg_matrix <- data.frame(g1,g2,g3,g4)
cor_matrix <- cor(avg_matrix,method="p")
			  				
plot_cor <- cor_matrix %>%
			rownames_to_column(sample_type) %>%
			pivot_longer()
			  							  