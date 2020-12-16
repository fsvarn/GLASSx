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

#Read in IDHwt data for analysis
q <- "
SELECT sc.aliquot_barcode, sc.cell_state AS scgp_cell_state, sc.fraction AS scgp_fraction,
nf.cell_state AS neftel_cell_state, nf.fraction AS neftel_fraction,
idh_codel_subtype AS idh_status
--CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.cibersortx_scgp sc
JOIN analysis.cibersortx_neftel nf ON nf.aliquot_barcode = sc.aliquot_barcode
JOIN analysis.rna_silver_set ps ON sc.aliquot_barcode = ps.tumor_barcode_a OR sc.aliquot_barcode = ps.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con, q)

dat <- dat %>% 
		mutate(scgp_cell_state = recode(scgp_cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		filter(scgp_cell_state %in% c("Diff.-like", "Stem-like","Prolif. stem-like")) %>%
		filter(neftel_cell_state %in% c("AClike", "Mesenchymal","NPClike","OPClike")) %>%
		group_by(aliquot_barcode) %>%
		summarise(scgp_cell_state, scgp_fraction = scgp_fraction/sum(scgp_fraction), neftel_cell_state, neftel_fraction = neftel_fraction/sum(neftel_fraction), idh_status) %>%
		as.data.frame()

idhwt_res <- dat %>% 
		mutate(scgp_cell_state = as_factor(scgp_cell_state)) %>%
		mutate(scgp_cell_state = fct_relevel(scgp_cell_state, "Prolif. stem-like","Stem-like","Diff.-like")) %>%
		filter(idh_status == "IDHwt") %>%
		group_by(scgp_cell_state, neftel_cell_state, idh_status) %>% 
		summarise(cor = cor(scgp_fraction, neftel_fraction, method="p")) %>%
		as.data.frame()
		
#Read in IDHmut data for analysis
q <- "
SELECT sc.aliquot_barcode, sc.cell_state AS scgp_cell_state, sc.fraction AS scgp_fraction,
nf.cell_state AS venteicher_cell_state, nf.fraction AS venteicher_fraction,
idh_codel_subtype AS idh_status
--CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.cibersortx_scgp sc
JOIN analysis.cibersortx_venteicher nf ON nf.aliquot_barcode = sc.aliquot_barcode
JOIN analysis.rna_silver_set ps ON sc.aliquot_barcode = ps.tumor_barcode_a OR sc.aliquot_barcode = ps.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con, q)

dat <- dat %>% 
		mutate(scgp_cell_state = recode(scgp_cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		filter(scgp_cell_state %in% c("Diff.-like", "Stem-like","Prolif. stem-like")) %>%
		filter(venteicher_cell_state %in% c("Astro", "Oligo","Undifferentiated")) %>%
		group_by(aliquot_barcode) %>%
		summarise(scgp_cell_state, scgp_fraction = scgp_fraction/sum(scgp_fraction), venteicher_cell_state, venteicher_fraction = venteicher_fraction/sum(venteicher_fraction), idh_status) %>%
		as.data.frame()

idhmut_res <- dat %>% 
		mutate(scgp_cell_state = as_factor(scgp_cell_state)) %>%
		mutate(scgp_cell_state = fct_relevel(scgp_cell_state, "Prolif. stem-like","Stem-like","Diff.-like")) %>%
		filter(idh_status != "IDHwt") %>%
		group_by(scgp_cell_state, venteicher_cell_state, idh_status) %>% 
		summarise(cor = cor(scgp_fraction, venteicher_fraction, method="p")) %>%
		as.data.frame()

colnames(idhwt_res) <- c("scgp_cell_state", "suva_state", "idh_status","cor")
colnames(idhmut_res) <- c("scgp_cell_state", "suva_state", "idh_status","cor")
full_res <- rbind(idhwt_res, idhmut_res) %>%
		mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut-noncodel","IDHmut-codel"))

				
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_neftel_venteicher_cor.pdf", width=2.75, height = 1.7) #,width=2.7,height=3)
ggplot(full_res, aes(x = suva_state, y = scgp_cell_state, fill=cor)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
facet_grid(. ~ idh_status,scales="free_x",space = "free_x") +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") 
dev.off()
