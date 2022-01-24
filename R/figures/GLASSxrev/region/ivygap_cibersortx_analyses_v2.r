###################################################
# Compare CIBERSORTx Ivy GAP profiles in initial vs recurrent tumors, transcriptional subtypes, and subtype switches
# Author: Frederick Varn
# Date: 2021.12.01
# Figures 2B, 2C, 2D
##################################################

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

# Compare initial/recurrent and transcriptional subtype

#Read in data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ag1.signature_name AS subtype_a, 
ag2.signature_name AS subtype_b,
ci1.cell_state AS region,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_ivygap ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
"

dat <- dbGetQuery(con,q)

grp_b <- dat %>% 
		 select("tumor_barcode_b","subtype_b", "region", "fraction_b") %>%
		 pivot_longer(-c("tumor_barcode_b","subtype_b", "region"), values_to = "fraction") %>%
		 rename(aliquot_barcode = tumor_barcode_b, subtype = subtype_b)
		 
ts_fract <- dat %>% 
		 select("tumor_barcode_a","subtype_a", "region", "fraction_a") %>%
		 pivot_longer(-c("tumor_barcode_a","subtype_a", "region"), values_to = "fraction") %>%
		 rename(aliquot_barcode = tumor_barcode_a, subtype = subtype_a) %>%
		 bind_rows(grp_b) %>%
		 filter(!grepl(",", subtype)) %>%	 
		 group_by(subtype, region) %>%
		 summarise(fraction = mean(fraction)) %>%
		 mutate(region = recode(region, "CT" = "Cellular tumor", "CTmvp" = "Microvascular proliferation", "CTpan" = "Pseudopalisading cells around necrosis", "LE" = "Leading edge")) %>%
		 mutate(region = as_factor(region)) %>%
		 mutate(region = fct_relevel(region, 
		 "Leading edge", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation"))

ts_fract <- ts_fract %>% mutate(fraction = fraction*100)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_stacked_barplot_ts_subtype.pdf",width=1,height=1.82)
ggplot(ts_fract, aes(fill=region, y=fraction, x=subtype)) + 
geom_bar(position="stack", stat="identity") +
scale_fill_manual(values=c("Leading edge" = "#009999", "Cellular tumor" = "#01b050", "Pseudopalisading cells around necrosis"="#02ffcc", "Microvascular proliferation"="#c00000")) +
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

idh_fract <- dat %>% 
			  group_by(idh_status, region) %>%
			  summarise(fraction_a = mean(fraction_a), fraction_b = mean(fraction_b)) %>%
			  pivot_longer(-c("idh_status","region"), values_to = "fraction") %>%
			  mutate(name = recode(name, "fraction_a" = "Initial", "fraction_b" = "Recurrent"),
			  region = recode(region, "CT" = "Cellular tumor", "CTmvp" = "Microvascular proliferation", "CTpan" = "Pseudopalisading cells around necrosis", "LE" = "Leading edge")) %>%
			  mutate(region = as_factor(region)) %>%
			  mutate(region = fct_relevel(region, 
			  "Leading edge", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation")) %>%
			  mutate(idh_status = as_factor(idh_status)) %>%
			  mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>%
			  mutate(fraction = fraction*100)

# Test changes
test <- dat %>%
group_by(region, idh_status) %>%
summarise(pval = wilcox.test(fraction_a, fraction_b, paired=TRUE)$p.value, eff = median(fraction_b - fraction_a))

dat %>%
group_by(region, idh_status) %>%
summarise(pval = t.test(fraction_a, fraction_b, paired=TRUE)$p.value, eff = mean(fraction_b - fraction_a))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_stacked_barplot_idh_subtype.pdf",width=1.286,height=1.9127)
ggplot(idh_fract, aes(fill=region, y=fraction, x=name)) + 
geom_bar(position="stack", stat="identity") +
scale_fill_manual(values=c("Leading edge" = "#009999", "Cellular tumor" = "#01b050", "Pseudopalisading cells around necrosis"="#02ffcc", "Microvascular proliferation"="#c00000")) +
labs(y = "Proportion") +
facet_grid(.~idh_status) +
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


ss_res <- dat %>%
		  filter(idh_status == "IDHwt") %>%
		  group_by(subtype_a, subtype_b, region) %>%
		  summarise(pval = t.test(fraction_a, fraction_b, paired=TRUE)$p.value, eff = mean(fraction_b - fraction_a)) %>%
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
		  mutate(region = fct_relevel(region, "Leading edge", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation")) %>%
		  as.data.frame()


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_subtype_switch_v2.pdf",width=3.4,height=1.5)
ggplot(ss_res, aes(x = subtype_a, y = subtype_b, fill=eff)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0,  limits = c(max(ss_res$eff)*-1, max(ss_res$eff))) +
facet_grid(. ~ region) +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") 
dev.off()		  

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_subtype_switch_legend.pdf",width=3,height=1.45)
ggplot(ss_res, aes(x = subtype_a, y = subtype_b, fill=eff)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0,  limits = c(max(ss_res$eff)*-1, max(ss_res$eff))) +
facet_grid(. ~ region) +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white")) 
dev.off()	


# Look at changes in cell state composition and region composition

#Read in IDHwt data for analysis
q <- "
SELECT ps.*, sc1.cell_state, sc2.fraction - sc1.fraction AS fraction_change, iv1.cell_state AS region, iv2.fraction - iv1.fraction AS ivy_change,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ps
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ps.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN analysis.cibersortx_ivygap iv1 ON iv1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.cibersortx_ivygap iv2 ON iv2.aliquot_barcode = ps.tumor_barcode_b AND iv2.cell_state = iv1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con, q)
	 
	 
res <- dat %>% 
	group_by(cell_state, region, idh_status) %>% 
	summarise(cor = cor(fraction_change, ivy_change, method="p"), pval = cor.test(fraction_change, ivy_change, method="p")$p.value) %>%
	ungroup() %>%
	mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>%
	mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligo.",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
	mutate(cell_state = fct_relevel(cell_state, rev(c("B cell", "Granulocyte","T cell", "Dendritic cell","Myeloid",
												"Oligo.", "Endothelial", "Pericyte", "Fibroblast",
												"Diff.-like","Stem-like","Prolif. stem-like")))) %>%
	mutate(idh_status = fct_relevel(idh_status, "IDHwt","IDHmut")) %>%
	mutate(region = recode(region, 
	"CT" = "Cellular tumor", "CTmvp" = "Microvascular proliferation", "CTpan" = "Pseudopalisading cells around necrosis", "LE" = "Leading edge")) %>%
	mutate(region = fct_relevel(region, "Microvascular proliferation", "Pseudopalisading cells around necrosis", "Cellular tumor","Leading edge")) %>%
	as.data.frame() 
	
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_ivygap_change_cor.pdf", width=6.49, height = 1.743) #,width=2.7,height=3)
ggplot(res, aes(x = cell_state, y = region, fill=cor)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
facet_grid(. ~ idh_status) +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") 
dev.off()

	
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_ivygap_change_cor_legend.pdf", width=6, height = 1.7) #,width=2.7,height=3)
ggplot(res, aes(x = cell_state, y = region, fill=cor)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
facet_grid(. ~ idh_status) +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"))
dev.off()