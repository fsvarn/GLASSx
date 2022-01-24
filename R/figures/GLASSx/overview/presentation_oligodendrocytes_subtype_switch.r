
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
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
"

dat <- dbGetQuery(con,q)

		  
# Remove the proneural to classical transition (n of 1)
ss_res <- dat %>%
		  filter(idh_status == "IDHwt" & !(subtype_a == "Proneural" & subtype_b =="Classical")) %>%
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
		  filter(region %in% c("oligodendrocyte","differentiated_tumor","stemcell_tumor","prolif_stemcell_tumor"))
		  
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_subtype_switch_v2.pdf",width=3.4,height=1.5)
ggplot(ss_res, aes(x = subtype_a, y = subtype_b, fill=logp)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0) +
facet_grid(. ~ region) +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") 
dev.off()		  
