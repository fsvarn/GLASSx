###################################################
# Compare CIBERSORTx profiles between initial and recurrent tumors with different extents of resection
# Author: Frederick Varn
# Date: 2021.12.29
# Figure S1F
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Figure 1: Rates of sub-total/total resection at each time point
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
su1.surgery_extent_of_resection AS extent_a, 
su2.surgery_extent_of_resection AS extent_b
FROM analysis.rna_silver_set ss
JOIN clinical.surgeries su1 ON su1.sample_barcode = substring(ss.tumor_barcode_a, 1, 15) 
JOIN clinical.surgeries su2 ON su2.sample_barcode = substring(ss.tumor_barcode_b, 1, 15) 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

c1 <- sum(dat$extent_a == "Subtotal", na.rm=TRUE)
c2 <- sum(dat$extent_a == "Total", na.rm=TRUE)
c3 <- sum(dat$extent_b == "Subtotal", na.rm=TRUE)
c4 <- sum(dat$extent_b == "Total", na.rm=TRUE)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))

fraction <- c(c1/sum(c1,c2), c2/sum(c1,c2),c3/sum(c3,c4), c4/sum(c3,c4))
status <- c("Initial", "Initial", "Recurrent","Recurrent")
resection <- c("Subtotal", "Total", "Subtotal", "Total")
plot_dat <- data.frame(status, resection, fraction)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/resection_over_time.pdf",width=1,height=1.7)
ggplot(plot_dat, aes(fill=resection, y=fraction*100, x=status)) + 
geom_bar(position="stack", stat="identity") +
scale_fill_manual(values=c("Total" = "#8B8B8B", "Subtotal" = "#ECECEC")) +
labs(y = "Proportion total resection (%)") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "none")
dev.off()



#Read in data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ci1.cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
su1.surgery_extent_of_resection AS extent_a, 
su2.surgery_extent_of_resection AS extent_b
FROM analysis.rna_silver_set ss
JOIN clinical.surgeries su1 ON su1.sample_barcode = substring(ss.tumor_barcode_a, 1, 15) 
JOIN clinical.surgeries su2 ON su2.sample_barcode = substring(ss.tumor_barcode_b, 1, 15) 
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt' --AND su1.surgery_extent_of_resection = 'Total' AND su2.surgery_extent_of_resection = 'Total'
"

dat <- dbGetQuery(con,q)

status <- rep("", nrow(dat))
status[which(dat$extent_a == "Total" & dat$extent_b == "Total")] <- "TT"
status[which(dat$extent_a == "Subtotal" & dat$extent_b == "Total")] <- "ST"
status[which(dat$extent_a == "Total" & dat$extent_b == "Subtotal")] <- "TS"
status[which(dat$extent_a == "Subtotal" & dat$extent_b == "Subtotal")] <- "SS"
dat$status <- status

dat %>%
filter(status != "") %>%
group_by(cell_state, status) %>%
summarise(pval = wilcox.test(fraction_a, fraction_b, paired=TRUE)$p.value,
		  eff = median(fraction_b - fraction_a)) %>%
data.frame()

plot_fract <- dat %>% 
			  filter(status != "") %>%
			  group_by(cell_state, status) %>%
			  summarise(fraction_a = mean(fraction_a), fraction_b = mean(fraction_b)) %>%
			  pivot_longer(-c("status", "cell_state"), values_to = "fraction") %>%
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
			  data.frame()


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/total_resection_change.pdf",width=1.58,height=1.913)
ggplot(plot_fract %>% filter(status == "TT" | status == "ST"), aes(fill=cell_state, y=fraction, x=name)) + 
geom_bar(position="stack", stat="identity") +
facet_grid(.~status) +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
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

