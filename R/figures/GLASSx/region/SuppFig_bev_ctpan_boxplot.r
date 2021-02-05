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

#Read in IDHwt data for analysis
q <- "
SELECT cc.*, iv1.cell_state, iv1.fraction AS fraction_a, iv2.fraction AS fraction_b
FROM analysis.tumor_rna_clinical_comparison cc
JOIN analysis.cibersortx_ivygap iv1 ON iv1.aliquot_barcode = cc.tumor_barcode_a
JOIN analysis.cibersortx_ivygap iv2 ON iv2.aliquot_barcode = cc.tumor_barcode_b AND iv1.cell_state = iv2.cell_state
WHERE received_bev --AND idh_codel_subtype = 'IDHwt'

"

dat <- dbGetQuery(con, q)

res <- dat %>% 
		mutate(cell_state = recode(cell_state, 
		"CT" = "Cellular tumor", "CTmvp" = "Microvascular proliferation", "CTpan" = "Pseudopalisading cells around necrosis", "LE" = "Leading edge")) %>%
		group_by(cell_state) %>%
		summarise(p.val = wilcox.test(fraction_a, fraction_b, paired=TRUE)$p.value, 
				  eff = median(fraction_b - fraction_a))

plot_dat <- dat %>% 
			dplyr::select(case_barcode, cell_state, fraction_a, fraction_b, idh_codel_subtype) %>%
			pivot_longer(c(fraction_a, fraction_b), names_to = "timepoint", values_to = "fraction") %>%
			mutate(timepoint=recode(timepoint, "fraction_a" = "Initial", "fraction_b" = "Recurrent"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/bev_prepost_mvp.pdf", width=1.25,height=1.75)
ggplot(data = plot_dat %>% filter(cell_state == "CTpan"), aes(x = timepoint, y = fraction*100)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size = 1, aes(colour = idh_codel_subtype)) +
scale_colour_manual(values=c("#00BA38", "#619CFF")) +
labs(y = "Proportion (%)", title = "PAN") +
theme_classic() +
theme(plot.title = element_text(size= 7, hjust = 0.5),
	axis.text.x = element_text(size=7, angle=45, hjust = 1),
	axis.text.y = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim=c(0,100))
dev.off()
