##################################################
# Compare subtype switching in samples with lost EGFR/PDGFRA amplifications
# Author: Frederick Varn
# Date: 2021.10.28
# Figure S3F
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(ggalluvial)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
WITH full_data AS
(
	SELECT ps.*, ds.idh_codel_subtype, cnv_driver_change_a, cnv_driver_change_b, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
	FROM analysis.platinum_set ps
	JOIN analysis.driver_status_cnv ds ON ds.tumor_barcode_a = ps.dna_barcode_a AND ds.tumor_barcode_b = ps.dna_barcode_b
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode =  ps.rna_barcode_a 
	JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode =  ps.rna_barcode_b 
	WHERE ds.idh_codel_subtype LIKE 'IDHwt' AND ((cnv_driver_change_b LIKE '%EGFR%' AND (cnv_driver_change_a NOT LIKE '%EGFR%' OR cnv_driver_change_a IS NULL)) OR
	(cnv_driver_change_b LIKE '%PDGFRA%' AND (cnv_driver_change_a NOT LIKE '%PDGFRA%' OR cnv_driver_change_a IS NULL)))
)
SELECT subtype_a, subtype_b, COUNT(*)
FROM full_data
GROUP BY subtype_a, subtype_b
"

dat <- dbGetQuery(con,q)
dat[,"count"] <- as.numeric(dat[,"count"])

# Fisher's exact tests
c1 <- sum(dat %>% filter(subtype_b == "Mesenchymal") %>% .$count)
c2 <- sum(dat %>% filter(subtype_b != "Mesenchymal") %>% .$count)
c3 <- sum(dat %>% filter(subtype_a == "Mesenchymal") %>% .$count)
c4 <- sum(dat %>% filter(subtype_a != "Mesenchymal") %>% .$count)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))

dat <- dat %>% 
	   mutate(subtype_a = recode(subtype_a, "Mesenchymal" = "Mes.", "Proneural" = "Pro.", "Classical" = "Class.")) %>%
	   mutate(subtype_b = recode(subtype_b, "Mesenchymal" = "Mes.", "Proneural" = "Pro.", "Classical" = "Class."))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/amp_loss_subtype_sankey.pdf",width=2,height=1.6)
ggplot(dat, aes(axis1 = subtype_a, axis2 = subtype_b, y=count)) +
geom_alluvium(aes(fill = subtype_a), width = 1/12, knot.pos = 0.4) +
geom_stratum(width = 1/4, fill = "white", color = "grey") +
scale_x_discrete(limits = c("Initial", "Recurrent"),expand = c(0,0)) +
scale_y_continuous(expand = c(0,0)) +
scale_fill_manual(values=c("#008A22","#8A0000","#00458A")) +
geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
theme_bw() +
theme(axis.text.x= element_text(size=7), axis.text.y=element_blank(),
axis.ticks = element_blank(),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust = 0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 
dev.off()
