###################################################
# Correlate CCF of some drivers with different CIBERSORTx components
# Updated: 2020.06.17
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

##################################################
rm(list=ls())
# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
SELECT ps.*, 
dr1.driver AS driver_a, dr1.driver_ccf AS driver_ccf_a, 
dr2.driver AS driver_b, dr2.driver_ccf AS driver_ccf_b, 
dr2.driver_ccf::numeric - dr1.driver_ccf::numeric AS driver_ccf_dif,
ci1.cell_state, ci1.fraction AS fraction_a, ci2.fraction AS fraction_b,
ci2.fraction - ci1.fraction AS fraction_dif
FROM analysis.platinum_set ps 
JOIN analysis.driver_status_snv_long dr1 ON ps.dna_barcode_a = dr1.aliquot_barcode
JOIN analysis.driver_status_snv_long dr2 ON ps.dna_barcode_b = dr2.aliquot_barcode AND dr1.driver = dr2.driver
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ps.rna_barcode_b AND ci2.cell_state = ci1.cell_state
"

dat <- dbGetQuery(con,q)

q <- "
SELECT * FROM ref.driver_genes WHERE has_mut
"
driver_table <- dbGetQuery(con,q)

drivers <- driver_table[-1,1]
cells <- unique(dat[,"cell_state"])

ccf_cor <- matrix(NA, nrow = length(cells), ncol = length(drivers))
rownames(ccf_cor) <- cells
colnames(ccf_cor) <- drivers
for(i in 1:length(cells))
{
	for(j in 1:length(drivers))
	{
		sub_dat <- dat %>%
			filter(cell_state == cells[i], grepl(drivers[j], driver_a))
		
		skip <- sum(!is.na(sub_dat[,"driver_ccf_dif"]))
		if(skip < 5){
			next}
			
		ccf_cor[i,j] <- cor(sub_dat[,"driver_ccf_dif"], sub_dat[,"fraction_dif"], method = "s", use="complete.obs")
	}
}
na_check <- apply(ccf_cor, 2, function(x)sum(is.na(x)))
ccf_cor <- ccf_cor[,-which(na_check > 0)]
ccf_cor <- as.data.frame(ccf_cor)
ccf_cor[,"cell_states"] <- rownames(ccf_cor)
ccf_cor <- ccf_cor %>% as.data.frame() %>% pivot_longer(-cell_states, names_to = "driver", values_to = "cor")

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/snv_ccf_cibersortx_heatmap.pdf",width=4,height=2)
ggplot(ccf_cor, aes(x=cell_states,y=driver,fill=cor)) +
geom_tile() +
scale_fill_gradient2(low="royalblue4", mid="#ffffff", high="tomato3",midpoint=0) +
theme_bw() +
theme(axis.title=element_blank(),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=5),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "right")
dev.off()



