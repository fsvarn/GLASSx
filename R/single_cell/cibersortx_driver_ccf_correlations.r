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
dr1.snv_driver AS snv_driver_a, dr1.snv_driver_ccf AS snv_driver_ccf_a, 
dr2.snv_driver AS snv_driver_b, dr2.snv_driver_ccf AS snv_driver_ccf_b, 
ci1.cell_state, ci1.fraction AS fraction_a, ci2.fraction AS fraction_b
FROM analysis.platinum_set ps 
JOIN analysis.drivers_status_snv_long dr1 ON ps.dna_barcode_a = dr1.aliquot_barcode
JOIN analysis.drivers_status_snv_long dr2 ON ps.dna_barcode_b = dr2.aliquot_barcode
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

for(i in 1:length(cells))
{
	for(j in 1:length(drivers))
	{
		sub_dat <- dat %>%
			filter(cell_state == cells[i], grepl(drivers[j], snv_driver_a))
			
		gene_ind_a <- sapply(strsplit(sub_dat[,"snv_driver_a"],","),function(x)grep(drivers[j],x))
		gene_ccf_a <- as.numeric(sapply(strsplit(sub_dat[,"snv_driver_ccf_a"],","),function(x)x[gene_ind_a]))
		
		gene_ind_b <- sapply(strsplit(sub_dat[,"snv_driver_b"],","),function(x)grep(drivers[j],x))
		gene_ccf_b <- as.numeric(sapply(strsplit(sub_dat[,"snv_driver_ccf_b"],","),function(x)x[gene_ind_a]))
		
		if(sum(!is.na(gene_ccf_a) < 5) | sum(!is.na(gene_ccf_a) < 5)){
			next}
			
		
	}
}



