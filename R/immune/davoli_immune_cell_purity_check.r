library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT pa.pair_barcode, im.signature_name, im.enrichment_score, tp.purity
FROM variants.titan_params tp
JOIN analysis.pairs pa ON pa.pair_barcode = tp.pair_barcode
JOIN analysis.analyte_sets an ON an.dna_barcode = pa.tumor_barcode
JOIN analysis.davoli_immune_score im ON im.aliquot_barcode = an.rna_barcode
"

dat <- dbGetQuery(con,q)

mycells <- unique(dat[,"signature_name"])

mycor <- rep(0,length(mycells))
names(mycor) <- mycells
for(i in 1:length(mycells))
{
	sub_dat <- dat[grep(mycells[i],dat[,"signature_name"]),]
	mycor[i] <- cor(sub_dat[,"enrichment_score"],sub_dat[,"purity"],method="s")
}

#            CD4.mature          CD8.effector              NK.cells 
#           -0.10032724           -0.10411075           -0.07383572 
#               B.cells                 T.reg             Dendritic 
#           -0.13181065           -0.10267669           -0.19782280 
# CD8.effector.NK.cells           Macrophages        Macrophages.M2 
#           -0.11801027           -0.08030337           -0.12844719 
#        Macrophages.M1 
#           -0.10434187 


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT pa.pair_barcode, im.signature_name, im.enrichment_score, sp.cellularity
FROM variants.seqz_params sp
JOIN analysis.pairs pa ON pa.pair_barcode = sp.pair_barcode
JOIN analysis.analyte_sets an ON an.dna_barcode = pa.tumor_barcode
JOIN analysis.davoli_immune_score im ON im.aliquot_barcode = an.rna_barcode
"

dat <- dbGetQuery(con,q)

mycells <- unique(dat[,"signature_name"])

mycor <- rep(0,length(mycells))
names(mycor) <- mycells
for(i in 1:length(mycells))
{
	sub_dat <- dat[grep(mycells[i],dat[,"signature_name"]),]
	mycor[i] <- cor(sub_dat[,"enrichment_score"],sub_dat[,"cellularity"],method="s")
}

#            CD4.mature          CD8.effector              NK.cells 
#            -0.2647990            -0.3116131            -0.2947136 
#               B.cells                 T.reg             Dendritic 
#            -0.3017697            -0.1395366            -0.4076268 
# CD8.effector.NK.cells           Macrophages        Macrophages.M2 
#            -0.3582346            -0.2447665            -0.3645780 
#        Macrophages.M1 
#            -0.3023455 