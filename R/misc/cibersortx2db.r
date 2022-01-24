###################################################
# Upload CIBERSORTx results calculated using CIBERSORTx results to db 
# Previous uploads are commented out in order
# Updated: 2020.08.05
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Final data freeze data 
#myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/CIBERSORTxGEP_GLASS_Fractions-Adjusted.txt"
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS_neftel_freeze/CIBERSORTxGEP_GLASS_Fractions-Adjusted.txt"
#myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS_venteicher_freeze/CIBERSORTxGEP_GLASS_Fractions-Adjusted.txt"
#myinf1 <- "/projects/verhaak-lab/GLASS-III/data/cibersortx/ivyGAP/CIBERSORTx_GLASS_20210923_ivyGAP_fractions.txt"

# Original GLASS Neftel data from August
#myinf1 <-  "/projects/verhaak-lab/GLASS-III/data/res/neftel_GLASS_cibersortx_prop_08052020.txt"
# Original GLASS SCGP data from June
#myinf1 <- "/projects/verhaak-lab/GLASS-III/data/res/10x_GLASS_cibersortx_prop_06112020.txt"

prop <-read.delim(myinf1)
prop <- prop[,1:(ncol(prop)-3)]

prop[,"Mixture"] <- gsub("\\.","-",prop[,"Mixture"])

db_table <- prop %>% pivot_longer(-Mixture, values_to = "fraction", names_to = "cell_state")
colnames(db_table)[1] <- c("aliquot_barcode")

#dbWriteTable(con, Id(schema="analysis", table="cibersortx_ivygap"), db_table, overwrite=TRUE)
dbWriteTable(con, Id(schema="analysis", table="cibersortx_neftel"), db_table, overwrite=TRUE)
#dbWriteTable(con, Id(schema="analysis", table="cibersortx_scgp"), db_table, overwrite=TRUE)

