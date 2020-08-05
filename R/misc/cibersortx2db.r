###################################################
# Upload CIBERSORTx results calculated using 10x SCGP sigs to db
# Updated: 2020.06.12
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

myinf1 <- "/projects/verhaak-lab/GLASS-III/data/res/10x_GLASS_cibersortx_prop_06112020.txt"

prop <-read.delim(myinf1)
prop <- prop[,1:13]

prop[,"Mixture"] <- gsub("\\.","-",prop[,"Mixture"])

db_table <- prop %>% pivot_longer(-Mixture, values_to = "fraction", names_to = "cell_state")
colnames(db_table)[1] <- c("aliquot_barcode")

dbWriteTable(con, Id(schema="analysis", table="cibersortx_scgp"), db_table, overwrite=TRUE)

