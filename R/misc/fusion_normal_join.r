###################################################
# Test how fusion counts associate with immune levels
# Updated: 2020.04.22
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)


#######################################################
rm(list=ls())

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Load previous PRADA fusion run on GLASS
q <- "
SELECT *
FROM variants.prada_fusions
"

dat <- dbGetQuery(con,q)

# Load fusions called in normal RNA from TCGA (Hu et al, NAR 2018)
myinf1 <- "/Users/varnf/Documents/Data/TCGA/normal_fusions_hu_nar_2018.txt"
norm_fusions <- read.delim(myinf1)

called_in_normal <- dat[,"junction"] %in% norm_fusions[,"Junction"]

dat <- cbind(dat, called_in_normal)

dbWriteTable(con, Id(schema="variants",table="prada_fusions"), dat, overwrite=TRUE)
