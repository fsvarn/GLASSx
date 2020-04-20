#######################################################
# Necessary packages:
library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

# Establish connection with db
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

myDir1 <- "results/prada/"
aliquot <- dir(myDir1)
myinf1 <- paste(myDir1, aliquot, "/fusion/prada.fus.summary.txt",sep="")

fusion_list <- list()
for(i in 1:length(myinf1))
{
	if(!file.exists(myinf1[i])){
		next}
		
	fus <- read.delim(myinf1[i])
	colnames(fus) <- tolower(colnames(fus))
	colnames(fus)[1] <- "gene_symbol_a"
	colnames(fus)[2] <- "gene_symbol_b"
	colnames(fus)[3] <- "chr_a"
	colnames(fus)[4] <- "chr_b"

	aliquot_barcode <- rep(aliquot[i], nrow(fus))

	fus <- cbind(aliquot_barcode, fus)

	fusion_list[[i]] <- fus
}

fusion_table <- do.call(rbind, fusion_list)

dbWriteTable(con, Id(schema="variants",table="prada_fusions"), fusion_table, overwrite=TRUE)