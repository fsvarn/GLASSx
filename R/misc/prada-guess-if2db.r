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
myinf1 <- paste(myDir1, aliquot, "/guess-if/EGFR/EGFR.GUESS-IF.summary.txt",sep="")

guess_list <- list()
for(i in 1:length(myinf1))
{
	if(!file.exists(myinf1[i])){
		next}
		
	guess <- readLines(myinf1[i])
	mystart <- grep(">junction", guess) + 1
	myend <- grep(">summary", guess)- 2
	
	if(mystart > myend){
		next
#		Currently excluding samples with 0 fusion reads from the table
# 		chr_a <- pos_a <- chr_b <- pos_b <- NA
# 		reads <- guess[grep("Number of Fusion Reads", guess)]
# 		reads <- sapply(strsplit(reads, " = "),function(x)as.numeric(x[2]))
	} else{
		intra_fus <- guess[mystart:myend]

		intra_fus <- strsplit(intra_fus,"\t")
	
		junction <- sapply(intra_fus, function(x)x[1])
		junction <- strsplit(junction, "\\.")
	
		chr_a <- sapply(junction, function(x) as.numeric(gsub("chr","", x[1])))
		pos_a <- sapply(junction, function(x) as.numeric(x[2]))
		chr_b <- sapply(junction, function(x) as.numeric(gsub("chr","", x[3])))
		pos_b <- sapply(junction, function(x) as.numeric(x[4]))
	
		reads <- as.numeric(sapply(intra_fus, function(x)x[2]))
	}
	
	aliquot_barcode <- rep(aliquot[i], length(reads))
	gene_symbol <- rep("EGFR", length(reads))
	guess_list[[i]] <- data.frame(aliquot_barcode, gene_symbol, chr_a, pos_a, chr_b, pos_b, reads)
}

guess_if_table <- do.call(rbind, guess_list)

dbWriteTable(con, Id(schema="variants",table="prada_guess_if"), guess_if_table, overwrite=T)