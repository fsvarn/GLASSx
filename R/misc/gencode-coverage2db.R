#######################################################
# Necessary packages:
library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

# Establish connection with db
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")


myDir1 <- "results/align/gencode-coverage/"
myinf1 <- dir(myDir1,)

#Specify cohort of interest if adding one specific cohort
#myinf1 <- myinf1[grep("GLSS-SN-",myinf1)]

aliquot_barcode <- gsub(".gencode-coverage.tsv","",myinf1)
myinf1 <- paste(myDir1,myinf1,sep="")

dat <- list()
for(i in 1:length(myinf1))
{
	aliquot_gencode <- read.delim(myinf1[i],header=FALSE)
	aliquot_gencode <- cbind(aliquot_barcode[i], aliquot_gencode)
	dat[[i]] <- aliquot_gencode
}

new_gencode <- do.call(rbind,dat)

#Remove sum of exon sizes column:
final_gencode <- new_gencode[,c(1,2,4)]
colnames(final_gencode) <- c("aliquot_barcode","ensembl_gene_id","gene_coverage")

dbWriteTable(con, Id(schema="analysis",table="gencode_coverage"), final_gencode, append=T)