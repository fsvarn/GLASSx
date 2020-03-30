#!/usr/bin/env Rscript

#######################################################
library(odbc)
library(DBI)
#######################################################

rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

geno <- read.delim(args[1],header=FALSE,stringsAsFactor=FALSE)
geno[which(geno[,3]=='X'),3] <- 23
geno[,3] <- as.numeric(geno[,3])

q <- "
SELECT * 
FROM variants.anno
WHERE chrom = %i AND pos = int4range(%s,'[]') AND alt = '%s'"

variant_id <- rep(0,nrow(geno))
for(i in 1:nrow(geno))
{
	cat(i,"\r")
	mychrom <- geno[i,3]
	
	mypos <- geno[i,4]
	mypos <- gsub("[[]","",mypos)
	mypos <- gsub("]","",mypos)
	
	myalt <- geno[i,5]
	
	myquery <- sprintf(q,mychrom,mypos,myalt)
	tab <- dbGetQuery(con,myquery)
	myvid <- as.numeric(tab[,1])
	variant_id[i] <- myvid
}

new_geno <- cbind(geno,variant_id)
new_geno <- new_geno[,c(1,ncol(new_geno),2:(ncol(new_geno)-1))]
colnames(new_geno) <- c("aliquot_barcode","variant_id","case_barcode","chrom","pos","alt","ad_ref","ad_alt","af","ssm2_pass_call","ssm2_saaf_none")

dbWriteTable(con, Id(schema="variants", table="geno"), new_geno, append=TRUE)
