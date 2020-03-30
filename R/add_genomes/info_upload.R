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

info <- read.delim(args[1],header=FALSE,stringsAsFactor=FALSE)
info[which(info[,2]=='X'),2] <- 23
info[,2] <- as.numeric(info[,2])

q <- "
SELECT * 
FROM variants.anno
WHERE chrom = %i AND pos = int4range(%s,'[]') AND alt = '%s'"

variant_id <- rep(0,nrow(info))
for(i in 1:nrow(info))
{
	cat(i,"\r")
	mychrom <- info[i,2]
	
	mypos <- info[i,3]
	mypos <- gsub("[[]","",mypos)
	mypos <- gsub("]","",mypos)
	
	myalt <- info[i,4]
	
	myquery <- sprintf(q,mychrom,mypos,myalt)
	tab <- dbGetQuery(con,myquery)
	myvid <- as.numeric(tab[,1])
	variant_id[i] <- myvid
}

new_info <- cbind(info,variant_id)
new_info <- new_info[,c(1,ncol(new_info),2:(ncol(new_info)-1))]
colnames(new_info) <- c("case_barcode","variant_id","chrom","pos","alt","filter","nlod","tlod")

dbWriteTable(con, Id(schema="variants", table="info"), new_info, append=TRUE)
