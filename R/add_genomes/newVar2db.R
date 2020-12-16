#!/usr/bin/env Rscript

library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

dat <- read.delim("/projects/verhaak-lab/GLASS-III/results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.tsv", colClasses = c("character","character",NULL,"character",rep(NULL,14)),header=FALSE)
dat[which(dat[,1]=='X'),1] <- 23
dat[,1] <- as.numeric(dat[,1])

q <- "
SELECT * 
FROM variants.anno
WHERE chrom = %i AND pos = int4range(%s,'[]') AND alt = '%s'"

myinds <- c()
for(i in 1:nrow(dat))
{
	cat(i,"\r")
	mychrom <- dat[i,1]
	
	mypos <- dat[i,2]
	mypos <- gsub("[[]","",mypos)
	mypos <- gsub("]","",mypos)
	
	myalt <- dat[i,4]
	
	myquery <- sprintf(q,mychrom,mypos,myalt)
	
	tab <- dbGetQuery(con,myquery)
	
	if(nrow(tab)==0){
		myinds <- c(myinds,i)}
}

new_dat <- dat[myinds,]
colnames(new_dat) <- c("chrom", "pos", "ref","alt","gene_symbol","variant_classification","secondary_variant_classification","variant_type","genome_change","transcript","transcript_strand","transcript_exon","transcript_position","cdna_change","cds_change","protein_change","gc_content","reference_context")

myoutf <- "/projects/verhaak-lab/GLASS-III/results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.GLSS-H2.tsv"
write.table(new_dat,myoutf,sep="\t",quote=F,row.names=F)

dbWriteTable(con, Id(schema="variants", table="anno"), new_dat, append=TRUE)

EOF