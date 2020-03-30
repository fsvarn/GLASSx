#PBS -l walltime=72:00:00,mem=32gb
/projects/varnf/SofWar/anaconda3/bin/R --vanilla <<EOF

library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

dat <- read.delim("/projects/varnf/GLASS-III/GLASS-III/results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.tsv", colClasses = c("character","character",NULL,"character",rep(NULL,14)),header=FALSE)
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

myoutf <- "/projects/varnf/GLASS-III/GLASS-III/results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.new.tsv"
write.table(new_dat,myoutf,sep="\t",quote=F,row.names=F)

dbWriteTable(con, Id(schema="variants", table="anno"), new_dat, append=TRUE)

EOF