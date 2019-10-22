library(ggplot2)
library(odbc)
library(DBI)
library(rjson)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

myDir1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/aliquot"
myinf1 <- paste(myDir1,dir(myDir1),sep="/")
myinf1 <- sapply(myinf1,function(x)paste(x,dir(x)[3],sep="/"),USE.NAMES=FALSE)

jsons <- lapply(myinf1,function(x)fromJSON(file=x))
names(jsons) <- dir(myDir1)

qc_table <- do.call(rbind,jsons)

qc_table <- as.data.frame(apply(qc_table,2,unlist))
qc_table[,1:7] <- apply(qc_table[,1:7],2,as.numeric)
range(qc_table[,"p_pseudoaligned"])
mean(qc_table[,"p_pseudoaligned"])
median(qc_table[,"p_pseudoaligned"])

qc_table[which(qc_table[,"p_pseudoaligned"] < 50),4:6]

qc_table <- qc_table[,1:7]
qc_table <- cbind(rownames(qc_table),qc_table)
rownames(qc_table) <- NULL
colnames(qc_table)[1] <- "aliquot_barcode"
dbWriteTable(con, Id(schema="analysis", table="kallisto_qc"), qc_table, overwrite=TRUE)