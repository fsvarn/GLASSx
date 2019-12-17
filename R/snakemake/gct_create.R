#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

matrix_path <- args[1]
gct_out <- args[2]

#Create gct file of gene expression matrix
data <- read.delim(matrix_path)
Description <- rep(NA,nrow(data))
data2 <- cbind(data[,1],Description,data[,2:ncol(data)])
colnames(data2)[1] <- "NAME"
write.table(data2, gct_out, sep="\t", quote=FALSE, row.names=FALSE)

conIn <- file(gct_out, "r")
rawfile = readLines(conIn)
close(conIn)

mytext <- c("#1.2", paste(nrow(data2),"\t",(ncol(data)-1),sep=""),rawfile)
conOut = file(gct_out, "w")
writeLines(mytext, conOut)
close(conOut)