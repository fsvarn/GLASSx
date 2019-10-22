#!/usr/bin/env Rscript

#######################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Please provide an input and a mapping file", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}

transcript_tpm <- args[1]
transcript_to_genes <- args[2]
output <- args[3]

cat("Reading in inputs...")
ttpm <- read.delim(transcript_tpm,stringsAsFactors=FALSE)

g2t <- read.delim(transcript_to_genes,header=FALSE,stringsAsFactors=FALSE)
genes <- unique(g2t[,3])
res <- matrix(0,nrow=length(genes),ncol=(ncol(ttpm)-2))
rownames(res) <- genes
colnames(res) <- colnames(ttpm[3:ncol(ttpm)])

cat("\nCollapsing transcripts into genes...\n")
for(i in 1:length(genes))
{
	cat("\r",trunc(i/length(genes)*100),"%")
	tmp_transcripts <- g2t[which(g2t[,3]==genes[i]),]
	tmp_tpm <- ttpm[which(ttpm[,"target_id"]%in%tmp_transcripts[,1]),3:ncol(ttpm)]

	res[i,] <- apply(tmp_tpm,2,sum)
}

cat("\nWriting output...\n")
res <- cbind(genes, res)
colnames(res)[1] <- "Gene_symbol"
write.table(res,output,,sep="\t",quote=FALSE,row.names=FALSE)
#dbWriteTable(con, Id(schema="analysis", table="gene_tpm"), res, overwrite=TRUE)

