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
  stop("Please provide an input", call.=FALSE)}

mixcr <- args[1]
receptor <- args[2]

cat("Reading in inputs...")
mxcr <- read.delim(mixcr,stringsAsFactors=FALSE,
colClasses=c(aliquot_barcode="character", cloneId="numeric",cloneCount="numeric",
				 cloneFraction="numeric", targetSequences="character",
				 targetQualities="character", allVHitsWithScore="character",
				 allDHitsWithScore="character",allJHitsWithScore="character",
				 allCHitsWithScore="character",allVAlignments="character",
				 allDAlignments="character",allJAlignments="character",
				 allCAlignments="character",nSeqFR1="character",
				 minQualFR1="numeric",nSeqCDR1="character",
				 minQualCDR1="numeric",nSeqFR2="character",
				 minQualFR2="numeric",nSeqCDR2="character",
				 minQualCDR2="numeric",nSeqFR3="character",
				 minQualFR3="numeric",nSeqCDR3="character",
				 minQualCDR3="numeric",nSeqFR4="character",
				 minQualFR4="numeric",aaSeqFR1="character",
				 aaSeqCDR1="character",aaSeqFR2="character",
				 aaSeqCDR2="character",aaSeqFR3="character",
				 aaSeqCDR3="character",aaSeqFR4="character",
				 refPoints="character"))
colnames(mxcr) <- tolower(colnames(mxcr))

table_out <- paste("mixcr_",receptor,sep="")

#dbWriteTable(con, Id(schema="analysis", table=table_out), mxcr, overwrite=TRUE)

