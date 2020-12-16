#######################################################
# Add ssm2 coverage counts to db
# Date: 2019.12.20 
# Author: Frederick Varn
#######################################################

## Parse snakemake
if(exists("snakemake")) {
  files = snakemake@input[["metrics"]]
  outfn = snakemake@output[["tsv"]]
} else {
  files = list.files("results/mutect2/geno2db", recursive = T, pattern = ".mfcov.tsv", full.names = T) # list("results/align/wgsmetrics/GLSS-DK-0012-NB-01D-WXS-ABCB18.WgsMetrics.txt", "results/align/wgsmetrics/GLSS-DK-0003-TP-01D-WXS-E43D26.WgsMetrics.txt")
}

outfn <- "results/align/GLSS-H2.mfcov.merged.tsv"
files <- files[grep("GLSS-HF-", files)]
files <- files[-grep("GLSS-HF-2", files)]
files <- files[-grep("GLSS-HF-3", files)]

# Necessary packages:
library(parallel)
library(tidyverse)
library(data.table)
library(DBI)
library(odbc)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

for(i in 1:length(files))
{
	mfcov <- read.delim(files[i],header=FALSE,stringsAsFactor=FALSE)
	colnames(mfcov) <- c("aliquot_barcode","ad_depth","ssm2_call_count")
	dbWriteTable(con, Id(schema="variants", table="ssm2_count"), mfcov, append=TRUE)
}
