#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

hla_file <- args[1]
hla_out <- args[2]

# Interactive run
# hla_file <- "/projects/varnf/GLASS-III/GLASS-III/results/optitype/HLA_calls/GLSS-CU-R010-NB-01D-WXS-0V2U96/GLSS-CU-R010-NB-01D-WXS-0V2U96_extended_result.tsv"
# hla_out <- "/projects/varnf/GLASS-III/analysis/lohhla_test/GLSS-CU-R010-TP-01-R1-01D-WXS/hlas"

hla_table <- read.delim(hla_file,stringsAsFactor=FALSE)

hlas <- hla_table[,2:(ncol(hla_table)-2)]
hlas <- as.character(hlas)

#Reformat for LOHHLA input
hlas <- tolower(hlas)
hlas <- gsub("\\*","_",hlas)
hlas <- gsub("\\:","_",hlas)
hlas <- paste("hla_",hlas,sep="")

conOut = file(hla_out, "w")
writeLines(hlas, conOut)
close(conOut)
