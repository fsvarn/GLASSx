#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

titan_file <- args[1]
aliquot_barcode <- args[2]
solutions_out <- args[3]

# Interactive run
#titan_file <- "/projects/varnf/GLASS-III/GLASS-III/results/cnv/titanfinal/params/GLSS-CU-R005-TP-01-NB-01D-WXS.params.txt"

conIn <- file(titan_file, "r")
titan_params <- readLines(conIn)

normal_contam <- titan_params[grep("Normal contamination",titan_params)]
normal_contam <- as.numeric(sapply(strsplit(normal_contam,"\t"),function(x)x[2]))
tumorPurity <- 1 - normal_contam

tumorPloidy <- titan_params[grep("Average tumour ploidy",titan_params)]
tumorPloidy <- as.numeric(sapply(strsplit(tumorPloidy,"\t"),function(x)x[2]))

Ploidy <- 2		#I believe this is normal ploidy
 
res <- data.frame(Ploidy, tumorPurity, tumorPloidy)
rownames(res) <- aliquot_barcode

write.table(res, solutions_out, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)