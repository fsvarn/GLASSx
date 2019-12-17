#!/usr/bin/env Rscript

#######################################################
library(odbc)
library(DBI)
library(ssgsea.GBM.classification)
#######################################################

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

gct_path <- args[1]
class_out <- args[2]

#Install Qianghu's classifier
#install.packages(pkgs="/projects/varnf/SofWar/R/ssgsea.GBM.classification/",repos=NULL)

#Run Qianghu's transcriptional classifier
runSsGSEAwithPermutation(gct_path,100)

#Read in results from classifier and upload to db

#Open db connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

res <- read.delim(class_out)
aliquot_barcode <- rep(rownames(res),3)
aliquot_barcode <- gsub("\\.","-",aliquot_barcode)
signature_name <- c(rep("Proneural",nrow(res)),rep("Classical",nrow(res)),rep("Mesenchymal",nrow(res)))
enrichment_score <- c(res[,"Proneural"],res[,"Classical"],res[,"Mesenchymal"])
p_value <- c(res[,"Proneural_pval"],res[,"Classical_pval"],res[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(aliquot_barcode,signature_name,enrichment_score,p_value)
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"aliquot_barcode"]),]

dbWriteTable(con, Id(schema="analysis", table="transcriptional_subtype"), transcriptional_subtype, overwrite=TRUE, row.names=FALSE)
