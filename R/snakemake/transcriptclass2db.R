#!/usr/bin/env Rscript

install.packages(pkgs="/projects/varnf/SofWar/R/ssgsea.GBM.classification/",repos=NULL)

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

#Run Qianghu's transcriptional classifier
runSsGSEAwithPermutation(gct_path,100)

#Read in results from classifier and convert to long format

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

#Calculate simplicity scores
aliquots <- unique(as.character(transcriptional_subtype[,"aliquot_barcode"]))

simplicity_score <- rep(0,length(aliquots))
for(i in 1:length(aliquots))
{
	sub_dat <- transcriptional_subtype[which(transcriptional_subtype[,"aliquot_barcode"] == aliquots[i]),]
	sub_dat[,"p_rank"] <- rank(sub_dat[,"p_value",],ties.method="min")
	r0 <- sub_dat[which(sub_dat[,"p_rank"] ==1),"p_value"][1]
	ri <- sub_dat[which(sub_dat[,"p_rank"] > 1),"p_value"]
	ri <- ri[order(ri)]
	
	adds <- sum(ri - r0)
	
	d <- abs(outer(ri,ri,"-"))
	diag(d) <- NA
	d[lower.tri(d)] <- NA
	adns <- sum(d,na.rm=TRUE)
	
	rn1 <- sub_dat[which(sub_dat[,"p_rank"] == max(sub_dat[,"p_rank"])),"p_value"][1]
	n1 <- 2	#Number of unique subtypes - 1
	simplicity_score[i] <- (adds - adns) * (rn1 - r0)/n1
}
simplicity_score <- data.frame(aliquots,simplicity_score)
colnames(simplicity_score) <- c("aliquot_barcode","simplicity_score")

#Upload to db
dbWriteTable(con, Id(schema="analysis", table="transcriptional_subtype"), transcriptional_subtype, overwrite=TRUE, row.names=FALSE)
dbWriteTable(con, Id(schema="analysis", table="simplicity_score"), simplicity_score, overwrite=TRUE, row.names=FALSE)
