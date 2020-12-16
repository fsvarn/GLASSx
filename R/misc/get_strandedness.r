library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

# Read in each the STAR ReadsPerGene output for each sample

myDir1 <- "results/rnafingerprint/star/"
mytag <- dir(myDir1)
myDir1 <- paste(myDir1, mytag, sep="")
myinf1 <- sapply(myDir1, function(x)paste(x, dir(x), sep = "/"))
myinf1 <- unlist(myinf1)
myinf1 <- myinf1[grep("ReadsPerGene.out.tab", myinf1)]

fraction3 <- fraction4 <- feat_diff <- rep(0, length(myinf1))
for(i in 1:length(myinf1))
{
	cat("\r",i)
	rpg <- read.delim(myinf1[i], header=FALSE)
	
	N_noFeature3 <- rpg[3,3]
	N_noFeature4 <- rpg[3,4]
	feat_diff[i] <- N_noFeature3 - N_noFeature4
	rpg <- rpg[5:nrow(rpg),]
	
	sums <- apply(rpg[,2:4],2,sum)
	fraction3[i] <- sums[2]/sums[1]
	fraction4[i] <- sums[3]/sums[1]
}
file_stem <- sapply(strsplit(myinf1,"/"),function(x)x[5])
aliquot_barcode <- sapply(strsplit(file_stem,"\\."),function(x)x[1])
readgroup_idtag <- sapply(strsplit(file_stem,"\\."),function(x)paste(x[2:3],collapse="."))
names(aliquot_barcode) <- names(readgroup_idtag) <- NULL

# Test alternative strand-checking approach
# test <- rep("", length(file_stem))
# test[which(abs(feat_diff) < 1000000)] <- "unstranded"
# test[which(feat_diff > 1000000)] <- "rf-stranded"
# test[which(feat_diff < -1000000)] <- "fr-stranded"

# Set rules for stranded/unstranded 
library_prep <- rep("", length(file_stem))
#library_prep[which(fraction3 > 0.6 | fraction4 < 0.4)] <- "fr-stranded"
#library_prep[which(fraction3 < 0.4 | fraction4 > 0.6)] <- "rf-stranded"
library_prep[which(fraction3 > 0.2 & fraction4 > 0.2)] <- "unstranded"
library_prep[which(fraction3 > 0.8)] <- "fr-stranded"
library_prep[which(fraction4 > 0.8)] <- "rf-stranded"

res <- data.frame(aliquot_barcode, readgroup_idtag, library_prep)
colnames(res) <- c("aliquot_barcode", "readgroup_idtag", "library_prep")

#Manually adjust strandedness based on prior knowledge (applies to LU FFPE RNAseq samples)
res[grep("GLSS-LU-",res[,"aliquot_barcode"]),"library_prep"] <- "rf-stranded"

write.table(res, "results/rnafingerprint/library_prep_strandedness.txt", sep = "\t", quote=FALSE, row.names=FALSE)

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Write results to db
dbWriteTable(con, Id(schema="analysis",table="rna_library_prep"), res, overwrite = TRUE)

