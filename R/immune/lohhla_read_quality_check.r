library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
myDir1 <- "results/lohhla/final/"
myDir1 <- paste(myDir1,dir(myDir1),sep="")
myinf1 <- unlist(sapply(myDir1,function(x)paste(x,dir(x),sep="/")))
myinf1 <- myinf1[-grep(".txt",myinf1)]

#Using 20 as the coverage filter, choose the 20 coverage filter run here
myinf1 <- myinf1[grep("/10$",myinf1)]

myinf1 <- unlist(sapply(myinf1,function(x)paste(x,dir(x),sep="/")))
myinf1 <- myinf1[-grep(".bed",myinf1)]
myinf1 <- myinf1[-grep(".xls",myinf1)]
myinf1 <- myinf1[-grep(".nix",myinf1)]
myinf1 <- myinf1[-grep(".fa",myinf1)]
myinf1 <- myinf1[-grep("flagstat",myinf1)]
myinf1 <- myinf1[-grep(".txt",myinf1)]
myinf1 <- myinf1[-grep("Figures",myinf1)]
myinf1 <- myinf1[-grep("50mer_uniq",myinf1)]
myinf1 <- myinf1[-grep("TCGA-19-0957",myinf1)]

myinf1 <- unlist(sapply(myinf1,function(x)paste(x,dir(x),sep="/")))
myinf1 <- myinf1[grep("chr6region.patient.reference.hlas.metrics",myinf1)]

aliquot_barcode <- sapply(strsplit(myinf1,"/"),function(x)x[6])
hash_length <- step_size  <- paired_reads <- proper_pairs <- read_sequences <- unique_alignment <- multi_mapped <- no_mapping_found <- read_length <- rep("",length(myinf1))
for(i in 1:length(myinf1))
{
	myfile <- readLines(myinf1[i])
	hl <- myfile[grep("Hash length",myfile)]
	hash_length[i] <- sapply(strsplit(hl," "),function(x)x[length(x)])
	
	ss <- myfile[grep("Step size",myfile)]
	step_size[i] <- sapply(strsplit(ss," "),function(x)x[length(x)])

	pr <- myfile[grep("Paired Reads",myfile)]
	paired_reads[i] <- sapply(strsplit(pr," "),function(x)x[length(x)])

	pp <- myfile[grep("Proper Pairs",myfile)]
	proper_pairs[i] <- sapply(strsplit(pp," "),function(x)x[length(x)-1])

	rs <- myfile[grep("Read Sequences",myfile)]
	read_sequences[i] <- sapply(strsplit(rs," "),function(x)x[length(x)])

	ua <- myfile[grep("Unique Alignment",myfile)]
	unique_alignment[i] <- sapply(strsplit(ua," "),function(x)x[length(x)-1])

# 	mm <- myfile[grep("Multi Mapped",myfile)]
# 	mm <- strsplit(mm," ")[[1]]
# 	mm <- mm[-grep("\\(",mm)]
# 	multi_mapped[i] <- mm[length(mm)-1]
# 
# 	nm <- myfile[grep("No Mapping Found",myfile)]
# 	nm <- strsplit(nm," ")[[1]]
# 	nm <- nm[-grep("\\(",nm)]
# 	no_mapping_found[i] <-  nm[length(nm)-1]

	rl <- myfile[grep("Mean",myfile)]
	rl <- sapply(strsplit(rl,","),function(x)x[1])
	rl <- gsub("# Mean","", rl)
	read_length[i] <- gsub(" ","", rl)
}

res <- data.frame(aliquot_barcode, hash_length, step_size, paired_reads, proper_pairs, read_sequences, unique_alignment, read_length)
res[,2:ncol(res)] <- apply(res[,2:ncol(res)],2,as.numeric)
res[,1] <- as.character(res[,1])
rownames(res) <- NULL
res <- res[-which(duplicated(res)),]


#Get tumor pairs table
mypairs <- dbReadTable(con,Id(schema="analysis",table="pairs"))
mypairs <- mypairs[which(mypairs[,"tumor_barcode"] %in% res[,"aliquot_barcode"] & (mypairs[,"normal_barcode"] %in% res[,"aliquot_barcode"])),]

#Set up tumor-normal table
comp <- list()
for(i in 1:nrow(mypairs))
{
	mytumor <- mypairs[i,2]
	mynormal <- mypairs[i,3]
	
	c1 <- res[which(res[,"aliquot_barcode"] == mytumor),c(1,4,5,7)]
	c2 <- res[which(res[,"aliquot_barcode"] == mynormal),c(1,4,5,7)]
	colnames(c1) <- paste("tumor_",colnames(c1),sep="")
	colnames(c1)[1] <- "tumor_barcode"
	colnames(c2) <- paste("normal_",colnames(c2),sep="")
	colnames(c2)[1] <- "normal_barcode"
	
	if(nrow(c2) > 1){cat(mynormal,"\n")}

	mycomp <- cbind(c1,c2)
	mycomp <- mycomp[,c(1,5,2,3,4,6,7,8)]
	
	mean_paired_reads <- mean(c(mycomp[,"tumor_paired_reads"],mycomp[,"normal_paired_reads"]))
	mean_proper_pairs <- mean(c(mycomp[,"tumor_proper_pairs"],mycomp[,"normal_proper_pairs"]))
	mean_unique_alignment <- mean(c(mycomp[,"tumor_unique_alignment"],mycomp[,"normal_unique_alignment"]))
	
	mycomp <- cbind(mycomp,mean_paired_reads,mean_proper_pairs,mean_unique_alignment)
	
	comp[[i]] <- mycomp
}

comp <- do.call(rbind, comp)


#Get the average stats for the tumor pair as a whole (normal is counted twice since it is used in both primary and recurrent)
mytumorpairs <- dbReadTable(con,Id(schema="analysis",table="tumor_pairs"))
mytumorpairs <- mytumorpairs[which(mytumorpairs[,"tumor_barcode_a"] %in% comp[,"tumor_barcode"] & (mytumorpairs[,"tumor_barcode_b"] %in% comp[,"tumor_barcode"])),]
tumor_comp <- list()

for(i in 1:nrow(mytumorpairs))
{
	mytumor_a <- mytumorpairs[i,3]
	mytumor_b <- mytumorpairs[i,4]
	c1 <- comp[which(comp[,"tumor_barcode"] == mytumor_a),c(1,9,10,11)]
	c2 <- comp[which(comp[,"tumor_barcode"] == mytumor_b),c(1,9,10,11)]
	colnames(c1) <- paste(colnames(c1),"_a",sep="")
	colnames(c2) <- paste(colnames(c2),"_b",sep="")
	
	mycomp <- cbind(c1,c2)
	mycomp <- mycomp[,c(1,5,2,3,4,6,7,8)]
	
	mean_paired_reads <- mean(c(mycomp[,"mean_paired_reads_a"],mycomp[,"mean_paired_reads_b"]))
	mean_proper_pairs <- mean(c(mycomp[,"mean_proper_pairs_a"],mycomp[,"mean_proper_pairs_b"]))
	mean_unique_alignment <- mean(c(mycomp[,"mean_unique_alignment_a"],mycomp[,"mean_unique_alignment_b"]))
	
	mycomp <- cbind(mycomp,mean_paired_reads,mean_proper_pairs,mean_unique_alignment)
	
	tumor_comp[[i]] <- mycomp
	
}

tumor_comp <- do.call(rbind, tumor_comp)

#For all aliquots that have matching WXS, WGS, identify the better sample based on average normal/tumor read count
portion_barcode <- sapply(strsplit(tumor_comp[,"tumor_barcode_a"],"-"),function(x)paste(x[1:5],collapse='-'))
portion_barcode <- unique(portion_barcode)
keep <- list()

for(i in 1:length(portion_barcode))
{
	myportion <- tumor_comp[grep(portion_barcode[i],tumor_comp[,"tumor_barcode_a"]),]
	keep[[i]] <- myportion[which(myportion[,"mean_proper_pairs"]==max(myportion[,"mean_proper_pairs"])),]
}

keep_list <- do.call(rbind, keep)
keep_barcodes <- keep_list[,c(1,2)]

portion_barcode <- sapply(strsplit(keep_barcodes[,"tumor_barcode_a"],"-"),function(x)paste(x[1:5],collapse='-'))
keep_barcodes <- cbind(portion_barcode, keep_barcodes)

#Intersect with the gold set to get the full "lohhla_set"

diamond_set <- dbReadTable(con,Id(schema="analysis",table="diamond_set"))
gold_portion <-  sapply(strsplit(diamond_set[,"tumor_barcode_a"],"-"),function(x)paste(x[1:5],collapse='-'))

lohhla_set <- keep_barcodes[which(keep_barcodes[,"portion_barcode"] %in% intersect(keep_barcodes[,"portion_barcode"],gold_portion)),]

tumor_pair_barcode <- mytumorpairs[match(lohhla_set[,"tumor_barcode_a"],mytumorpairs[,"tumor_barcode_a"]),"tumor_pair_barcode"]
case_barcode <- mytumorpairs[match(lohhla_set[,"tumor_barcode_a"],mytumorpairs[,"tumor_barcode_a"]),"case_barcode"]

lohhla_set <- cbind(tumor_pair_barcode, case_barcode, lohhla_set)

dbWriteTable(con, Id(schema="analysis", table="lohhla_set"), lohhla_set, overwrite=TRUE,row.names=FALSE)
