#######################################################
# Necessary packages:
library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

# Read in selected transcripts file
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

selected_transcripts <- dbReadTable(con, Id(schema="ref", table="ensembl_prada_transcripts"))

# Read in exon reference file
ref_file <- "/projects/varnf/SofWar/pyPRADA_1.2/ref/PRADA-reference-hg19/Homo_sapiens.GRCh37.64.gtf"
exon_gtf <- read.delim(ref_file, header=FALSE, stringsAsFactor=FALSE)
exon_gtf <- exon_gtf[which(exon_gtf[,1] %in% c(as.character(1:22),"X")),]
exon_gtf <- exon_gtf[which(exon_gtf[,3] == "CDS"),]

transcript_info <- strsplit(exon_gtf[,9],";")
ensembl_transcript_id <- sapply(transcript_info,function(x)x[2])
ensembl_transcript_id <- gsub(" transcript_id ","",ensembl_transcript_id)

exon_number <- sapply(transcript_info,function(x)x[3])
exon_number <- as.numeric(gsub(" exon_number ","",exon_number))
# Zero indexing for downstream analyses
exon_number <- exon_number - 1

gene_name <- sapply(transcript_info,function(x)x[4])
gene_name <- gsub(" gene_name ","",gene_name)

exon_id <- paste(ensembl_transcript_id, exon_number, sep = "_")
chrom <- exon_gtf[,1]
chrom[which(chrom=="X")] <- "23"
chrom <- as.numeric(chrom)
start <- exon_gtf[,4]
end <- exon_gtf[,5]
exon_link <- data.frame(exon_id, ensembl_transcript_id, gene_name, chrom, start, end)

#exon_link <- exon_link[which(exon_link[,"ensembl_transcript_id"] %in% selected_transcripts[,"ensembl_transcript_id"]),]

# Create preprocessing file list
myDir1 <- "results/prada/"
aliquot <- dir(myDir1)
myinf1 <- paste(myDir1, aliquot, "/preprocessing/",aliquot,"/", aliquot, "/", aliquot, ".metrics.tmp.txt.intronReport.txt_exonOnly.txt", sep="")

# Read in fusions and create fusion table
myDir2 <- "results/prada/"
aliquot <- dir(myDir2)
myinf2 <- paste(myDir2, aliquot, "/fusion/prada.fus.summary.txt",sep="")

fusion_list <- list()
for(i in 1:length(myinf2))
{
	if(!file.exists(myinf2[i])){
		next}
		
	fus <- read.delim(myinf2[i])
	colnames(fus) <- tolower(colnames(fus))
	colnames(fus)[1] <- "gene_symbol_a"
	colnames(fus)[2] <- "gene_symbol_b"
	colnames(fus)[3] <- "chr_a"
	colnames(fus)[4] <- "chr_b"

	aliquot_barcode <- rep(aliquot[i], nrow(fus))

	fus <- cbind(aliquot_barcode, fus)

	fusion_list[[i]] <- fus
}

fusion_table <- do.call(rbind, fusion_list)

fusion_table[which(fusion_table[,"chr_a"]=="X"),"chr_a"] <- 23
fusion_table[which(fusion_table[,"chr_b"]=="X"),"chr_b"] <- 23

fusion_table[,"chr_a"] <- as.integer(fusion_table[,"chr_a"])
fusion_table[,"chr_b"] <- as.integer(fusion_table[,"chr_b"])


# Calculate transcript allele fraction:

final_taf <- rep("", nrow(fusion_table))
ct <- 1
for(i in 1:length(myinf1))
{
	exon_lines <- readLines(myinf1[i])
	exon_head <- exon_lines[1]
	exon_head <- unlist(strsplit(exon_head,"\t"))
	
	# Additional for loop for each unique fusion within a patient
	sub_fusion <- fusion_table[which(fusion_table[,1] == aliquot[i]),]
	
	for(j in 1:nrow(sub_fusion))
	{
		myjunction <- as.character(sub_fusion[j,"junction"])
		myjunction <- unlist(strsplit(myjunction,"\\|"))
	
		taf_a <- taf_b <- rep("", length(myjunction))
		for(k in 1:length(myjunction))
		{
			cat("\r",ct)
			sub_junction <- myjunction[k]
			sub_junction <- unlist(strsplit(sub_junction,","))
			numerator <- as.numeric(sub_junction[2])
		
			join_a <- unlist(strsplit(sub_junction[1],"_"))[1]
			join_b <- unlist(strsplit(sub_junction[1],"_"))[2]
		
			gene_a <- unlist(strsplit(join_a,"\\:"))[1]
			chr_a <- as.numeric(unlist(strsplit(join_a,"\\:"))[2])
			pos_a <- as.numeric(unlist(strsplit(join_a,"\\:"))[3])
		
			index_a <- exon_link[which(exon_link[,"chrom"] == chr_a &
							 exon_link[,"start"] <= pos_a & exon_link[,"end"] >= pos_a),]
			exon_a <- as.character(index_a[,"exon_id"])
		
			exon_data_a <- exon_lines[grep(paste(exon_a,"\t",sep=""), exon_lines)]
			exon_data_a <- unlist(strsplit(exon_data_a,"\t"))
			names(exon_data_a) <- exon_head
		
			denom_a <- as.numeric(exon_data_a["Exon_Reads"])

			gene_b <- unlist(strsplit(join_b,"\\:"))[1]
			chr_b <- as.numeric(unlist(strsplit(join_b,"\\:"))[2])
			pos_b <- as.numeric(unlist(strsplit(join_b,"\\:"))[3])
		
			index_b <- exon_link[which(exon_link[,"chrom"] == chr_b &
							 exon_link[,"start"] <= pos_b & exon_link[,"end"] >= pos_b),]
			exon_b <- as.character(index_b[,"exon_id"])

			exon_data_b <- exon_lines[grep(paste(exon_b,"\t",sep=""), fixed=TRUE, exon_lines)]
			exon_data_b <- unlist(strsplit(exon_data_b,"\t"))
			names(exon_data_b) <- exon_head
		
			denom_b <- as.numeric(exon_data_b["Exon_Reads"])
		
			taf_a[j] <- paste(gene_a, round(numerator/denom_a,3) ,sep=",")
			taf_b[j] <- paste(gene_b, round(numerator/denom_b,3) ,sep=",")
		}
		taf <- paste(taf_a, taf_b, sep="_")
		final_taf[ct] <- paste(taf, collapse="|")
		ct <- ct + 1
	}
}


fusion_table <- data.frame(fusion_table, final_taf)

# Establish connection with db
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

dbWriteTable(con, Id(schema="variants",table="prada_fusions"), fusion_table, overwrite=TRUE)