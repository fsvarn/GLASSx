#######################################################
library(odbc)
library(DBI)
library(tidyverse)
library(stringi)
#######################################################

rm(list=ls())

file_sizes <- read.delim("/Users/varnf/Documents/Projects/GLASS/CU_anti-PD1/aligned_bams/CU_aPD1_bam_file_size.txt",header=FALSE,row.names=NULL,stringsAsFactor=FALSE)
md5sums <- read.delim("/Users/varnf/Documents/Projects/GLASS/CU_anti-PD1/aligned_bams/CU_aPD1_bam_md5sums.txt",sep="",header=FALSE,row.names=NULL,stringsAsFactor=FALSE)

colnames(file_sizes) <- c("file_size","file_path")
colnames(md5sums) <- c("file_md5sum","file_path")

#build file_additions table
file_additions <- merge(file_sizes,md5sums,by="file_path")
file_additions[,"file_name"] <- sapply(strsplit(file_additions[,"file_path"],"/"),function(x)x[length(x)])  
#file_additions[,"file_path"] <- sapply(strsplit(file_additions[,"file_path"],"/"),function(x)paste(x[1:length(x)-1],collapse="/"))  
file_additions[,"file_format"] <- rep("aligned BAM",nrow(file_additions))

#Extract aliquot barcode from file name
aliquot_barcode <- sapply(strsplit(file_additions[,"file_path"],"/"),function(x)x[length(x)])  
aliquot_barcode <- gsub(".realn.mdup.bqsr.bam","",aliquot_barcode)
file_additions[,"aliquot_barcode"] <- aliquot_barcode

#reorder columns in files table
file_additions <- file_additions[,c("aliquot_barcode","file_name","file_size","file_md5sum","file_format","file_path")]

write.table(files_readgroups,"/Users/varnf/Documents/Projects/GLASS/CU_anti-PD1/aligned_bams/files.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

#Upload
#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

### Files ###
files_master <- file_additions %>% mutate(aliquot_barcode = aliquot_barcode,
                                 file_name = file_name,
                                 file_size = file_size,
                                 file_md5sum = file_md5sum,
                                 file_format = file_format)

#Write to database
dbWriteTable(con, Id(schema="analysis", table="files"), files_master, append=TRUE)


