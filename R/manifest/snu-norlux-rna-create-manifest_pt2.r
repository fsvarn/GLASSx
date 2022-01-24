library(odbc)
library(DBI)
library(tidyverse)
library(stringr)
library(stringi)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/")
set.seed(11)

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

cases <- read.delim("file_metadata/NorLux/norlux-cases-jax-database-20210907-mRNA_exists_pt2.txt",stringsAsFactors=FALSE)
samples <- read.delim("file_metadata/NorLux/norlux-samples-jax-database-plus-linker-20210907-mRNA_exists_pt2.txt",stringsAsFactors=FALSE)

# Cases table
#dbWriteTable(con, Id(schema="clinical",table = "cases"), cases, append = TRUE)

# Samples table
samples_upload <- samples[,c("case_barcode", "sample_barcode")]
sample_type <- substring(samples_upload$sample_barcode, 14, 16)
samples_upload <- data.frame(samples_upload, sample_type)

#dbWriteTable(con, Id(schema="biospecimen",table = "samples"), samples_upload, append = TRUE)

#-------------------------
# Create aliquots table

# Use this script to read in the samples table for a few samples we are adding reads to:
#sample_upload <- dbReadTable(con, Id(schema = "biospecimen", table = "samples"))
#samps2add <- c("GLSS-SN-0009-R1", "GLSS-SN-0013-R2", "GLSS-LX-0357-R3", "GLSS-LX-0561-R1")
#aliquots2add <- c("GLSS-SN-0009-R1", "GLSS-SN-0013-R2")
#samples_upload <- sample_upload %>% filter(sample_barcode %in% aliquots2add)

norlux_sample_barcode <- samples_upload$sample_barcode

sample_barcode <- norlux_sample_barcode

aliquot_analyte_type <- rep("R",length(sample_barcode))
aliquot_analysis_type <- rep("RNA",length(sample_barcode))
aliquot_portion <- rep(1,length(sample_barcode))
#aliquot_batch <- paste(substring(sample_barcode,1,8),"RNA",sep="")
aliquot_batch <- rep("GLSS-SN-RNA",length(sample_barcode))  # Using this as the batch because all samples were processed together

# Read in all previous UUIDs to ensure that they do not conflict
old_table <- read.delim("file_metadata/SNU/synapse_aliquots_release_20210823.txt",stringsAsFactors=FALSE)
old_uuids <- old_table$aliquot_uuid_short
db_table <- dbReadTable(con, Id(schema = "biospecimen", table="aliquots"))
db_uuids <- db_table$aliquot_uuid_short
old_uuids <- union(old_uuids, db_uuids)

#Add a short 6-character UUID to the aliquot barcodes that need it
#This snippit of code was pulled from the GLASS Github (GLASS/R/Manifest/life-history-barcode-generation.R)
aliquot_uuid <- stri_rand_strings((length(old_uuids) + length(sample_barcode) + 17), 6, "[A-Z0-9]") # Adding 17 because several uuids were removed post-hoc
aliquot_uuid_short <- aliquot_uuid[(length(aliquot_uuid)-(length(sample_barcode)-1)):length(aliquot_uuid)]

#Check to make sure there is no overlap between new UUIDs and current ones
ifelse(sum(old_uuids%in%aliquot_uuid_short) == 0, 
       message("UUIDs do not overlap"), 
       message("WARNING! Overlap"))

#Check to make sure each UUID is unique
ifelse(n_distinct(aliquot_uuid_short)==length(aliquot_uuid_short)[1], 
       message("UUIDs are unique."), 
       message("WARNING! Not unique"))

aliquot_barcode <- paste(sample_barcode, "-", str_pad(aliquot_portion, 2, pad = "0"), aliquot_analyte_type, "-", aliquot_analysis_type, "-", aliquot_uuid_short,sep="")

aliquots <- data.frame(aliquot_barcode, sample_barcode, aliquot_uuid_short, aliquot_analyte_type, aliquot_analysis_type, aliquot_portion, aliquot_batch)

#dbWriteTable(con, Id(schema="biospecimen",table = "aliquots"), aliquots, append = TRUE)

#-------------------------
# Create readgroups table

#conIn <- file("file_metadata/NorLux/rna_fastq_headers.txt","r")
conIn <- file("file_metadata/NorLux/rna_fastq_headers_pt3.txt","r")
fastq <- readLines(conIn)
close(conIn)

file_name <- fastq[seq(from = 1, to = length(fastq), by=2)]
file_path <- paste("/fastscratch/varnf/tmp_storage/SNU_NorLux/",file_name,sep="")
fastq_header <- fastq[seq(from = 2, to = length(fastq), by=2)]

file_info <- data.frame(file_path, file_name, fastq_header,stringsAsFactors = FALSE)

# Extract legacy sample ID to map fastq to GLASS barcode
legacy_sample_id <- sapply(strsplit(file_name,"_",),function(x)x[1])

# Extract readgroup_idtag info
flowcell_id <- sapply(strsplit(fastq_header,":"),function(x)x[3])
lane <- sapply(strsplit(fastq_header,":"),function(x)x[4])
readgroup_idtag <- paste(substring(flowcell_id,1, 5), ".", lane,sep="")

# Create idtag_legacy and platform columns (same for all samples)
readgroup_idtag_legacy <- rep(NA, nrow(file_info))
readgroup_platform <- rep("Illumina", nrow(file_info))

# Create readgroup_library and readgroup_platform_unit ({flowcell_barcode}.{lane}.{library_id})
readgroup_library <- sapply(strsplit(file_name, "_"),function(x)x[2])
readgroup_platform_unit <- paste(readgroup_idtag, readgroup_library,sep=".")

# Create readgroup_center (all samples done at JAX)
readgroup_center <- rep("JAX", nrow(file_info))
readgroup_timestamp <- rep(NA, nrow(file_info))

# Get aliquot_barcode and readgroup_sample_id (same thing for this project) by using the linker table
# Add the NorLux table
norlux_linker <- samples[,c("sampleId","sample_barcode")]
norlux_linker <- norlux_linker %>%
  inner_join(aliquots, by="sample_barcode") %>%
  select(sampleId, sample_barcode, aliquot_barcode)
colnames(norlux_linker) <- c("frozen_id", "glass_id", "aliquot_barcode")
frozen_glass <- norlux_linker

sub_readgroups <- data.frame(legacy_sample_id, readgroup_idtag, readgroup_idtag_legacy, readgroup_platform, readgroup_library, readgroup_platform_unit, readgroup_center, readgroup_timestamp)

# Manual addition for some late reads
###
#aliquot_barcode <- c("GLSS-SN-0009-R1-01R-RNA-T7P20X", "GLSS-SN-0009-R1-01R-RNA-T7P20X", "GLSS-SN-0013-R2-01R-RNA-A9XVF7", "GLSS-SN-0013-R2-01R-RNA-A9XVF7",
#                     "GLSS-LX-0357-R3-01R-RNA-X524U5", "GLSS-LX-0357-R3-01R-RNA-X524U5", "GLSS-LX-0561-R1-01R-RNA-JXGAE0", "GLSS-LX-0561-R1-01R-RNA-JXGAE0")
#sub_readgroups[,"aliquot_barcode"] <- aliquot_barcode
###
readgroups <- sub_readgroups %>% inner_join(frozen_glass, by = c("legacy_sample_id" = "frozen_id"))
readgroups <- readgroups[,c("aliquot_barcode", "readgroup_idtag", "readgroup_idtag_legacy", "readgroup_platform", "readgroup_platform_unit", "readgroup_library", "readgroup_center", "aliquot_barcode", "readgroup_timestamp")]
colnames(readgroups) <- c("aliquot_barcode", "readgroup_idtag", "readgroup_idtag_legacy", "readgroup_platform", "readgroup_platform_unit", "readgroup_library", "readgroup_center", "readgroup_sample_id", "readgroup_timestamp")

# Remove duplicates (paired-end reads so 2 entries per file)
readgroups <- readgroups[-which(duplicated(readgroups)),]
#dbWriteTable(con, Id(schema="biospecimen",table = "readgroups"), readgroups, append = TRUE)

# Create analysis.files table
file_sizes <- read.delim("file_metadata/NorLux/snu_norlux_size2.txt",sep="\t", header=FALSE)
md5sum <- read.delim("file_metadata/NorLux/snu_norlux_md5sum2.txt",sep=" ",header=FALSE)
file_info <- file_sizes %>%
  inner_join(md5sum, by = c("V2" = "V3"))
colnames(file_info) <- c("file_size", "file_path", "file_md5sum", "NA")
file_info <- file_info[,-which(colnames(file_info) == "NA")]
file_info$file_path <- paste("/fastscratch/varnf/tmp_storage/SNU_NorLux/",file_info$file_path,sep="")

file_info[,"file_name"] <- sapply(strsplit(file_info$file_path,"/"),function(x)x[length(x)])
file_info[,"file_format"] <- "FASTQ"
file_info <- file_info[,c("file_name","file_size", "file_md5sum", "file_format", "file_path")]

# Manual addition of late reads
#aliquot_barcode <- c("GLSS-SN-0009-R1-01R-RNA-T7P20X", "GLSS-SN-0009-R1-01R-RNA-T7P20X", "GLSS-SN-0013-R2-01R-RNA-A9XVF7", "GLSS-SN-0013-R2-01R-RNA-A9XVF7",
#                     "GLSS-LX-0357-R3-01R-RNA-X524U5", "GLSS-LX-0357-R3-01R-RNA-X524U5", "GLSS-LX-0561-R1-01R-RNA-JXGAE0", "GLSS-LX-0561-R1-01R-RNA-JXGAE0")
#file_info[,"aliquot_barcode"] <- aliquot_barcode

# Link aliquot barcode to files table
file_info[,"legacy_id"] <- sapply(strsplit(file_info$file_name,"_"), function(x)x[1])
file_info <- file_info %>%
  inner_join(frozen_glass, by=c("legacy_id" = "frozen_id"))

files <- file_info[,c("aliquot_barcode", "file_name", "file_size", "file_md5sum", "file_format", "file_path")]
#dbWriteTable(con, Id(schema="analysis",table = "files"), files, append = TRUE)


# Create files_readgroups
files_readgroups <- readgroups %>%
  inner_join(files, by= c("readgroup_sample_id" = "aliquot_barcode")) %>%
  select(file_name, readgroup_idtag, readgroup_sample_id)

#dbWriteTable(con, Id(schema="analysis",table = "files_readgroups"), files_readgroups, append = TRUE)

