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

cases <- read.delim("file_metadata/SNU/snu_cases_20210820.txt",stringsAsFactors=FALSE)
samples <- read.delim("file_metadata/SNU/snu_samples_20210823.txt",stringsAsFactors=FALSE)

#dbWriteTable(con, Id(schema="clinical",table = "cases"), cases, append = TRUE)
#dbWriteTable(con, Id(schema="biospecimen",table = "samples"), samples, append = TRUE)

#-------------------------

# Create aliquots table

sample_barcode <- samples$sample_barcode
aliquot_analyte_type <- rep("D",length(sample_barcode))
aliquot_analysis_type <- rep("WGS",length(sample_barcode))
aliquot_portion <- rep(1,length(sample_barcode))
aliquot_batch <- rep("GLSS-SN-WGS",length(sample_barcode))

# Read in all previous UUIDs to ensure that they do not conflict
old_table <- read.delim("file_metadata/SNU/synapse_aliquots_release_20210823.txt",stringsAsFactors=FALSE)
old_uuids <- old_table$aliquot_uuid_short

#Add a short 6-character UUID to the aliquot barcodes that need it
#This snippit of code was pulled from the GLASS Github (GLASS/R/Manifest/life-history-barcode-generation.R)
aliquot_uuid <- stri_rand_strings((length(old_uuids) + length(sample_barcode)), 6, "[A-Z0-9]")
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

conIn <- file("file_metadata/SNU/snu_fastq_header.txt","r")
fastq <- readLines(conIn)
close(conIn)

file_path <- fastq[seq(from = 1, to = length(fastq), by=2)]
fastq_header <- fastq[seq(from = 2, to = length(fastq), by=2)]
file_name <- sapply(strsplit(file_path,"/"),function(x)x[2])

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
# Read in linker table
linker <- read.delim("/Users/varnf/Documents/Projects/IMC_Glioma/SNU/sample_storage/snu_cohort_mapping_table.txt",stringsAsFactors=FALSE)
frozen_id <- c(linker$Frozen.banking.number..initial., linker$Frozen.banking.number..recurrent.,linker$Frozen.banking.number..normal)
glass_id <- c(linker$sample_barcode..initial., linker$sample_barcode..recurrent.,linker$sample_barcode..normal.)
frozen_glass <- data.frame(frozen_id, glass_id)
frozen_glass <- frozen_glass %>%
inner_join(aliquots %>% select(aliquot_barcode, sample_barcode), by = c("glass_id" = "sample_barcode")) %>%
mutate(frozen_id = recode(frozen_id, "4655/4656" = "4655", "4111-_/4111-_/4111-_" = "4111", "4045/4046/4047" = "4045", "4792/4793/4794" = "4792"))

sub_readgroups <- data.frame(legacy_sample_id, readgroup_idtag, readgroup_idtag_legacy, readgroup_platform, readgroup_library, readgroup_platform_unit, readgroup_center, readgroup_timestamp)

readgroups <- sub_readgroups %>% inner_join(frozen_glass, by = c("legacy_sample_id" = "frozen_id"))
readgroups <- readgroups[,c("aliquot_barcode", "readgroup_idtag", "readgroup_idtag_legacy", "readgroup_platform", "readgroup_platform_unit", "readgroup_library", "readgroup_center", "aliquot_barcode", "readgroup_timestamp")]
colnames(readgroups) <- c("aliquot_barcode", "readgroup_idtag", "readgroup_idtag_legacy", "readgroup_platform", "readgroup_platform_unit", "readgroup_library", "readgroup_center", "readgroup_sample_id", "readgroup_timestamp")

# Remove duplicates (paired-end reads so 2 entries per file)
readgroups <- readgroups[-which(duplicated(readgroups)),]
#dbWriteTable(con, Id(schema="biospecimen",table = "readgroups"), readgroups, append = TRUE)

# Create analysis.files table
file_sizes <- read.delim("file_metadata/SNU/SNU_file_size.txt",sep="\t", header=FALSE)
md5sum <- read.delim("file_metadata/SNU/snu_md5sum.txt",sep=" ",header=FALSE)
file_info <- file_sizes %>%
         inner_join(md5sum, by = c("V2" = "V3"))
colnames(file_info) <- c("file_size", "file_path", "file_md5sum", "NA")
file_info <- file_info[,-which(colnames(file_info) == "NA")]
file_info$file_path <- paste("/projects/verhaak-lab/GLASS-SNU/wgs/",file_info$file_path,sep="")

file_info[,"file_name"] <- sapply(strsplit(file_info$file_path,"/"),function(x)x[length(x)])
file_info[,"file_format"] <- "FASTQ"
file_info <- file_info[,c("file_name","file_size", "file_md5sum", "file_format", "file_path")]

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

# Manually remove the incorrect pairings
ind1 <- which(grepl("4804_GT21-14705_GATTCTGC-CTCTCGTC_S114_L004_R", files_readgroups$file_name) & files_readgroups$readgroup_idtag =="HGWCK.1")
ind2 <- which(grepl("4804_GT21-14705_GATTCTGC-CTCTCGTC_S2_L001_R", files_readgroups$file_name) & files_readgroups$readgroup_idtag =="HKJWJ.4")
ind3 <- which(grepl("4906_GT21-14719_GCGCTCTA-GTCGGAGC_S118_L004_R", files_readgroups$file_name) & files_readgroups$readgroup_idtag =="HGWCK.4")
ind4 <- which(grepl("4906_GT21-14719_GCGCTCTA-GTCGGAGC_S26_L004_R", files_readgroups$file_name) & files_readgroups$readgroup_idtag =="HKJWJ.4")

files_readgroups <- files_readgroups[-c(ind1,ind2,ind3,ind4),]

#dbWriteTable(con, Id(schema="analysis",table = "files_readgroups"), files_readgroups, append = TRUE)


# Create analysis.files table for aligned bams
file_sizes <- read.delim("file_metadata/SNU/GLSS-SN-bam_size.txt",sep="\t", header=FALSE)
md5sum <- read.delim("file_metadata/SNU/GLSS-SN-bam_md5.txt",sep=" ",header=FALSE)
file_info <- file_sizes %>%
  inner_join(md5sum, by = c("V2" = "V3"))
colnames(file_info) <- c("file_size", "file_path", "file_md5sum", "NA")
file_info <- file_info[,-which(colnames(file_info) == "NA")]
file_info$file_path <- paste("/projects/verhaak-lab/GLASS-III/results/align/bqsr/",file_info$file_path,sep="")

file_info[,"file_name"] <- sapply(strsplit(file_info$file_path,"/"),function(x)x[length(x)])
file_info[,"file_format"] <- "aligned BAM"
file_info <- file_info[,c("file_name","file_size", "file_md5sum", "file_format", "file_path")]
file_info[,"aliquot_barcode"] <- substring(file_info[,"file_name"], 1, 30)

files <- file_info[,c("aliquot_barcode", "file_name", "file_size", "file_md5sum", "file_format", "file_path")]

#dbWriteTable(con, Id(schema="analysis",table = "files"), files, append = TRUE)


# Upload clinical.surgeries table
surgeries <- read.delim("file_metadata/SNU/snu_surgeries_20210913_upload_FINAL.txt",stringsAsFactors=FALSE)

surgeries$idh_codel_subtype <- paste(surgeries$idh_status, "-", surgeries$codel_status, sep="")
surgeries[grep("NA-", surgeries$idh_codel_subtype),"idh_codel_subtype"] <- NA

dbWriteTable(con, Id(schema="clinical",table = "surgeries"), surgeries, append = TRUE)

# Upload DNA blocklist table
blocklist <- read.delim("file_metadata/SNU/snu-dna-blocklist-20210913.txt",stringsAsFactors=FALSE)

dbWriteTable(con, Id(schema="analysis",table = "blocklist"), blocklist, append = TRUE)


# Upload RNA blocklist table
rna_blocklist <- read.delim("file_metadata/Norlux/snu_norlux-rna-blocklist-20210913.txt",stringsAsFactors=FALSE)

#dbWriteTable(con, Id(schema="analysis",table = "rna_blocklist"), rna_blocklist, append = TRUE)

