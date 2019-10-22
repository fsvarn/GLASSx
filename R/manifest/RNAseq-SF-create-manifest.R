#######################################################
# Create manifest for SF RNAseq samples
# Date: 2019.10.22
# Author: Fred V.
#######################################################
#Local Directory for github repo
mybasedir = "/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/"
setwd(mybasedir)

#Files with information about fastq information and barcodes
case_file <- "data/metadata/sf_cases.txt"
sample_file <- "data/metadata/sf_samples.txt"
aliquot_file <- "data/metadata/sf_aliquots.txt"
readgroup_file <- "data/metadata/sf_readgroups.txt"
file_file <- "data/metadata/sf_files.txt"
files_readgroups_file <- "data/metadata/sf_files_readgroups.txt"

#######################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# We need to generate the following fields required by the SNV snakemake pipeline:
# aliquots, files, cases, samples, pairs, and readgroups.
### Cases ###
cases <- read.delim(case_file, stringsAsFactor=FALSE)
cases_master <- cases %>% mutate(case_barcode = case_barcode,
                                 case_project = case_project,
                                 case_source = case_source,
                                 case_sex = recode(case_sex,"M"="male","F"="female"),
                                 case_age_diagnosis_years = as.numeric(case_age_diagnosis_years),
                                 case_vital_status = case_vital_status,
                                 case_overall_survival_mo = as.integer(case_overall_survival_mo))

#Write to database
dbWriteTable(con, Id(schema="clinical", table="cases"), cases_master, append=TRUE)

### Samples ###
samples <- read.delim(sample_file, stringsAsFactor=FALSE)
samples_master <- samples %>% mutate(case_barcode = case_barcode,
                                     sample_barcode = sample_barcode,
                                     sample_type = sample_type)

#Write to database
dbWriteTable(con, Id(schema="biospecimen", table="samples"), samples_master, append=TRUE)

### Aliquots ###
aliquots <- read.delim(aliquot_file, stringsAsFactor=FALSE)
aliquots_master <- aliquots %>% mutate(aliquot_barcode = aliquot_barcode,
                                       aliquot_sample = aliquot_sample,
                                       aliquot_uuid_short = aliquot_uuid_short,
                                       aliquot_analyte_type = aliquot_analyte_type,
                                       aliquot_analysis_type = aliquot_analysis_type,
                                       aliquot_batch = aliquot_batch) %>%
  rename(sample_barcode = aliquot_sample)

#Write to database
dbWriteTable(con, Id(schema="biospecimen", table="aliquots"), aliquots_master, append=TRUE)

### Readgroups ###
readgroups <- read.delim(readgroup_file, stringsAsFactor=FALSE)
readgroups_master <- readgroups %>% mutate(aliquot_barcode = aliquot_barcode,
                                           readgroup_idtag = readgroup_idtag,
                                           readgroup_idtag_legacy = readgroup_idtag_legacy,
                                           readgroup_platform = readgroup_platform,
                                           readgroup_platform_unit = readgroup_platform_unit,
                                           readgroup_library = readgroup_library,
                                           readgroup_center = readgroup_center,
                                           readgroup_sample_id = readgroup_sample_id,
                                           readgroup_timestamp = readgroup_timestamp)

#Write to database
dbWriteTable(con, Id(schema="biospecimen", table="readgroups"), readgroups_master, append=TRUE)

### Files ###
files <- read.delim(file_file, stringsAsFactor=FALSE)
files_master <- files %>% mutate(aliquot_barcode = aliquot_barcode,
                                 file_name = file_name,
                                 file_size = file_size,
                                 file_md5sum = file_md5sum,
                                 file_format = file_format)

#Write to database
dbWriteTable(con, Id(schema="analysis", table="files"), files_master, append=TRUE)

### files_readgroups ###
files_readgroups <- read.delim(files_readgroups_file, stringsAsFactor=FALSE)
files_readgroups_master <- files_readgroups %>% mutate(file_name = file_name,
                                                       readgroup_idtag = readgroup_idtag,
                                                       readgroup_sample_id = readgroup_sample_id)

#Write to database
dbWriteTable(con, Id(schema="analysis", table="files_readgroups"), files_readgroups_master, append=TRUE)
