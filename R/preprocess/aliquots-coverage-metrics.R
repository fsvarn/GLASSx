#######################################################
# Enumerate cumulative coverage per aliquot for WGS/WXS
# Date: 2018.11.06 
# Author: Kevin J., modified by Fred V on 10/14/2020 for GLASSv3 
#######################################################

# Directory for GLASS analysis.
datadir  = 'results/align/wgsmetrics/'
pattern   = '.WgsMetrics.txt$'

#######################################################

# Necessary packages:
library(parallel)
library(tidyverse)
library(data.table)
library(DBI)

#######################################################
# Establish connection with the database.
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

## Read in an example "*.WgsMetrics.txt" file to test the calling.
files = list.files(datadir, full.names = T, pattern = pattern, recursive=T)

# If it is desirable to include the sample names.
samples = data.frame(sample_id=gsub(".WgsMetrics.txt", "", basename(files)), library_type = substring(basename(files), 21, 23))

# The first 10 rows of each file represent a header of additional information.
cov_dat = mclapply(files, function(f){
  dat = tryCatch(read.delim(f,as.is=T, header=T, row.names = NULL, skip = 10), error=function(e) e)
  if(inherits(dat,'error')) {
    message(f, '\n', dat, '\n')
    return()
  }
  # Truncate the file name to just the sample_id.
  dat = dat %>%
    mutate(sample_id = gsub(".WgsMetrics.txt", "", basename(f))) # %>%  
#    filter(coverage!="0") # Filter out those bases with `0` coverage.
  
  return(dat)
  
}, mc.cores=20)

## Combine all the samples from the GLASS cohort.
glass_cov = data.table::rbindlist(cov_dat)

# Cumulatively add the number of bases at each level:
glass_samples_cumulative_cov = glass_cov %>% 
  group_by(sample_id) %>% 
  mutate(cumulative_coverage = rev(cumsum(rev(high_quality_coverage_count)))) %>% 
  # Make sure colnames are formatting right.
  select(aliquot_barcode = sample_id, coverage, high_quality_coverage_count, cumulative_coverage) 
  

# Total number should be 1166 (2019.03.08).
n_distinct(glass_samples_cumulative_cov$aliquot_barcode)

# Pull out H2 (2020.10.14)
hf_cumulative_cov <- glass_samples_cumulative_cov %>%
					 filter(grepl("GLSS-HF-", aliquot_barcode)) %>%
					 filter(!grepl("GLSS-HF-2", aliquot_barcode)) %>%
					 filter(!grepl("GLSS-HF-3", aliquot_barcode)) 

# Write output as one table or a table for each file:
# write.table(glass_samples_cumulative_cov, file = "/Users/johnsk/Documents/Life-History/GLASS-WG/data/ref/glass-cumulative-coverage.txt", sep="\t", row.names = F, col.names = T, quote = F)

# Write to cumulative coverage files to database.
dbWriteTable(con, Id(schema="analysis",table="coverage"), hf_cumulative_cov, append=T)
