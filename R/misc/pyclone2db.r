### push pyclone results to db

library(tidyverse)
library(DBI)
library(odbc)
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

myDir <- 'results/pyclone/run/'
files <- dir(myDir)
files <- paste(myDir, files, sep = '')

# Select H2, PD, or SN cohorts:
# files <- files[grep("GLSS-CU-P|GLSS-HF", files)]
# files <- files[-grep("GLSS-HF-2", files)]
# files <- files[-grep("GLSS-HF-3", files)]
#files <- files[grep("GLSS-SN-", files)]

clust <- paste(files, '/tables/cluster.tsv', sep='')
loci <- paste(files, '/tables/loci.tsv', sep='')

clust <- clust[which(file.exists(clust))]
loci <- loci[which(file.exists(loci))]

lapply(clust, function(f){
  message(f)
  dat <- read.delim(f, as.is=T, header=T, row.names = NULL)
  df <- dat %>%
    transmute(aliquot_barcode = sample_id,
              cluster_id = cluster_id,
              size = size,
              mean = mean,
              std = std)

    #dbWriteTable(con, Id(schema="variants",table="pyclone_cluster"), df, append=T)
    Sys.sleep(1)
})

lapply(loci, function(f){
  message(f)
  dat <- read.delim(f, as.is=T, header=T, row.names = NULL)
  df <- dat %>%
    transmute(variant_id = mutation_id,
    		  aliquot_barcode = sample_id,
              cluster_id = cluster_id,
              cellular_prevalence = cellular_prevalence,
              cellular_prevalence_sd = cellular_prevalence_std,
              variant_allele_frequency = variant_allele_frequency)

    #dbWriteTable(con, Id(schema="variants",table="pyclone_loci"), df, append=T)
    Sys.sleep(1)
})