#######################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

transcript_tpm <- "results/kallisto/kallisto/final/transcript_tpms_all_samples.tsv"

cat("Reading in inputs...")
ttpm <- read.delim(transcript_tpm,stringsAsFactors=FALSE)

dbWriteTable(con, Id(schema="analysis", table="transcript_tpm"), ttpm, overwrite=TRUE)

