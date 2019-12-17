#!/usr/bin/env Rscript

#######################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please provide an input", call.=FALSE)}

transcript_tpm <- args[1]

cat("Reading in inputs...")
ttpm <- read.delim(transcript_tpm,stringsAsFactors=FALSE)

dbWriteTable(con, Id(schema="analysis", table="transcript_tpm"), ttpm, overwrite=TRUE)

