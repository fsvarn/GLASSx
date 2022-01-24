#!/usr/bin/env Rscript

library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

new_dat <- read.delim("/projects/verhaak-lab/GLASS-III/results/mutect2/consensusvcf/consensus.normalized.sorted.funcotated.GLSS-SN.tsv",stringsAsFactor=FALSE)
dbWriteTable(con, Id(schema="variants", table="anno"), new_dat, append=TRUE)

EOF