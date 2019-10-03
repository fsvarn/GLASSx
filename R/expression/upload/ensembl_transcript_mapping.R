#######################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

myinf1 <- "/projects/varnf/SofWar/kallisto/Ensembl_v96_homo_sapiens/homo_sapiens/transcripts_to_genes.txt"

data <- read.delim(myinf1,header=FALSE)
colnames(data) <- c("ensembl_transcript_id", "ensembl_gene_id", "gene_symbol")

dbWriteTable(con, Id(schema="ref", table="ensembl_transcript_mapping"), data, overwrite=TRUE, row.names=FALSE)

