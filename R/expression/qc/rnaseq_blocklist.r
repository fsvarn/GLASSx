library(ggplot2)
library(odbc)
library(DBI)
library(rjson)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT kq.*
FROM analysis.kallisto_qc kq"
ka <- dbGetQuery(con,q)

unique_thr <- mean(ka[,"p_unique"]) - 2*sd(ka[,"p_unique"])

complexity_exclusion <- ifelse(ka[,"p_unique"] < mean(ka[,"p_unique"])-2*sd(ka[,"p_unique"]),"block","allow")
complexity_exclusion_reason <- ifelse(ka[,"p_unique"] < mean(ka[,"p_unique"])-2*sd(ka[,"p_unique"]),"low_complexity",NA)

res <- data.frame(ka[,"aliquot_barcode"],complexity_exclusion,complexity_exclusion_reason)
colnames(res) <- c("aliquot_barcode","complexity_exclusion","complexity_exclusion_reason")

dbWriteTable(con, Id(schema="analysis", table="rna_blocklist"), res, overwrite=TRUE)
