### push titan seg into db

library(tidyverse)
library(DBI)
library(odbc)
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

segfiles <- list.files('results/cnv/titanfinal/seg', full.names = TRUE)

q <- "SELECT * FROM analysis.pairs"
pair <- dbGetQuery(con, q)

lapply(segfiles, function(f){
  message(f)
  dat <- read.delim(f, as.is=T, header=T, row.names = NULL)
  df <- dat %>%
    transmute(pair_barcode = Sample,
              chrom = Chromosome,
              pos = sprintf("[%s,%s]",Start_Position.bp.,End_Position.bp.),
              num_snp = Length.snp.,
              median_ratio = Median_Ratio,
              median_logr = Median_logR,
              titan_state = TITAN_state,
              titan_call = TITAN_call,
              copy_number = Copy_Number,
              major_cn = MajorCN,
              minor_cn = MinorCN,
              clonal_cluster = Clonal_Cluster,
              cellular_prevalence = Cellular_Prevalence,
              logr_copy_number = logR_Copy_Number,
              corrected_copy_number = Corrected_Copy_Number,
              corrected_call = Corrected_Call)
	
	df[which(df[,"chrom"]=="X"),"chrom"] <- 23
	df[which(df[,"chrom"]=="Y"),"chrom"] <- 24

    dbWriteTable(con, Id(schema="variants",table="titan_seg"), df, append=T)
    Sys.sleep(1)
})