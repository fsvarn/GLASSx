library(odbc)
library(DBI)
library(tidyverse)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in aggregate pathologist annotations
myinf1 <- "data/pathology_annot/Histo_review_aggregate_v2.txt"

dat <- read.delim(myinf1)

dat_long <- dat %>%
  select(-Notes) %>%
  pivot_longer(-c(CASES,Pathologist),names_to ="feature", values_to = "metric")

feat_class <- rep("Primary",nrow(dat_long))
feat_class[grepl("NecTreatEff|Fibrosis|DepopTum", dat_long$feature)] <- "Recurrent"
dat_long$feat_class <- feat_class

metric_class <- rep("Percent",nrow(dat_long))
metric_class[grepl("_score", dat_long$feature)] <- "Score"
dat_long$metric_class <- metric_class

dat_long <- dat_long %>%
  mutate(feature = fct_relevel(feature, unique(dat_long$feature)))

# Transcriptomic analysis

# Read in IvyGAP table
ivy <- dbReadTable(con, Id(schema="analysis",table="cibersortx_ivygap"))

ivy <- ivy %>%
  filter(grepl("GLSS-HF|GLSS-LU", aliquot_barcode)) 
ivy$aliquot_barcode <- sapply(strsplit(ivy$aliquot_barcode,"-"),function(x)paste(x[1:4],collapse="-"))
ivy$fraction <- ivy$fraction * 100
ivy <- ivy %>%
       pivot_wider(names_from = "cell_state", values_from = "fraction")

pathologist <- unique(dat$Pathologist)
init_le <- init_ct <- init_nec <- init_mvp <- rec_le <- rec_ct <- rec_nec <- rec_mvp <- rep(0, length(pathologist))
for(i in 1:length(pathologist))
{
    sub_dat <- dat_long %>% 
               filter(Pathologist == pathologist[i]) %>%
               select(CASES, feature, metric) %>%
               pivot_wider(names_from = feature, values_from = metric) %>%
               inner_join(ivy, by = c("CASES" = "aliquot_barcode"))
    init_cor <- sub_dat %>% 
               filter(grepl("-TP",CASES)) %>%
               column_to_rownames("CASES") %>%
               cor(., method = "p")
    rec_cor <- sub_dat %>% 
               filter(grepl("-R1|-R2",CASES)) %>%
               column_to_rownames("CASES") %>%
               cor(., method = "p") 
    init_le[i] <- init_cor["LE_percent","LE"]
    init_ct[i] <- init_cor["CT_percent","CT"]
    init_nec[i] <- init_cor["Necrosis_percent","CTpan"]
    init_mvp[i] <- init_cor["Necrosis_percent","CTmvp"]
    rec_le[i] <- rec_cor["LE_percent","LE"]
    rec_ct[i] <- rec_cor["CT_percent","CT"]
    rec_nec[i] <- rec_cor["Necrosis_percent","CTpan"]
    rec_mvp[i] <- rec_cor["Necrosis_percent","CTmvp"]
    
}  
init_res <- data.frame(pathologist, init_le, init_ct, init_nec, init_mvp)
rec_res <- data.frame(pathologist, rec_le, rec_ct, rec_nec, rec_mvp)

init_res
rec_res