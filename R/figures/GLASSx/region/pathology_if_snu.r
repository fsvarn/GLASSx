library(odbc)
library(DBI)
library(tidyverse)
library(corrr)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in aggregate pathologist annotations
myinf1 <- "data/pathology_annot/neuropath_template_final_aggregate.txt"

dat <- read.delim(myinf1)


sn_dat <- dat %>%
filter(folder=="SN") %>%
group_by(sample_barcode) %>%
select(-c(timepoint, folder, annotations, pathologist)) %>%
summarise_if(is_numeric, mean) %>%
mutate(sample_barcode = recode(sample_barcode, "Pt10_primary" = "S13 18774 4", "Pt10_recurrence" = "S14 14567 2",
                                              "Pt16_primary" = "S17 40599 2", "Pt16_recurrence" = "S18 14642 2",
                                              "Pt17_primary" = "S16 24169 1", "Pt17_recurrence" = "S19 71211 6",
                                              "Pt6_primary" = "S12 47771 4", "Pt6_recurrence" = "S13 69431 1",
                                              "Pt8_primary" = "S12 63073", "Pt8_recurrence" = "S13 38572 12",
))

sn_if <- read.delim("/Users/varnf/Documents/Projects/IMC_Glioma/SNU/histocytometry/Glioma composition.txt")

sn_comp <- sn_dat %>%
           inner_join(sn_if, "sample_barcode")

# No IDH mutant
pri_comp <- sn_comp %>%
            filter(sample_barcode %in% c("S13 18774 4","S17 40599 2","S12 47771 4","S12 63073"))
rec_comp <- sn_comp %>%
  filter((sample_barcode %in% c("S14 14567 2","S18 14642 2","S13 69431 1","S13 38572 12")))

#pri_comp <- sn_comp %>%
#            filter(sample_barcode %in% c("S13 18774 4","S17 40599 2","S12 47771 4","S12 63073","S16 24169 1"))
#rec_comp <- sn_comp %>%
#  filter((sample_barcode %in% c("S14 14567 2","S18 14642 2","S13 69431 1","S13 38572 12","S19 71211 6")))

cor(pri_comp[,4:11], pri_comp[,12:17])
cor(rec_comp[,4:11], rec_comp[,12:17])

idhwt_comp <- sn_comp %>%
  filter(!(sample_barcode %in% c("S16 24169 1", "S19 71211 6")))
cor(idhwt_comp[,4:11], idhwt_comp[,12:17])

pri_comp$ct_percent <- pri_comp$cellular_tumor/(pri_comp$cellular_tumor+pri_comp$leading_edge)
pri_comp$le_percent <- pri_comp$leading_edge/(pri_comp$cellular_tumor+pri_comp$leading_edge)
cor(pri_comp[,c("le_percent","ct_percent")], pri_comp[,12:17],method="p")
