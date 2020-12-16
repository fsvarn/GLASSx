###################################################
# Create stacked barplots (transcriptional classifier/simplicity score/CIBERSORTx) of each GLASS sample
# Updated: 2020.07.06
# Author: Frederick Varn
##################################################
library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(survival)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
dat <- dbReadTable(con, Id(schema="analysis", table="tumor_rna_clinical_comparison"))

# Test whether proneural to mesenchymal transition is associated with treatment in IDHwt

idhwt <- dat %>% filter(idh_codel_subtype == 'IDHwt' & !is.na(received_treatment))

idhwt[,"PMT"] <- idhwt[,"subtype_a"] == "Proneural" & idhwt[,"subtype_b"] == "Mesenchymal"
g1 <- idhwt %>% filter(PMT & received_rt == 1) %>% nrow()
g2 <- idhwt %>% filter(PMT & received_rt == 0) %>% nrow()
g3 <- idhwt %>% filter(!PMT & received_rt == 1) %>% nrow()
g4 <- idhwt %>% filter(!PMT & received_rt == 0) %>% nrow()

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)

g1 <- idhwt %>% filter(PMT & received_tmz == 1) %>% nrow()
g2 <- idhwt %>% filter(PMT & received_tmz == 0) %>% nrow()
g3 <- idhwt %>% filter(!PMT & received_tmz == 1) %>% nrow()
g4 <- idhwt %>% filter(!PMT & received_tmz == 0) %>% nrow()

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)


g1 <- idhwt %>% filter(PMT & received_bev == 1) %>% nrow()
g2 <- idhwt %>% filter(PMT & received_bev == 0) %>% nrow()
g3 <- idhwt %>% filter(!PMT & received_bev == 1) %>% nrow()
g4 <- idhwt %>% filter(!PMT & received_bev == 0) %>% nrow()

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)

# Compare PMT to samples that started proneural and didn't change
idhwt[,"PP"] <- idhwt[,"subtype_a"] == "Proneural" & idhwt[,"subtype_b"] == "Proneural"
g1 <- idhwt %>% filter(PMT & received_tmz == 1) %>% nrow()
g2 <- idhwt %>% filter(PMT & received_tmz == 0) %>% nrow()
g3 <- idhwt %>% filter(PP & received_tmz == 1) %>% nrow()
g4 <- idhwt %>% filter(PP & received_tmz == 0) %>% nrow()

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)


# Test whether proneural to mesenchymal transition is associated with treatment in IDHmut

idhmut <- dat %>% filter(idh_codel_subtype != 'IDHwt' & !is.na(received_tmz) & !is.na(received_rt))

idhmut[,"PMT"] <- idhmut[,"subtype_a"] == "Proneural" & idhmut[,"subtype_b"] == "Mesenchymal"
g1 <- idhmut %>% filter(PMT & received_rt == 1) %>% nrow()
g2 <- idhmut %>% filter(PMT & received_rt == 0) %>% nrow()
g3 <- idhmut %>% filter(!PMT & received_rt == 1) %>% nrow()
g4 <- idhmut %>% filter(!PMT & received_rt == 0) %>% nrow()

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)

g1 <- idhmut %>% filter(PMT & received_tmz == 1) %>% nrow()
g2 <- idhmut %>% filter(PMT & received_tmz == 0) %>% nrow()
g3 <- idhmut %>% filter(!PMT & received_tmz == 1) %>% nrow()
g4 <- idhmut %>% filter(!PMT & received_tmz == 0) %>% nrow()

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)

# Compare PMT to samples that started proneural and didn't change
idhmut[,"PP"] <- idhmut[,"subtype_a"] == "Proneural" & idhmut[,"subtype_b"] == "Proneural"
g1 <- idhmut %>% filter(PMT & received_tmz == 1) %>% nrow()
g2 <- idhmut %>% filter(PMT & received_tmz == 0) %>% nrow()
g3 <- idhmut %>% filter(PP & received_tmz == 1) %>% nrow()
g4 <- idhmut %>% filter(PP & received_tmz == 0) %>% nrow()

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)


# Patient prognosis
#######################################################

dat[,"PMT"] <- dat[,"subtype_a"] == "Proneural" & dat[,"subtype_b"] == "Mesenchymal"
dat[,"recur_status"] <- 1
dat[,"into_mes"] <- dat[,"subtype_a"] != "Mesenchymal" & dat[,"subtype_b"] == "Mesenchymal"
mycox <- coxph(Surv(surgical_interval, recur_status) ~ PMT + idh_codel_subtype, data=dat)
mycox <- summary(mycox)

mycox <- coxph(Surv(surgical_interval, recur_status) ~ into_mes + idh_codel_subtype, data=dat)
mycox <- summary(mycox)


