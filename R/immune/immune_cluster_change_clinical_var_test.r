library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT cc.*, ic.change
FROM analysis.tumor_clinical_comparison cc
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = cc.tumor_barcode_a AND gs.tumor_barcode_b = cc.tumor_barcode_b
JOIN analysis.rna_dna_pairs rd ON rd.dna_barcode_a = cc.tumor_barcode_a AND rd.dna_barcode_b = cc.tumor_barcode_b
JOIN analysis.immune_cluster_change ic ON ic.case_barcode = cc.case_barcode
"
dat <- dbGetQuery(con,q)

mytreatmentvar <- c("received_tmz","received_rt","received_alk","received_treatment","hypermutator_status","alk_assoc_hypermutator_status")

#Treatment check/hypermutator check
#---------------------------------------------
inc_res <- dec_res <- matrix(0,nrow=length(mytreatmentvar),ncol=3)
colnames(inc_res) <- colnames(dec_res) <- c("change_fraction","no_change_fraction","p.val")
rownames(inc_res) <- rownames(dec_res) <- mytreatmentvar
for(i in 1:length(mytreatmentvar))
{
	#Immune increase test
	g1 <- nrow(dat[which(dat[,"change"]=="Increase" & dat[,mytreatmentvar[i]]==1),])
	g2 <- nrow(dat[which(dat[,"change"]!="Increase" & dat[,mytreatmentvar[i]]==1),])
	g3 <- nrow(dat[which(dat[,"change"]=="Increase" & dat[,mytreatmentvar[i]]==0),])
	g4 <- nrow(dat[which(dat[,"change"]!="Increase" & dat[,mytreatmentvar[i]]==0),])
	
	inc_res[mytreatmentvar[i],"change_fraction"] <- g1/sum(g1,g3)
	inc_res[mytreatmentvar[i],"no_change_fraction"] <- g2/sum(g2,g4)
	inc_res[mytreatmentvar[i],"p.val"] <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value

	#Immune decrease test
	g1 <- nrow(dat[which(dat[,"change"]=="Decrease" & dat[,mytreatmentvar[i]]==1),])
	g2 <- nrow(dat[which(dat[,"change"]!="Decrease" & dat[,mytreatmentvar[i]]==1),])
	g3 <- nrow(dat[which(dat[,"change"]=="Decrease" & dat[,mytreatmentvar[i]]==0),])
	g4 <- nrow(dat[which(dat[,"change"]!="Decrease" & dat[,mytreatmentvar[i]]==0),])
	
	dec_res[mytreatmentvar[i],"change_fraction"] <- g1/sum(g1,g3)
	dec_res[mytreatmentvar[i],"no_change_fraction"] <- g2/sum(g2,g4)
	dec_res[mytreatmentvar[i],"p.val"] <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value
}

#Recurrence location
#---------------------------------------------

#Immune increase test
g1 <- nrow(dat[which(dat[,"change"]=="Increase" & dat[,"recurrence_location"]=="Local"),])
g2 <- nrow(dat[which(dat[,"change"]!="Increase" & dat[,"recurrence_location"]=="Local"),])
g3 <- nrow(dat[which(dat[,"change"]=="Increase" & dat[,"recurrence_location"]=="Distal"),])
g4 <- nrow(dat[which(dat[,"change"]!="Increase" & dat[,"recurrence_location"]=="Distal"),])

recloc_inc_fraction <- g1/sum(g1,g3)
recloc_no_inc_fraction <- g2/sum(g2,g4)
recloc_inc_p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value

#Immune decrease test
g1 <- nrow(dat[which(dat[,"change"]=="Decrease" & dat[,"recurrence_location"]=="Local"),])
g2 <- nrow(dat[which(dat[,"change"]!="Decrease" & dat[,"recurrence_location"]=="Local"),])
g3 <- nrow(dat[which(dat[,"change"]=="Decrease" & dat[,"recurrence_location"]=="Distal"),])
g4 <- nrow(dat[which(dat[,"change"]!="Decrease" & dat[,"recurrence_location"]=="Distal"),])

recloc_dec_fraction <- g1/sum(g1,g3)
recloc_no_dec_fraction <- g2/sum(g2,g4)
recloc_dec_p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value

#Grade change
#---------------------------------------------

#Immune increase test
g1 <- nrow(dat[which(dat[,"change"]=="Increase" & dat[,"grade_change"]=="Grade up"),])
g2 <- nrow(dat[which(dat[,"change"]!="Increase" & dat[,"grade_change"]=="Grade up"),])
g3 <- nrow(dat[which(dat[,"change"]=="Increase" & dat[,"grade_change"]=="Grade stable"),])
g4 <- nrow(dat[which(dat[,"change"]!="Increase" & dat[,"grade_change"]=="Grade stable"),])

grade_inc_fraction <- g1/sum(g1,g3)
grade_no_inc_fraction <- g2/sum(g2,g4)
grade_inc_p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value

#Immune decrease test
g1 <- nrow(dat[which(dat[,"change"]=="Decrease" & dat[,"grade_change"]=="Grade up"),])
g2 <- nrow(dat[which(dat[,"change"]!="Decrease" & dat[,"grade_change"]=="Grade up"),])
g3 <- nrow(dat[which(dat[,"change"]=="Decrease" & dat[,"grade_change"]=="Grade stable"),])
g4 <- nrow(dat[which(dat[,"change"]!="Decrease" & dat[,"grade_change"]=="Grade stable"),])

grade_dec_fraction <- g1/sum(g1,g3)
grade_no_dec_fraction <- g2/sum(g2,g4)
grade_dec_p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value


#Surgical interval
g1 <- dat[which(dat[,"change"]=="Increase"),"surgical_interval"]
g2 <- dat[which(dat[,"change"]!="Increase"),"surgical_interval"]

g3 <- dat[which(dat[,"change"]=="Decrease"),"surgical_interval"]
g4 <- dat[which(dat[,"change"]!="Increase"),"surgical_interval"]

wilcox.test(g1,g2)
wilcox.test(g3,g4)

wilcox.test(g1,g3)

