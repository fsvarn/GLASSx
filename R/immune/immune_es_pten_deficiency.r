#######################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
dat <- dbGetQuery(con,read_file("/projects/varnf/GLASS-III/GLASS-III/sql/expression/pten_deficiency_immune.sql"))

cells <- unique(dat[,"signature_name"])
pval <- eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub.dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	g1 <- sub.dat[which(sub.dat[,"mutect2_call_a"]=="1"),"enrichment_score_a"]
	g2 <- sub.dat[which(sub.dat[,"mutect2_call_a"]=="0" | is.na(sub.dat[,"mutect2_call_a"])),"enrichment_score_a"]
	pval[i] <- wilcox.test(g1,g2)$p.value
	eff[i] <- median(g1)-median(g2)
}
res <- data.frame(cells,eff,pval)
res <- res[order(res[,"pval"]),]

#                    cells          eff         pval
# 8         Macrophages.M2  0.057480941 0.0006949576
# 5              Dendritic  0.067316641 0.0028898990
# 6            Macrophages  0.041817352 0.0045841298
# 10                 T.reg  0.064158398 0.0071330187
# 2             CD4.mature  0.058826797 0.0328509752
# 3           CD8.effector  0.031469847 0.0373099250
# 4  CD8.effector.NK.cells  0.042623444 0.0412387370
# 7         Macrophages.M1  0.036963970 0.0696268809
# 1                B.cells -0.004506207 0.6353573083
# 9               NK.cells  0.017734294 0.7807553873


#Examine acquired PTEN mutation effecst on immune infiltrate

pair.dat <- dat[which(dat[,"mutect2_call_a"]=="0" & dat[,"mutect2_call_b"]=="1"),]
pval <- eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub.dat <- pair.dat[which(pair.dat[,"signature_name"]==cells[i]),]
	g1 <- sub.dat[,"enrichment_score_a"]
	g2 <- sub.dat[,"enrichment_score_b"]
	pval[i] <- wilcox.test(g1,g2,paired=TRUE)$p.value
	eff[i] <- median(g1)-median(g2)
}
res2 <- data.frame(cells,eff,pval)
res2 <- res2[order(res2[,"pval"]),]

#                    cells          eff      pval
# 5                B.cells -0.048670765 0.1075309
# 2             CD4.mature  0.043409549 0.1694820
# 8  CD8.effector.NK.cells  0.008897839 0.4557790
# 4               NK.cells -0.005831893 0.5559944
# 3           CD8.effector  0.003740786 0.6098451
# 1         Macrophages.M1  0.104172192 0.7238779
# 9            Macrophages  0.048745798 0.7238779
# 10        Macrophages.M2  0.041798752 0.7238779
# 6                  T.reg -0.015036524 0.8444011
# 7              Dendritic  0.032512119 0.8444011
