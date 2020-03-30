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
SELECT * FROM variants.hla_loh_change
"

dat <- dbGetQuery(con,q)

frame <- dat[,c("case_barcode","tumor_barcode_a","tumor_barcode_b","idh_codel_subtype")]
frame <- frame[-which(duplicated(frame)),]

hla_a_tab <- dat[grep("HLA-A",dat[,"hla_type1"]),]
hla_b_tab <- dat[grep("HLA-B",dat[,"hla_type1"]),]
hla_c_tab <- dat[grep("HLA-C",dat[,"hla_type1"]),]

hla_a <- hla_a_tab[,"hla_loh_change"]
names(hla_a) <- hla_a_tab[,"case_barcode"]
hla_a <- hla_a[frame[,"case_barcode"]]

hla_b <- hla_b_tab[,"hla_loh_change"]
names(hla_b) <- hla_b_tab[,"case_barcode"]
hla_b <- hla_b[frame[,"case_barcode"]]

hla_c <- hla_c_tab[,"hla_loh_change"]
names(hla_c) <- hla_c_tab[,"case_barcode"]
hla_c <- hla_c[frame[,"case_barcode"]]

mat <- cbind(frame,hla_a,hla_b,hla_c)
mat[,"any_gain"] <- apply(mat[,c("hla_a","hla_b","hla_c")],1,function(x)sum(x=="gain",na.rm=TRUE)>0)
mat[,"any_loss"] <- apply(mat[,c("hla_a","hla_b","hla_c")],1,function(x)sum(x=="loss",na.rm=TRUE)>0)
mat[,"any_stable"] <- apply(mat[,c("hla_a","hla_b","hla_c")],1,function(x)sum(x=="stable",na.rm=TRUE)>0)
mat[,"none"] <- apply(mat[,c("hla_a","hla_b","hla_c")],1,function(x)sum(x=="none",na.rm=TRUE)==(sum(!is.na(x))))

subtypes <- unique(mat[,"idh_codel_subtype"])

res <- matrix(0,nrow = 4, ncol = length(subtypes)+1)
rownames(res) <- c("prop_gain","prop_loss","prop_stable","prop_none")
colnames(res) <- c(subtypes,"total")
for(i in 1:length(subtypes))
{
	sub_mat <- mat[which(mat[,"idh_codel_subtype"]==subtypes[i]),]
	
	prop_gain <- sum(sub_mat[,"any_gain"])/nrow(sub_mat)
	prop_loss <- sum(sub_mat[,"any_loss"])/nrow(sub_mat)
	prop_stable <- sum(sub_mat[,"any_stable"])/nrow(sub_mat)
	prop_none <- sum(sub_mat[,"none"])/nrow(sub_mat)
	
	res[,i] <- c(prop_gain, prop_loss, prop_stable, prop_none)
}

prop_gain <- sum(mat[,"any_gain"])/nrow(mat)
prop_loss <- sum(mat[,"any_loss"])/nrow(mat)
prop_stable <- sum(mat[,"any_stable"])/nrow(mat)
prop_none <- sum(mat[,"none"])/nrow(mat)

res[,4] <- c(prop_gain, prop_loss, prop_stable, prop_none)
