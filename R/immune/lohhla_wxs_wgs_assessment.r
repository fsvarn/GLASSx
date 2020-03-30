library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

dat <- dbReadTable(con,Id(schema="variants",table="lohhla"))

pair_id <- substr(dat[,1],1,25)
analysis_id <- substr(dat[,1],27,29)

uni_pairs <- unique(pair_id)

both <- c()
for(i in 1:length(uni_pairs))
{
	mypair <- pair_id[which(pair_id==uni_pairs[i])]
	myanalysis <- analysis_id[which(pair_id==uni_pairs[i])]
	
	if(sum("WXS" %in% myanalysis) > 0 & sum("WGS" %in% myanalysis)){
		both <- c(both, unique(mypair))}
}

sub_dat <- dat[which(substr(dat[,1],1,25) %in%  both),]

wxs_dat <- sub_dat[grep("-WXS",sub_dat[,"pair_barcode"]),]
wgs_dat <- sub_dat[grep("-WGS",sub_dat[,"pair_barcode"]),]

#wxs_dat <- wxs_dat[-which(duplicated(wxs_dat)),]
#wgs_dat <- wgs_dat[-which(duplicated(wgs_dat)),]

mylist <- list()
ct <- 1
for(i in 1:nrow(wxs_dat))
{
	myrow <- wxs_dat[i,]
	matchrow <- wgs_dat[which(substr(wgs_dat[,1],1,25) == substr(myrow[,1],1,25) &
		wgs_dat[,"coverage_filter"] == myrow[,"coverage_filter"] &
		wgs_dat[,"hla_type1"] == myrow[,"hla_type1"]),]
	
	if(nrow(matchrow)==0){
		next}
	
	subrow1 <- myrow[,c("pair_barcode","coverage_filter","hla_type1","hla_type2","pval","loss_allele","num_mismatch_sites_cov","prop_supportive_sites")]
	subrow2 <- matchrow[,c("pval","loss_allele","num_mismatch_sites_cov","prop_supportive_sites")]
	colnames(subrow2)	<- paste(colnames(subrow2),"_wgs",sep="")
	newrow <- cbind(subrow1,subrow2)
	
	mylist[[ct]] <- newrow
	ct <- ct + 1
}

res <- do.call(rbind,mylist)

res <- res[,c(1,2,3,4,5,9,6,10,7,11,8,12)]

coverage_filter <- c(5,10,20,30)
na_wxs <- na_wgs <- tot <- wxs_pval_conc <- wgs_pval_conc <- wxs_wgs_disc <- sig_wxs <- sig_wgs <- rep(0,length(coverage_filter))
for(i in 1:length(coverage_filter))
{
	sub_res <- res[which(res[,"coverage_filter"]==coverage_filter[i]),]
	
	tot[i] <- nrow(sub_res)
	
	na_wxs[i] <- sum(is.na(sub_res[,"pval"]))
	na_wgs[i] <- sum(is.na(sub_res[,"pval_wgs"]))
	
	sig_wxs[i] <- sum(sub_res[,"pval"] < 0.05)
	sig_wgs[i] <- sum(sub_res[,"pval_wgs"] < 0.05,na.rm=TRUE)
	
	wxs_pval_conc[i] <- sum(sub_res[,"pval"]<0.05 & sub_res[,"pval_wgs"] <0.05 & sub_res[,"loss_allele"] == sub_res[,"loss_allele_wgs"],na.rm=TRUE)/sum(sub_res[,"pval"] < 0.05)
	wgs_pval_conc[i] <- sum(sub_res[,"pval"]<0.05 & sub_res[,"pval_wgs"] <0.05 & sub_res[,"loss_allele"] == sub_res[,"loss_allele_wgs"],na.rm=TRUE)/sum(sub_res[,"pval_wgs"] < 0.05,na.rm=TRUE)

	wxs_wgs_disc[i] <- sum(sub_res[,"pval"]<0.05 & sub_res[,"pval_wgs"] <0.05 & sub_res[,"loss_allele"] != sub_res[,"loss_allele_wgs"],na.rm=TRUE)
}

data.frame(coverage_filter, tot, na_wxs, na_wgs, sig_wxs, sig_wgs, wxs_pval_conc, wgs_pval_conc,wxs_wgs_disc)



