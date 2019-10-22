library(ggplot2)
library(odbc)
library(DBI)
library(rjson)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

myDir1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/fastp"
myinf1 <- paste(myDir1,dir(myDir1),sep="/")
myinf1 <- sapply(myinf1,function(x)paste(x,dir(x)[2],sep="/"),USE.NAMES=FALSE)

#Load clinical data
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT al.aliquot_barcode, sa.sample_type, su.*, aliquot_batch
FROM biospecimen.aliquots al 
LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode
LEFT JOIN biospecimen.samples sa ON al.sample_barcode = sa.sample_barcode
WHERE aliquot_analyte_type = 'R'"

clin <- dbGetQuery(con,q)

mybatch <- unique(clin[,"aliquot_batch"])

myrate <- rep(0,length(myinf1))
for(i in 1:length(myinf1))
{
	result <- fromJSON(file = myinf1[i])
	myrate[i] <- result$duplication$rate
	
}

names(myrate) <- dir(myDir1)

aliquot_batch <- gsub("-RNA","",clin[match(names(myrate),clin[,"aliquot_barcode"]),"aliquot_batch"])
plot_res <- data.frame(myrate,aliquot_batch)


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/dup_by_batch.pdf",width=7,height=4)
p1 <- ggplot(plot_res,aes(y = myrate, x = aliquot_batch)) +
	geom_boxplot() +
	labs(x = "Aliquot batch", y = "PCR duplication rate") +
	theme_bw() +
	theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
	axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position="none") +
	coord_cartesian(ylim=c(0,1))
p1
dev.off()