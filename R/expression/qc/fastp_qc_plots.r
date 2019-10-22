#Code to create a qc table from the fastp output. Table contains the following information for each fastq:
#	--total reads
#	--total bases
#	--q20/q30 bases
#	--q20/q30 rate
#	--mean length of read 1
#	--mean length of read 2
#	--gc_content
#	--duplication rate
#-----------------------------------------------------

library(ggplot2)
library(odbc)
library(DBI)
library(rjson)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

myDir1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/fastp"
myinf1 <- paste(myDir1,dir(myDir1),sep="/")
myinf1 <- unlist(sapply(myinf1,function(x)paste(x,"/",dir(x)[grep(".json",dir(x))],sep="")))
names(myinf1) <- gsub("/projects/varnf/GLASS-III/GLASS-III/results/kallisto/fastp/","",names(myinf1))

beforelist <- afterlist <- list()
duplication_rate <- rep(0,length(myinf1))
for(i in 1:length(myinf1))
{
	result <- fromJSON(file = myinf1[i])
	mysummary <- result$summary
	before <- mysummary$before_filtering
	after <- mysummary$after_filtering
	
	mycommand <- result$command
	before$file_path_R1 <- after$file_path_R1 <- sapply(strsplit(mycommand," "),function(x)x[3])
	before$file_path_R2 <- after$file_path_R2 <- sapply(strsplit(mycommand," "),function(x)x[5])
	
	beforelist[[i]] <- before
	afterlist[[i]] <- after
	
	result <- fromJSON(file = myinf1[i])
	duplication_rate[i] <- result$duplication$rate

}
names(beforelist) <- names(afterlist) <- names(myinf1)

before_table <- do.call(rbind,beforelist)
after_table <- do.call(rbind,afterlist)

before_table <- as.data.frame(apply(before_table,2,unlist))
before_table[,1:(ncol(before_table)-2)] <- apply(before_table[,1:(ncol(before_table)-2)],2,as.numeric)
before_table[,"aliquot_barcode"] <- substr(rownames(before_table),1,30)
rownames(before_table) <- NULL
before_table <- before_table[,c(ncol(before_table),1:(ncol(before_table)-1))]
before_table <- cbind(before_table,duplication_rate)

after_table <- as.data.frame(apply(after_table,2,unlist))
after_table[,1:(ncol(after_table)-2)] <- apply(after_table[,1:(ncol(after_table)-2)],2,as.numeric)
after_table[,"aliquot_barcode"] <- substr(rownames(after_table),1,30)
rownames(after_table) <- NULL
after_table <- after_table[,c(ncol(after_table),1:(ncol(after_table)-1))]
after_table <- cbind(after_table,duplication_rate)

fastp_applied <- as.factor(c(rep(0,nrow(before_table)),rep(1,nrow(after_table))))
final_table <- rbind(before_table,after_table)
final_table <- data.frame(final_table,fastp_applied)
final_table <- final_table[,c(1:10,13,11,12,14)]
final_table[,"fastp_applied"] <- as.numeric(as.character(final_table[,"fastp_applied"]))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/read_count.pdf")
ggplot(final_table, aes(x=total_reads,colour=fastp_applied)) + 
geom_density()
dev.off()

dbWriteTable(con, Id(schema="analysis", table="fastp_metrics"), final_table, overwrite=TRUE)
