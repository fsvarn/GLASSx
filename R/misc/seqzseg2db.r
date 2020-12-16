library(DBI)
library(odbc)
library(tidyverse)

rm(list=ls())

myDir1 <- "results/sequenza/seqzR/"
mytag <- dir(myDir1)
myDir1 <- paste(myDir1, mytag, sep="")
myinf1 <- sapply(myDir1, function(x)paste(x, dir(x), sep = "/"))
myinf1 <- unlist(myinf1)
myinf1 <- myinf1[grep("_segments.txt", myinf1)]

full_table <- do.call(read.delim, myinf1)

tabs <- list()
for(i in 1:length(myinf1))
{
	cat("\r",i)
	seqz <- read.delim(myinf1[i], header=TRUE)
	pair_barcode <- mytag[i]
	
	tabs[[i]] <- data.frame(pair_barcode, seqz)
}

full_table <- do.call(rbind, tabs)

df <- full_table %>%
   	  transmute(pair_barcode = pair_barcode,
              chrom = chromosome,
              pos = sprintf("[%s,%s]", start.pos, end.pos),
              baf = Bf,
              baf_n = N.BAF,
              baf_sd = sd.BAF,
              ratio = depth.ratio,
              ratio_n = N.ratio,
              ratio_sd = sd.ratio,
              copy_number = CNt,
              major_cn = A,
              minor_cn = B,
              log_posterior_proba = LPP)
              
# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

old <- dbReadTable(con, Id(schema="variants",table="seqz_seg"))
old_pairs <- unique(old[,"pair_barcode"])

df_new <- df[-which(df[,"pair_barcode"] %in% old_pairs),]

# Limit to new cohorts (GLSS-PD-WXS, GLSS-H2-WXS)
df_new <- df_new[grep("GLSS-HF-|GLSS-CU-P",df_new[,"pair_barcode"]),]
df_new[,"chrom"] <- as.character(df_new[,"chrom"])
df_new[which(df_new[,"chrom"]=='X'),"chrom"] <- 23
df_new[,"chrom"] <- as.numeric(df_new[,"chrom"])
# Write new results to db
#dbWriteTable(con, Id(schema="variants",table="seqz_seg"), df_new, append = TRUE)

# Upload seq_params
#-------------------------------------------
# Bash script to build seqz_params:
# echo -e "pair_barcode\tcellularity\tploidy\tslpp" > seqz_params.txt
# for i in seqzR/*/*_alternative_solutions.txt
# do 
# 	V1=$(echo $i | cut -c 7-35)
# 	V2=$(head -2 "$i" | tail -1)
# 	echo -e "$V1\t$V2" >> seqz_params.txt
# done 

# Add parameters to db
dat <- read.delim("results/sequenza/seqz_params.txt")

old <- dbReadTable(con, Id(schema="variants",table="seqz_params"))
old_pairs <- unique(old[,"pair_barcode"])

dat_new <- dat[-which(dat[,"pair_barcode"] %in% old_pairs),]

# Limit to new cohorts (GLSS-PD-WXS, GLSS-H2-WXS)
dat_new <- dat_new[grep("GLSS-HF-|GLSS-CU-P",dat_new[,"pair_barcode"]),]

# Write new results to db
dbWriteTable(con, Id(schema="variants",table="seqz_params"), dat_new, append = TRUE)
