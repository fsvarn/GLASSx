library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)

#######################################################
rm(list=ls())

myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/fpkm_table.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/rows-genes.csv"
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/columns-samples.csv"

gct_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/fpkm_table.gct"

fpkm <- read.csv(myinf1,header=TRUE)
fpkm <- fpkm[,-1]

genes <- read.csv(myinf2,stringsAsFactor=FALSE)
samps <- read.csv(myinf3,stringsAsFactor=FALSE)

rownames(fpkm) <- genes$gene_symbol

# ids <- paste(samps$tumor_id,"__",samps$structure_id,sep="")
# names(ids) <- samps$rna_well_id
# ids <- ids[colnames(fpkm)]

#colnames(fpkm) <- ids

fpkm <- data.frame(rownames(fpkm), fpkm)
colnames(fpkm) <- gsub("^X","",colnames(fpkm))

#Create gct file of gene expression matrix
Description <- rep(NA,nrow(fpkm))
fpkm2 <- cbind(fpkm[,1],Description,fpkm[,2:ncol(fpkm)])
colnames(fpkm2)[1] <- "NAME"
write.table(fpkm2, gct_path, sep="\t", quote=FALSE, row.names=FALSE)

conIn <- file(gct_path, "r")
rawfile = readLines(conIn)
close(conIn)

mytext <- c("#1.2", paste(nrow(fpkm2),"\t",(ncol(fpkm)-1),sep=""),rawfile)
conOut = file(gct_path, "w")
writeLines(mytext, conOut)
close(conOut)


# Run Qianghu's transcriptional classifier (only needed once)
runSsGSEAwithPermutation(gct_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/p_result_fpkm_table.gct.txt", stringsAsFactors=FALSE)
rownames(subtype_ssgsea) <- gsub("^X","",rownames(subtype_ssgsea))

subtype_ssgsea <- data.frame(rownames(subtype_ssgsea),subtype_ssgsea,stringsAsFactors=FALSE)

res <- subtype_ssgsea
rna_well_id <- rep(rownames(res),3)
signature_name <- c(rep("Proneural",nrow(res)),rep("Classical",nrow(res)),rep("Mesenchymal",nrow(res)))
enrichment_score <- c(res[,"Proneural"],res[,"Classical"],res[,"Mesenchymal"])
p_value <- c(res[,"Proneural_pval"],res[,"Classical_pval"],res[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(rna_well_id,signature_name,enrichment_score,p_value,stringsAsFactors=FALSE)
transcriptional_subtype[,"rna_well_id"] <- as.integer(transcriptional_subtype[,"rna_well_id"])
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"rna_well_id"]),]


full <- transcriptional_subtype %>%
inner_join(samps,by="rna_well_id") %>%
filter(!tumor_name %in% c("W4-1-1", "W3-1-1", "W10-1-1", "W31-1-1","W35-1-1","W48-1-1","W53-1-1")) %>% #Remove IDHmut
filter(p_value < 0.05) %>%
group_by(rna_well_id) %>%
arrange(signature_name) %>%
summarise(signature_name = paste(signature_name, collapse=","),tumor_id=first(tumor_id)) # ,tumor_name,structure_name) 

patient <- full %>%
group_by(tumor_id) %>%
arrange(signature_name) %>%
summarise(signature_name = paste(unique(signature_name), collapse=",")) %>%
ungroup() %>%
mutate(signature_name = unlist(sapply(strsplit(signature_name, ","),function(x)paste(unique(x[order(x)]),collapse=",")))) %>%
mutate(signature_name = unlist(sapply(strsplit(signature_name, ","),function(x)paste(unique(x[order(x)]),collapse=","))))
 
counts <- table(patient$signature_name)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- 
"WITH full_data AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	CONCAT(ts1.signature_name, ',', ts2.signature_name) AS subtype_switch
	FROM analysis.rna_silver_set ss
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b 
	JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
	WHERE idh_codel_subtype = 'IDHwt'
)
SELECT subtype_switch, COUNT(*)
FROM full_data
GROUP BY subtype_switch
"
dat <- dbGetQuery(con,q)
dat$count <- as.numeric(dat$count)

ct <- matrix(c(9,sum(dat$count)-9,7,sum(counts)-7),nrow=2,ncol=2)
fisher.test(ct)

ct <- matrix(c(20,sum(dat$count)-20,12,sum(counts)-12),nrow=2,ncol=2)
fisher.test(ct)
