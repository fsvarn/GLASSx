#Code to compare the PRADA RPKM results from Wang et al Cancer Cell 2015 to the kallisto results now
#Performs correlation analyses of RPKM (PRADA) and TPM (kallisto) for each aliquot present in both datasets
#-----------------------------------------------------

library(DBI)
library(odbc)
library(ggplot2)
library(reshape)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT r1.sample_barcode AS rna_barcode, tpm.aliquot_barcode AS tpm_barcode, r1.gene_symbol, rpkm, tpm
FROM analysis.rnaseq r1
JOIN analysis.gene_tpm tpm ON substr(tpm.aliquot_barcode,1,15) = r1.sample_barcode AND tpm.gene_symbol = r1.gene_symbol"

expr <- dbGetQuery(con,q)

aliquots <- unique(expr[,"tpm_barcode"])
mycor <- rep(0,length(aliquots))
for(i in 1:length(aliquots))
{
	cat("\r",i)
	sub_expr <- expr[which(expr[,"tpm_barcode"]==aliquots[i]),]
	mycor[i] <- cor(sub_expr[,"rpkm"],sub_expr[,"tpm"],method="s",use="complete.obs")
}
res <- data.frame(aliquots,mycor)
range(res[,"mycor"])
mean(res[,"mycor"])
median(res[,"mycor"])

# > range(res[,"mycor"])
# [1] 0.8941404 0.9551054
# > mean(res[,"mycor"])
# [1] 0.9383039
# > median(res[,"mycor"])
# [1] 0.9439455

#Check inter-sample TPM vs RPKM
myres <- matrix(0,nrow=length(aliquots),ncol=length(aliquots))
rownames(myres) <- colnames(myres) <- aliquots
for(i in 1:length(aliquots))
{
	for(j in 1:length(aliquots))
	{
		cat("\r",i,"-->",j)
		sub_expr1 <- expr[which(expr[,"tpm_barcode"]==aliquots[i]),]
		sub_expr2 <- expr[which(expr[,"tpm_barcode"]==aliquots[j]),]
		myres[i,j] <- cor(sub_expr1[,"tpm"],sub_expr2[,"rpkm"],method="s",use="complete.obs")
	}
}
range(diag(myres))
mean(diag(myres))
median(diag(myres))

# > range(diag(myres))
# [1] 0.8941404 0.9551054
# > mean(diag(myres))
# [1] 0.9383039
# > median(diag(myres))
# [1] 0.9439455

range(myres[upper.tri(myres)])
mean(myres[upper.tri(myres)])
median(myres[upper.tri(myres)])

# > range(myres[upper.tri(myres)])
# [1] -0.02629383  0.94715435
# > mean(myres[upper.tri(myres)])
# [1] 0.01414496
# > median(myres[upper.tri(myres)])

#Check inter-sample TPM vs TPM
myres <- matrix(0,nrow=length(aliquots),ncol=length(aliquots))
rownames(myres) <- colnames(myres) <- aliquots
for(i in 1:length(aliquots))
{
	for(j in 1:length(aliquots))
	{
		cat("\r",i,"-->",j)
		sub_expr1 <- expr[which(expr[,"tpm_barcode"]==aliquots[i]),]
		sub_expr2 <- expr[which(expr[,"tpm_barcode"]==aliquots[j]),]
		myres[i,j] <- cor(sub_expr1[,"tpm"],sub_expr2[,"tpm"],method="s",use="complete.obs")
	}
}
range(diag(myres))
mean(diag(myres))
median(diag(myres))

# > range(diag(myres))
# [1] 1 1
# > mean(diag(myres))
# [1] 1
# > median(diag(myres))
# [1] 1

range(myres[upper.tri(myres)])
mean(myres[upper.tri(myres)])
median(myres[upper.tri(myres)])

# > range(myres[upper.tri(myres)])
# [1] -0.02657313  0.98609227
# > mean(myres[upper.tri(myres)])
# [1] 0.01468341
# > median(myres[upper.tri(myres)])
# [1] 0.001844946

#---------------------------------

#Compare TPM matrix to TPM long format (FIXED)

library(DBI)
library(odbc)

gene_mat <- read.delim("/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv")
rownames(gene_mat) <- gene_mat[,1]
gene_mat <- gene_mat[,-1]
colnames(gene_mat) <- gsub("\\.","-",colnames(gene_mat))

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT aliquot_barcode, gene_symbol, tpm
FROM analysis.gene_tpm"

tpm <- dbGetQuery(con,q)
tpm <- tpm[-which(is.na(tpm[,"gene_symbol"])),]

aliquots <- colnames(gene_mat)

mycor <- rep(0,length(aliquots))
for(i in 1:length(aliquots))
{
	mysamp <- tpm[which(tpm[,"aliquot_barcode"]==aliquots[i]),]
	mymat <- gene_mat[mysamp[,"gene_symbol"],aliquots[i]]
	
	mycor[i] <- cor(mysamp[,"tpm"],mymat,method="s")
}

range(mycor)
mean(mycor)
median(mycor)

# > range(mycor)
# [1] 1 1
# > mean(mycor)
# [1] 1
# > median(mycor)
# [1] 1

#---------------------------------
#Compare RPKM matrix to TPM matrix (FIXED)

expr <- read.delim("/projects/varnf/GLASS/data/CIBERSORT/synapse_final_RNAseq.txt",sep="\t",header=T,stringsAsFactor=F)
expr <- expr[,c(1,20:ncol(expr))]
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- t(expr)

gene_mat <- read.delim("/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv")
rownames(gene_mat) <- gene_mat[,1]
gene_mat <- gene_mat[,-1]

colnames(gene_mat) <- substr(colnames(gene_mat),1,15)
colnames(gene_mat) <- gsub("\\.","-",colnames(gene_mat))

comxx <- intersect(rownames(expr),rownames(gene_mat))
expr <- expr[comxx,]
gene_mat <- gene_mat[comxx,]

comxx <- intersect(colnames(expr),colnames(gene_mat))
expr <- expr[,comxx]
gene_mat <- gene_mat[,comxx]

mycor <- rep(0,ncol(expr))
for(i in 1:ncol(expr))
{
	mycor[i] <- cor(expr[,i],gene_mat[,i],method="s",use="complete")
}

range(mycor)
mean(mycor)
median(mycor)

# > range(mycor)
# [1] 0.8900890 0.9461537
# > mean(mycor)
# [1] 0.9307283
# > median(mycor)
# [1] 0.9323431
