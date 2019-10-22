library(estimate)
library(DBI)
library(odbc)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())

#Load and prepare gene expression matrix

filterCommonGenes("/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv", "/projects/varnf/GLASS-III/GLASS-III/data/sigs/gene_tpm_matrix_all_samples_estimate_filtered.tsv",id="GeneSymbol")
estimateScore("/projects/varnf/GLASS-III/GLASS-III/data/sigs/gene_tpm_matrix_all_samples_estimate_filtered.tsv","/projects/varnf/GLASS-III/GLASS-III/data/sigs/gene_tpm_matrix_all_samples_estimate_scores.txt","affymetrix")

estimate <- t(read.delim("/projects/varnf/GLASS-III/GLASS-III/data/sigs/gene_tpm_matrix_all_samples_estimate_scores.txt"))
colnames(estimate) <- estimate[2,]
estimate <- estimate[-c(1,2),]
aliquot_barcode <- gsub("\\.","-",estimate[,2])
estimate <- estimate[,3:5]
estimate <- apply(estimate,2,as.numeric)
est_purity <- cos(0.6049872018+0.0001467884 * estimate[,"ESTIMATEScore"])
estimate <- data.frame(aliquot_barcode,estimate,est_purity)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Gold set
q <- "SELECT al.sample_barcode, tumor_barcode, an.rna_barcode, tp.purity AS titan_purity, sp.cellularity AS seqz_cellularity, al.aliquot_batch
FROM variants.titan_params tp
JOIN variants.seqz_params sp ON sp.pair_barcode = tp.pair_barcode
JOIN analysis.pairs pa ON pa.pair_barcode = tp.pair_barcode
JOIN analysis.analyte_sets an ON pa.tumor_barcode = an.dna_barcode
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = an.dna_barcode OR gs.tumor_barcode_b = an.dna_barcode
JOIN biospecimen.aliquots al ON an.rna_barcode = al.aliquot_barcode"

pur <- dbGetQuery(con,q)

all_purity <- merge(pur,estimate,by.x="rna_barcode",by.y="aliquot_barcode")

wxs_purity <- all_purity[grep("-WXS-",all_purity[,"tumor_barcode"]),]

mycors <- cor(wxs_purity[,c(4,5,10)],method="s")
mycorp <- cor(wxs_purity[,c(4,5,10)],method="p")
titan_cor <- round(mycorp[1,3],2)
seqz_cor <- round(mycorp[2,3],2)

titan_lm <- lm(wxs_purity[,"titan_purity"]~wxs_purity[,"est_purity"])
seqz_lm <- lm(wxs_purity[,"seqz_cellularity"]~wxs_purity[,"est_purity"])


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/purity_plots.pdf",width=7,height=7)
ggplot(data = wxs_purity, aes(x = est_purity, y = titan_purity,colour=aliquot_batch)) +
  geom_point() +
  ggtitle("ESTIMATE purity versus titan purity") +
  scale_colour_manual(values=brewer.pal(11,"Set3")) +
  annotate("text", x = 0.1, y = 0, label = paste("R = ",titan_cor,sep="")) +
  theme_classic()+
  geom_abline(intercept = titan_lm$coefficients[1], slope = titan_lm$coefficients[2]) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))


ggplot(data = wxs_purity, aes(x = est_purity, y = seqz_cellularity,colour=aliquot_batch)) +
  geom_point() +
  ggtitle("ESTIMATE purity versus sequenza cellularity") +
  scale_colour_manual(values=brewer.pal(11,"Set3")) +
  annotate("text", x = 0.1, y = 0, label = paste("R = ",seqz_cor,sep="")) +
  theme_classic() +
  geom_abline(intercept = seqz_lm$coefficients[1], slope = seqz_lm$coefficients[2]) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))

dev.off()


#All samples
q <- "SELECT al.sample_barcode, tumor_barcode, an.rna_barcode, tp.purity AS titan_purity, sp.cellularity AS seqz_cellularity, al.aliquot_batch
FROM variants.titan_params tp
JOIN variants.seqz_params sp ON sp.pair_barcode = tp.pair_barcode
JOIN analysis.pairs pa ON pa.pair_barcode = tp.pair_barcode
JOIN analysis.analyte_sets an ON pa.tumor_barcode = an.dna_barcode
JOIN biospecimen.aliquots al ON an.rna_barcode = al.aliquot_barcode
"

pur <- dbGetQuery(con,q)

all_purity <- merge(pur,estimate,by.x="rna_barcode",by.y="aliquot_barcode")

wxs_purity <- all_purity[grep("-WXS-",all_purity[,"tumor_barcode"]),]

cor(wxs_purity[,c(4,5,10)],method="s")

#---------------------------------

#Check TCGA purity with previous values

tcga_pur <- read.delim("/projects/varnf/Data/TCGA/TCGA_purity.txt")

glass_tcga <- all_purity[grep("TCGA-",all_purity[,"sample_barcode"]),]
glass_tcga[,"sample_barcode"] <- gsub("-TP$","-01A", glass_tcga[,"sample_barcode"])
glass_tcga[,"sample_barcode"] <- gsub("-R1$","-02A", glass_tcga[,"sample_barcode"])
glass_tcga[,"sample_barcode"] <- gsub("-R2$","-02B", glass_tcga[,"sample_barcode"])

#Manual change
glass_tcga[13,"sample_barcode"] <- "TCGA-14-1034-02B"
glass_tcga[32,"sample_barcode"] <- "TCGA-14-1034-01B"

tcga_all <- merge(glass_tcga,tcga_pur,by.x="sample_barcode",by.y="Sample.ID")
cor(tcga_all[,"est_purity"],tcga_all[,"ESTIMATE"],method="s",use="complete")

#---------------------------------
#Compare RPKM ESTIMATE purity with TPM ESTIMATE purity

expr <- read.delim("/projects/varnf/GLASS/data/CIBERSORT/synapse_final_RNAseq.txt",sep="\t",header=T,stringsAsFactor=F)
expr <- expr[,c(1,20:ncol(expr))]
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- t(expr)
gene_symbol <- rownames(expr)
expr <- data.frame(gene_symbol,expr)
write.table(expr,"/projects/varnf/GLASS/data/CIBERSORT/RPKM_matrix.txt",sep="\t",row.names=FALSE,quote=FALSE)

filterCommonGenes("/projects/varnf/GLASS/data/CIBERSORT/RPKM_matrix.txt", "/projects/varnf/GLASS/data/CIBERSORT/RPKM_matrix_filtered.txt",id="GeneSymbol")
estimateScore("/projects/varnf/GLASS/data/CIBERSORT/RPKM_matrix_filtered.txt","/projects/varnf/GLASS/data/CIBERSORT/RPKM_matrix_filtered_estimate.txt","affymetrix")


estimate <- t(read.delim("/projects/varnf/GLASS-III/GLASS-III/data/sigs/gene_tpm_matrix_all_samples_estimate_scores.txt"))
colnames(estimate) <- estimate[2,]
estimate <- estimate[-c(1,2),]
aliquot_barcode <- gsub("\\.","-",estimate[,2])
estimate <- estimate[,3:5]
estimate <- apply(estimate,2,as.numeric)
est_purity <- cos(0.6049872018+0.0001467884 * estimate[,"ESTIMATEScore"])
estimate <- data.frame(aliquot_barcode,estimate,est_purity)
tpm_estimate <- estimate

estimate <- t(read.delim("/projects/varnf/GLASS/data/CIBERSORT/RPKM_matrix_filtered_estimate.txt"))
colnames(estimate) <- estimate[2,]
estimate <- estimate[-c(1,2),]
aliquot_barcode <- gsub("\\.","-",estimate[,2])
estimate <- estimate[,3:5]
estimate <- apply(estimate,2,as.numeric)
est_purity <- cos(0.6049872018+0.0001467884 * estimate[,"ESTIMATEScore"])
estimate <- data.frame(aliquot_barcode,estimate,est_purity)
rpkm_estimate <- estimate

tpm_estimate[,"sample_barcode"] <- substr(tpm_estimate[,"aliquot_barcode"],1,15)

full_estimate <- merge(tpm_estimate,rpkm_estimate,by.x="sample_barcode",by.y="aliquot_barcode")

cor(full_estimate[,"est_purity.x"],full_estimate[,"est_purity.y"],method="s")

