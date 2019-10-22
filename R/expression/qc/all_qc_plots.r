library(ggplot2)
library(odbc)
library(DBI)
library(rjson)

rm(list=ls())

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT fq.aliquot_barcode, kq.n_processed, kq.n_pseudoaligned, kq.n_unique, kq.p_pseudoaligned, kq.p_unique, fq.total_reads, fq.q20_rate, fq.q30_rate, fq.gc_content,fq.duplication_rate,al.aliquot_batch 
FROM analysis.fastp_metrics fq 
LEFT JOIN analysis.kallisto_qc kq ON kq.aliquot_barcode = fq.aliquot_barcode
LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = fq.aliquot_barcode
WHERE fastp_applied = 1 
ORDER BY  kq.p_pseudoaligned"

qc <- dbGetQuery(con,q)

qc[,"aliquot_batch"] <- gsub("-RNA","",qc[,"aliquot_batch"])

q <- "SELECT kq.*,al.aliquot_batch
FROM analysis.kallisto_qc kq 
LEFT JOIN biospecimen.aliquots al ON al.aliquot_barcode = kq.aliquot_barcode
ORDER BY  kq.p_pseudoaligned"

ka <- dbGetQuery(con,q)

qc[,"aliquot_batch"] <- gsub("-RNA","",qc[,"aliquot_batch"])
ka[,"aliquot_batch"] <- gsub("-RNA","",ka[,"aliquot_batch"])


unique_thr <- mean(ka[,"p_unique"]) - 2*sd(ka[,"p_unique"])
qc[,"premium"] <- as.factor(as.numeric(qc[,"p_unique"] < unique_thr))
ka[,"premium"] <- as.factor(as.numeric(ka[,"p_unique"] < unique_thr))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/qc_metrics_by_batch.pdf",width=7,height=4)

#total reads 
ggplot(qc,aes(y = total_reads, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "Total reads", x = "Aliquot batch", y = "Total reads") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")

#q20_rate
ggplot(qc,aes(y = q20_rate, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "q20 rate", x = "Aliquot batch", y = "q20 rate") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,1))

#q30_rate
ggplot(qc,aes(y = q30_rate, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "q30 rate", x = "Aliquot batch", y = "q30 rate") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,1))

#gc_content
ggplot(qc,aes(y = gc_content, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "GC content", x = "Aliquot batch", y = "GC content") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,1))

#duplication_rate
ggplot(qc,aes(y = duplication_rate, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "Duplication rate", x = "Aliquot batch", y = "Duplication rate") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,1))

#n_processed 
ggplot(ka,aes(y = n_processed, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "Transcripts processed", x = "Aliquot batch", y = "Transcripts processed") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")

#Proportion pseudoaligned
ggplot(ka,aes(y = p_pseudoaligned, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "Proportion pseudoaligned", x = "Aliquot batch", y = "Proportion pseudoaligned") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,100))

#Proportion unique
ggplot(ka,aes(y = p_unique, x = aliquot_batch)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(width=0.3,size=0.8, aes(colour=premium)) +
scale_colour_manual(values=c("black","red")) +
labs(title = "Proportion unique", x = "Aliquot batch", y = "Proportion unique") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,100))

dev.off()