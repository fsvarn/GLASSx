###################################################
# Compare TPM for NUPR1a tpm between primary and recurrent samples
# Collaboration with H.K. Ng group
# Updated: 2020.09.23
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)
library(survival)

#######################################################
rm(list=ls())

#Establish connection
myinf1 <- "results/kallisto/kallisto/final/transcript_tpm_matrix_all_samples.tsv"

# Read in expression matrix
orig_dat <- read.delim(myinf1,row.names=1)
colnames(orig_dat) <- gsub("\\.","-",colnames(orig_dat))

# Subset to publicly available GLASS data:
dat <- orig_dat[,c("GLSS-SM-R060-TP-01R-RNA-ZI9YVM", "GLSS-SM-R064-TP-01R-RNA-W6PIW2", "GLSS-SM-R060-R1-01R-RNA-6SGM43", "GLSS-SM-R056-R2-01R-RNA-DHLJ2T", "GLSS-SM-R064-R1-01R-RNA-GLLBN3", "GLSS-HF-2548-TP-01R-RNA-4R6KR6", "GLSS-HF-2829-TP-01R-RNA-90I02Y", "GLSS-HF-2869-TP-01R-RNA-2H6MOF", "GLSS-HF-2548-R1-01R-RNA-F7XJ9B", "GLSS-HF-2919-TP-01R-RNA-3AV46F", "GLSS-HF-2934-TP-01R-RNA-WUQKHB", "GLSS-HF-2829-R1-01R-RNA-D71HI8", "GLSS-HF-2869-R3-01R-RNA-TF6NT8", "GLSS-HF-2998-TP-01R-RNA-5H7ODC", "GLSS-HF-2998-R1-01R-RNA-2HE04E", "GLSS-HF-2934-R1-01R-RNA-IBE7TN", "GLSS-HF-3050-TP-01R-RNA-DR0HT8", "GLSS-HF-2919-R1-01R-RNA-R6L8SO", "GLSS-HF-3081-TP-01R-RNA-EE79OQ", "GLSS-HF-3162-TP-01R-RNA-K5P3FT", "GLSS-HF-3081-R2-01R-RNA-939287", "GLSS-HF-3162-R1-01R-RNA-JUILFB", "GLSS-HF-3050-R1-01R-RNA-7FET5L", "GLSS-SM-R072-R1-01R-RNA-1DQIND", "GLSS-SM-R099-TP-01R-RNA-YGXA72", "GLSS-SM-R100-TP-01R-RNA-EUI7AZ", "GLSS-SM-R106-TP-01R-RNA-P5JI35", "GLSS-SM-R104-TP-01R-RNA-G52UET", "GLSS-SM-R109-TP-01R-RNA-I4HTPP", "GLSS-SM-R099-R1-01R-RNA-MNTPMI", "GLSS-SM-R100-R1-01R-RNA-46UW5U", "GLSS-SM-R101-TP-02R-RNA-V4TRVD", "GLSS-SM-R106-R1-01R-RNA-MFHIBQ", "GLSS-SM-R107-TP-01R-RNA-ZSYP37", "GLSS-SM-R108-TP-01R-RNA-KUU8QP", "GLSS-SM-R104-R1-01R-RNA-YF3DT7", "GLSS-SM-R109-R1-01R-RNA-R5I1C7", "GLSS-SM-R110-TP-02R-RNA-A85T37", "GLSS-SM-R101-R1-02R-RNA-BIXXBZ", "GLSS-SM-R110-R1-01R-RNA-F7JINL", "GLSS-SM-R108-R1-01R-RNA-NICL65", "GLSS-SM-R107-R1-02R-RNA-1K9XUJ", "GLSS-SM-R103-TP-01R-RNA-7YEFHB", "GLSS-SM-R103-R1-01R-RNA-TSX1I7", "GLSS-SM-R061-TP-01R-RNA-AZSPBW", "GLSS-SM-R061-R1-01R-RNA-0NG539", "GLSS-SM-R095-TP-01R-RNA-5T2HT4", "GLSS-SM-R095-R1-01R-RNA-BQ7VMI", "GLSS-SM-R072-TP-01R-RNA-K18OXT", "GLSS-SM-R056-TP-01R-RNA-JWKQQG", "GLSS-SM-R067-TP-01R-RNA-JYXPKW", "GLSS-SM-R067-R1-01R-RNA-O82WJY", "TCGA-06-0125-TP-01R-RNA-RCQ5QS", "TCGA-06-0125-R1-11R-RNA-EKH06U", "TCGA-06-0190-TP-01R-RNA-S2074S", "TCGA-06-0190-R1-01R-RNA-HWTLG4", "TCGA-06-0210-TP-01R-RNA-JMBENW", "TCGA-06-0210-R1-01R-RNA-WB10MK", "TCGA-06-0211-TP-01R-RNA-Y9DKHR", "TCGA-06-0211-R1-02R-RNA-3YQ1G4", "TCGA-14-1034-TP-01R-RNA-S6SBZH", "TCGA-14-1034-R1-01R-RNA-GSNR8T", "TCGA-DH-A669-TP-12R-RNA-XG7CTS", "TCGA-DH-A669-R1-11R-RNA-KJQZGQ", "TCGA-DU-5870-TP-11R-RNA-CSQWQZ", "TCGA-DU-5870-R1-12R-RNA-GIZUWN", "TCGA-DU-5872-TP-11R-RNA-HS9DNI", "TCGA-DU-5872-R1-21R-RNA-0EW2ZH", "TCGA-DU-6397-TP-11R-RNA-NYNT65", "TCGA-DU-6397-R1-12R-RNA-N760C1", "TCGA-DU-6404-TP-11R-RNA-A7XMXJ", "TCGA-DU-6404-R1-21R-RNA-LJZDP3", "TCGA-DU-6407-TP-13R-RNA-6NHIOE", "TCGA-DU-6407-R1-12R-RNA-T3UMCW", "TCGA-DU-7304-TP-12R-RNA-23D1EM", "TCGA-DU-7304-R1-12R-RNA-OMXRQS", "TCGA-FG-5963-TP-11R-RNA-2Y6I4P", "TCGA-FG-5963-R1-12R-RNA-1EJ1SV", "TCGA-FG-5965-TP-11R-RNA-0XBPNY", "TCGA-FG-5965-R1-11R-RNA-QSA7Y3", "TCGA-FG-A4MT-TP-11R-RNA-K4BQG0", "TCGA-FG-A4MT-R1-11R-RNA-9OU9UC", "TCGA-TM-A7CF-TP-11R-RNA-YLP5PU", "TCGA-TM-A7CF-R1-11R-RNA-EZWWA2", "TCGA-TQ-A7RK-TP-11R-RNA-LRGF1D", "TCGA-TQ-A7RK-R1-11R-RNA-PKFXB3", "TCGA-TQ-A7RV-TP-21R-RNA-6F23HF", "TCGA-TQ-A7RV-R1-11R-RNA-TIRW8I", "TCGA-TQ-A8XE-TP-11R-RNA-0YPLAZ", "TCGA-TQ-A8XE-R1-11R-RNA-47VZAQ", "GLSS-SM-R071-TP-01R-RNA-OAXGI8", "GLSS-SM-R071-R1-01R-RNA-7AZ6G2", "GLSS-SM-R068-TP-01R-RNA-0UPMYO", "GLSS-SM-R068-R1-01R-RNA-7I5H9P", "GLSS-SM-R065-TP-01R-RNA-OM7S2C", "GLSS-SM-R065-R1-01R-RNA-V35ETA", "GLSS-SM-R070-TP-01R-RNA-QWG60O", "GLSS-SM-R070-R1-01R-RNA-PXOL5Z", "GLSS-SM-R063-TP-01R-RNA-WIVUSX", "GLSS-SM-R063-R1-01R-RNA-YQOR0G", "GLSS-SM-R066-TP-01R-RNA-G9A6CK", "GLSS-SM-R066-R1-01R-RNA-MKHTJ4")]
#dat <- orig_dat
dat <- log10(dat+1)
#pub_inds <- c(grep("GLSS-HF-2",colnames(dat)), grep("GLSS-HF-3",colnames(dat)),grep("GLSS-SM",colnames(dat)), grep("TCGA-",colnames(dat)))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in clinical data
q <- "
SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
--WHERE ts1.signature_name != 'Mesenchymal' AND ts2.signature_name = 'Mesenchymal'
"

info <- dbGetQuery(con, q)
info <- info[which(info[,"tumor_barcode_a"] %in% colnames(dat) & info[,"tumor_barcode_b"] %in% colnames(dat)),]

# Test in IDHwt tumors only:
g1 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]=="IDHwt"),"tumor_barcode_a"]])
g2 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]=="IDHwt"),"tumor_barcode_b"]])
idhwt_pval <- formatC(wilcox.test(g1,g2,paired=TRUE)$p.value, format = "e", digits = 0)
mean(g1)
sd(g1)/sqrt(length(g1))
mean(g2)
sd(g2)/sqrt(length(g2))


g1 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]!="IDHwt"),"tumor_barcode_a"]])
g2 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]!="IDHwt"),"tumor_barcode_b"]])
idhmut_pval <- formatC(wilcox.test(g1,g2,paired=TRUE)$p.value, format = "e", digits = 0)
mean(g1)
sd(g1)/sqrt(length(g1))
mean(g2)
sd(g2)/sqrt(length(g2))

# Make a plot
aliquot_barcode <- c(info[,"tumor_barcode_a"], info[,"tumor_barcode_b"])
case_barcode <- substr(aliquot_barcode,1,12)
idh_codel_subtype <- rep(info[,"idh_codel_subtype"],2)
idhmut <- as.character(idh_codel_subtype != "IDHwt")  %>% recode("TRUE" = "IDHmut","FALSE" = "IDHwt")
timepoint <- c(rep("Initial",nrow(info)),rep("Recurrent",nrow(info)))
tpm <- t(dat["ENST00000395641",aliquot_barcode])
plot_dat <- data.frame(aliquot_barcode, case_barcode, tpm, timepoint, idh_codel_subtype,idhmut)
rownames(plot_dat) <- NULL
colnames(plot_dat) <- c("aliquot_barcode","case_barcode","tpm","timepoint","idh_codel_subtype","idhmut")

p_val <- c(deparse((bquote(italic("P") ~" = " ~ .(idhwt_pval)))),
		   deparse((bquote(italic("P") ~" = " ~ .(idhmut_pval)))))
annotation_text <- data.frame(idhmut = factor(c("IDHwt","IDHmut")),
							  timepoint = 1.1,
							  tpm = 2.0,
							  p_val)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ENST00000395641_longitudinal_boxplot.pdf",width=4,height=2)
#pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ENST00000395641_longitudinal_boxplot_all_samps.pdf",width=4,height=2)
p1 <- ggplot(plot_dat, aes(x = timepoint, y = tpm)) +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode,colour= idh_codel_subtype)) +
geom_point(size=1) +
facet_grid(.~idhmut) +
labs(main = "ENST00000395641", y = "log10(TPM)") +
geom_text(data=annotation_text,label=p_val, size=2.5, parse=TRUE) +
scale_fill_manual(values=c("white","white")) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7))
p1
dev.off()


# Read in clinical data. Only choose tumors that acquire mesenchymal
q <- "
SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE ts1.signature_name != 'Mesenchymal' AND ts2.signature_name = 'Mesenchymal'
"

info <- dbGetQuery(con, q)
info <- info[which(info[,"tumor_barcode_a"] %in% colnames(dat) & info[,"tumor_barcode_b"] %in% colnames(dat)),]

# Test in IDHwt tumors only:
g1 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]=="IDHwt"),"tumor_barcode_a"]])
g2 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]=="IDHwt"),"tumor_barcode_b"]])
#idhwt_pval <- formatC(wilcox.test(g1,g2,paired=TRUE)$p.value, format = "e", digits = 0)
idhwt_pval <- round(wilcox.test(g1,g2,paired=TRUE)$p.value, 2)

mean(g1)
sd(g1)/sqrt(length(g1))
mean(g2)
sd(g2)/sqrt(length(g2))


g1 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]!="IDHwt"),"tumor_barcode_a"]])
g2 <- as.numeric(dat["ENST00000395641",info[which(info[,"idh_codel_subtype"]!="IDHwt"),"tumor_barcode_b"]])
idhmut_pval <- formatC(wilcox.test(g1,g2,paired=TRUE)$p.value, format = "e", digits = 0)



# Make a plot
aliquot_barcode <- c(info[,"tumor_barcode_a"], info[,"tumor_barcode_b"])
case_barcode <- substr(aliquot_barcode,1,12)
idh_codel_subtype <- rep(info[,"idh_codel_subtype"],2)
idhmut <- as.character(idh_codel_subtype != "IDHwt")  %>% recode("TRUE" = "IDHmut","FALSE" = "IDHwt")
timepoint <- c(rep("Initial",nrow(info)),rep("Recurrent",nrow(info)))
tpm <- t(dat["ENST00000395641",aliquot_barcode])
plot_dat <- data.frame(aliquot_barcode, case_barcode, tpm, timepoint, idh_codel_subtype,idhmut)
rownames(plot_dat) <- NULL
colnames(plot_dat) <- c("aliquot_barcode","case_barcode","tpm","timepoint","idh_codel_subtype","idhmut")

p_val <- c(deparse((bquote(italic("P") ~" = " ~ .(idhwt_pval)))))
annotation_text <- data.frame(timepoint = 1.5,
							  tpm = 2.0,
							  p_val)

plot_dat <- plot_dat %>% filter(idhmut=="IDHwt")
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ENST00000395641_longitudinal_mesenchymal_boxplot.pdf",width=1.5,height=2)
#pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ENST00000395641_longitudinal_mesenchymal_boxplot_all.pdf",width=1.5,height=2)
p1 <- ggplot(plot_dat, aes(x = timepoint, y = tpm)) +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "#619CFF") +
geom_point(size=1) +
labs(main = "ENST00000395641", y = "log10(TPM)") +
geom_text(data=annotation_text,label=p_val, size=2.5, parse=TRUE) +
scale_fill_manual(values=c("white","white")) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7))
p1
dev.off()



# Survival analyses
q <- "
SELECT * FROM analysis.tumor_rna_clinical_comparison WHERE idh_codel_subtype = 'IDHwt'
"
info <- dbGetQuery(con, q)
info <- info[which(info[,"tumor_barcode_a"] %in% colnames(dat) & info[,"tumor_barcode_b"] %in% colnames(dat)),]

diff <- dat["ENST00000395641",info[,"tumor_barcode_b"]] - dat["ENST00000395641",info[,"tumor_barcode_a"]]

g1 <- names(diff[which(diff > 0)])
g2 <- names(diff[which(diff < 0)])

surg_int <- info[,"surgical_interval"]
names(surg_int) <- info[,"tumor_barcode_b"]

reg <- c(rep("up", length(g1)), rep("dn", length(g2)))
event <- rep(1, length(reg))
int <- surg_int[c(g1, g2)]
survframe <- data.frame(reg, event, int)
mycox <- coxph(Surv(int, event) ~ reg, data = survframe)
summary(mycox)
survdiff(Surv(int, event) ~ reg, data = survframe)




# Quick query to look at average single gene expression levels across subtypes:
SELECT idh_codel_subtype, AVG(gt1.tpm), AVG(gt2.tpm)
FROM analysis.rna_silver_set ss
JOIN analysis.gene_tpm gt1 ON gt1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.gene_tpm gt2 ON gt2.aliquot_barcode = ss.tumor_barcode_b AND gt2.gene_symbol = gt1.gene_symbol
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
WHERE gt1.gene_symbol = 'NUPR1'
GROUP BY cs.idh_codel_subtype
