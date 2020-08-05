###################################################
# Compare TPM for RP11-139K4.2 lncRNA between IDHmut and IDHwt samples
# Collaboration with Sunit Das group
# Updated: 2020.06.09
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)


#######################################################
rm(list=ls())

#Establish connection
myinf1 <- "results/kallisto/noncoding/final/noncoding_tpm_matrix_all_samples.tsv"

# Read in expression matrix
dat <- read.delim(myinf1,row.names=1)
colnames(dat) <- gsub("\\.","-",colnames(dat))

# Subset to publicly available GLASS data:
dat <- dat[,which(colnames(dat) )]
dat <- dat[,c("GLSS-SM-R060-TP-01R-RNA-ZI9YVM", "GLSS-SM-R064-TP-01R-RNA-W6PIW2", "GLSS-SM-R060-R1-01R-RNA-6SGM43", "GLSS-SM-R056-R2-01R-RNA-DHLJ2T", "GLSS-SM-R064-R1-01R-RNA-GLLBN3", "GLSS-HF-2548-TP-01R-RNA-4R6KR6", "GLSS-HF-2829-TP-01R-RNA-90I02Y", "GLSS-HF-2869-TP-01R-RNA-2H6MOF", "GLSS-HF-2548-R1-01R-RNA-F7XJ9B", "GLSS-HF-2919-TP-01R-RNA-3AV46F", "GLSS-HF-2934-TP-01R-RNA-WUQKHB", "GLSS-HF-2829-R1-01R-RNA-D71HI8", "GLSS-HF-2869-R3-01R-RNA-TF6NT8", "GLSS-HF-2998-TP-01R-RNA-5H7ODC", "GLSS-HF-2998-R1-01R-RNA-2HE04E", "GLSS-HF-2934-R1-01R-RNA-IBE7TN", "GLSS-HF-3050-TP-01R-RNA-DR0HT8", "GLSS-HF-2919-R1-01R-RNA-R6L8SO", "GLSS-HF-3081-TP-01R-RNA-EE79OQ", "GLSS-HF-3162-TP-01R-RNA-K5P3FT", "GLSS-HF-3081-R2-01R-RNA-939287", "GLSS-HF-3162-R1-01R-RNA-JUILFB", "GLSS-HF-3050-R1-01R-RNA-7FET5L", "GLSS-SM-R072-R1-01R-RNA-1DQIND", "GLSS-SM-R099-TP-01R-RNA-YGXA72", "GLSS-SM-R100-TP-01R-RNA-EUI7AZ", "GLSS-SM-R106-TP-01R-RNA-P5JI35", "GLSS-SM-R104-TP-01R-RNA-G52UET", "GLSS-SM-R109-TP-01R-RNA-I4HTPP", "GLSS-SM-R099-R1-01R-RNA-MNTPMI", "GLSS-SM-R100-R1-01R-RNA-46UW5U", "GLSS-SM-R101-TP-02R-RNA-V4TRVD", "GLSS-SM-R106-R1-01R-RNA-MFHIBQ", "GLSS-SM-R107-TP-01R-RNA-ZSYP37", "GLSS-SM-R108-TP-01R-RNA-KUU8QP", "GLSS-SM-R104-R1-01R-RNA-YF3DT7", "GLSS-SM-R109-R1-01R-RNA-R5I1C7", "GLSS-SM-R110-TP-02R-RNA-A85T37", "GLSS-SM-R101-R1-02R-RNA-BIXXBZ", "GLSS-SM-R110-R1-01R-RNA-F7JINL", "GLSS-SM-R108-R1-01R-RNA-NICL65", "GLSS-SM-R107-R1-02R-RNA-1K9XUJ", "GLSS-SM-R103-TP-01R-RNA-7YEFHB", "GLSS-SM-R103-R1-01R-RNA-TSX1I7", "GLSS-SM-R061-TP-01R-RNA-AZSPBW", "GLSS-SM-R061-R1-01R-RNA-0NG539", "GLSS-SM-R095-TP-01R-RNA-5T2HT4", "GLSS-SM-R095-R1-01R-RNA-BQ7VMI", "GLSS-SM-R072-TP-01R-RNA-K18OXT", "GLSS-SM-R056-TP-01R-RNA-JWKQQG", "GLSS-SM-R067-TP-01R-RNA-JYXPKW", "GLSS-SM-R067-R1-01R-RNA-O82WJY", "TCGA-06-0125-TP-01R-RNA-RCQ5QS", "TCGA-06-0125-R1-11R-RNA-EKH06U", "TCGA-06-0190-TP-01R-RNA-S2074S", "TCGA-06-0190-R1-01R-RNA-HWTLG4", "TCGA-06-0210-TP-01R-RNA-JMBENW", "TCGA-06-0210-R1-01R-RNA-WB10MK", "TCGA-06-0211-TP-01R-RNA-Y9DKHR", "TCGA-06-0211-R1-02R-RNA-3YQ1G4", "TCGA-14-1034-TP-01R-RNA-S6SBZH", "TCGA-14-1034-R1-01R-RNA-GSNR8T", "TCGA-DH-A669-TP-12R-RNA-XG7CTS", "TCGA-DH-A669-R1-11R-RNA-KJQZGQ", "TCGA-DU-5870-TP-11R-RNA-CSQWQZ", "TCGA-DU-5870-R1-12R-RNA-GIZUWN", "TCGA-DU-5872-TP-11R-RNA-HS9DNI", "TCGA-DU-5872-R1-21R-RNA-0EW2ZH", "TCGA-DU-6397-TP-11R-RNA-NYNT65", "TCGA-DU-6397-R1-12R-RNA-N760C1", "TCGA-DU-6404-TP-11R-RNA-A7XMXJ", "TCGA-DU-6404-R1-21R-RNA-LJZDP3", "TCGA-DU-6407-TP-13R-RNA-6NHIOE", "TCGA-DU-6407-R1-12R-RNA-T3UMCW", "TCGA-DU-7304-TP-12R-RNA-23D1EM", "TCGA-DU-7304-R1-12R-RNA-OMXRQS", "TCGA-FG-5963-TP-11R-RNA-2Y6I4P", "TCGA-FG-5963-R1-12R-RNA-1EJ1SV", "TCGA-FG-5965-TP-11R-RNA-0XBPNY", "TCGA-FG-5965-R1-11R-RNA-QSA7Y3", "TCGA-FG-A4MT-TP-11R-RNA-K4BQG0", "TCGA-FG-A4MT-R1-11R-RNA-9OU9UC", "TCGA-TM-A7CF-TP-11R-RNA-YLP5PU", "TCGA-TM-A7CF-R1-11R-RNA-EZWWA2", "TCGA-TQ-A7RK-TP-11R-RNA-LRGF1D", "TCGA-TQ-A7RK-R1-11R-RNA-PKFXB3", "TCGA-TQ-A7RV-TP-21R-RNA-6F23HF", "TCGA-TQ-A7RV-R1-11R-RNA-TIRW8I", "TCGA-TQ-A8XE-TP-11R-RNA-0YPLAZ", "TCGA-TQ-A8XE-R1-11R-RNA-47VZAQ", "GLSS-SM-R071-TP-01R-RNA-OAXGI8", "GLSS-SM-R071-R1-01R-RNA-7AZ6G2", "GLSS-SM-R068-TP-01R-RNA-0UPMYO", "GLSS-SM-R068-R1-01R-RNA-7I5H9P", "GLSS-SM-R065-TP-01R-RNA-OM7S2C", "GLSS-SM-R065-R1-01R-RNA-V35ETA", "GLSS-SM-R070-TP-01R-RNA-QWG60O", "GLSS-SM-R070-R1-01R-RNA-PXOL5Z", "GLSS-SM-R063-TP-01R-RNA-WIVUSX", "GLSS-SM-R063-R1-01R-RNA-YQOR0G", "GLSS-SM-R066-TP-01R-RNA-G9A6CK", "GLSS-SM-R066-R1-01R-RNA-MKHTJ4")]
dat <- log10(dat+1)
#pub_inds <- c(grep("GLSS-HF",colnames(dat)), grep("GLSS-SM",colnames(dat)), grep("TCGA-",colnames(dat)))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in clinical data
q <- "
SELECT ps.*, cs.idh_codel_subtype
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

info <- dbGetQuery(con, q)
info <- info[which(info[,"tumor_barcode_a"] %in% colnames(dat) & info[,"tumor_barcode_b"] %in% colnames(dat)),]

# Test in initial tumors only:
g1 <- as.numeric(dat["RP11-139K4.2",info[which(info[,"idh_codel_subtype"]=="IDHwt"),"tumor_barcode_a"]])
g2 <- as.numeric(dat["RP11-139K4.2",info[which(info[,"idh_codel_subtype"]!="IDHwt"),"tumor_barcode_a"]])
init_pval <- formatC(wilcox.test(g1,g2)$p.value, format = "e", digits = 0)


g1 <- as.numeric(dat["RP11-139K4.2",info[which(info[,"idh_codel_subtype"]=="IDHwt"),"tumor_barcode_b"]])
g2 <- as.numeric(dat["RP11-139K4.2",info[which(info[,"idh_codel_subtype"]!="IDHwt"),"tumor_barcode_b"]])
rec_pval <- formatC(wilcox.test(g1,g2)$p.value, format = "e", digits = 0)

# Make a plot
aliquot_barcode <- c(info[,"tumor_barcode_a"], info[,"tumor_barcode_b"])
idh_codel_subtype <- rep(info[,"idh_codel_subtype"],2)
idhmut <- as.character(idh_codel_subtype != "IDHwt")  %>% recode("TRUE" = "IDHmut","FALSE" = "IDHwt")
timepoint <- c(rep("Initial",nrow(info)),rep("Recurrent",nrow(info)))
tpm <- t(dat["RP11-139K4.2",aliquot_barcode])
plot_dat <- data.frame(aliquot_barcode, tpm, timepoint, idh_codel_subtype,idhmut)
rownames(plot_dat) <- NULL
colnames(plot_dat) <- c("aliquot_barcode","tpm","timepoint","idh_codel_subtype","idhmut")

p_val <- c(deparse((bquote(italic("P") ~" = " ~ .(init_pval)))),
		   deparse((bquote(italic("P") ~" = " ~ .(rec_pval)))))
annotation_text <- data.frame(timepoint = factor(c("Initial","Recurrent")),
							  idhmut = 1.1,
							  tpm = 0.95,
							  p_val)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/RP11-139K4.2_idhmut_boxplot.pdf",width=4,height=2)
p1 <- ggplot(plot_dat, aes(x = idhmut, y = tpm)) +
geom_boxplot(outlier.shape=NA) +
geom_jitter(size=0.4, aes(colour = idh_codel_subtype)) +
facet_grid(.~timepoint) +
labs(main = "RP11-139K4.2", y = "log10(TPM)") +
geom_text(data=annotation_text,label=p_val, size=2.5, parse=TRUE) +
scale_fill_manual(values=c("white","white")) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7)) +
coord_cartesian(ylim=c(0,1))
p1
dev.off()



