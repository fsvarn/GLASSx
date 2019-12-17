library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.shannon AS tcr_shannon_a, 
im2.shannon AS tcr_shannon_b, 
im1.evenness AS tcr_evenness_a,
im2.evenness AS tcr_evenness_b,
im1.richness AS tcr_richness_a,
im2.richness AS tcr_richness_b,
im1.total_tcr AS total_tcr_a,
im2.total_tcr AS total_tcr_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b
FROM analysis.rna_silver_set ps
JOIN analysis.tcr_stats im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.tcr_stats im2 ON im2.aliquot_barcode = ps.tumor_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
WHERE im1.total_tcr > 0 AND im2.total_tcr > 0 --compare samples with TCR expression
ORDER BY 1, 2, 6"

dat <- dbGetQuery(con, q)

wilcox.test(dat[,"total_tcr_a"], dat[,"total_tcr_b"],paired=TRUE)
wilcox.test(dat[,"tcr_evenness_a"], dat[,"tcr_evenness_b"],paired=TRUE)
wilcox.test(dat[,"tcr_richness_a"], dat[,"tcr_richness_b"],paired=TRUE)
wilcox.test(dat[,"tcr_shannon_a"], dat[,"tcr_shannon_b"],paired=TRUE)

aliquot_barcode <- c(dat[,"tumor_barcode_a"],dat[,"tumor_barcode_b"])
sample_type <- c(dat[,"sample_type_a"],dat[,"sample_type_b"])
total_tcr <- c(dat[,"total_tcr_a"],dat[,"total_tcr_b"])
richness <- c(dat[,"tcr_richness_a"],dat[,"tcr_richness_b"])
shannon <- c(dat[,"tcr_shannon_a"],dat[,"tcr_shannon_b"])
evenness <- c(dat[,"tcr_evenness_a"],dat[,"tcr_evenness_b"])

plot_dat <- data.frame(aliquot_barcode,sample_type,total_tcr,richness,shannon,evenness)

#Boxplots: CURRENTLY UNFINISHED

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/tcr_boxplots_initial_recurrent.pdf",width=2,height=2)
ggplot(plot_dat,aes(x = sample_type, y = total_tcr,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "log TCR abundance (initial)", y = "log TCR abundance (recurrent)", title="TCR abundance over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0,log(300)),ylim=c(0,log(300)))

ggplot(plot_dat,aes(x = sample_type, y = richness,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "TCR Shannon entropy (initial)", y = "TCR Shannon entropy (recurrent)", title="TCR Shannon entropy over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0,6),ylim=c(0,6))

ggplot(plot_dat,aes(x = sample_type, y = shannon,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "log TCR richness (initial)", y = "log TCR richness (recurrent)", title="TCR richness over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0,log(250)),ylim=c(0,log(250)))

ggplot(plot_dat,aes(x = sample_type, y = evenness, colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "TCR evenness (initial)", y = "TCR evenness (recurrent)", title="TCR evenness over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0.7,1),ylim=c(0.7,1))
dev.off()




#######################################################

#Read in data for plotting (logs)
q <- "SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.shannon AS tcr_shannon_a, 
im2.shannon AS tcr_shannon_b, 
im1.evenness AS tcr_evenness_a,
im2.evenness AS tcr_evenness_b,
COALESCE(ln(NULLIF(im1.richness,0)),0) AS tcr_richness_a,
COALESCE(ln(NULLIF(im2.richness,0)),0) AS tcr_richness_b,
COALESCE(ln(NULLIF(im1.total_tcr,0)),0) AS total_tcr_a,
COALESCE(ln(NULLIF(im2.total_tcr,0)),0) AS total_tcr_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b
FROM analysis.rna_silver_set ps
JOIN analysis.tcr_stats im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.tcr_stats im2 ON im2.aliquot_barcode = ps.tumor_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
ORDER BY 1, 2, 6"

dat <- dbGetQuery(con, q)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/tcr_scatterplots_initial_recurrent.pdf",width=2,height=2)
ggplot(dat,aes(x = total_tcr_a, y = total_tcr_b,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "log TCR abundance (initial)", y = "log TCR abundance (recurrent)", title="TCR abundance over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0,log(300)),ylim=c(0,log(300)))

ggplot(dat,aes(x = tcr_shannon_a, y = tcr_shannon_b,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "TCR Shannon entropy (initial)", y = "TCR Shannon entropy (recurrent)", title="TCR Shannon entropy over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0,6),ylim=c(0,6))

ggplot(dat,aes(x = tcr_richness_a, y = tcr_richness_b,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "log TCR richness (initial)", y = "log TCR richness (recurrent)", title="TCR richness over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0,log(250)),ylim=c(0,log(250)))

ggplot(dat,aes(x = tcr_evenness_a, y = tcr_evenness_b,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
labs(x = "TCR evenness (initial)", y = "TCR evenness (recurrent)", title="TCR evenness over time") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0.7,1),ylim=c(0.7,1))
dev.off()


#######################################################
#Stratifying by subtypes

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.shannon AS tcr_shannon_a, 
im2.shannon AS tcr_shannon_b, 
im1.evenness AS tcr_evenness_a,
im2.evenness AS tcr_evenness_b,
im1.richness AS tcr_richness_a,
im2.richness AS tcr_richness_b,
im1.total_tcr AS total_tcr_a,
im2.total_tcr AS total_tcr_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b
FROM analysis.rna_silver_set ps
JOIN analysis.tcr_stats im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.tcr_stats im2 ON im2.aliquot_barcode = ps.tumor_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
ORDER BY 1, 2, 6"

dat <- dbGetQuery(con, q)

subtype <- unique(dat[,"subtype_a"])

total.tcr.pval <-  richness.pval <- shannon.pval <- evenness.pval <- total.tcr.eff <-  richness.eff <- shannon.eff <- evenness.eff <- rep(0,length(subtype))
for(i in 1:length(subtype))
{
	sub.dat <- dat[which(dat[,"subtype_a"]==subtype[i]),]
	
	total.tcr.pval[i] <- wilcox.test(sub.dat[,"total_tcr_b"],sub.dat[,"total_tcr_a"],paired=TRUE)$p.value
	total.tcr.eff[i] <- median(sub.dat[,"total_tcr_b"]) - median(sub.dat[,"total_tcr_a"])
	
	richness.pval[i] <- wilcox.test(sub.dat[,"tcr_richness_b"],sub.dat[,"tcr_richness_a"],paired=TRUE)$p.value
	richness.eff[i] <- median(sub.dat[,"tcr_richness_b"]) - median(sub.dat[,"tcr_richness_a"])
	
	shannon.pval[i] <- wilcox.test(sub.dat[,"tcr_shannon_b"],sub.dat[,"tcr_shannon_a"],paired=TRUE)$p.value
	shannon.eff[i] <- median(sub.dat[,"tcr_shannon_b"]) - median(sub.dat[,"tcr_shannon_a"])
	
	evenness.pval[i] <- wilcox.test(sub.dat[,"tcr_evenness_b"],sub.dat[,"tcr_evenness_a"],paired=TRUE)$p.value
	evenness.eff[i] <- median(sub.dat[,"tcr_evenness_b"]) - median(sub.dat[,"tcr_evenness_a"])
	
}

data.frame(subtype, total.tcr.pval, total.tcr.eff, richness.pval, richness.eff, shannon.pval, shannon.eff, evenness.pval, evenness.eff)

#           subtype total.tcr.pval total.tcr.eff richness.pval richness.eff
# 1           IDHwt      0.2069364           1.0     0.1761123         -0.5
# 2 IDHmut-noncodel      0.6134133           1.0     0.5469475          1.0
# 3    IDHmut-codel      0.2692941          -0.5     0.2692941          0.0
#   shannon.pval shannon.eff evenness.pval evenness.eff
# 1    0.5238074 -0.07557516     0.1893492           NA
# 2    0.7668982  0.00000000     0.2500000           NA
# 3    0.7892680 -0.02831651     1.0000000           NA