###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS/CIBERSORTxHiRes_GLASS_myeloid_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
#geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE cs.idh_codel_subtype = 'IDHwt' AND ts1.signature_name != 'Mesenchymal' AND ts2.signature_name = 'Mesenchymal'"

dat <- dbGetQuery(con,q)

g1 <- geps[,dat[,"tumor_barcode_a"]]
g2 <- geps[,dat[,"tumor_barcode_b"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(i in 1:nrow(geps))
{
	group1 <- as.numeric(g1[i,])
	group2 <- as.numeric(g2[i,])
	
	p.val[i] <- wilcox.test(group1,group2,paired=TRUE)$p.value
	eff[i] <- log2(median(group2/group1))
	#eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
idhwt_res <- data.frame(p.val, q.val, eff)
idhwt_res <- idhwt_res[order(p.val),]

idhwt_res[which(idhwt_res[,"q.val"] < 0.2 & idhwt_res[,"eff"] > 0),]

idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.2

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_myeloid_volcano_nonmes_mes.pdf",width=3,height=2)  
ggplot(idhwt_res, aes(x=eff, y=logp)) + 
geom_point(size=0.5,aes(colour = sig)) +
labs(x = "log2 fold-change", y = "-log10(p-value)") +
scale_colour_manual(values = c("royalblue4", "tomato3")) +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none") 
dev.off()







q <- "SELECT *
FROM analysis.platinum_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE cs.idh_codel_subtype LIKE 'IDHmut%'"

dat <- dbGetQuery(con,q)

g1 <- geps[,dat[,"rna_barcode_a"]]
g2 <- geps[,dat[,"rna_barcode_b"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(i in 1:nrow(geps))
{
	group1 <- as.numeric(g1[i,])
	group2 <- as.numeric(g2[i,])
	
	p.val[i] <- wilcox.test(group1,group2,paired=TRUE)$p.value
	eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
idhmut_res <- data.frame(p.val, q.val, eff)
idhmut_res <- idhmut_res[order(p.val),]
idhmut_res <- idhmut_res

idhmut_res[which(idhmut_res[,"q.val"] < 0.2 & idhmut_res[,"eff"] > 0),]

