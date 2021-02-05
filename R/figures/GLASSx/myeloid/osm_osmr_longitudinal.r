
library(tidyverse)
library(odbc)
library(DBI)
library(grid)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name, ts1.enrichment_score AS score_a, ts2.enrichment_score AS score_b, ts2.enrichment_score - ts1.enrichment_score AS change
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE cs.idh_codel_subtype = 'IDHwt' AND ts1.signature_name ='Mesenchymal' AND ts2.signature_name = 'Mesenchymal'
ORDER BY 1"

dat <- dbGetQuery(con,q)


geps <- read.delim(myinf1[2], row.names=1)
#geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))



g1 <- as.numeric(geps["OSM",dat[,"tumor_barcode_a"]])
g2 <- as.numeric(geps["OSM",dat[,"tumor_barcode_b"]])
gene_diff <- g2 - g1

cor(gene_diff, dat$change, method="p")

tpm_change <- gene_diff
mes_change <- dat$change
gene <- rep("Myeloid: OSM", length(gene_diff))

geps <- read.delim(myinf1[1], row.names=1)
#geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))


g1 <- as.numeric(geps["OSMR",dat[,"tumor_barcode_a"]])
g2 <- as.numeric(geps["OSMR",dat[,"tumor_barcode_b"]])
gene_diff <- g2 - g1

cor(gene_diff, dat$change, method="p")

tpm_change <- c(tpm_change, gene_diff)
mes_change <- c(mes_change, dat$change)
gene <- c(gene, rep("Diff.-like: OSMR", length(gene_diff)))

res <- data.frame(tpm_change, mes_change, gene)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/glass_osm_osmr_bulk_mesenchymal_change.pdf",width=2.5,height=1.5)
ggplot(res, aes(tpm_change, mes_change/1000)) + 
geom_point()  +
geom_smooth(method="lm", se = FALSE, fullrange=TRUE) +
facet_wrap(vars(gene),scales="free_x",nrow=1) +
theme_bw() +
labs(x="Change in inferred expression", y = "Change in mes. score") +
#scale_x_continuous(breaks = breaks_fun, limits = c(0, NA)) +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") 
dev.off()
