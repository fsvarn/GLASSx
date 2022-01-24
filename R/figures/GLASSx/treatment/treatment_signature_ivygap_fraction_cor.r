library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)

#  Post-treatment signature
myinf2 <- paste("data/res/CIBERSORTx/analysis/GLASS_idhwt_", cell_state,"_postreatment_result.txt",sep="")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in IDHwt data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ag1.signature_name AS subtype_a, 
ag2.signature_name AS subtype_b,
ci1.cell_state AS region,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_ivygap ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con, q)


sig_score <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	geps <- log10(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]	
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]
	colnames(geps) <- gsub("\\.","-",colnames(geps))


	mysig <- read.delim(myinf2[i])
	mysig <- rownames(mysig %>% filter(sig, eff > 0))
		
	sig_score[[i]] <- apply(geps[mysig,], 2, mean)
}

aliquot_barcode <- names(sig_score[[1]])
sig_mat <- data.frame(sig_score[[1]], sig_score[[2]], sig_score[[3]])
colnames(sig_mat) <- cell_state

# Make a table with these scores
#write.table(sig_mat, "data/res/CIBERSORTx/analysis/idhwt_signature_scores.txt", sep = "\t", row.names=TRUE, quote=FALSE)

sig_mat_a <- sig_mat[dat[,"tumor_barcode_a"],]
colnames(sig_mat_a) <- paste(colnames(sig_mat_a), "_a", sep="")
dat <- cbind(dat,sig_mat_a)

sig_mat_b <- sig_mat[dat[,"tumor_barcode_b"],]
colnames(sig_mat_b) <- paste(colnames(sig_mat_b), "_b", sep="")
dat <- cbind(dat,sig_mat_b)

rownames(dat) <- NULL

# Initial
res1 <- dat %>%
group_by(region) %>%
summarise(
differentiated_tumor = cor(fraction_a, differentiated_tumor_a),
stemcell_tumor = cor(fraction_a, stemcell_tumor_a),
prolif_stemcell_tumor = cor(fraction_a, prolif_stemcell_tumor_a)) %>%
data.frame()

# Recurrent
res2 <- dat %>%
group_by(region) %>%
summarise(
differentiated_tumor = cor(fraction_b, differentiated_tumor_b),
stemcell_tumor = cor(fraction_b, stemcell_tumor_b),
prolif_stemcell_tumor = cor(fraction_b, prolif_stemcell_tumor_b)) %>%
data.frame()

le_dat <- dat %>% filter(region == "LE")
score <- c(le_dat$differentiated_tumor_a, le_dat$prolif_stemcell_tumor_a, le_dat$stemcell_tumor_a, le_dat$differentiated_tumor_b, le_dat$prolif_stemcell_tumor_b, le_dat$stemcell_tumor_b)
cell_state <- rep(c(rep("Diff.-like", nrow(le_dat)), rep("Prolif. stem-like", nrow(le_dat)), rep("Stem-like", nrow(le_dat))), 2)
fraction <- c(le_dat$fraction_a, le_dat$fraction_a, le_dat$fraction_a, le_dat$fraction_b, le_dat$fraction_b, le_dat$fraction_b)
status <- c(rep("Initial", nrow(le_dat)*3), rep("Recurrent", nrow(le_dat)*3))
case_barcode <- rep(le_dat$case_barcode, 6)
plot_res <- data.frame(case_barcode, status, cell_state, fraction, score) %>%
		    mutate(cell_state = fct_relevel(cell_state, "Diff.-like","Stem-like","Prolif. stem-like"))

res1[,2:4] <- round(res1[,2:4],2)
res2[,2:4] <- round(res2[,2:4],2)

pcor <- c(deparse((bquote(italic("R") ~" = " ~ .(as.character(res1[4,2]))))),
		  deparse((bquote(italic("R") ~" = " ~ .(as.character(res1[4,3]))))),
		  deparse((bquote(italic("R") ~" = " ~ .(as.character(res1[4,4]))))),
		  deparse((bquote(italic("R") ~" = " ~ .(as.character(res2[4,2]))))),
		  deparse((bquote(italic("R") ~" = " ~ .(as.character(res2[4,3]))))),
		  deparse((bquote(italic("R") ~" = " ~ .(as.character(res2[4,4]))))))
		   
annotation_text <- data.frame(status = factor(rep(c("Initial","Recurrent"),each=3)),
							  cell_state = rep(c("Diff.-like","Stem-like","Prolif. stem-like"), 2), 
							  score = 0.8,
							  fraction = .7,
							  pcor)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_leading_edge_sig_cors_v2.pdf", width=3.3,height=2.2)
ggplot(plot_res, aes(x = fraction * 100, y = score)) + 
	geom_point() +
	geom_smooth(method = "lm", se=FALSE, aes(colour=cell_state)) +
	geom_text(data=annotation_text,label=pcor, size=2.5, parse=TRUE) +
	scale_colour_manual(values = c("Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15", "Stem-like" = "#fb6a4a")) + 
	facet_grid(status ~ cell_state) + 
	labs(x="Leading edge fraction (%)", y="Signature score") +
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.text = element_text(size=7),
	axis.title = element_text(size=7),
	strip.background = element_blank(),
	strip.text = element_text(size=7),
	legend.position="none")
dev.off()

stem_res <- plot_res %>% filter(cell_state == "Stem-like")
stem_text = annotation_text %>% filter(cell_state == "Stem-like")
stem_text[,"score"] <- c(1.0, 1.0)
stemcor <- pcor[c(2,5)]

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_leading_edge_sig_cors_stem.pdf", width=1.55,height=2.2)
ggplot(stem_res, aes(x = fraction * 100, y = score)) + 
	geom_point() +
	geom_smooth(method = "lm", se=FALSE, aes(colour=cell_state)) +
	geom_text(data=stem_text,label=stemcor, size=2.5, parse=TRUE) +
	scale_colour_manual(values = c("Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15", "Stem-like" = "#fb6a4a")) + 
	facet_grid(status ~ .) + 
	labs(x="Leading edge fraction (%)", y="Signature score") +
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.text = element_text(size=7),
	axis.title = element_text(size=7),
	strip.background = element_blank(),
	strip.text = element_text(size=7),
	legend.position="none")
dev.off()
