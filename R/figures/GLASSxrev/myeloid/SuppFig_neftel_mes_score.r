###################################################
# Examine MES-like signature expression of neoplastic cells undergoing a mesenchymal transition
# Author: Frederick Varn
# Date: 2022.01.06
# Figure S5C
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)
library(survival)
library(topGO)

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

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype = 'IDHwt' AND (subtype_a != 'Mesenchymal' AND subtype_b = 'Mesenchymal')"

dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))

# Read in mesenchymal signature
q <- "SELECT * FROM ref.neftel_sig WHERE cell_state = 'MES'"
sig <- dbGetQuery(con, q)
full_mes <- unique(sig[,1])

state <- c("Diff.-like", "Prolif. stem-like", "Stem-like")
res <- list()
p.val <- rep(0, length(myinf1))
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
	nrow(geps)

	inter_mes <- intersect(rownames(geps), full_mes)

	g1 <- geps[inter_mes,dat[,"tumor_barcode_a"]]
	g2 <- geps[inter_mes,dat[,"tumor_barcode_b"]]

	g1 <- apply(g1,2,mean)
	g2 <- apply(g2,2,mean)
	
	#median normalize
	med <- median(c(g1,g2))
	g1 <- g1 - med
	g2 <- g2 - med
	
	score <- c(g1, g2)
	aliquot_barcode <- c(names(g1), names(g2))
	case_barcode <- substring(aliquot_barcode, 1,12)
	timepoint <- c(rep("Initial", length(g1)), rep("Recurrent", length(g2)))
	sub_res <- data.frame(aliquot_barcode, case_barcode, timepoint, score)
	sub_res[,"state"] <- state[i]
	res[[i]] <- sub_res

	p.val[i] <- wilcox.test(g1,g2,paired=TRUE)$p.value
}

plot_res <- do.call(rbind, res)

plot_res <- plot_res %>% filter(state != "Prolif. stem-like")

# Create subtype vector for plotting
subtype <- c(dat$subtype_a, dat$subtype_b)
names(subtype) <- c(dat$tumor_barcode_a, dat$tumor_barcode_b)

plot_res$subtype <- subtype[as.character(plot_res$aliquot_barcode)]

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_MES-like_neftel_score.pdf",width=2,height=1.5)
ggplot(plot_res, aes(x=timepoint, y=score)) + 
geom_boxplot(outlier.shape = NA)  +
geom_line(size=0.5,alpha=0.4, aes(group= case_barcode)) +
geom_point(size=1,aes(colour=subtype)) +
facet_grid(.~state) +
theme_classic() +
labs(y = "MES-like signature score") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position = "none") +
scale_colour_manual(values=c("#008A22", "#8A0000","#00458A")) +
coord_cartesian(ylim = c(-0.32, 0.32))
dev.off()