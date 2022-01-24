###################################################
# Compare myeloid cell expression between initial and recurrent IDH-mutant tumors
# Author: Frederick Varn
# Date: 2022.01.07
# Figures S5G
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
myinf1 <- myinf1[grep("myeloid",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

sigs <- p <- se <- list()

# IDHwt analysis
#---------------------------

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype = 'IDHwt' AND received_treatment AND subtype_a = subtype_b"

dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))

geps <- read.delim(myinf1[1], row.names=1)
#geps <- log10(geps+1)		No need because we are taking log2 fold-change later
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]	
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))
nrow(geps)

g1 <- geps[,dat[,"tumor_barcode_a"]]
g2 <- geps[,dat[,"tumor_barcode_b"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(j in 1:nrow(geps))
{
	group1 <- as.numeric(g1[j,])
	group2 <- as.numeric(g2[j,])

	p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
	eff[j] <- log2(mean(group2)/mean(group1))
	#eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
idhwt_res <- data.frame(p.val, q.val, eff)
idhwt_res <- idhwt_res[order(eff),]

idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.05

idhwt_res <- idhwt_res[order(idhwt_res$p.val),]

#myoutf <- paste("data/res/CIBERSORTx/analysis/GLASS_idhwt_", cell_state[i],"_postreatment_result.txt",sep="")

#write.table(idhwt_res, myoutf, sep="\t",quote=FALSE,row.names=TRUE) 

sigs[[1]] <- idhwt_res

# Plot heatmaps
mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
sig_matrix <- geps[mygenes,c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"])]

sig_matrix <- t(apply(sig_matrix, 1, function(x) log10(x) - median(log10(x))))

grp1 <- sig_matrix[,dat[,"tumor_barcode_a"]]
grp1 <- apply(grp1, 1, mean)
grp2 <- sig_matrix[,dat[,"tumor_barcode_b"]]
grp2 <- apply(grp2, 1, mean)

expr <- c(grp1, grp2)
# Rescale for color
#expr[which(expr > quantile(expr,.99))] <- quantile(expr,.99)
expr[which(expr > 0.17)] <- 0.17 # Plug the number in from the whole dataset (obtained below) back in for scaling
gene_symbol <- c(names(grp1), names(grp2))
timepoint <- c(rep("Initial", length(grp1)), rep("Recurrent", length(grp2)))

plot_hm <- data.frame(gene_symbol, expr, timepoint)
lev <- names(grp1)[order(grp2 - grp1)]
plot_hm$gene_symbol <- factor(plot_hm$gene_symbol, levels = lev)
plot_hm[,"cell"] <-"myeloid"

p[[1]] <- plot_hm

se[[1]] <- ggplot(data = plot_hm, aes(x = timepoint, y = gene_symbol)) +
geom_tile(aes(fill=expr)) +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
space = "Lab", na.value = "grey50", guide = "colourbar",
aesthetics = "fill", limits = c(-0.17, 0.17)) + 
theme_void() +
theme(axis.text = element_blank(),
axis.title= element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
legend.position = "none")	


# IDHmut analysis
#---------------------------

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype LIKE 'IDHmut%' AND received_treatment AND subtype_a = subtype_b"

dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))


geps <- read.delim(myinf1, row.names=1)
#geps <- log10(geps+1)		No need because we are taking log2 fold-change later
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]	
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))
nrow(geps)

g1 <- geps[,dat[,"tumor_barcode_a"]]
g2 <- geps[,dat[,"tumor_barcode_b"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(j in 1:nrow(geps))
{
	group1 <- as.numeric(g1[j,])
	group2 <- as.numeric(g2[j,])

	p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
	eff[j] <- log2(mean(group2)/mean(group1))
	#eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
idhmut_res <- data.frame(p.val, q.val, eff)
idhmut_res <- idhmut_res[order(eff),]

idhmut_res[,"logp"] <- -log10(idhmut_res[,"p.val"])
idhmut_res[,"sig"] <- idhmut_res[,"q.val"] < 0.05

idhmut_res <- idhmut_res[order(idhmut_res$p.val),]

#myoutf <- paste("data/res/CIBERSORTx/analysis/GLASS_idhmut_", cell_state[i],"_postreatment_result.txt",sep="")

#write.table(idhmut_res, myoutf, sep="\t",quote=FALSE,row.names=TRUE) 

sigs[[2]] <- idhmut_res

# Plot heatmaps
mygenes <- rownames(idhmut_res[which(idhmut_res[,"sig"]),])
sig_matrix <- geps[mygenes,c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"])]

sig_matrix <- t(apply(sig_matrix, 1, function(x)log10(x) - median(log10(x))))

grp1 <- sig_matrix[,dat[,"tumor_barcode_a"]]
grp1 <- apply(grp1, 1, mean)
grp2 <- sig_matrix[,dat[,"tumor_barcode_b"]]
grp2 <- apply(grp2, 1, mean)

expr <- c(grp1, grp2)
# Rescale for color
#expr[which(expr > quantile(expr,.99))] <- quantile(expr,.99)
expr[which(expr > 0.17)] <- 0.17 # Plug the number in from the whole dataset (obtained below) back in for scaling
gene_symbol <- c(names(grp1), names(grp2))
timepoint <- c(rep("Initial", length(grp1)), rep("Recurrent", length(grp2)))

plot_hm <- data.frame(gene_symbol, expr, timepoint)
lev <- names(grp1)[order(grp2 - grp1)]
plot_hm$gene_symbol <- factor(plot_hm$gene_symbol, levels = lev)
plot_hm[,"cell"] <- "myeloid"

p[[2]] <- plot_hm

se[[2]] <- ggplot(data = plot_hm %>% filter(), aes(x = timepoint, y = gene_symbol)) +
geom_tile(aes(fill=expr)) +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
space = "Lab", na.value = "grey50", guide = "colourbar",
aesthetics = "fill", limits = c(-0.17, 0.17)) + 
theme_void() +
theme(axis.text = element_blank(),
axis.title= element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
legend.position = "none")	

plot_res <- do.call(rbind, p)
# Identify the 99th percentile globally
quantile(plot_res$expr,.99) #0.17

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/treatment_pre_post_myeloid_csx_heatmaps_scale.pdf",width=7, height =3)
grid.arrange(se[[1]],se[[2]],nrow=1)
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/treatment_pre_post_myeloid_csx_legends.pdf",width=7, height =3)
ggplot(data = plot_res, aes(x = timepoint, y = gene_symbol)) +
	geom_tile(aes(fill=expr)) +
	scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
	space = "Lab", na.value = "grey50", guide = "colourbar",
	aesthetics = "fill", limits = c(-0.17, 0.17)) + 
	theme_void() +
	theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank())	
dev.off()