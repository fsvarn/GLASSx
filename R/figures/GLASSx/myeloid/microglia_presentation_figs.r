library(tidyverse)
library(odbc)
library(DBI)
library(gridExtra)
library(ssgsea.GBM.classification)
library(umap)

# TCGA IDHwt vs IDHmut microglia signature

#######################################################
rm(list=ls())
set.seed(11)
##################################################
# Step 1: Subtype each TCGA sample
##################################################

gct_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct"
clin_path <- "/projects/verhaak-lab/USERS/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt"

info <- read.delim(clin_path)

# Run Qianghu's transcriptional classifier (only needed once)
#runSsGSEAwithPermutation(gct_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/p_result_TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct.txt",stringsAsFactor=FALSE)
rownames(subtype_ssgsea) <- gsub("\\.","-",rownames(subtype_ssgsea))

aliquot_barcode <- rep(rownames(subtype_ssgsea),3)
signature_name <- c(rep("Proneural",nrow(subtype_ssgsea)),rep("Classical",nrow(subtype_ssgsea)),rep("Mesenchymal",nrow(subtype_ssgsea)))
enrichment_score <- c(subtype_ssgsea[,"Proneural"],subtype_ssgsea[,"Classical"],subtype_ssgsea[,"Mesenchymal"])
p_value <- c(subtype_ssgsea[,"Proneural_pval"],subtype_ssgsea[,"Classical_pval"],subtype_ssgsea[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(aliquot_barcode,signature_name,enrichment_score,p_value)
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"aliquot_barcode"]),]

transcriptional_subtype[,"signif"] <- transcriptional_subtype[,"p_value"] < 0.05

sig_sub <- transcriptional_subtype %>%
		   group_by(aliquot_barcode, signif) %>%
		   summarise(signature_name = paste(signature_name, collapse=",")) %>%
		   data.frame()
# Pull out significant subtypes first
sig_sub1 <- sig_sub %>%
		    filter(signif)
signif_ali <- sig_sub1[,"aliquot_barcode"]

sig_sub2 <- sig_sub %>%
			filter(!(aliquot_barcode %in% signif_ali) & !signif)
sig_sub <- rbind(sig_sub1, sig_sub2)
  
sig_sub[which(sig_sub[,"signature_name"] == "Proneural,Classical,Mesenchymal"),"signature_name"] <- "Mixed"

sig_sub[,"case_barcode"] <- sapply(strsplit(as.character(sig_sub[,"aliquot_barcode"]),"-"),function(x)paste(x[1:3],collapse="-"))
sig_sub <- sig_sub %>%
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case"))

##################################################
# Step 2: UMAP of each CIBERSORTx GEP
##################################################

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated|t_cell|oligodendrocyte",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_TCGA_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

diff_list <- embedding <- gep_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	colnames(geps) <- paste(mytag[i], colnames(geps), sep="__")
	geps <- log10(geps+1)
	gep_list[[i]] <- t(geps)

	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]


	embedding[[i]] <- umap(t(geps))$layout

}

full_embedding <- do.call(rbind, embedding)

plot_res <- full_embedding
aliquot_barcode <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[2])
cell_state <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[1])
plot_res <- data.frame(aliquot_barcode, plot_res, cell_state)
colnames(plot_res) <- c("aliquot_barcode","UMAP1","UMAP2", "cell_state")
plot_res <- plot_res %>%
	inner_join(sig_sub, by = "aliquot_barcode")
rownames(plot_res) <- paste(plot_res[,"cell_state"], plot_res[,"aliquot_barcode"], sep = "__")

# UMAP plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_facet_umap.pdf",width=7,height=3.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
geom_point()  +
theme_bw() +
facet_wrap(~ cell_state, scales = 'free', nrow = 2) +
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
#scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_umap.pdf",width=4,height=2.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid"), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
geom_point()  +
theme_bw() +
labs(title = "Myeloid") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.text = element_text(size=7),
legend.title = element_text(size=7)) +
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_idhwt_myeloid_umap.pdf",width=4,height=2.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid", IDH.codel.subtype=="IDHwt"), aes(UMAP1, UMAP2, color = signature_name)) + 
geom_point(shape = "square")  +
theme_bw() +
labs(title = "Myeloid") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.text = element_text(size=7),
legend.title = element_text(size=7)) +
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B")) + 
coord_cartesian()
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_umap_legend.pdf",width=4,height=8)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid"), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
geom_point()  +
theme_bw() +
labs(title = "Myeloid") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.text = element_text(size=7),
legend.title = element_text(size=7)) +
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
dev.off()

##################################################
# Step 3: Calcualte macrophage and microglia score
##################################################

# Calculate average score for macrophages and microglia 

# Establish connection to db and get macrophage and microglia score
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
sigs <- sigs %>%
		filter(signature_set == "Muller")

myeloid_gep <- gep_list[[2]]
case_barcode <- sapply(strsplit(sapply(strsplit(rownames(myeloid_gep),"__"),function(x)x[2]),"-"),function(x)paste(x[1:3],collapse="-"))
myeloid_gep <- data.frame(case_barcode, myeloid_gep)

#MDM score
mdm_sig <- sigs %>%
		   filter(signature_name == "Macrophages") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mdm_gep <- myeloid_gep[,c("case_barcode",mdm_sig)]
mdm_rem <- apply(mdm_gep, 2, function(x)sum(is.na(x)))
mdm_gep <- mdm_gep[,-which(mdm_rem == nrow(mdm_gep))]
mdm_score <- apply(mdm_gep[,2:ncol(mdm_gep)], 1, mean)
mdm_score <- data.frame(case_barcode, mdm_score)

#MG score
mg_sig <- sigs %>%
		   filter(signature_name == "Microglia") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mg_gep <- myeloid_gep[,c("case_barcode",mg_sig)]
mg_rem <- apply(mg_gep, 2, function(x)sum(is.na(x)))
mg_gep <- mg_gep[,-which(mg_rem == nrow(mg_gep))]
mg_score <- apply(mg_gep[,2:ncol(mg_gep)], 1, mean)
mg_score <- data.frame(case_barcode, mg_score)

plot_mdm <- plot_res %>%
			filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid") %>%
		    inner_join(mdm_score, by = "case_barcode") %>%
		    rename(myeloid_score = mdm_score) %>%
		    add_column(cell_type = "Macrophages")
plot_mg <- plot_res %>%
			filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid") %>%
		    inner_join(mg_score, by = "case_barcode") %>%
		    rename(myeloid_score = mg_score) %>%
		    add_column(cell_type = "Microglia")
plot_score <- rbind(plot_mdm, plot_mg)

colors <- c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B")
names(colors) <- levels(factor(plot_score$signature_name))

##################################################
# Step 4: Create boxplot showing this more explicitly
##################################################

plot_score <- plot_score %>%
mutate(signature_name = fct_relevel(signature_name, "Proneural", "Proneural,Classical","Mixed","Classical","Proneural,Mesenchymal","Classical,Mesenchymal","Mesenchymal")) %>%
mutate(IDH.status = fct_relevel(IDH.status, "WT", "Mutant"))
reord_colors <- colors[levels(plot_score$signature_name)]


p2 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Microglia"), aes(IDH.status, myeloid_score, fill = IDH.status)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
labs(title = "Microglia") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=c("#619CFF","#00BA38")) +
coord_cartesian(ylim = c(3.07, 3.45))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_mg_idh_status.pdf",width=1.3, height=1.6)
p2
dev.off()

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHwt" & cell_type == "Microglia") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype != "IDHwt" & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)


# IDHmut grade change microglia score

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

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype LIKE 'IDHmut%' AND subtype_a = subtype_b"


dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))

geps <- read.delim(myinf1, row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]	
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))
nrow(geps)

# Read in microglia signature
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
mg_sig <- sigs %>% filter(signature_set == "Muller", signature_name == "Microglia") %>% .$gene_symbol
mdm_sig <- sigs %>% filter(signature_set == "Muller", signature_name == "Macrophages") %>% .$gene_symbol

mg_sig <- intersect(mg_sig, rownames(geps))
mdm_sig <- intersect(mdm_sig, rownames(geps))

mg_score <- apply(geps[mg_sig,], 2, mean)
mdm_score <- apply(geps[mdm_sig,], 2, mean)

up_dat <- dat %>% filter(grade_change == "Grade up")
g1 <- up_dat$tumor_barcode_a
g2 <- up_dat$tumor_barcode_b

wilcox.test(mg_score[g1],mg_score[g2], paired=TRUE)
wilcox.test(mdm_score[g1],mdm_score[g2], paired=TRUE)


st_dat <- dat %>% filter(grade_change == "Grade stable")
g1 <- st_dat$tumor_barcode_a
g2 <- st_dat$tumor_barcode_b

wilcox.test(mg_score[g1],mg_score[g2], paired=TRUE)
wilcox.test(mdm_score[g1],mdm_score[g2], paired=TRUE)

case_barcode <- rep(dat$case_barcode, 4)
aliquot_barcode <- rep(c(dat$tumor_barcode_a, dat$tumor_barcode_b), 2)
grade_change <- rep(dat$grade_change, 4)
idh_codel_subtype <- rep(dat$idh_codel_subtype, 4) 
status <- rep(rep(c("Initial", "Recurrent"), each = nrow(dat)), 2)
score <- c(mg_score[dat$tumor_barcode_a], mg_score[dat$tumor_barcode_b], mdm_score[dat$tumor_barcode_a], mdm_score[dat$tumor_barcode_b])
cell_state <- c(rep("Microglia", nrow(dat)*2), rep("Macrophages", nrow(dat)*2))
plot_dat <- data.frame(case_barcode, aliquot_barcode, status, grade_change, idh_codel_subtype, cell_state, score)

plot_dat <- plot_dat %>%
			mutate(grade_change = fct_relevel(grade_change, "Grade up", "Grade stable"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_glass_idhmut_mdm_score.pdf",width=2.5,height=1.4)
ggplot(plot_dat %>% filter(cell_state == "Macrophages"), aes(x=status, y=score)) + 
geom_boxplot(outlier.shape = NA)  +
geom_line(size=0.5,alpha=0.4, aes(group= case_barcode)) +
geom_point(size=1,aes(colour=idh_codel_subtype)) +
facet_wrap(. ~ grade_change) + 
scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks = c(1.6, 1.7, 1.8, 1.9)) +
theme_bw() +
labs(y = "Mean expression") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7,),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") +
scale_colour_manual(values=c("#F8766D", "#00BA38"))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_glass_idhmut_mg_score_pres.pdf",width=2.5,height=1.65)
ggplot(plot_dat %>% filter(cell_state == "Microglia"), aes(x=status, y=score)) + 
geom_boxplot(outlier.shape = NA)  +
geom_line(size=0.5,alpha=0.4, aes(group= case_barcode)) +
geom_point(size=1,aes(colour=idh_codel_subtype)) +
facet_wrap(. ~ grade_change) + 
scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
theme_bw() +
labs(y = "Mean expression") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7,),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") +
scale_colour_manual(values=c("#00BA38", "#00BA38"))
dev.off()
