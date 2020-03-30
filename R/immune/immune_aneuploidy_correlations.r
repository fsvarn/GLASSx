library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)

#######################################################
rm(list=ls())
#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
SELECT ps.*, an1.prop_aneuploidy AS prop_an1, an2.prop_aneuploidy AS prop_an2, (an2.prop_aneuploidy-an1.prop_aneuploidy) AS prop_dif, an1.aneuploidy_score AS an_score1, an2.aneuploidy_score AS an_score2, (an2.aneuploidy_score - an1.aneuploidy_score) AS score_dif, 
im1.signature_name, im1.enrichment_score AS es_a, im2.enrichment_score AS es_b, (im2.enrichment_score - im1.enrichment_score) AS es_dif, cs.idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im1.signature_name = im2.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
--WHERE cs.idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)


cells <- c("CD4.mature","CD8.effector","NK.cells","B.cells","Dendritic","Macrophages")

score_a <- score_b <- prop_a <- prop_b <- score_dif <- prop_dif <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	score_a[i] <- cor(sub_dat[,"es_a"],sub_dat[,"an_score1"],method="s")
	score_b[i] <- cor(sub_dat[,"es_b"],sub_dat[,"an_score2"],method="s")
	
	prop_a[i] <- cor(sub_dat[,"es_a"],sub_dat[,"prop_an1"],method="s")
	prop_b[i] <- cor(sub_dat[,"es_b"],sub_dat[,"prop_an2"],method="s")
	
	score_dif[i] <- cor(sub_dat[,"score_dif"],sub_dat[,"es_dif"],method="s")
	prop_dif[i] <- cor(sub_dat[,"prop_dif"],sub_dat[,"es_dif"],method="s")
}

res <- data.frame(cells,score_a,score_b,prop_a,prop_b,score_dif,prop_dif)

sub_dat <- dat %>% filter(signature_name == "CD8.effector")
sub_dat[,"score_dif"] <- as.numeric(sub_dat[,"score_dif"])
cor1 <- round(cor(sub_dat[,"prop_dif"],sub_dat[,"es_dif"],method="p"),2)

se1 <- ggplot(data = sub_dat, aes(x = prop_dif, y = es_dif, colour= idh_codel_subtype)) +
geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
geom_point(size=1) +
labs(x = "Change in chromosomal instability", y = "Immune change", title = "CD8+ T cells") +
annotate(geom="text", size=2.5, 
	x=max(sub_dat[,"prop_dif"]) - (0.80*(max(sub_dat[,"prop_dif"]) - min(sub_dat[,"prop_dif"]))), 
	y=max(sub_dat[,"es_dif"]) - (0.97*(max(sub_dat[,"es_dif"]) - min(sub_dat[,"es_dif"]))), 
	label=deparse(bquote(italic("R")~" = "~.(cor1))),parse=TRUE) +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")

sub_dat <- dat %>% filter(signature_name == "Macrophages")
sub_dat[,"score_dif"] <- as.numeric(sub_dat[,"score_dif"])
cor2 <- round(cor(sub_dat[,"prop_dif"],sub_dat[,"es_dif"],method="p"),2)

se2 <- ggplot(data = sub_dat, aes(x = prop_dif, y = es_dif, colour= idh_codel_subtype)) +
geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
geom_point(size=1) +
labs(x = "Change in chromosomal instability", y = "Immune change", title = "Macrophages") +
annotate(geom="text", size=2.5, 
	x=max(sub_dat[,"prop_dif"]) - (0.80*(max(sub_dat[,"prop_dif"]) - min(sub_dat[,"prop_dif"]))), 
	y=max(sub_dat[,"es_dif"]) - (0.97*(max(sub_dat[,"es_dif"]) - min(sub_dat[,"es_dif"]))), 
	label=deparse(bquote(italic("R")~" = "~.(cor2))),parse=TRUE) +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")
	
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/copy_number_immune_change_scatter.pdf",width=2,height=4)
grid.arrange(se1,se2,nrow=2)
dev.off()

# 
# dat[,"aneu_status"] <- as.factor(dat[,"prop_dif"] >= 0)
# 
# sub_dat <- dat %>% filter(signature_name == "CD8.effector") %>%
# 	mutate(aneu_status=recode(aneu_status, "FALSE" = "Loss", "TRUE" = "Gain"))
# 	
# se1 <- ggplot(data = sub_dat, aes(x = aneu_status, y = es_dif, fill= idh_codel_subtype)) +
# geom_boxplot() +
# labs(y = "Proportion of genome altered",title="CD8+ T cells") +
# theme_bw() +
# theme(axis.text = element_text(size=7),
# 	axis.title.x= element_blank(),
# 	axis.title.y= element_text(size=7),
# 	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
# 	strip.text = element_text(size=7),
# 	strip.background = element_rect(colour="white",fill="white"),
# 	legend.position = "none")
# 	
# sub_dat <- dat %>% filter(signature_name == "Macrophages") %>%
# 	mutate(aneu_status=recode(aneu_status, "FALSE" = "Loss", "TRUE" = "Gain"))
# 	
# se2 <- ggplot(data = sub_dat, aes(x = aneu_status, y = es_dif, fill= idh_codel_subtype)) +
# geom_boxplot() +
# labs(y = "Proportion of genome altered",title="Macrophages") +
# theme_bw() +
# theme(axis.text = element_text(size=7),
# 	axis.title.x= element_blank(),
# 	axis.title.y= element_text(size=7),
# 	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
# 	strip.text = element_text(size=7),
# 	strip.background = element_rect(colour="white",fill="white"),
# 	legend.position = "none")
# 		
# pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/aneuploidy_immune_change_box.pdf",width=3,height=2)
# grid.arrange(se1,se2,nrow=1)
# dev.off()
# 
