library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q1 <- "
WITH subtype_label AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	ts1.signature_name, 
	ts1.p_value AS pval_a,
	ts2.p_value AS pval_b,
	CASE WHEN ts1.p_value < 0.05 THEN ts1.signature_name ELSE NULL END AS subtype_a,
	CASE WHEN ts2.p_value < 0.05 THEN ts2.signature_name ELSE NULL END AS subtype_b
	FROM analysis.rna_silver_set ss
	JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
	WHERE ts1.p_value < 0.05 OR ts2.p_value < 0.05
),
subtype_switch AS
(
	SELECT sl.tumor_pair_barcode,
	sl.case_barcode,
	string_agg(subtype_a,',') AS subtype_a,
	string_agg(subtype_b,',') AS subtype_b
	FROM subtype_label sl
	GROUP BY tumor_pair_barcode, sl.case_barcode
	ORDER BY 1
)
SELECT ss.*, 
si1.simplicity_score AS ss_a,
si2.simplicity_score AS ss_b,
CASE WHEN sw.subtype_a IS NULL THEN 'Mixed' ELSE sw.subtype_a END AS subtype_a,
CASE WHEN sw.subtype_b IS NULL THEN 'Mixed' ELSE sw.subtype_b END AS subtype_b,
CASE WHEN sw.subtype_a = sw.subtype_b THEN 'None' ELSE 'Switch' END AS switch,
su.idh_codel_subtype,
cn.cnv_driver_change_a,
cn.cnv_driver_change_b,
sn.snv_driver_change_a,
sn.snv_driver_change_b
FROM analysis.rna_silver_set ss
JOIN analysis.rna_dna_pairs rd ON ss.tumor_barcode_a = rd.rna_barcode_a AND ss.tumor_barcode_b = rd.rna_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.surgeries su ON su.sample_barcode = al1.sample_barcode
JOIN analysis.simplicity_score si1 ON ss.tumor_barcode_a = si1.aliquot_barcode
JOIN analysis.simplicity_score si2 ON ss.tumor_barcode_b = si2.aliquot_barcode
JOIN subtype_switch sw ON ss.tumor_pair_barcode = sw.tumor_pair_barcode
JOIN analysis.driver_status_cnv cn ON cn.tumor_pair_barcode = rd.dna_pair_barcode
JOIN analysis.driver_status_snv sn ON sn.tumor_pair_barcode = rd.dna_pair_barcode
WHERE su.idh_codel_subtype IS NOT NULL -- ts1.p_value < 0.05 OR ts2.p_value < 0.05
ORDER BY 1, 2
"

dat <- dbGetQuery(con,q1)

genes_cnv <- unlist(strsplit(c(dat[,"cnv_driver_change_a"],dat[,"cnv_driver_change_b"]),","))
genes_cnv <- gsub("-","",genes_cnv)
genes_cnv <- gsub("\\+","",genes_cnv)
#genes_cnv <- gsub("amp","",genes_cnv)
#genes_cnv <- gsub("del","",genes_cnv)
#genes_cnv <- gsub(" ","",genes_cnv)
genes_cnv <- trimws(genes_cnv,which="both")
genes_cnv <- genes_cnv[-which(is.na(genes_cnv))]
genes_cnv <- unique(genes_cnv)

genes_snv <- unlist(strsplit(c(dat[,"snv_driver_change_a"],dat[,"snv_driver_change_b"]),","))
genes_snv <- gsub("p..*","",genes_snv)
genes_snv <- gsub("-","",genes_snv)
genes_snv <- gsub("\\+","",genes_snv)
#genes_snv <- gsub(" ","",genes_snv)
genes_snv <- genes_snv[-which(is.na(genes_snv))]
genes_snv <- trimws(genes_snv,which="both")
genes_snv <- unique(genes_snv)

mysubtypes <- c("Proneural","Classical","Mesenchymal")

cn_pvals <- sn_pvals <- list()
for(i in 1:length(mysubtypes))
{
	#Check for alterations associated with subtype gains and losses
	myind_a <- grep(mysubtypes[i],dat[,"subtype_a"])
	myind_b <- grep(mysubtypes[i],dat[,"subtype_b"])
	
	#Check for alterations associated with subtype gains and losses	
	sub_gain <- dat[ myind_b[-which(myind_b %in% myind_a)],]
	sub_loss <- dat[ myind_a[-which(myind_a %in% myind_b)],]
	
	#Copy number
	gain_amp_pval <- loss_amp_pval <- gain_del_pval <- loss_del_pval <- rep(0,length(genes_cnv))
	for(j in 1:length(genes_cnv))
	{
		#Copy number gain associations with subtype gain
		g1 <- nrow(sub_gain[grep(genes_cnv[j],sub_gain[,"cnv_driver_change_a"]),])	
		g2 <- nrow(dat[grep(genes_cnv[j],dat[,"cnv_driver_change_a"]),]) - g1
		g3 <- nrow(sub_gain) - g1
		g4 <- nrow(dat) - nrow(sub_gain) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		gain_amp_pval[j] <- fisher.test(ct)$p.value *eff
		
		#Copy number gain associations with subtype loss
		g1 <- nrow(sub_loss[grep(genes_cnv[j],sub_loss[,"cnv_driver_change_a"]),])	
		g2 <- nrow(dat[grep(genes_cnv[j],dat[,"cnv_driver_change_a"]),]) - g1
		g3 <- nrow(sub_loss) - g1
		g4 <- nrow(dat) - nrow(sub_loss) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		loss_amp_pval[j] <- fisher.test(ct)$p.value * eff
		
		#Copy number loss associations with subtype gain
		g1 <- nrow(sub_gain[grep(genes_cnv[j],sub_gain[,"cnv_driver_change_b"]),])	
		g2 <- nrow(dat[grep(genes_cnv[j],dat[,"cnv_driver_change_b"]),]) - g1
		g3 <- nrow(sub_gain) - g1
		g4 <- nrow(dat) - nrow(sub_gain) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		gain_del_pval[j] <- fisher.test(ct)$p.value * eff
		
		#Copy number loss associations with subtype loss
		g1 <- nrow(sub_loss[grep(genes_cnv[j],sub_loss[,"cnv_driver_change_b"]),])	
		g2 <- nrow(dat[grep(genes_cnv[j],dat[,"cnv_driver_change_b"]),]) - g1
		g3 <- nrow(sub_loss) - g1
		g4 <- nrow(dat) - nrow(sub_loss) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		loss_del_pval[j] <- fisher.test(ct)$p.value * eff
	}
	
		cn_res <- data.frame(gain_amp_pval,loss_amp_pval,gain_del_pval,loss_del_pval)
		#sum_signif <- apply(cn_res,2,function(x)sum(abs(x)<0.05))
		#cn_res <- rbind(cn_res,sum_signif)
		cn_res <- data.frame(genes_cnv,cn_res)
		cn_pvals[[i]] <- cn_res
}		
names(cn_pvals) <- mysubtypes
	
for(i in 1:length(mysubtypes))
{	
	#Check for alterations associated with subtype gains and losses
	myind_a <- grep(mysubtypes[i],dat[,"subtype_a"])
	myind_b <- grep(mysubtypes[i],dat[,"subtype_b"])
	
	#Check for alterations associated with subtype gains and losses	
	sub_gain <- dat[ myind_b[-which(myind_b %in% myind_a)],]
	sub_loss <- dat[ myind_a[-which(myind_a %in% myind_b)],]
			
	#Somatic alterations
	gain_amp_pval <- loss_amp_pval <- gain_del_pval <- loss_del_pval <- rep(0,length(genes_snv))
	for(j in 1:length(genes_snv))
	{
		#Mutation gain associations with subtype gain
		g1 <- nrow(sub_gain[grep(genes_snv[j],sub_gain[,"snv_driver_change_a"]),])	
		g2 <- nrow(dat[grep(genes_snv[j],dat[,"snv_driver_change_a"]),]) - g1
		g3 <- nrow(sub_gain) - g1
		g4 <- nrow(dat) - nrow(sub_gain) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		gain_amp_pval[j] <- fisher.test(ct)$p.value *eff
		
		#Mutation gain associations with subtype loss
		g1 <- nrow(sub_loss[grep(genes_snv[j],sub_loss[,"snv_driver_change_a"]),])	
		g2 <- nrow(dat[grep(genes_snv[j],dat[,"snv_driver_change_a"]),]) - g1
		g3 <- nrow(sub_loss) - g1
		g4 <- nrow(dat) - nrow(sub_loss) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		loss_amp_pval[j] <- fisher.test(ct)$p.value * eff
		
		#Mutation loss associations with subtype gain
		g1 <- nrow(sub_gain[grep(genes_snv[j],sub_gain[,"snv_driver_change_b"]),])	
		g2 <- nrow(dat[grep(genes_snv[j],dat[,"snv_driver_change_b"]),]) - g1
		g3 <- nrow(sub_gain) - g1
		g4 <- nrow(dat) - nrow(sub_gain) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		gain_del_pval[j] <- fisher.test(ct)$p.value * eff
		#Mutation gain associations with subtype loss
		g1 <- nrow(sub_loss[grep(genes_snv[j],sub_loss[,"snv_driver_change_b"]),])	
		g2 <- nrow(dat[grep(genes_snv[j],dat[,"snv_driver_change_b"]),]) - g1
		g3 <- nrow(sub_loss) - g1
		g4 <- nrow(dat) - nrow(sub_loss) - g2
		
		ct <- matrix(c(g1,g2,g3,g4),nrow=2,byrow=2)
		eff <- ifelse(g1/(g1+g3) > g2/(g2+g4),1,0)
		loss_del_pval[j] <- fisher.test(ct)$p.value * eff 
	}
	
		sn_res <- data.frame(gain_amp_pval,loss_amp_pval,gain_del_pval,loss_del_pval)
		#sum_signif <- apply(sn_res,2,function(x)sum(abs(x)<0.05))
		#sn_res <- rbind(sn_res,sum_signif)
		sn_res <- data.frame(genes_snv,sn_res)
		sn_pvals[[i]] <- sn_res
}
names(sn_pvals) <- mysubtypes

######################## 
## Common plotting elements
########################

plot_grid     <- facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme    <- theme_bw(base_size = 10) + theme(axis.title = element_text(size = 10),
                                                  axis.text = element_text(size=10),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank())
null_legend   <- theme(legend.position = 'none')
null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_text(size=7),
                       axis.ticks.y=element_blank())
bottom_x      <- theme(axis.text.x=element_blank())
null_facet    <- theme(strip.background = element_blank(),
                       strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))

########################
## View test plot
########################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  else
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

nolabels <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + null_x + null_y + null_facet + null_legend
  else
    gg + plot_theme + null_x + null_facet + null_legend
} 

gg_cbind <- function(..., widths = NULL) {
  if(length(match.call()) - 2 != length(widths))
    message("Number of widths does not match number of columns")
  gg <- gtable_cbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$widths[panels] <- unit(widths, "null")
  return(gg)
}


gg_rbind <- function(..., heights = NULL, ncol = 2) {
  if(length(match.call()) - 3 != length(heights))
    message("Number of heights does not match number of rows")
  gg <- gtable_rbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$heights[panels] <- unit(rep(heights,each = ncol), "null")
  return(gg)
}

## Extract legend
gg_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

########################
## Blank plot
########################

#Copy number plot

cnv_melt <- melt(cn_pvals)
label <- rep("+",nrow(cnv_melt))
label[grep("_del",cnv_melt[,"variable"])] <- "-"
cnv_melt[,"label"] <- paste(label, cnv_melt[,"genes_cnv"],sep=" ")

cnv_melt_gain <- cnv_melt[grep("gain_",cnv_melt[,"variable"]),]
cnv_melt_loss <- cnv_melt[grep("loss_",cnv_melt[,"variable"]),]
gain_pvalue <- cnv_melt_gain[,"value"]
gain_pvalue[which(gain_pvalue == 0)] <- 1
loss_pvalue <- cnv_melt_loss[,"value"]
loss_pvalue[which(loss_pvalue == 0)] <- 1

cnv_plot <- cnv_melt[,c("L1","label")]
cnv_plot <- cnv_plot[-which(duplicated(cnv_plot)),]
cnv_plot <- data.frame(cnv_plot,gain_pvalue,loss_pvalue)
mymin <- apply(cnv_plot[,c("gain_pvalue","loss_pvalue")],1,min)
cnv_plot <- cbind(cnv_plot,mymin)

colours <- rep("none",nrow(cnv_plot))
colours[which(cnv_plot[,"gain_pvalue"] <= cnv_plot[,"loss_pvalue"])] <- "pos"
colours[which(cnv_plot[,"gain_pvalue"] > cnv_plot[,"loss_pvalue"])] <- "neg"
colours[which(cnv_plot[,"mymin"] > 0.05)] <- "none"
cnv_plot[,"colour"] <- colours

colnames(cnv_plot) <- c("subtype","alteration","gain_pvalue","loss_pvalue","mymin","colour")#,"subtype_change")

cnv_plot[,"colour"] <- factor(cnv_plot[,"colour"],levels=c("none","neg","pos"))
cnv_plot[,"subtype"] <- factor(cnv_plot[,"subtype"],levels=c("Mesenchymal","Classical","Proneural"))
cnv_plot[,"alteration"] <- factor(cnv_plot[,"alteration"],levels=unique(cnv_plot[order(cnv_plot[,"alteration"]),"alteration"]))

#SNV plot

snv_melt <- melt(sn_pvals)
label <- rep("+",nrow(snv_melt))
label[grep("_del",snv_melt[,"variable"])] <- "-"
snv_melt[,"label"] <- paste(label, snv_melt[,"genes_snv"],sep=" ")

snv_melt_gain <- snv_melt[grep("gain_",snv_melt[,"variable"]),]
snv_melt_loss <- snv_melt[grep("loss_",snv_melt[,"variable"]),]
gain_pvalue <- snv_melt_gain[,"value"]
gain_pvalue[which(gain_pvalue == 0)] <- 1
loss_pvalue <- snv_melt_loss[,"value"]
loss_pvalue[which(loss_pvalue == 0)] <- 1

snv_plot <- snv_melt[,c("L1","label")]
snv_plot <- snv_plot[-which(duplicated(snv_plot)),]
snv_plot <- data.frame(snv_plot,gain_pvalue,loss_pvalue)
mymin <- apply(snv_plot[,c("gain_pvalue","loss_pvalue")],1,min)
snv_plot <- cbind(snv_plot,mymin)

colours <- rep("none",nrow(snv_plot))
colours[which(snv_plot[,"gain_pvalue"] <= snv_plot[,"loss_pvalue"])] <- "pos"
colours[which(snv_plot[,"gain_pvalue"] > snv_plot[,"loss_pvalue"])] <- "neg"
colours[which(snv_plot[,"mymin"] > 0.05)] <- "none"
snv_plot[,"colour"] <- colours

snv_plot[,"label"] <- paste(snv_plot[,"label"]," mut")

colnames(snv_plot) <- c("subtype","alteration","gain_pvalue","loss_pvalue","mymin","colour")#,"subtype_change")

snv_plot[,"colour"] <- factor(snv_plot[,"colour"],levels=c("none","neg","pos"))
snv_plot[,"subtype"] <- factor(snv_plot[,"subtype"],levels=c("Mesenchymal","Classical","Proneural"))
snv_plot[,"alteration"] <- factor(snv_plot[,"alteration"],levels=unique(snv_plot[order(snv_plot[,"alteration"]),"alteration"]))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/molecular_alterations.pdf",width=5.1,height=2.8)
gg_cn <-
  ggplot(cnv_plot, aes(x=alteration,y=subtype,fill=colour)) +
  geom_tile() +
  theme_bw() +
  theme(axis.title=element_blank(),
  axis.text.x=element_text(size=7,angle=45,hjust=1),
  axis.text.y=element_text(size=7),
  legend.position="none") +
  scale_fill_manual(values=c("#FFFFFF","royalblue4","tomato3"))
  
gg_sn <-
  ggplot(snv_plot, aes(x=alteration,y=subtype,fill=colour)) +
  geom_tile() +
  theme_bw() +
  theme(axis.title=element_blank(),
  axis.text.x=element_text(size=7,angle=45,hjust=1),
  axis.text.y=element_text(size=7),
  legend.position="none") +
  scale_fill_manual(values=c("#FFFFFF","royalblue4","tomato3"))
  
grid.arrange(gg_cn,gg_sn,nrow=2)
dev.off()