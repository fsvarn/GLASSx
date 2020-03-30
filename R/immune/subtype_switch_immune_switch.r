library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
WITH subtype_label AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	ts1.signature_name, 
	ts1.p_value AS pval_a,
	ts2.p_value AS pval_b,
	CASE WHEN ts1.p_value < 0.05 THEN ts1.signature_name ELSE NULL END AS subtype_a,
	CASE WHEN ts2.p_value < 0.05 THEN ts2.signature_name ELSE NULL END AS subtype_b,
	su.idh_codel_subtype
	FROM analysis.rna_silver_set ss
	JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
	JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
	JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
	JOIN clinical.surgeries su ON su.sample_barcode = al1.sample_barcode
	WHERE idh_codel_subtype IS NOT NULL AND
	ts1.p_value < 0.05 OR ts2.p_value < 0.05
),
subtype_collapse AS
(
	SELECT sl.tumor_pair_barcode,
	sl.case_barcode,
	string_agg(subtype_a,',') AS subtype_a,
	string_agg(subtype_b,',') AS subtype_b,
	idh_codel_subtype
	FROM subtype_label sl
	GROUP BY tumor_pair_barcode, sl.case_barcode, idh_codel_subtype
	ORDER BY 1
)
SELECT sc.*,
ic.change AS immune_change
FROM subtype_collapse sc
JOIN analysis.immune_cluster_change ic ON ic.case_barcode = sc.case_barcode

"
dat <- dbGetQuery(con, q)
dat[which(is.na(dat[,"subtype_b"])),"subtype_b"] <- "Mixed"
dat[which(is.na(dat[,"subtype_a"])),"subtype_a"] <- "Mixed"

mesenchymal_gain <-	!(str_detect("Mesenchymal",dat[,"subtype_a"])) & str_detect("Mesenchymal",dat[,"subtype_b"]) |
					str_detect(",",dat[,"subtype_a"]) & str_detect("Mesenchymal",dat[,"subtype_b"])
mesenchymal_loss <-	str_detect("Mesenchymal",dat[,"subtype_a"]) & !(str_detect("Mesenchymal",dat[,"subtype_b"]))

#Check for relationship between immune increase and mesenchymal gain
g1 <- nrow(dat[which(dat[,"immune_change"]=="Increase" & mesenchymal_gain),])
g2 <- nrow(dat[which(dat[,"immune_change"]!="Increase" & mesenchymal_gain),])
g3 <- nrow(dat[which(dat[,"immune_change"]=="Increase" & !mesenchymal_gain),])
g4 <- nrow(dat[which(dat[,"immune_change"]!="Increase" & !mesenchymal_gain),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2,byrow=TRUE)
fisher.test(ct)

counts <- c(g1,g2,g3,g4)
subtype_switch <- c("Mesenchymal gain","Mesenchymal gain","No mesenchymal gain","No mesenchymal gain")
immune_switch <- c("Increase","No increase","Increase","No increase")
plot_gain <- data.frame(subtype_switch,immune_switch,counts)

plot_gain[,"percent"] <- c(counts[1:2]/sum(counts[1:2]), counts[3:4]/sum(counts[3:4]))

#Check for relationship between immune decrease and mesenchymal loss
g1 <- nrow(dat[which(dat[,"immune_change"]=="Decrease" & mesenchymal_loss),])
g2 <- nrow(dat[which(dat[,"immune_change"]!="Decrease" & mesenchymal_loss),])
g3 <- nrow(dat[which(dat[,"immune_change"]=="Decrease" & !mesenchymal_loss),])
g4 <- nrow(dat[which(dat[,"immune_change"]!="Decrease" & !mesenchymal_loss),])

ct <- matrix(c(g1,g2,g3,g4),nrow=2,ncol=2,byrow=TRUE)
fisher.test(ct)

counts <- c(g1,g2,g3,g4)
subtype_switch <- c("Mesenchymal loss","Mesenchymal loss","No mesenchymal loss","No mesenchymal loss")
immune_switch <- c("Decrease","No decrease","Decrease","No decrease")
plot_loss <- data.frame(subtype_switch,immune_switch,counts)

plot_loss[,"percent"] <- c(counts[1:2]/sum(counts[1:2]), counts[3:4]/sum(counts[3:4]))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_subtype_switch_bar.pdf",width=2,height=1.2)
p1 <- ggplot(data = plot_gain, aes(x = subtype_switch, y = percent)) +
geom_bar(aes(fill=immune_switch),stat="identity") +
scale_fill_manual(values=c("#BD1E2D","gray50")) +
labs(y = "Percent") +
theme_bw() +
theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
	axis.text = element_text(size=7),
	axis.title.x = element_blank(), axis.title.y=element_text(size=7),
	legend.position = "none")
	
p2 <- ggplot(data = plot_loss, aes(x = subtype_switch, y = percent)) +
geom_bar(aes(fill=immune_switch),stat="identity") +
scale_fill_manual(values=c("#26ABE2","gray50")) +
labs(y = "Percent") +
theme_bw() +
theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
	axis.text = element_text(size=7),
	axis.title.x = element_blank(), axis.title.y=element_text(size=7),
	legend.position = "none")
grid.arrange(p1,p2,ncol=2)
dev.off()

