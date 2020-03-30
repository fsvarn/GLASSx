library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(ggalluvial)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
WITH subtype_rank AS
(
	SELECT *,
	RANK() OVER (PARTITION BY aliquot_barcode ORDER BY p_value ASC) AS p_rank
	FROM analysis.transcriptional_subtype
),
top_rank AS
(
	SELECT *
	FROM subtype_rank
	WHERE p_rank = 1
),
agg AS
(
	SELECT aliquot_barcode, 
	string_agg(signature_name,',') AS subtype 
	FROM top_rank
	GROUP BY aliquot_barcode
)
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ag1.subtype AS subtype_a, 
ag2.subtype AS subtype_b
FROM analysis.rna_silver_set ss
JOIN agg ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN agg ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
"

dat <- dbGetQuery(con,q)

dat[which(dat[,"subtype_a"]=="Mesenchymal,Classical"),"subtype_a"]  <- "Classical,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Classical"),"subtype_b"]  <- "Classical,Mesenchymal"

dat[which(dat[,"subtype_a"]=="Mesenchymal,Proneural"),"subtype_a"]  <- "Proneural,Mesenchymal"
dat[which(dat[,"subtype_b"]=="Mesenchymal,Proneural"),"subtype_b"]  <- "Proneural,Mesenchymal"

#Convert for the purposes of plotting
dat[which(dat[,"subtype_a"]=="Proneural,Classical"),"subtype_a"]  <- "Proneural"
dat[which(dat[,"subtype_b"]=="Proneural,Classical"),"subtype_b"]  <- "Proneural"

dat[which(dat[,"subtype_a"]=="Proneural,Mesenchymal"),"subtype_a"]  <- "Mesenchymal"
dat[which(dat[,"subtype_b"]=="Proneural,Mesenchymal"),"subtype_b"]  <- "Mesenchymal"

dat[which(dat[,"subtype_a"]=="Classical,Mesenchymal"),"subtype_a"]  <- "Classical"
dat[which(dat[,"subtype_b"]=="Classical,Mesenchymal"),"subtype_b"]  <- "Classical"

#dat[,"subtype_a"] <- recode(dat[,"subtype_a"], "Proneural" = "Pro.","Classical" = "Cla.","Mesenchymal" = "Mes.")
#dat[,"subtype_b"] <- recode(dat[,"subtype_b"], "Proneural" = "Pro.","Classical" = "Cla.","Mesenchymal" = "Mes.")

#Sankey plot for transcriptional subtype
dat[,"subtype_a"] <- factor(dat[,"subtype_a"],levels=c("Proneural","Classical","Mesenchymal"))
dat[,"subtype_b"] <- factor(dat[,"subtype_b"],levels=c("Proneural","Classical","Mesenchymal"))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/subtype_sankey_plot.pdf",width=3,height=2)
ggplot(dat, aes(axis1 = subtype_a, axis2 = subtype_b, y=case_barcode)) +
geom_alluvium(aes(fill = subtype_a), width = 1/12, knot.pos = 0.4) +
geom_stratum(width = 1/3, fill = "white", color = "grey") +
scale_x_discrete(limits = c("Initial", "Recurrent"),expand = c(0,0)) +
# scale_fill_manual(values=c("#00458a","#008a22","#8a0000","#8a7900","#5a008a","#008a66")) +
scale_fill_manual(values=c("#00458a","#008a22","#8a0000")) +
geom_text(stat = "stratum", infer.label=TRUE, reverse=TRUE,size=2.5) +
theme_bw() +
theme(axis.text.x= element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 
dev.off()
