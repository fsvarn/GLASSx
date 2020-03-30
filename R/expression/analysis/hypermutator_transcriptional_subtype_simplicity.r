library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(lme4)
library(car)
library(reshape)
library(grid)

#######################################################
rm(list=ls())
# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
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
SELECT ps.*, 
ag1.subtype AS subtype_a, 
ag2.subtype AS subtype_b,
CASE WHEN ag1.subtype = ag2.subtype THEN 'None' ELSE 'Switch' END AS switch,
si1.simplicity_score AS ss_a,
si2.simplicity_score AS ss_b,
tc.grade_change,
tc.received_tmz,
tc.received_rt,
tc.alk_assoc_hypermutator_status,
tc.hypermutator_status,
cs.idh_codel_subtype
FROM analysis.platinum_set ps
JOIN agg ag1 ON ag1.aliquot_barcode = ps.rna_barcode_a
JOIN agg ag2 ON ag2.aliquot_barcode = ps.rna_barcode_b
JOIN analysis.tumor_clinical_comparison tc ON tc.tumor_barcode_a = ps.dna_barcode_a AND tc.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.simplicity_score si1 ON ps.rna_barcode_a = si1.aliquot_barcode
JOIN analysis.simplicity_score si2 ON ps.rna_barcode_b = si2.aliquot_barcode
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE cs.idh_codel_subtype IS NOT NULL 
ORDER BY 1, 2"

dat <- dbGetQuery(con,q)


###################################################
# Step 1: Examine whether hypermutation is associated with subtype switching
##################################################

g1 <- nrow(dat[which(dat[,"switch"]=="Switch" & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g2 <- nrow(dat[which(dat[,"switch"]=="None" & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g3 <- nrow(dat[which(dat[,"switch"]=="Switch" & dat[,"alk_assoc_hypermutator_status"]=="1"),])
g4 <- nrow(dat[which(dat[,"switch"]=="None" & dat[,"alk_assoc_hypermutator_status"]=="1"),])

fisher.test(matrix(c(g1,g2,g3,g4),nrow=2))	# P = 0.75

###################################################
# Step 2: Examine whether hypermutation is associated with a specific subtype
##################################################

g1 <- nrow(dat[which(grepl("Mesenchymal",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g2 <- nrow(dat[which(!grepl("Mesenchymal",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g3 <- nrow(dat[which(grepl("Mesenchymal",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="1"),])
g4 <- nrow(dat[which(!grepl("Mesenchymal",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="1"),])

fisher.test(matrix(c(g1,g2,g3,g4),nrow=2))	# P = 1

g1 <- nrow(dat[which(grepl("Proneural",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g2 <- nrow(dat[which(!grepl("Proneural",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g3 <- nrow(dat[which(grepl("Proneural",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="1"),])
g4 <- nrow(dat[which(!grepl("Proneural",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="1"),])

fisher.test(matrix(c(g1,g2,g3,g4),nrow=2))	# P = 0.2; likely due to higher rate of IDHmut in hypermutator cases

g1 <- nrow(dat[which(grepl("Classical",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g2 <- nrow(dat[which(!grepl("Classical",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="0"),])
g3 <- nrow(dat[which(grepl("Classical",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="1"),])
g4 <- nrow(dat[which(!grepl("Classical",dat[,"subtype_b"]) & dat[,"alk_assoc_hypermutator_status"]=="1"),])

fisher.test(matrix(c(g1,g2,g3,g4),nrow=2))	# P = 0.31

###################################################
# Step 3: Compare the simplicity score between initial and recurrent hypermutated tumors
##################################################

# Non-hypermutators
nhm_dat <- dat[which(dat[,"hypermutator_status"]=="0"),]
pval1 <- wilcox.test(nhm_dat[,"ss_a"],nhm_dat[,"ss_b"],paired=TRUE)$p.value	#P = 0.62
pval1 <- round(pval1, 2)

# Only hypermutators
hm_dat <- dat[which(dat[,"hypermutator_status"]=="1"),]
pval2 <- wilcox.test(hm_dat[,"ss_a"],hm_dat[,"ss_b"],paired=TRUE)$p.value	#P = 0.06
pval2 <- round(pval2, 2)

# Samples treated by alkylating agents that are not hypermutators (not plotted)
alk_dat <- dat[which(dat[,"hypermutator_status"]=="0" & dat[,"received_tmz"]==1),]
pval3 <- wilcox.test(alk_dat[,"ss_a"],alk_dat[,"ss_b"],paired=TRUE)$p.value	#P = 0.16
pval3 <- round(pval1, 2)


###################################################
# Step 4: Plot hypermutator simplicity association
##################################################

case_barcode <- rep(dat[,"case_barcode"],2)
aliquot_barcode <- c(dat[,"rna_barcode_a"], dat[,"rna_barcode_b"])
simplicity_score <- c(dat[,"ss_a"],dat[,"ss_b"])
subtype <- rep(dat[,"idh_codel_subtype"],2)
timepoint <- c(rep("Initial",nrow(dat)),rep("Recurrent",nrow(dat)))
hypermutator <- factor(rep(dat[,"hypermutator_status"],2),levels=c(0,1),labels=c("Non-hypermutator","Hypermutator"))
plot_res <- data.frame(case_barcode, aliquot_barcode, simplicity_score, subtype, timepoint, hypermutator)

p_value <- c(deparse((bquote(italic("P") ~" = " ~ .(pval1)))),
			 deparse((bquote(italic("P") ~" = " ~ .(pval2)))))
annotation_text <- data.frame(hypermutator = factor(c("Non-hypermutator","Hypermutator")),
							  timepoint = 1.1,
							  simplicity_score = min(plot_res[,"simplicity_score"])-(0.15 * (max(plot_res[,"simplicity_score"]) - min(plot_res[,"simplicity_score"]))) + 
							  0.04 * abs(min(plot_res[,"simplicity_score"])-(0.15 * (max(plot_res[,"simplicity_score"]) - min(plot_res[,"simplicity_score"])))),
							  p_value)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/hypermutator_simplicity.pdf",width=2.5,height=2)
ggplot(plot_res, aes(x = timepoint, y = simplicity_score)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode,colour=subtype)) +
geom_point(size=1,colour="black") +
facet_grid(.~hypermutator, scales = "free_x") +
labs(y = "Simplicity score") +
geom_text(data=annotation_text,label=p_value, size=2.5, parse=TRUE) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(min(plot_res[,"simplicity_score"])-(0.15 * (max(plot_res[,"simplicity_score"]) - min(plot_res[,"simplicity_score"]))),
max(plot_res[,"simplicity_score"])))
dev.off()

