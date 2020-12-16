library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(qusage)

rm(list=ls())

myinf1 <- "/projects/verhaak-lab/GLASS-III/data/res/verhaak_prelim/msigdb_hallmark_enrichment_idhwt.txt"

es <- read.delim(myinf1)
colnames(es) <- gsub("\\.","-",colnames(es))
es <- t(es)


#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT * 
FROM analysis.tumor_rna_clinical_comparison
WHERE idh_codel_subtype = 'IDHwt' AND received_rt AND tumor_barcode_a LIKE '%-TP-%' AND tumor_barcode_b LIKE '%-R1-%'
"

info <- dbGetQuery(con, q)


# Define responders as patients whose tumors do not come back within a year
resp_a <- info %>% filter(surgical_interval >= 12) %>% .$tumor_barcode_a
resi_a <- info %>% filter(surgical_interval < 12) %>% .$tumor_barcode_a
resp_b <- info %>% filter(surgical_interval >= 12) %>% .$tumor_barcode_b
resi_b <- info %>% filter(surgical_interval < 12) %>% .$tumor_barcode_b

resp_p <- resi_p <- resp_eff <- resi_eff <- rep(0, ncol(es))
for(i in 1:ncol(es))
{
	resp_p[i] <- wilcox.test(es[resp_a,i], es[resp_b,i], paired=TRUE)$p.value
	resi_p[i] <- wilcox.test(es[resi_a,i], es[resi_b,i], paired=TRUE)$p.value
	
	resp_eff[i] <- median(es[resp_b,i] - es[resp_a,i])
	resi_eff[i] <- median(es[resi_b,i] - es[resi_a,i])
}
gene_set <- colnames(es)
res <- data.frame(gene_set, resp_p, resp_eff, resi_p, resi_eff)
res[order(res[,"resi_p"]),]

# Create a long table
gene_set <- rep(gene_set,2)
logp <- c(-log10(resi_p), -log10(resp_p))
eff <- c(resi_eff < 0, resp_eff < 0)
eff <- ifelse(eff, "Negative","Positive")
eff[which(logp < -log10(0.05))] <- "Non"
response <- c(rep("Recurrence < 1 year", length(resi_p)), rep("Recurrence > 1 year", length(resp_p)))
gene_set <- factor(gene_set, levels = rev(gene_set[order(resi_p)]))

plot_res <- data.frame(gene_set, logp, eff, response)

# Plot the relevant result here
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/hallmark_pathway_response_resistance.pdf", height=5, width=7)
ggplot(plot_res, aes(x = gene_set, y = logp, fill=eff)) +
geom_bar(stat="identity") +
geom_hline(yintercept=-log10(0.05),linetype="longdash") +
labs(y = "-log(p-value)") +
facet_grid(.~response) +
scale_fill_manual(values=c("royalblue4","gray","tomato3")) +
theme_bw() +
theme(axis.text.x= element_text(size=7), axis.text.y=element_text(size=7),
axis.title.x = element_text(size=7), axis.title.y = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_flip()
dev.off()
