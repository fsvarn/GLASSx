library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggbeeswarm)


#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT im.*, aliquot_batch
FROM biospecimen.aliquots al
JOIN analysis.davoli_immune_score im ON im.aliquot_barcode = al.aliquot_barcode
ORDER BY 1, 2"

res <- dbGetQuery(con, q)

aliquot_barcode <- unique(res[,"aliquot_barcode"])
cells <- unique(res[,"signature_name"])

mat <- matrix(0,nrow=length(unique(res[,"signature_name"])),ncol=length(aliquot_barcode))
rownames(mat) <- cells
colnames(mat) <- aliquot_barcode
batch <- rep("",length(aliquot_barcode))
names(batch) <- aliquot_barcode
for(i in 1:length(aliquot_barcode))
{
	mat[,i] <- res[which(res[,"aliquot_barcode"]==aliquot_barcode[i]),"enrichment_score"]
	batch[i] <- res[which(res[,"aliquot_barcode"]==aliquot_barcode[i]),"aliquot_batch"][1]
}

mypca <- prcomp(t(mat),scale.=TRUE)
scores <- as.data.frame(mypca$x)
scores[,"aliquot_batch"] <- as.factor(batch)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/pca_immune.pdf",width=7,height=7)
ggplot(data = scores, aes(x = PC1, y = PC2, colour=aliquot_batch)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  scale_colour_manual(values=brewer.pal(12,"Set3")) +
  ggtitle("PCA plot of GLASS immune scores by aliquot") +
  theme_classic()
dev.off()

#NO BATCH EFFECT <3