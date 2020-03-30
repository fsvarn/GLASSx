library(tidyverse)
library(GSVA)
library(DBI)
library(odbc)
library(reshape)

rm(list=ls())

#Load and prepare gene expression matrix
myinf1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"

expr <- read.delim(myinf1,row.names=1)
colnames(expr) <- gsub("\\.","-",colnames(expr))
expr <- as.matrix(expr)

#only use genes expressed in at least half of the cohort
sums <- apply(expr,1,function(x)sum(x>0))
#means <- apply(expr,1,mean)
expr <- expr[which(sums>(0.5*max(sums))),]

#Load Davoli gene sets
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT * FROM ref.immune_signatures WHERE signature_set = 'Davoli'"

imm_sigs <- dbGetQuery(con, q)

#Convert table to signature list
immune_cell <- unique(imm_sigs[,"signature_name"])
gene_list <- list()
for(i in 1:length(immune_cell))
{
	sub_imm_sig <- imm_sigs[which(imm_sigs[,"signature_name"] == immune_cell[i]),]
	sig_genes <- sub_imm_sig[,"gene_symbol"]
	gene_list[[i]] <- sig_genes
}
names(gene_list) <- immune_cell

res <- gsva(expr, gene_list, method="ssgsea")

cor(t(res),method="s")

#                       CD4.mature CD8.effector  NK.cells   B.cells     T.reg
# CD4.mature             1.0000000    0.4703657 0.5410595 0.5385371 0.4697823
# CD8.effector           0.4703657    1.0000000 0.7592786 0.5473699 0.5223386
# NK.cells               0.5410595    0.7592786 1.0000000 0.5361608 0.4717043
# B.cells                0.5385371    0.5473699 0.5361608 1.0000000 0.4095615
# T.reg                  0.4697823    0.5223386 0.4717043 0.4095615 1.0000000
# Dendritic              0.5395383    0.5848674 0.4695568 0.5534723 0.6556318
# CD8.effector.NK.cells  0.5050089    0.8519542 0.7838420 0.5519157 0.4863984
# Macrophages            0.4785151    0.4246606 0.3376630 0.4031107 0.5732111
# Macrophages.M2         0.4952228    0.4328153 0.3392596 0.4148362 0.6093651
# Macrophages.M1         0.5875026    0.5214634 0.4728293 0.4805586 0.6373450
#                       Dendritic CD8.effector.NK.cells Macrophages
# CD4.mature            0.5395383             0.5050089   0.4785151
# CD8.effector          0.5848674             0.8519542   0.4246606
# NK.cells              0.4695568             0.7838420   0.3376630
# B.cells               0.5534723             0.5519157   0.4031107
# T.reg                 0.6556318             0.4863984   0.5732111
# Dendritic             1.0000000             0.5939385   0.7255236
# CD8.effector.NK.cells 0.5939385             1.0000000   0.4166863
# Macrophages           0.7255236             0.4166863   1.0000000
# Macrophages.M2        0.7480854             0.4432382   0.9685218
# Macrophages.M1        0.7242990             0.5620241   0.9016750
#                       Macrophages.M2 Macrophages.M1
# CD4.mature                 0.4952228      0.5875026
# CD8.effector               0.4328153      0.5214634
# NK.cells                   0.3392596      0.4728293
# B.cells                    0.4148362      0.4805586
# T.reg                      0.6093651      0.6373450
# Dendritic                  0.7480854      0.7242990
# CD8.effector.NK.cells      0.4432382      0.5620241
# Macrophages                0.9685218      0.9016750
# Macrophages.M2             1.0000000      0.9003250
# Macrophages.M1             0.9003250      1.0000000



#PCA plot (to confirm no batch effect)

res <- as.data.frame(res)
res[,"signature_name"] <- rownames(res)
res <- melt(res, id="signature_name")

res <- res[,c(2,1,3)]
colnames(res) <- c("aliquot_barcode","signature_name","enrichment_score")
head(res)

dbWriteTable(con, Id(schema="analysis", table="davoli_immune_score"), res, overwrite=TRUE)

