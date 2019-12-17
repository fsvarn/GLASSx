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
# CD4.mature             1.0000000    0.4994420 0.3041680 0.5081072 0.4851560
# CD8.effector           0.4994420    1.0000000 0.6972606 0.4736223 0.5243411
# NK.cells               0.3041680    0.6972606 1.0000000 0.3891978 0.2342843
# B.cells                0.5081072    0.4736223 0.3891978 1.0000000 0.3825280
# T.reg                  0.4851560    0.5243411 0.2342843 0.3825280 1.0000000
# Dendritic              0.5418421    0.5790097 0.2647034 0.4530233 0.6393420
# CD8.effector.NK.cells  0.5426428    0.8563622 0.7714938 0.5169264 0.5072840
# Macrophages            0.4714316    0.4505275 0.1414959 0.4260011 0.5674009
# Macrophages.M2         0.4869551    0.4509262 0.1155329 0.4233100 0.6017204
# Macrophages.M1         0.5775275    0.5484313 0.2687855 0.4853260 0.6458493
#                       Dendritic CD8.effector.NK.cells Macrophages
# CD4.mature            0.5418421             0.5426428   0.4714316
# CD8.effector          0.5790097             0.8563622   0.4505275
# NK.cells              0.2647034             0.7714938   0.1414959
# B.cells               0.4530233             0.5169264   0.4260011
# T.reg                 0.6393420             0.5072840   0.5674009
# Dendritic             1.0000000             0.6017597   0.7470737
# CD8.effector.NK.cells 0.6017597             1.0000000   0.4627895
# Macrophages           0.7470737             0.4627895   1.0000000
# Macrophages.M2        0.7688800             0.4707170   0.9713643
# Macrophages.M1        0.7442438             0.5985528   0.9091088
#                       Macrophages.M2 Macrophages.M1
# CD4.mature                 0.4869551      0.5775275
# CD8.effector               0.4509262      0.5484313
# NK.cells                   0.1155329      0.2687855
# B.cells                    0.4233100      0.4853260
# T.reg                      0.6017204      0.6458493
# Dendritic                  0.7688800      0.7442438
# CD8.effector.NK.cells      0.4707170      0.5985528
# Macrophages                0.9713643      0.9091088
# Macrophages.M2             1.0000000      0.9033194
# Macrophages.M1             0.9033194      1.0000000



#PCA plot (to confirm no batch effect)

res <- as.data.frame(res)
res[,"signature_name"] <- rownames(res)
res <- melt(res, id="signature_name")

res <- res[,c(2,1,3)]
colnames(res) <- c("aliquot_barcode","signature_name","enrichment_score")
head(res)

dbWriteTable(con, Id(schema="analysis", table="davoli_immune_score"), res, overwrite=TRUE)

