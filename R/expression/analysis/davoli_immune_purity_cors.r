library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
SELECT an.rna_barcode, an.dna_barcode, im.signature_name, im.enrichment_score, sq.cellularity
FROM analysis.davoli_immune_score im
JOIN analysis.platinum_set ps ON im.aliquot_barcode = ps.rna_barcode_a OR im.aliquot_barcode = ps.rna_barcode_b
JOIN analysis.analyte_sets an ON an.rna_barcode = im.aliquot_barcode
JOIN analysis.pairs pa ON pa.tumor_barcode = an.dna_barcode
JOIN variants.seqz_params sq ON pa.pair_barcode = sq.pair_barcode
ORDER BY 1,3
"

imm_pur <- dbGetQuery(con,q)

cells <- unique(imm_pur[,"signature_name"])

spearman <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_imm_pur <- imm_pur[which(imm_pur[,"signature_name"]==cells[i]),]
	spearman[i] <- cor(sub_imm_pur[,"enrichment_score"],sub_imm_pur[,"cellularity"],method="s")
	p.value[i] <- cor.test(sub_imm_pur[,"enrichment_score"],sub_imm_pur[,"cellularity"],method="s")$p.value
}
res <- data.frame(cells,spearman,p.value)

#                    cells   spearman      p.value
# 1                B.cells -0.2194587 1.569112e-03
# 2             CD4.mature -0.2350478 6.927121e-04
# 3           CD8.effector -0.2964095 1.593412e-05
# 4  CD8.effector.NK.cells -0.3139115 4.573983e-06
# 5              Dendritic -0.3996962 2.900755e-09
# 6            Macrophages -0.3770809 2.502978e-08
# 7         Macrophages.M1 -0.3209851 2.699252e-06
# 8         Macrophages.M2 -0.3838556 1.334979e-08
# 9               NK.cells -0.1907395 6.153797e-03
# 10                 T.reg -0.1062457 1.294689e-01