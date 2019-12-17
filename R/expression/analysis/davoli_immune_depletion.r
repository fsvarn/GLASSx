library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
SELECT an.rna_barcode, an.dna_barcode, im.signature_name, im.enrichment_score, nd.rneo
FROM analysis.davoli_immune_score im
JOIN analysis.platinum_set ps ON im.aliquot_barcode = ps.rna_barcode_a OR im.aliquot_barcode = ps.rna_barcode_b
JOIN analysis.analyte_sets an ON an.rna_barcode = im.aliquot_barcode
JOIN analysis.neoantigen_depletion nd ON nd.aliquot_barcode = an.dna_barcode
WHERE nd.nobs >= 3
ORDER BY 1,3
"

imm_rneo <- dbGetQuery(con,q)

cells <- unique(imm_rneo[,"signature_name"])

spearman <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_imm_rneo <- imm_rneo[which(imm_rneo[,"signature_name"]==cells[i]),]
	spearman[i] <- cor(sub_imm_rneo[,"enrichment_score"],sub_imm_rneo[,"rneo"],method="s")
	p.value[i] <- cor.test(sub_imm_rneo[,"enrichment_score"],sub_imm_rneo[,"rneo"],method="s")$p.value
}
res <- data.frame(cells,spearman,p.value)