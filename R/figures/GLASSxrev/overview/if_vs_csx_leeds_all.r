library(odbc)
library(DBI)
library(tidyverse)
library(corrr)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
q <- "SELECT * 
      FROM analysis.cibersortx_scgp 
      WHERE aliquot_barcode LIKE 'GLSS-LU-0%B%' ORDER BY 1, 2"
dat <- dbGetQuery(con, q)

# Read in aggregate pathologist annotations
myinf1 <- dir("data/if/")
mytag <- myinf1
myinf1 <- paste("data/if/",mytag,sep="")
mytag <- mytag[grep(".csv",myinf1)]
myinf1 <- myinf1[grep(".csv", myinf1)]

thr <- c(29, 24, 13, 20, 30, 52)
names(thr) <- c("Ki67", "Olig2", "Sox2", "CD14", "CD44", "EGFR")

KI67 <- OLIG2 <- SOX2 <- CD14 <- CD44 <- EGFR <- stem <- prolifstem <- diff <- diffsox <- soxoligneg <- myeloid <- tumor_stem <- tumor_prolifstem <- tumor_diff <- tumor_ac <- tumor_mes <- tumor_oligneg <- rep(0, length(myinf1))
counts <- total <- rep(0, length(myinf1))
for(i in 1:length(myinf1))
{
  cat("\r", i)
  imm <- read.csv(myinf1[i])
  
  denom <- (sum((imm[,"Sox2"] >= thr["Sox2"] &
                #imm[,"CD14"] < thr["CD14"] & 
                ((imm[,"CD44"] >= thr["CD44"] | imm[,"EGFR"] >= thr["EGFR"]) | imm[,"Olig2"] >= thr["Olig2"] | imm[,"Ki67"] >= thr["Ki67"])) |
            (imm[,"CD14"] >= thr["CD14"]) & imm[,"Sox2"] < thr["Sox2"]))
  
  diff_count <- sum(imm[,"Sox2"] >= thr["Sox2"] & imm[,"Ki67"] < thr["Ki67"] & (imm[,"CD44"] >= thr["CD44"] | imm[,"EGFR"] >= thr["EGFR"]))
  stem_count <- sum(imm[,"Sox2"] >= thr["Sox2"] & imm[,"Ki67"] < thr["Ki67"] & imm[,"CD44"] < thr["CD44"] & imm[,"EGFR"] < thr["EGFR"] & imm[,"Olig2"] >= thr["Olig2"])
  prolifstem_count <- sum(imm[,"Sox2"] >= thr["Sox2"] &  imm[,"Ki67"] >= thr["Ki67"])
  myeloid_count <- sum( imm[,"CD14"] >= thr["CD14"] & imm[,"Sox2"] < thr["Sox2"])
  tumor_diff[i] <- diff_count/denom
  tumor_stem[i] <- stem_count/denom
  tumor_prolifstem[i] <- prolifstem_count/denom
  CD14[i] <-  myeloid_count/denom
  
  counts[i] <- sum(diff_count, stem_count, prolifstem_count,myeloid_count)
  total[i] <- denom
  
  # Option 2 (SOX2/CD14 is myeloid)
  #cat("\r", i)
  #imm <- read.csv(myinf1[i])
  #
  #denom <- (sum((imm[,"Sox2"] > thr["Sox2"] &
  #                 imm[,"CD14"] < thr["CD14"] & 
  #                 ((imm[,"CD44"] > thr["CD44"] | imm[,"EGFR"] > thr["EGFR"]) | imm[,"Olig2"] > thr["Olig2"] | imm[,"Ki67"] > thr["Ki67"])) |
  #                (imm[,"CD14"] > thr["CD14"])))
  #tumor_diff[i] <- sum(imm[,"Sox2"] > thr["Sox2"] & imm[,"CD14"] < thr["CD14"] & imm[,"Ki67"] < thr["Ki67"] & (imm[,"CD44"] > thr["CD44"] | imm[,"EGFR"] > thr["EGFR"]))/denom
  #tumor_stem[i] <- sum(imm[,"Sox2"] > thr["Sox2"] & imm[,"CD14"] < thr["CD14"] & imm[,"Ki67"] < thr["Ki67"] & imm[,"CD44"] < thr["CD44"] & imm[,"EGFR"] < thr["EGFR"] & imm[,"Olig2"] > thr["Olig2"])/denom
  #tumor_prolifstem[i] <- sum(imm[,"Sox2"] > thr["Sox2"] & imm[,"CD14"] < thr["CD14"] & imm[,"Ki67"] > thr["Ki67"])/denom
  #CD14[i] <-  sum( imm[,"CD14"] > thr["CD14"])/denom
}
sum(counts == total) #8, all cells accounted for

if_pct <- data.frame(mytag, tumor_diff, tumor_stem, tumor_prolifstem, CD14)
apply(if_pct[,2:5], 1, sum)

if_pct <- if_pct %>%
  mutate(mytag = recode(mytag, "B10_prim_gated.csv" = "GLSS-LU-0B10-TP-01R-RNA-LNPBXB", "B10_rec.csv" = "GLSS-LU-0B10-R1-01R-RNA-L1NER4",
                        "B12_prim_gated.csv" = "GLSS-LU-0B12-TP-01R-RNA-KNI3OS", "B12_rec_gated.csv" = "GLSS-LU-0B12-R1-01R-RNA-3P44CC",
                        "B13_prim_gated.csv" = "GLSS-LU-0B13-TP-01R-RNA-96LHFD", "B13_rec_gated.csv" = "GLSS-LU-0B13-R1-01R-RNA-CZ1AHH",
                        "B9_prim_gated.csv" = "GLSS-LU-00B9-TP-01R-RNA-S9XA8D", "B9_rec_gated.csv" = "GLSS-LU-00B9-R1-01R-RNA-1IVJNJ"))

tumor_dat <- dat %>%
  filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor", "prolif_stemcell_tumor", "myeloid","granulocyte","dendritic_cell")) %>%
  group_by(aliquot_barcode) %>%
  summarise(cell_state, fraction = fraction/sum(fraction)) %>%
  pivot_wider(names_from = cell_state, values_from = fraction) %>%
  mutate(myeloid = sum(myeloid, granulocyte, dendritic_cell)) %>%
  select(aliquot_barcode, differentiated_tumor, stemcell_tumor, prolif_stemcell_tumor, myeloid) %>%
  inner_join(if_pct, by = c("aliquot_barcode" = "mytag")) %>%
  filter(!(aliquot_barcode %in% c("GLSS-LU-0B12-R1-01R-RNA-3P44CC", "GLSS-LU-00B9-R1-01R-RNA-1IVJNJ")))

cor(tumor_dat[,2:5], tumor_dat[,6:ncol(tumor_dat)],method="p")
cor(tumor_dat[,2:5], tumor_dat[,6:ncol(tumor_dat)],method="s")

cor.test(tumor_dat$prolif_stemcell_tumor, tumor_dat$tumor_prolifstem, method="p")
cor.test(tumor_dat$stemcell_tumor, tumor_dat$tumor_stem, method="p")
cor.test(tumor_dat$differentiated_tumor, tumor_dat$tumor_diff,method="p")
cor.test(tumor_dat$myeloid, tumor_dat$CD14,method="p")


#plot(tumor_dat$differentiated_tumor, tumor_dat$CD44)
cor(c(tumor_dat$tumor_stem, tumor_dat$tumor_prolifstem, tumor_dat$tumor_diff), c(tumor_dat$stemcell_tumor, tumor_dat$prolif_stemcell_tumor, tumor_dat$differentiated_tumor), method="s")
cor(c(tumor_dat$tumor_stem, tumor_dat$tumor_prolifstem, tumor_dat$tumor_diff), c(tumor_dat$stemcell_tumor, tumor_dat$prolif_stemcell_tumor, tumor_dat$differentiated_tumor), method="p")
cor.test(c(tumor_dat$tumor_stem, tumor_dat$tumor_prolifstem, tumor_dat$tumor_diff), c(tumor_dat$stemcell_tumor, tumor_dat$prolif_stemcell_tumor, tumor_dat$differentiated_tumor), method="p")$p.value

aliquot_barcode <- rep(tumor_dat$aliquot_barcode, 4)
if_dat <- c(tumor_dat$tumor_stem, tumor_dat$tumor_prolifstem, tumor_dat$tumor_diff, tumor_dat$CD14)
csx_dat <- c(tumor_dat$stemcell_tumor, tumor_dat$prolif_stemcell_tumor, tumor_dat$differentiated_tumor, tumor_dat$myeloid)
cell_state <- rep(c("Stem-like", "Prolif. stem-like","Diff.-like", "Myeloid"), each = nrow(tumor_dat))

plot_dat <- data.frame(aliquot_barcode, if_dat, csx_dat, cell_state)

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions/all_cor.pdf", height = 3, width = 3)
ggplot(plot_dat, aes(x = if_dat * 100, y = csx_dat * 100)) +
  geom_point(aes(colour=cell_state)) +
  geom_smooth(method = lm, se = FALSE) + 
  scale_colour_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                               "Oligodendrocyte" = "#2ca25f",
                               "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                               "Fibroblast" = "#feb24c",
                               "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  theme_classic() +
  labs(x = "Immunofluorescence (%)", y = "CIBERSORTx (%)") +
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        legend.position="none") 
dev.off()

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions/all_cor_split.pdf",width=4,height=1.3)
ggplot(plot_dat, aes(x=if_dat*100,y=csx_dat*100)) +
  geom_point(colour="black",size=1) +
  geom_smooth(method = lm, se = FALSE, aes(colour = cell_state)) + 
  facet_wrap(.~cell_state, scales = "free", ncol = 4) + 
  scale_colour_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                               "Oligodendrocyte" = "#2ca25f",
                               "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                               "Fibroblast" = "#feb24c",
                               "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  theme_classic() +
  labs(x = "Immunofluorescence (%)", y = "CIBERSORTx (%)") +
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=7), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=7),
        legend.position="none") 
#scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) #+
#scale_y_continuous(limits = c(0, NA)) +
#scale_x_continuous(limits = c(0, NA))
dev.off()

