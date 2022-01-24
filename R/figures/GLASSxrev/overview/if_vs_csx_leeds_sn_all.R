###################################################
# Final immunofluorescence versus CIBERSORTx analysis; includes all samples (Leeds + IDH mutant)
# Author: Frederick Varn
# Date: 2021.12.30
# Figures S1C, 3B, 3D
##################################################

library(odbc)
library(DBI)
library(tidyverse)
library(corrr)
library(gridExtra)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
q <- "SELECT * 
      FROM analysis.cibersortx_scgp 
      WHERE aliquot_barcode LIKE 'GLSS-LU-0%B%' OR aliquot_barcode LIKE 'GLSS-SN-0017-%' ORDER BY 1, 2"
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
  
  diff_count <- sum(imm[,"Sox2"] >= thr["Sox2"] & imm[,"Ki67"] < thr["Ki67"] & (imm[,"CD44"] >= thr["CD44"] | imm[,"EGFR"] >= thr["EGFR"]) & imm[,"Olig2"] < thr["Olig2"])
  stem_count <- sum(imm[,"Sox2"] >= thr["Sox2"] & imm[,"Ki67"] < thr["Ki67"] & imm[,"Olig2"] >= thr["Olig2"])
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
                        "B9_prim_gated.csv" = "GLSS-LU-00B9-TP-01R-RNA-S9XA8D", "B9_rec_gated.csv" = "GLSS-LU-00B9-R1-01R-RNA-1IVJNJ",
                        "S16_IDH_mut.csv" = "GLSS-SN-0017-TP-01R-RNA-FX5HYW", "S19_IDH_mut.csv" = "GLSS-SN-0017-R1-01R-RNA-OEK8D6"))

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
plot_dat <- plot_dat %>% 
mutate(cell_state = fct_relevel(cell_state, "Prolif. stem-like", "Stem-like", "Diff.-like", "Myeloid"))

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions/all_cor_idhmut.pdf", height = 3, width = 3)
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

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions_final/Fig1/all_cor_split_idhmut_facet.pdf",width=3.86,height=1.3)
ggplot(plot_dat, aes(x=if_dat*100,y=csx_dat*100)) +
  geom_point(colour="black",size=1) +
  geom_smooth(method = lm, se = FALSE, aes(colour = cell_state),fullrange=TRUE) + 
  #geom_abline(aes(colour = cell_state, intercept = 0, slope = 1)) +
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
        legend.position="none") +
scale_y_continuous(labels = scales::number_format(accuracy = 1), limits = c(0, NA)) +
scale_x_continuous(limits = c(0, NA))
#scale_y_continuous(labels = scales::number_format(accuracy = 1), limits = c(0, 100)) +
#scale_x_continuous(limits = c(0, 100))
dev.off()

state <- unique(plot_dat$cell_state)
se <- list()
for(i in 1:length(state))
{
  #tmp <- plot_dat %>% filter(cell_state == state[i])
  #xlim = range(tmp$if_dat) * 100
  #ylim = range(tmp$csx_dat) * 100
  #lims <- c(0, max(c(xlim, ylim) + 0.05*max(c(xlim,ylim))))
  
  se[[i]] <- ggplot(plot_dat %>% filter(cell_state == state[i]), aes(x=if_dat*100,y=csx_dat*100)) +
    geom_point(colour="black",size=1) +
    geom_smooth(method = lm, se = FALSE, aes(colour = cell_state),fullrange=TRUE) + 
    #geom_abline(aes(colour = cell_state, intercept = 0, slope = 1)) +
    #facet_wrap(.~cell_state, scales = "free", ncol = 4) + 
    scale_colour_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                                 "Oligodendrocyte" = "#2ca25f",
                                 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                                 "Fibroblast" = "#feb24c",
                                 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
    theme_classic() +
    labs(title = as.character(state[i]), x = "Immunofluorescence (%)", y = "CIBERSORTx (%)") +
    theme(axis.text = element_text(size=7),
          axis.title = element_text(size=7), 
          plot.title = element_text(size=7,hjust = 0.5),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size=7),
          legend.position="none") +
    scale_y_continuous(labels = scales::number_format(accuracy = 1), limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA))
}

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions_final/Fig1/all_cor_split_idhmut.pdf",width=4.4,height=1.255)
grid.arrange(se[[2]], se[[1]], se[[3]], se[[4]],nrow = 1)
dev.off()


# Quantify the Ki67 expression in SOX2 cells in IDHmut and B13 (hypermutant)

# IDHmut
imm1 <-   read.csv("data/if/S16_IDH_mut.csv")
imm2 <-   read.csv("data/if/S19_IDH_mut.csv")

c1 <- nrow(imm1 %>% filter(Sox2 >= thr["Sox2"], Ki67 >= thr["Ki67"]))
c2 <- nrow(imm1 %>% filter(Sox2 >= thr["Sox2"], Ki67 < thr["Ki67"]))
c3 <- nrow(imm2 %>% filter(Sox2 >= thr["Sox2"], Ki67 >= thr["Ki67"]))
c4 <- nrow(imm2 %>% filter(Sox2 >= thr["Sox2"], Ki67 < thr["Ki67"]))                      

mut_ct <- matrix(c(c1,c2,c3,c4),nrow=2)
fisher.test(mut_ct)$p.value

fraction <- c(c1/(c1+c2),c2/(c1+c2),c3/(c3+c4),c4/(c3+c4))
timepoint <- c("Init.","Init.","Rec.","Rec.")
status <- c("Ki67+", "Ki67-","Ki67+","Ki67-")
plot_mut <- data.frame(fraction,timepoint,status)

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions_final/Fig3/idhmut_ki67.pdf",width = 1,height=1.6)
ggplot(plot_mut, aes(fill=status, y=fraction*100, x=timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("Ki67+" = "#8B8B8B", "Ki67-" = "#ECECEC")) +
  labs(y = "Proportion of SOX2+ cells (%)") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7,hjust=0.5),
        axis.text.y=element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "none")
dev.off()

# IDHwt
imm1 <-   read.csv("data/if/B13_prim_gated.csv")
imm2 <-   read.csv("data/if/B13_rec_gated.csv")

c1 <- nrow(imm1 %>% filter(Sox2 >= thr["Sox2"], Ki67 >= thr["Ki67"]))
c2 <- nrow(imm1 %>% filter(Sox2 >= thr["Sox2"], Ki67 < thr["Ki67"]))
c3 <- nrow(imm2 %>% filter(Sox2 >= thr["Sox2"], Ki67 >= thr["Ki67"]))
c4 <- nrow(imm2 %>% filter(Sox2 >= thr["Sox2"], Ki67 < thr["Ki67"]))                      

wt_ct <- matrix(c(c1,c2,c3,c4),nrow=2)
fisher.test(wt_ct)$p.value

fraction <- c(c1/(c1+c2),c2/(c1+c2),c3/(c3+c4),c4/(c3+c4))
timepoint <- c("Init.","Init.","Rec.","Rec.")
status <- c("Ki67+", "Ki67-","Ki67+","Ki67-")
plot_wt <- data.frame(fraction,timepoint,status)

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions_final/Fig3/idhwt_ki67.pdf",width = 1,height=1.6)
ggplot(plot_wt, aes(fill=status, y=fraction*100, x=timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("Ki67+" = "#8B8B8B", "Ki67-" = "#ECECEC")) +
  labs(y = "Proportion of SOX2+ cells (%)") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7,hjust=0.5),
        axis.text.y=element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "none")
dev.off()



# IDHmut
imm1 <-   read.csv("data/if/S16_IDH_mut.csv")
imm2 <-   read.csv("data/if/S19_IDH_mut.csv")

imm1_sox <- imm1 %>% filter(Sox2 >= thr["Sox2"])
imm2_sox <- imm2 %>% filter(Sox2 >= thr["Sox2"])

wilcox.test(imm1$Ki67, imm2$Ki67)

mut_init <- imm1_sox[,"Ki67"]
mut_rec <- imm2_sox[,"Ki67"]
mut_ki67 <- c(mut_init, mut_rec)
timepoint <- c(rep("Init.", length(imm1_sox[,"Ki67"])), rep("Rec.", length(imm2_sox[,"Ki67"])))
plot_mut <- data.frame(mut_ki67, timepoint)
  
pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions_final/Fig3/idhmut_ki67.pdf",width = 1,height=1.5)
ggplot(plot_mut, aes(x=timepoint,y=mut_ki67)) +
  geom_violin(colour="black",size=1) +
  geom_boxplot(colour="black",outlier.shape =NA, width = 0.2) +
  #geom_smooth(method = lm, se = FALSE, aes(colour = cell_state)) + 
  theme_classic() +
  labs(y = "Ki67 expression (SOX2+)") +
  theme(axis.text = element_text(size=7),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=7), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=7),
        legend.position="none") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1), limits = c(0, 50))
dev.off()

# IDHwt
imm1 <-   read.csv("data/if/B13_prim_gated.csv")
imm2 <-   read.csv("data/if/B13_rec_gated.csv")

imm1_sox <- imm1 %>% filter(Sox2 >= thr["Sox2"])
imm2_sox <- imm2 %>% filter(Sox2 >= thr["Sox2"])

wilcox.test(imm1$Ki67, imm2$Ki67)

wt_init <- imm1_sox[,"Ki67"]
wt_rec <- imm2_sox[,"Ki67"]
wt_ki67 <- c(wt_init, wt_rec)
timepoint <- c(rep("Init.", length(imm1_sox[,"Ki67"])), rep("Rec.", length(imm2_sox[,"Ki67"])))
plot_wt <- data.frame(wt_ki67, timepoint)

pdf("/Users/varnf/Desktop/Figs/GLASS-III/revisions_final/Fig3/idhwt_ki67.pdf",width = 1,height=1.5)
ggplot(plot_wt, aes(x=timepoint,y=wt_ki67)) +
  geom_violin(colour="black",size=1) +
  geom_boxplot(colour="black",outlier.shape =NA, width = 0.2) +
  #geom_smooth(method = lm, se = FALSE, aes(colour = cell_state)) + 
  theme_classic() +
  labs(y = "Ki67 expression (SOX2+)") +
  theme(axis.text = element_text(size=7),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=7), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=7),
        legend.position="none") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1), limits = c(0, 50))
dev.off()

