library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

myinf1 <- "data/immunofluorescence/snu_single_positive_percentages.txt"
myinf2 <- "data/immunofluorescence/snu_gated_percentages.txt"
myinf3 <- "results/cibersortx/hires/GLASS-freeze/CIBERSORTxGEP_GLASS_Fractions-Adjusted.txt"

sp <- read.delim(myinf1)
gp <- read.delim(myinf2)
frac <- read.delim(myinf3)

gp[,2:ncol(gp)] <- gp[,2:ncol(gp)]/100
sp[,2:ncol(sp)] <- sp[,2:ncol(sp)]/100

frac[,"Mixture"] <- gsub("\\.","-",frac[,"Mixture"])
frac[,"Mixture"] <- sapply(strsplit(frac[,"Mixture"],"-"),function(x)paste(x[1:4],collapse="-"))

frac <- frac %>%
		inner_join(sp, by = c("Mixture" = "sample_barcode")) %>%
		inner_join(gp, by = c("Mixture" = "sample_barcode")) %>%
		select(-c("P.value","Correlation","RMSE"))
	
########		
tmp <- frac %>% column_to_rownames("Mixture")
stemdif <- frac[,"stemcell_tumor"] - frac[,"stem_like_sox2_olig2"]
diffdif <- frac[,"differentiated_tumor"] - frac[,"Sox2_olig2neg"]
prolifdif <- frac[,"prolif_stemcell_tumor"] - (frac[,"prolif_stem_like_sox2_olig2_ki67"] * frac[,"stem_like_sox2_olig2"])
myeloid_dif <- frac[,"myeloid"] - frac[,"CD14"]

diff_frame <- data.frame(stemdif, diffdif, prolifdif, myeloid_dif)
rownames(diff_frame) <- rownames(tmp)

 avg_diff <- apply(diff_frame, 1, function(x) mean(abs(x)))
 avg_diff <- avg_diff[order(avg_diff)]

cor(frac)
#frac <- frac[-which(frac[,"Mixture"] == "GLSS-SN-0009-R1"),]
select <- c("GLSS-SN-0010-R1", "GLSS-SN-0016-R1", "GLSS-SN-0016-TP", "GLSS-SN-0017-TP", "GLSS-SN-0017-R1", "GLSS-SN-0015-R2","GLSS-SN-0009-R1")

cor(frac[,"differentiated_tumor"], frac[,"Sox2_olig2neg"])
cor(frac[,"stemcell_tumor"],frac[,"stem_like_sox2_olig2"])
cor(frac[,"prolif_stemcell_tumor"],frac[,"prolif_stem_like_sox2_olig2_ki67"] * frac[,"stem_like_sox2_olig2"])
cor(frac[,"myeloid"], frac[,"CD14"])
cor(frac[,"differentiated_tumor"] + frac[,"stemcell_tumor"] + frac[,"prolif_stemcell_tumor"], frac[,"Sox2"])

###########

renorm_tumor <- t(apply(frac[,c("stemcell_tumor", "differentiated_tumor", "prolif_stemcell_tumor")], 1, function(x) x/sum(x)))
rownames(renorm_tumor) <- rownames(tmp)
renorm_tumor <- data.frame(renorm_tumor, tmp[,13:ncol(tmp)])

cor(renorm_tumor[,"differentiated_tumor"], renorm_tumor[,"Sox2_olig2neg"])
cor(renorm_tumor[,"stemcell_tumor"],renorm_tumor[,"stem_like_sox2_olig2"])
cor(renorm_tumor[,"prolif_stemcell_tumor"],renorm_tumor[,"prolif_stem_like_sox2_olig2_ki67"] * renorm_tumor[,"stem_like_sox2_olig2"])

# tp <- frac[grep("-TP", frac$Mixture),]
# cor(tp[,2:ncol(tp)]
# 
# rec <- frac[-grep("-TP", frac$Mixture),]
# cor(rec[,2:ncol(rec)]

imm <- c(frac$Sox2_olig2neg, frac$stem_like_sox2_olig2, frac$prolif_stem_like_sox2_olig2_ki67 * frac[,"stem_like_sox2_olig2"], frac$CD14, frac$Sox2)
cibersortx <- c(frac$differentiated_tumor, frac$stemcell_tumor, frac$prolif_stemcell_tumor, frac$myeloid, frac$differentiated_tumor + frac$stemcell_tumor + frac$prolif_stemcell_tumor)
comp <- rep(c("Differentiated-like", "Stem-like", "Prolif. stem-like", "Myeloid", "Tumor"), each = nrow(frac))
sample_barcode <- rep(frac$Mixture, 5)
timepoint <- ifelse(grepl("-TP", sample_barcode), "Initial", "Recurrent")
plot_dat <- data.frame(sample_barcode, imm, cibersortx, comp, timepoint)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/snu_if_cibersortx.pdf", height = 3, width = 4)
ggplot(data = plot_dat, aes(x = imm*100, y = cibersortx*100)) +
  geom_point(aes(shape = timepoint)) +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Immunofluorescence (%)", y = "CIBERSORTx (%)") +
  facet_wrap(.~comp, ncol = 3, scales = "free") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position = "none")
dev.off()

sub_dat <- plot_dat %>% filter(comp != "Tumor")
cor(sub_dat$imm, sub_dat$cibersortx)


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/snu_if_together_cibersortx.pdf", height = 2, width = 2)
ggplot(data = sub_dat, aes(x = imm*100, y = cibersortx*100)) +
  geom_point(aes(fill = comp), colour="black",pch=21,size=2) +
  geom_smooth(method="lm",se = FALSE) +
  scale_fill_manual(values = c("Differentiated-like" = "#fcbba1", "Stem-like" = "#fb6a4a", "Prolif. stem-like" = "#a50f15","Myeloid" = "#08519c")) +
  labs(x = "Immunofluorescence (%)", y = "CIBERSORTx (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position = "none")
dev.off()


tmp <- data.frame(frac[,"differentiated_tumor"], frac[,"Sox2"] - frac[,"stem_like_sox2_olig2"] - frac[,"prolif_stem_like_sox2_olig2_ki67"])
colnames(tmp) <- c("V1","V2")
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/diff_olig2_subtrac_cibersortx.pdf", height = 2, width = 2)
ggplot(data = tmp, aes(x = V2, y = V1*100)) +
  geom_point() +
  geom_smooth(method="lm",se = FALSE) +
  #scale_fill_manual(values = c("Differentiated-like" = "#fcbba1", "Stem-like" = "#fb6a4a", "Prolif. stem-like" = "#a50f15","Myeloid" = "#08519c")) +
  labs(x = "Immunofluorescence (%)", y = "CIBERSORTx (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position = "none")
dev.off()


