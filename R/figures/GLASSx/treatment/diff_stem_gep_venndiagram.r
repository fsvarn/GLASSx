###################################################
# Identify overlap of gene signatures, Venn diagrams
# Updated: 2020.03.21
# Author: Frederick Varn
##################################################

library(tidyverse)
library(VennDiagram)

#######################################################
rm(list=ls())

myinf1 <- paste("data/res/CIBERSORTx/analysis/GLASS_idhwt_differentiated_tumor_postreatment_result.txt",sep="")
myinf2 <- paste("data/res/CIBERSORTx/analysis/GLASS_idhwt_stemcell_tumor_postreatment_result.txt",sep="")

dat1 <- read.delim(myinf1)
dat2 <- read.delim(myinf2)

diff_upsig <- rownames(dat1 %>% filter(sig, eff > 0))
diff_dnsig <- rownames(dat1 %>% filter(sig, eff < 0))
stem_upsig <- rownames(dat2 %>% filter(sig, eff > 0))
stem_dnsig <- rownames(dat2 %>% filter(sig, eff < 0))

#opc <- c("OLIG1", "OMG", "PLP1", "PLLP", "TNR", "ALCAM")

overlap_up <- intersect(diff_upsig, stem_upsig)
overlap_dn <- intersect(diff_dnsig, stem_dnsig)
overlap <- unique(c(overlap_up, overlap_dn))

diffsig <- c(diff_upsig, diff_dnsig)
diffsig <- diffsig[-which(diffsig %in% overlap)]
stemsig <- c(stem_upsig, stem_dnsig)
stemsig <- stemsig[-which(stemsig %in% overlap)]


myCol <- c("#fcbba1","#fb6a4a")

venn.diagram(
		x = list(diff_sig, stem_sig),
		category.names = c("Diff.-like" , "Stem-like "),
		filename = '/projects/verhaak-lab/GLASS-III/figures/analysis/diff_stem_venndiagram.tiff',
		output=TRUE,
		
        # Output features
        imagetype="tiff",
        height = 1.5 , 
        width = 1.5 , 
        resolution = 600,
        compression = "lzw",
        units = "in",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .6,
        fontfamily = "Helvetica",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
       # cat.pos = c(-27, 27, 135),
       # cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "Helvetica"
)


#######################################################
rm(list=ls())

myinf1 <- paste("data/res/CIBERSORTx/analysis/GLASS_idhmut_differentiated_tumor_postreatment_result.txt",sep="")
myinf2 <- paste("data/res/CIBERSORTx/analysis/GLASS_idhmut_stemcell_tumor_postreatment_result.txt",sep="")

dat1 <- read.delim(myinf1)
dat2 <- read.delim(myinf2)

diff_upsig <- rownames(dat1 %>% filter(sig, eff > 0))
diff_dnsig <- rownames(dat1 %>% filter(sig, eff < 0))
stem_upsig <- rownames(dat2 %>% filter(sig, eff > 0))
stem_dnsig <- rownames(dat2 %>% filter(sig, eff < 0))

overlap_up <- intersect(diff_upsig, stem_upsig)
overlap_dn <- intersect(diff_dnsig, stem_dnsig)
overlap <- unique(c(overlap_up, overlap_dn))

diffsig <- c(diff_upsig, diff_dnsig)
diffsig <- diffsig[-which(diffsig %in% overlap)]
stemsig <- c(stem_upsig, stem_dnsig)
stemsig <- stemsig[-which(stemsig %in% overlap)]


myCol <- c("#fcbba1","#fb6a4a")

venn.diagram(
		x = list(diff_sig, stem_sig),
		category.names = c("Diff.-like" , "Stem-like "),
		filename = '/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_diff_stem_venndiagram.tiff',
		output=TRUE,
		
        # Output features
        imagetype="tiff",
        height = 1.5 , 
        width = 1.5 , 
        resolution = 600,
        compression = "lzw",
        units = "in",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .6,
        fontfamily = "Helvetica",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
       # cat.pos = c(-27, 27, 135),
       # cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "Helvetica"
)
