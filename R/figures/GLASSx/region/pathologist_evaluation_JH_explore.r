library(odbc)
library(DBI)
library(tidyverse)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in aggregate pathologist annotations
myinf1 <- "data/pathology_annot/Histo_review_aggregate_v2.txt"

dat <- read.delim(myinf1)

dat_long <- dat %>%
  select(-Notes) %>%
  pivot_longer(-c(CASES,Pathologist),names_to ="feature", values_to = "metric")

feat_class <- rep("Primary",nrow(dat_long))
feat_class[grepl("NecTreatEff|Fibrosis|DepopTum", dat_long$feature)] <- "Recurrent"
dat_long$feat_class <- feat_class

metric_class <- rep("Percent",nrow(dat_long))
metric_class[grepl("_score", dat_long$feature)] <- "Score"
dat_long$metric_class <- metric_class

dat_long <- dat_long %>%
  mutate(feature = fct_relevel(feature, unique(dat_long$feature)))

dat_comp <- dat_long %>%
  select(-c(Pathologist, feat_class, metric_class)) %>%
  pivot_wider(names_from = feature, values_from = metric)

# Transcriptomic analysis

# Read in IvyGAP table
ivy <- dbReadTable(con, Id(schema="analysis",table="cibersortx_ivygap"))


ivy$aliquot_barcode <- sapply(strsplit(ivy$aliquot_barcode,"-"),function(x)paste(x[1:4],collapse="-"))
ivy <- ivy %>% filter(aliquot_barcode %in% unique(dat_long$CASES))

ivy$fraction <- ivy$fraction * 100

ivy_comp <- ivy %>%
  pivot_wider(names_from = cell_state, values_from = fraction)

# Compare both metrics
full_comp <- dat_comp %>%
  inner_join(ivy_comp, by = c("CASES" = "aliquot_barcode")) %>%
  column_to_rownames("CASES")
full_comp[,"all_CT"]  <- apply(full_comp[,c("CT","CTpan","CTmvp")], 1, sum)
cor(full_comp, method="s")
heatmap_plot <- cor(full_comp, method="p")

# Plot the full comp heatmap
heatmap_plot <- heatmap_plot %>% 
  data.frame() %>%
  rownames_to_column(var = "path") %>%
  select(path, CT, CTmvp, CTpan, LE) %>%
  pivot_longer(-path, names_to = "ivygap", values_to = "cor") %>% 
  filter(!path %in% c("CT","CTmvp","CTpan","LE")) %>%
  mutate(path = fct_relevel(path, rev(unique(path)))) %>%
  mutate(ivygap = fct_relevel(ivygap, "LE","CT","CTpan","CTmvp"))

pdf("figures/pathology/heatmap_cor_JH.pdf", height = 6, width = 6)
ggplot(data = heatmap_plot, aes(x = ivygap, y = path, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0) +
  labs(x = "Transcriptomic", y = "Pathology") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
dev.off()


#-------------------


# Compare both metrics
full_comp <- dat_comp %>%
  inner_join(ivy_comp, by = c("CASES" = "aliquot_barcode")) %>%
  #filter(grepl("-HF-",CASES)) %>%
  column_to_rownames("CASES")
full_comp[,"all_CT"]  <- apply(full_comp[,c("CT","CTpan","CTmvp")], 1, sum)
cor(full_comp, method="p")
heatmap_plot <- cor(full_comp, method="p")

# Plot the full comp heatmap
heatmap_plot <- heatmap_plot %>% 
  data.frame() %>%
  rownames_to_column(var = "path") %>%
  select(path, CT, CTmvp, CTpan, LE) %>%
  pivot_longer(-path, names_to = "ivygap", values_to = "cor") %>% 
  filter(!path %in% c("CT","CTmvp","CTpan","LE")) %>%
  mutate(path = fct_relevel(path, rev(unique(path)))) %>%
  mutate(ivygap = fct_relevel(ivygap, "LE","CT","CTpan","CTmvp"))

pdf("figures/pathology/heatmap_cor_HF_JH.pdf", height = 6, width = 6)
ggplot(data = heatmap_plot, aes(x = ivygap, y = path, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0) +
  labs(x = "Transcriptomic", y = "Pathology") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
dev.off()





# Composite correlation plot: PAN vs necrosis
ivy_frac <- c(full_comp[,"LE"], full_comp[,"CT"], full_comp[,"CTpan"])
path_frac <- c(full_comp[,"LE_percent"], full_comp[,"CT_percent"],full_comp[,"Necrosis_percent"])
feat <- rep(c("LE","CT","PAN"), each = nrow(full_comp))

cor_plot <- data.frame(feat, ivy_frac, path_frac)
cor(cor_plot[,3], cor_plot[,2])

pdf("figures/pathology/feature_cor_with_nec_JH.pdf", height = 3, width = 4)
ggplot(data = cor_plot, aes(x = path_frac, y = ivy_frac)) +
  geom_point(aes(colour = feat)) +
  geom_smooth(method="lm",) +
  labs(x = "Pathology (%)", y = "Transcriptomic (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
dev.off()


# Composite correlation plot: PAN vs palisading necrosis
ivy_frac <- c(full_comp[,"LE"], full_comp[,"CT"], full_comp[,"CTpan"])
path_frac <- c(full_comp[,"LE_percent"], full_comp[,"CT_percent"],full_comp[,"PN_percent"])
feat <- rep(c("LE","CT","PAN"), each = nrow(full_comp))

cor_plot <- data.frame(feat, ivy_frac, path_frac)
cor(cor_plot[,3], cor_plot[,2])
#pdf("figures/pathology/feature_cor_with_pn.pdf", height = 3, width = 4)
ggplot(data = cor_plot, aes(x = path_frac, y = ivy_frac)) +
  geom_point(aes(colour = feat)) +
  geom_smooth(method="lm",) +
  labs(x = "Pathology (%)", y = "Transcriptomic (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
#dev.off()

# Composite correlation plot: PAN vs total necrosis sum
ivy_frac <- c(full_comp[,"LE"], full_comp[,"CT"], full_comp[,"CTpan"])
necrosis_sum <- apply(full_comp[,c("PN_percent", "Necrosis_percent")], 1, sum)
path_frac <- c(full_comp[,"LE_percent"], full_comp[,"CT_percent"],necrosis_sum)
feat <- rep(c("LE","CT","PAN"), each = nrow(full_comp))

cor_plot <- data.frame(feat, ivy_frac, path_frac)
cor(cor_plot[,3], cor_plot[,2])
#pdf("figures/pathology/feature_cor_with_necsum.pdf", height = 3, width = 4)
ggplot(data = cor_plot, aes(x = path_frac, y = ivy_frac)) +
  geom_point(aes(colour = feat)) +
  geom_smooth(method="lm",) +
  labs(x = "Pathology (%)", y = "Transcriptomic (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
#dev.off()

# Split this out by primary vs recurrent
status <- sapply(strsplit(rownames(full_comp),"-"),function(x)x[4])
ivy_frac <- c(full_comp[,"LE"], full_comp[,"CT"], full_comp[,"CTpan"])
path_frac <- c(full_comp[,"LE_percent"], full_comp[,"CT_percent"],full_comp[,"Necrosis_percent"])
feat <- rep(c("LE","CT","PAN"), each = nrow(full_comp))
status <- rep(status, 3)

cor_plot <- data.frame(feat, ivy_frac, path_frac,status)
cor_plot <- cor_plot %>% 
  mutate(status = recode(status, "R1" = "Recurrent", "R2" = "Recurrent"))
  mutate(status = fct_relevel(status, "TP","Recurrent"))
pri_plot <- cor_plot %>% filter(status == "TP")
rec_plot <- cor_plot %>% filter(status == "Recurrent")
cor(pri_plot[,3], pri_plot[,2])
cor(rec_plot[,3], rec_plot[,2])

pdf("figures/pathology/feature_cor_with_necsum_by_status.pdf", height = 3, width = 8)
ggplot(data = cor_plot, aes(x = path_frac, y = ivy_frac)) +
  geom_point(aes(colour = feat)) +
  geom_smooth(method="lm",) +
  facet_grid(.~status) +
  labs(x = "Pathology (%)", y = "Transcriptomic (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
dev.off()

# Visualizations of each feature

comp_plot <- full_comp[,c("CT_percent","CT", "LE_percent","Fibrosis_percent", "LE", "PN_percent", "Necrosis_percent", "CTpan")]
comp_plot <- comp_plot %>%
  rownames_to_column("sample_barcode") %>%
  pivot_longer(-sample_barcode, names_to = "feature", values_to = "fraction")

feature_class <- rep("", nrow(comp_plot))
feature_class[grep("CT", comp_plot$feature)] <- "CT"
feature_class[grep("LE", comp_plot$feature)] <- "LE"
feature_class[grep("Fibrosis_percent", comp_plot$feature)] <- "LE"
feature_class[grep("CTpan", comp_plot$feature)] <- "Nec"
feature_class[grep("PN_percent", comp_plot$feature)] <- "Nec"
feature_class[grep("Necrosis_percent", comp_plot$feature)] <- "Nec"
comp_plot$feature_class <- feature_class

comp_plot$feature <- gsub("_percent", "_path", as.character(comp_plot$feature))
comp_plot$status <- sapply(strsplit(comp_plot$sample_barcode, "-"),function(x)x[4])
comp_plot <- comp_plot %>%
  mutate(status = recode(status, "R1" = "Recurrent", "R2" = "Recurrent"))

pdf("figures/pathology/feature_by_feature_tp.pdf", height = 5, width = 5)
ggplot(data = comp_plot %>% filter(status == "TP"), aes(x = feature, y = fraction)) +
  geom_boxplot(colour="black") +
  geom_point(size=1) +
  facet_grid(.~feature_class,scales = "free_x") +
  scale_colour_manual(values="#619CFF") +
  labs(y = "Proportion (%)") +
  theme_classic() +
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position = "none") +
  coord_cartesian(ylim=c(0,100))
dev.off()


pdf("figures/pathology/feature_by_feature_r1.pdf", height = 5, width = 5)
ggplot(data = comp_plot %>% filter(status == "R1"), aes(x = feature, y = fraction)) +
  geom_boxplot(colour="black") +
  geom_point(size=1) +
  facet_grid(.~feature_class,scales = "free_x") +
  scale_colour_manual(values="#619CFF") +
  labs(y = "Proportion (%)") +
  theme_classic() +
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title.x= element_blank(),
        axis.title.y= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position = "none") +
  coord_cartesian(ylim=c(0,100))
dev.off()

# Category plot for MVP
cat_plot <- full_comp[,c("CT_score","CT", "LE_score","LE", "Necrosis_score", "PN_score", "CTpan", "MVP_score", "HBV_score", "CTmvp")]
cat_plot <- cat_plot %>%
  rownames_to_column("sample_barcode") %>%
  pivot_longer(cols = ends_with("_score"), names_to = "feature", values_to = "score")

feature_class <- rep("", nrow(cat_plot))
feature_class[grep("CT", cat_plot$feature)] <- "CT"
feature_class[grep("LE", cat_plot$feature)] <- "LE"
feature_class[grep("CTpan", cat_plot$feature)] <- "Nec"
feature_class[grep("PN_", cat_plot$feature)] <- "Nec"
feature_class[grep("Necrosis_", cat_plot$feature)] <- "Nec"
feature_class[grep("CTmvp", cat_plot$feature)] <- "MVP"
feature_class[grep("MVP_", cat_plot$feature)] <- "MVP"
feature_class[grep("HBV_", cat_plot$feature)] <- "MVP"
cat_plot$feature_class <- feature_class

pdf("figures/pathology/mvp_score.pdf", height = 5, width = 5)
ggplot(data = cat_plot %>% filter(feature == "MVP_score"), aes(x = score, y = CTmvp)) +
  geom_point(colour="black") +
  geom_smooth(method="lm",) +
  facet_grid(.~feature_class,scales = "free_x") +
  labs(x = "Pathology (Score)", y = "Transcriptomic (%)") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position = "none") 
dev.off()


# Sandbox

jh <- dat_long %>% filter(Pathologist == "JH", feature == "HBV_score")
head(jh)
ivy

ivy %>%
  inner_join(jh, c("aliquot_barcode" = "CASES")) %>%
  filter(cell_state == "CTmvp")

full_comp[,"mes"] <- apply(full_comp[,c("CTpan","CTmvp")], 1, sum)

tp_comp <- cor(full_comp[grep("-TP",rownames(full_comp)),])
r1_comp <- cor(full_comp[grep("-R1",rownames(full_comp)),])

plot(full_comp[grep("-TP",rownames(full_comp)),c("LE_score","LE")])
plot(full_comp[grep("-TP",rownames(full_comp)),c("CT_score","all_CT")])
plot(full_comp[grep("-TP",rownames(full_comp)),c("Necrosis_score","mes")])


plot(full_comp[grep("-R1",rownames(full_comp)),c("LE_score","LE")])
plot(full_comp[grep("-R1",rownames(full_comp)),c("CT_score","all_CT")])
plot(full_comp[grep("-R1",rownames(full_comp)),c("Necrosis_score","mes")])


plot(full_comp[grep("-TP",rownames(full_comp)),c("LE_percent","LE")])
plot(full_comp[grep("-TP",rownames(full_comp)),c("CT_percent","all_CT")])
plot(full_comp[grep("-TP",rownames(full_comp)),c("Necrosis_percent","mes")])


plot(full_comp[grep("-R1",rownames(full_comp)),c("LE_percent","LE")])
plot(full_comp[grep("-R1",rownames(full_comp)),c("CT_percent","all_CT")])
plot(full_comp[grep("-R1",rownames(full_comp)),c("Necrosis_percent","mes")])

plot(full_comp[,c("LE_percent","LE")])
plot(full_comp[,c("CT_percent","all_CT")])
plot(full_comp[,c("Necrosis_percent","mes")])

plot(full_comp[,c("LE_score","LE")])
plot(full_comp[,c("CT_score","all_CT")])
plot(full_comp[,c("Necrosis_score","mes")])

thr1 <- median(full_comp$LE_score)
g1 <- full_comp[which(full_comp$LE_score > thr1),"LE"]
g2 <- full_comp[which(full_comp$LE_score <= thr1),"LE"]
wilcox.test(g1,g2)

thr2 <- median(full_comp$CT_score)
g1 <- full_comp[which(full_comp$CT_score > thr2),"all_CT"]
g2 <- full_comp[which(full_comp$CT_score <= thr2),"all_CT"]
wilcox.test(g1,g2)

thr3 <- median(full_comp$Necrosis_score)
g1 <- full_comp[which(full_comp$Necrosis_score > thr3),"CTpan"]
g2 <- full_comp[which(full_comp$Necrosis_score <= thr3),"CTpan"]
wilcox.test(g1,g2)



thr1 <- median(full_comp$LE_percent)
g1 <- full_comp[which(full_comp$LE_percent > thr1),"LE"]
g2 <- full_comp[which(full_comp$LE_percent <= thr1),"LE"]
wilcox.test(g1,g2)

thr2 <- median(full_comp$CT_percent)
g1 <- full_comp[which(full_comp$CT_percent > thr2),"CT"]
g2 <- full_comp[which(full_comp$CT_percent <= thr2),"CT"]
wilcox.test(g1,g2)

thr3 <- median(full_comp$Necrosis_percent)
g1 <- full_comp[which(full_comp$Necrosis_percent > thr3),"CTpan"]
g2 <- full_comp[which(full_comp$Necrosis_percent <= thr3),"CTpan"]
wilcox.test(g1,g2)

