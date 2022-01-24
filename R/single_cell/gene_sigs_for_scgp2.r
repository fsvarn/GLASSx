library(tidyverse)

rm(list=ls())

#Step 1: Identify neuron and astrocyte-specific genes in the Darmanis et al dataset

myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/darmanis_2017/GBM_normalized_gene_counts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/darmanis_2017/GBM_metadata.csv"

data <- read.delim(myinf1,sep=" ")
colnames(data) <- gsub("^X","",colnames(data))
info <- read.delim(myinf2,sep=" ")

# Define neuron signature
neurons <- rownames(info %>% filter(Selection == "Neurons(Thy1)"))
non_neurons <- colnames(data)[-which(colnames(data) %in% neurons)]
p.val <- apply(data, 1, function(x)t.test(x[neurons],x[non_neurons])$p.value)
eff <- apply(data, 1, function(x)mean(x[neurons])/median(x[non_neurons]))
q.val <- p.adjust(p.val, "BH")

neuron_genes <- data.frame(p.val, eff, q.val)
neuron_genes <- neuron_genes[order(neuron_genes$q.val),]

# Define astrocyte signature
astrocytes <- rownames(info %>% filter(Selection == "Astrocytes(HEPACAM)"))
non_astrocytes <- colnames(data)[-which(colnames(data) %in% astrocytes)]
p.val <- apply(data, 1, function(x)t.test(x[astrocytes],x[non_astrocytes])$p.value)
eff <- apply(data, 1, function(x)mean(x[astrocytes])/median(x[non_astrocytes]))
q.val <- p.adjust(p.val, "BH")

astrocyte_genes <- data.frame(p.val, eff, q.val)