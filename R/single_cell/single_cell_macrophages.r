###################################################
# Define GBM-specific macrophage programs from single cell data
# Reference for methods and dataset: An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma: Cell
# Updated: 2020.05.07
# Author: Frederick Varn
##################################################

library(amap)
library(tidyverse)
library(broom)
library(scran)

rm(list=ls())
myinf1 <- "/projects/varnf/Data/Neftel_IDHwt_GBM/IDHwtGBM.processed.SS2.logTPM.txt"
myinf2 <- "/projects/varnf/Data/Neftel_IDHwt_GBM/IDHwt.GBM.Metadata.SS2.txt"

# Loads expression matrix where values = log2((tpm/10) + 1)
tpm <- read.delim(myinf1,row.names=1)
colnames(tpm) <- gsub("\\.","-",colnames(tpm))

# Remove malignant cells, only interested in macrophages
info <- read.delim(myinf2)
mac_info <- info[which(info[,"CellAssignment"] == "Macrophage"),]

mac_tpm <- tpm[,mac_info[,"NAME"]]
patient <- unique(sapply(substr(colnames(mac_tpm),1,6),function(x)x[1]))
patient <- sub("-","",patient)

# Calculate aggregate gene expression to determine what genes to include (see methods)
agg_expr <- apply(mac_tpm,1,function(x)log2(mean(((2^x) -1) * 10) + 1))
mac_filt_tpm <- mac_tpm[names(agg_expr[which(agg_expr > 4)]),]

# Calculate relative expression by centering expression levels per gene about mean (no need)
mean_tpm <- apply(mac_filt_tpm, 1, mean)
mac_rel_tpm <- mac_filt_tpm - mean_tpm


##################################################
# Step 1: Cluster the single macrophage cells from each patient separately
##################################################
 
patient_clusters <- list()
for(i in 1:length(patient))
{
	cat("\r",i)
	if(length(grep(patient[i],colnames(mac_rel_tpm))) == 1){
		next
	}
	sub_mac <- mac_rel_tpm[,grep(patient[i],colnames(mac_rel_tpm))]
	
	# Average linkage clustering using 1 - Pearson distance (see methods)
	tree <- hclust(d = Dist(x = t(sub_mac), method = "pearson"), method = "average")
	patient_clusters[[i]] <- tree
}
names(patient_clusters) <- patient
patient_clusters[sapply(patient_clusters, is.null)] <- NULL

##################################################
# Step 2: Obtain all possible clusters
##################################################

cluster_id <- cell_id <- patient_id <- c()
for(i in 1:length(patient_clusters))
{
	mycluster <- patient_clusters[[i]]
	mypatient <- names(patient_clusters)[i]
	
	# Get total number of cells
	n_objects <- length(mycluster[["order"]])
	
	pat_clust_id <- 0
	
	# This loop makes a cut in the tree to get j number of clusters
	# The number of objects is equal to the total number of possible clusters
	for(j in 1:n_objects)
	{
		sub_cluster <- cutree(mycluster, k=j)
		cluster_counts <- table(sub_cluster)
		
		# Filter out all clusters that contain more than 80% of cells in a patient or fewer than 5 cells total (see methods)
		cluster_save <- cluster_counts[which(cluster_counts/n_objects <= 0.8 & cluster_counts >= 5)]
		if(length(cluster_save) == 0){
			next}
		
		# Vector of clusters/patients that qualify
		sub_cluster <- sub_cluster[which(sub_cluster %in% names(cluster_save))]
		
		# Create new cluster ID so that each cluster is treated separately		
		add_clust <- 1:length(cluster_save)
		new_clust_id <- pat_clust_id + add_clust
		pat_clust_id <- max(new_clust_id)

		# Create vector with new cluster ID by converting old cluster ID to factors and then reverting
		cluster_id <- c(cluster_id, as.numeric(as.character(factor(sub_cluster, labels =as.character(new_clust_id)))))
		cell_id <- c(cell_id, names(sub_cluster))
		patient_id <- c(patient_id,rep(mypatient,length(sub_cluster)))
	}

}

# Create a long format table with the cells in each patient and their respective cluster
cell_clusters <- data.frame(patient_id, cell_id, cluster_id, stringsAsFactors=FALSE)
cell_clusters <- cell_clusters[order(cell_clusters[,"patient_id"],cell_clusters[,"cluster_id"], cell_clusters[,"cell_id"]),]

##################################################
# Step 3: Perform differential expression analysis to identify signature genes for each cluster
##################################################

# Create an index for looping
index <- cell_clusters %>%
	select(patient_id, cluster_id) %>%
	distinct()

sig <- patient_id <- cluster_id <- c()
for(i in 1:nrow(index))
{	
	cat(i, "out of ", nrow(index),"\r")

	# Combine cell_clusters and macrophage tpm
	in_cells <- cell_clusters[which(
	cell_clusters[,"patient_id"] == index[i,"patient_id"] &
	cell_clusters[,"cluster_id"] == index[i,"cluster_id"]),"cell_id"]
	
	sub_mac <- data.matrix(mac_rel_tpm[,grep(index[i,"patient_id"],colnames(mac_rel_tpm))])
	group <- as.character(colnames(sub_mac) %in% in_cells) %>%
		recode("FALSE"="out","TRUE"="in")
	
	in_genes <- as.data.frame(findMarkers(sub_mac, group, test.type="t")[["in"]]) %>%
		rownames_to_column(var = "gene_symbol") %>%
		as_tibble() %>%
		filter(logFC.out > 3, FDR < 0.05)
	
	# Only include clusters that have:
	# > 50 genes with FDR < 0.05 + fold-change > 3 AND 
	# > 10 genes with FDR < 0.005 + fold-change > 3
	#----------------------------------------------	
	if(nrow(in_genes) > 50 & nrow(in_genes %>% filter(FDR < 0.005)) > 10){
			sig <- c(sig, pull(in_genes, gene_symbol))
			cluster_id <- c(cluster_id, rep(index[i,"cluster_id"], nrow(in_genes)))
			patient_id <- c(patient_id, rep(index[i,"patient_id"], nrow(in_genes)))
	}
}


cluster_genes <- data.frame(patient_id, cluster_id, sig,stringsAsFactors=FALSE)

##################################################
# Step 4: Calculate pairwise Jaccard index to remove redundant clusters
##################################################

# Calculate pairwise jaccard index

# Create an index for looping
index <- cluster_genes %>%
	select(patient_id, cluster_id) %>%
	distinct() 
	
index_names	<- index %>%
	unite(index, sep = "_") %>%
	pull(index)
	
jaccard <- matrix(NA, nrow = length(index_names), ncol = length(index_names))
rownames(jaccard) <- colnames(jaccard) <- index_names
for(i in 1:nrow(index))
{
	cluster1 <- cluster_genes %>%
	filter(patient_id == index[i,1], cluster_id == index[i,2]) %>%
	pull(sig) %>%
	as.character()
	
	for(j in 1:i)
	{
		cat(i, "-->", j,"\r")
				
		cluster2 <- cluster_genes %>%
			filter(patient_id == index[j,1], cluster_id == index[j,2]) %>%
			pull(sig) %>%
			as.character()
			
		jaccard[i,j] <- length(intersect(cluster1, cluster2))/length(union(cluster1, cluster2))
	}
}
	
# Identify signatures that had no overlap in other patients:
diag(jaccard) <- 0
sum_loss1 <- apply(jaccard, 1, function(x)sum(x>0.75,na.rm=TRUE))
sum_loss2 <- apply(jaccard, 2, function(x)sum(x>0.75,na.rm=TRUE))
solo_sigs <- data.frame(sum_loss1, sum_loss2) %>%
	rownames_to_column(var = "sig") %>%
	filter(sum_loss1 == 0 & sum_loss2 == 0) %>%
	select(sig) %>%
	separate(sig, "_", into = c("pat1", "sig1")) 

# Pick overlapping signatures: Pivot_long Jaccard matrix, subset to > 0.75, and eliminate all rows with losers (fewer sig1)
jaccard_long <- jaccard %>%
	as.data.frame() %>%
	rownames_to_column(var = "sig1") %>%
	pivot_longer(-sig1, names_to = "sig2", values_to = "jaccard") %>%
	filter(!is.na(jaccard), jaccard > 0.75, sig1 != sig2) %>%
	separate(sig1, "_", into = c("pat1", "sig1")) %>%
	separate(sig2, "_", into = c("pat2", "sig2")) 


# Make sure this loop is functioning the way I want it to!!!
jaccard_new <- jaccard_long
pat_winner <- sig_winner <- c()
for(i in 1:nrow(jaccard_long))
{
	cat("\r",i)
	mypat1 <- as.character(jaccard_long[i,"pat1"])
	mysig1 <- as.numeric(jaccard_long[i,"sig1"])
	
	n1 <- cluster_genes %>% 
		as_tibble() %>%
		filter(patient_id == mypat1, cluster_id == mysig1) %>%
		nrow()
		
	mypat2 <- as.character(jaccard_long[i,"pat2"])
	mysig2 <- as.numeric(jaccard_long[i,"sig2"])
	
	n2 <- cluster_genes %>% 
		as_tibble() %>%
		filter(patient_id == mypat2, cluster_id == mysig2) %>%
		nrow()
		
	pat_loser <- ifelse(n1 >= n2, mypat1, mypat2)
	sig_loser <- ifelse(n1 >= n2, mysig1, mysig2)
	
	jaccard_new <- jaccard_new %>%
		filter(pat1 != pat_loser | sig1 != sig_loser)
}

winners <- jaccard_new %>%
	select(pat1, sig1) %>%
	distinct() %>%
	bind_rows(solo_sigs)

winner_sigs <- cluster_genes %>%
	as_tibble() %>%
	mutate(cluster_id = as.character(cluster_id)) %>%
	inner_join(winners, by = c("patient_id" = "pat1", "cluster_id" = "sig1"))

save.image(file='/projects/varnf/GLASS-III/analysis/single_cell/neftel_macrophage_signatures.RData')
load(file='/projects/varnf/GLASS-III/analysis/single_cell/neftel_macrophage_signatures.RData')

# End
