# WARNING: BEFORE GOING FORWARD
# Make sure this loop is functioning the way I want it to!!!


	
	df <- cluster_genes %>% 
	  # Compute all pairwise co-occurrences
	  full_join(cluster_genes, by="sig") %>% 
	  group_by(tp.x, tp.y) %>% 
	  summarise(N_genes = length(unique(genes))) %>%
	  # Exclude self pairs (gene in h2 and h2)
	  filter(tp.x != tp.y) %>% 
	  # Check for duplicate pairs in reverse order (h2-h4 and h4-h2)
	  mutate(comb=ifelse(tp.x < tp.y, 
						 yes = paste(tp.x, tp.y), 
						 no  = paste(tp.y, tp.x)))
	# Remove the duplicate pairs
	df <- df[which(!duplicated(df$comb)),]

	
	# Get number of genes in each signature for ordering
	n_genes <- table(sub_genes[,"cluster_id"])
	
	sub_genes %>% pivot_wider(
	names_from = cluster_id, 
 	values_from = sig1,
  	values_fn = list(sig1 = intersect(sig1)))
	
	# Convert long table into Jaccard index matrix
	clusters <- unique(sub_genes[,"cluster_id"])
	for(j in 1:)
	jaccard_mat <- acast(sub_genes, cluster_id + sig1 ~ cluster_id + sig1, fun.aggregate = function(x) intersect()/union())

}



#Tidyverse technique (too slow)
sig <- patient_id <- cluster_id <- c()
for(i in 1:nrow(index))
{	
	cat(i, "out of ", nrow(index),"\r")

	# Combine cell_clusters and long_mac (tpm)
	sub_mac <- long_mac %>%
		filter(patient == index[i,"patient_id"]) %>%
		inner_join(cell_clusters, c("cell" = "cell_id", "patient" = "patient_id")) %>%
		group_by(gene_symbol)
		
	# Remove genes that have 0 TPM variance (makes t-test fail), add column indicating cluster of interest, perform t-test, multiple testing correction
	t_test_results <- sub_mac %>%
		summarise(variance = var(tpm)) %>%
		left_join(sub_mac, by = "gene_symbol") %>%
		filter(variance > 0) %>%
		mutate(group = if_else(cluster_id == index[i,"cluster_id"], 'in_clust', 'out_clust')) %>%
		select(gene_symbol, tpm, group) %>%
		group_by(gene_symbol, group) %>%
		nest() %>%
		spread(key = group, value = data) %>%
		mutate(
			t_test = map2(in_clust, out_clust, ~{t.test(.x, .y) %>% tidy()}),
			fold_change = map2(in_clust, out_clust, ~{mean(.x$tpm)/mean(.y$tpm)}),
			in_clust = map(in_clust, nrow),
			out_clust = map(out_clust, nrow)
		) %>% 
  		unnest(cols = c(in_clust, out_clust, fold_change, t_test)) %>%
  		ungroup() %>%
  		mutate(fdr = p.adjust(p.value, method="BH")) %>%
  		filter(fold_change > 3, fdr < 0.05)
  			
	# Only include clusters that have:
	# > 50 genes with FDR < 0.05 + fold-change > 3 AND 
	# > 10 genes with FDR < 0.005 + fold-change > 3
	#----------------------------------------------			
	if(nrow(t_test_results) > 50 & nrow(t_test_results %>% filter(fdr < 0.005)) > 10){
			sig <- c(sig, pull(t_test_results, gene_symbol))
			cluster_id <- c(cluster_id, rep(index[i,"cluster_id"], nrow(t_test_results)))
			patient_id <- c(patient_id, rep(index[i,"patient_id"], nrow(t_test_results)))
	}
}
		
	
#Old Base R (do not use)

	mypatient <- patient[i]
	pat_clusters <- cell_clusters[which(cell_clusters[,"patient_id"] == mypatient),]
	uniq_clusters <- unique(pat_clusters[,"cluster_id"])
	
	# This loop goes through the clusters in patient i and compares all the cells in cluster j to all the cells not in cluster j
	# Differential expression analyses are performed with a t.test and then filtered using the criteria in the Neftel manuscript
	for(j in 1:length(uniq_clusters))
	{
		cat(i,"-->",j,"\r")
		in_cells <- as.character(pat_clusters[which(pat_clusters[,"cluster_id"] == uniq_clusters[j]),"cell_id"])
		
		sub_mac <- mac_rel_tpm[,grep(mypatient,colnames(mac_rel_tpm))]
		sub_mac <- melt(sub_mac)
		
		
		
		# Remove genes with 0 variance (makes t-test fail)
		sub_mac <- sub_mac[which(apply(sub_mac,1,var) > 0),]
		
		# Get cells in cluster j
		ind1 <- which(colnames(sub_mac) %in% in_cells)
		# Get cells outside cluster j
		ind2 <- which(!(colnames(sub_mac) %in% in_cells))
		
		fold_change <- apply(sub_mac, 1, function(x) mean(x[ind1]) / mean(x[ind2]))	
		p_value <- apply(sub_mac, 1, function(x) t.test(x[ind1], x[ind2])$p.value)
		fdr <- p.adjust(p_value, "BH")
		
		# Only include clusters that have:
		# > 50 genes with FDR < 0.05 + fold-change > 3 AND 
		# > 10 genes with FDR < 0.005 + fold-change > 3
		#----------------------------------------------
		sig_res <- data.frame(fold_change, p_value, fdr)
		genes1 <- rownames(sig_res)[which(sig_res[,"fold_change"] > 3 & sig_res[,"fdr"] < 0.05)]
		genes2 <- rownames(sig_res)[which(sig_res[,"fold_change"] > 3 & sig_res[,"fdr"] < 0.005)]
		
		if(length(genes1) > 50 & length(genes2) > 10){
				sig <- c(sig, genes1)
				cluster_id <- c(cluster_id, rep(uniq_clusters[j], length(genes1)))
				patient_id <- c(patient_id, rep(mypatient, length(genes1)))
		}
	}
}

cluster_genes <- data.frame(patient_id, cluster_id, sig)


final_clusters <- cluster_genes %>% 
full_join(cluster_genes, by="patient_id") %>% 
group_by(cluster_id.x, cluster_id.y) %>% 
  summarise(test = length(intersect(sig1.x, sig1.y))/length(union(sig1.x, sig1.y)))


# Collapse signatures based on Jaccard index > 0.75

patient_label <- as.character(unique(cell_clusters[,"patient_id"]))

for(i in 1:patient_label)
{
	mypatient <- patient_label[i]
	sub_genes <- cluster_genes[which(cluster_genes[,"patient_id"] == mypatient),2:3]
	
	df <- sub_genes %>% 
	  # Compute all pairwise Jaccard indexes
	  full_join(sub_genes, by="sig1") %>% 
	  group_by(cluster_id.x, cluster_id.y) %>% 
	  summarise(test = length(intersect(sig1.x, sig1.y))/length(union(sig1.x, sig1.y)))

	# Remove the duplicate pairs
	df <- df[which(!duplicated(df$comb)),]
	
	# Get number of genes in each signature for ordering
	n_genes <- table(sub_genes[,"cluster_id"])
	
	sub_genes %>% pivot_wider(
	names_from = cluster_id, 
 	values_from = sig1,
  	values_fn = list(sig1 = intersect(sig1)))
	
	# Convert long table into Jaccard index matrix
	clusters <- unique(sub_genes[,"cluster_id"])
	for(j in 1:)
	jaccard_mat <- acast(sub_genes, cluster_id + sig1 ~ cluster_id + sig1, fun.aggregate = function(x) intersect()/union())

}

	}
	
	for(j in 1:)	
		if(j == 1){
			nsig1[[1]] <- sig1
			nsig2[[1]] <- sig2
		} else{
			int1 <- sapply(nsig1, function(x) length(intersect(x, sig1)))
			uni1 <- sapply(nsig1, function(x) length(union(x, sig1)))
			jac1 <- int1/uni1
			
			
			int2 <- sapply(nsig2, function(x), length(intersect(x, sig2)))
			uni2 <- sapply(nsig2, function(x), length(union(x, sig2)))
			jac2 <- int2/uni2

		}


### old method with relative TPM:

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

# Calculate relative expression by centering expression levels per gene about mean
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
	sub_mac <- mac_filt_tpm[,grep(patient[i],colnames(mac_filt_tpm))]
	
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
	
	sub_mac <- data.matrix(mac_filt_tpm[,grep(index[i,"patient_id"],colnames(mac_filt_tpm))])
	group <- as.character(colnames(sub_mac) %in% in_cells) %>%
		recode("FALSE"="out","TRUE"="in")
	
	in_genes2 <- as.data.frame(findMarkers(sub_mac, group, test.type="t")[["in"]]) %>%
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

# Pick overlapping signatures: Melt Jaccard matrix, subset to > 0.75, and eliminate all rows with losers (fewer sig1)
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

# End
