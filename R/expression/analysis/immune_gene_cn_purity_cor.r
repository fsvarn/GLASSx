#######################################################

library(tidyverse)
library(odbc)
library(DBI)
library(ggbeeswarm)
library(RColorBrewer)
library(gridExtra)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
expr <- dbGetQuery(con,read_file("/projects/varnf/GLASS-III/GLASS-III/sql/expression/immune_sig_purity_copy_num.sql"))


method <- unique(expr[,"signature_set"])
pur_cor <- cn_cor <- gene_symbol <- sig_set <- pur_pval <- cn_pval <- c()
for(i in 1:length(method))
{
	sub_meth <- expr[which(expr[,"signature_set"]==method[i]),]
	gene <- unique(sub_meth[,"gene_symbol"])
	for(j in 1:length(gene))
	{
		sub_gene <- sub_meth[which(sub_meth[,"gene_symbol"]==gene[j]),]
		pur_cor <- c(pur_cor,cor(sub_gene[,"tpm"],sub_gene[,"cellularity"],method="s"))
		pur_pval <- c(pur_pval,cor.test(sub_gene[,"tpm"],sub_gene[,"cellularity"],method="s")$p.value)
		cn_cor <- c(cn_cor,cor(sub_gene[,"tpm"],sub_gene[,"wcr"],method="s"))
		cn_pval <- c(cn_pval,cor.test(sub_gene[,"tpm"],sub_gene[,"wcr"],method="s")$p.value)
		gene_symbol <- c(gene_symbol, gene[j])
		sig_set <- c(sig_set, method[i])
	}
}

pur_sig <- cn_sig <- rep(0,length(pur_pval))

pur_sig[which(pur_pval < 0.001)] <- 1
pur_sig[which(pur_pval >= 0.001 & pur_pval < 0.01)] <- 2
pur_sig[which(pur_pval >= 0.01 & pur_pval < 0.05)] <- 3
pur_sig[which(pur_pval >= 0.05)] <- 4
cn_sig[which(cn_pval < 0.001)] <- 1
cn_sig[which(cn_pval >= 0.001 & cn_pval < 0.01)] <- 2
cn_sig[which(cn_pval >= 0.01 & cn_pval < 0.05)] <- 3
cn_sig[which(cn_pval >= 0.05)] <- 4

pur_sig <- as.factor(pur_sig)
cn_sig <- as.factor(cn_sig)

res <- data.frame(sig_set, gene_symbol, pur_cor, cn_cor,pur_sig,cn_sig)

n <- pur_pct <- cn_pct <- rep(0,length(method))
for(i in 1:length(method))
{
	sub_res <- res[which(res[,"sig_set"]==method[i]),]
	n[i] <- nrow(sub_res)
	pur_pct[i] <- sum(as.numeric(sub_res[,"pur_sig"])<4 & sub_res[,"pur_cor"] < 0)/nrow(sub_res)
	cn_pct[i] <- sum(as.numeric(sub_res[,"cn_sig"])<4 & sub_res[,"cn_cor"] > 0)/nrow(sub_res)
}
pct_res <- data.frame(method,n,pur_pct,cn_pct)

p1 <-ggplot(res,aes(x = sig_set, y = pur_cor,colour=pur_sig)) +
geom_beeswarm(size=0.3) +
stat_summary(fun.y=median,geom="line",colour="red") +
scale_colour_manual(values=rev(brewer.pal(4,"Blues"))) +
labs(title = "Purity", y = "Spearman's rho") +
geom_hline(yintercept=0,linetype="dashed") +
theme_bw() +
theme(axis.text.x=element_blank(),axis.text.y = element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")

p2 <-ggplot(pct_res,aes(x = method, y = pur_pct)) +
geom_bar(stat="identity") +
labs(y = "Genes correlated") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,1))

p3 <- ggplot(res,aes(x = sig_set, y = cn_cor,colour=cn_sig)) +
geom_beeswarm(size=0.3) +
stat_summary(fun.y=median,geom="line",colour="red") +
scale_colour_manual(values=rev(brewer.pal(4,"Blues"))) +
labs(title = "Purity", y = "Spearman's rho") +
geom_hline(yintercept=0,linetype="dashed") +
theme_bw() +
theme(axis.text.x=element_blank(),axis.text.y = element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")

p4 <-ggplot(pct_res,aes(x = method, y = cn_pct)) +
geom_bar(stat="identity") +
labs(y = "Genes correlated") +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,1))

# gb1 <- ggplot_build(p1)
# gb2 <- ggplot_build(p2)
# gb3 <- ggplot_build(p3)
# gb4 <- ggplot_build(p4)
# 
# n1 <- length(gb1$layout$panel_params[[1]]$y.labels)
# n2 <- length(gb2$layout$panel_params[[1]]$y.labels)
# n3 <- length(gb3$layout$panel_params[[1]]$y.labels)
# n4 <- length(gb4$layout$panel_params[[1]]$y.labels)
# 
# gA <- ggplot_gtable(gb1)
# gB <- ggplot_gtable(gb2)
# gC <- ggplot_gtable(gb3)
# gD <- ggplot_gtable(gb4)
# 
# g <- gtable:::rbind_gtable(gA, gB, "last")
# g <- gtable:::rbind_gtable(g, gC, "last")
# g <- gtable:::rbind_gtable(g, gD, "last")
# 
# panels <- g$layout$t[grep("panel", g$layout$name)]
# g$heights[panels[1]] <- unit(n1/1, "null") 
# g$heights[panels[2]] <- unit(n1/6,"null")
# g$heights[panels[3]] <- unit(n1/1,"null")
# g$heights[panels[4]] <- unit(n1/6,"null")
# 
# #width_reference <- gb1
# 
# grid.newpage()
# 
# #widths <-length(width_reference[[1]]$layout$panel_params[[1]]$x.labels)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/rnaseq_gene_benchmark.pdf",width=4,height=3)
p1
p2
p3
p4
dev.off()


#--------------------------------------

#Immune cell-specific analysis

#######################################################

library(tidyverse)
library(odbc)
library(DBI)
library(ggbeeswarm)
library(RColorBrewer)
library(gridExtra)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
expr <- dbGetQuery(con,read_file("/projects/varnf/GLASS-III/GLASS-III/sql/expression/immune_sig_purity_copy_num.sql"))

method <- paste(expr[,"signature_set"],expr[,"signature_name"],sep="__")
method <- unique(method)
method <- method[order(method)]
pur_cor <- cn_cor <- gene_symbol <- sig_set <- pur_pval <- cn_pval <- c()
for(i in 1:length(method))
{
	sub_meth <- expr[which(paste(expr[,"signature_set"],expr[,"signature_name"],sep="__")==method[i]),]
	gene <- unique(sub_meth[,"gene_symbol"])
	for(j in 1:length(gene))
	{
		sub_gene <- sub_meth[which(sub_meth[,"gene_symbol"]==gene[j]),]
		pur_cor <- c(pur_cor,cor(sub_gene[,"tpm"],sub_gene[,"cellularity"],method="s"))
		pur_pval <- c(pur_pval,cor.test(sub_gene[,"tpm"],sub_gene[,"cellularity"],method="s")$p.value)
		cn_cor <- c(cn_cor,cor(sub_gene[,"tpm"],sub_gene[,"wcr"],method="s"))
		cn_pval <- c(cn_pval,cor.test(sub_gene[,"tpm"],sub_gene[,"wcr"],method="s")$p.value)
		gene_symbol <- c(gene_symbol, gene[j])
		sig_set <- c(sig_set, method[i])
	}
}

pur_sig <- cn_sig <- rep(0,length(pur_pval))

pur_sig[which(pur_pval < 0.001)] <- 1
pur_sig[which(pur_pval >= 0.001 & pur_pval < 0.01)] <- 2
pur_sig[which(pur_pval >= 0.01 & pur_pval < 0.05)] <- 3
pur_sig[which(pur_pval >= 0.05)] <- 4
cn_sig[which(cn_pval < 0.001)] <- 1
cn_sig[which(cn_pval >= 0.001 & cn_pval < 0.01)] <- 2
cn_sig[which(cn_pval >= 0.01 & cn_pval < 0.05)] <- 3
cn_sig[which(cn_pval >= 0.05)] <- 4

pur_sig <- as.factor(pur_sig)
cn_sig <- as.factor(cn_sig)

res <- data.frame(sig_set, gene_symbol, pur_cor, cn_cor,pur_sig,cn_sig)

n <- pur_pct <- cn_pct <- rep(0,length(method))
for(i in 1:length(method))
{
	sub_res <- res[which(res[,"sig_set"]==method[i]),]
	n[i] <- nrow(sub_res)
	pur_pct[i] <- sum(as.numeric(sub_res[,"pur_sig"])<4 & sub_res[,"pur_cor"] < 0)/nrow(sub_res)
	cn_pct[i] <- sum(as.numeric(sub_res[,"cn_sig"])<4 & sub_res[,"cn_cor"] > 0)/nrow(sub_res)
}
pct_res <- data.frame(method,n,pur_pct,cn_pct)


