library(tidyverse)
library(odbc)
library(DBI)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
SELECT im.*,ts.total_tcr
FROM analysis.tcr_stats ts
LEFT JOIN analysis.davoli_immune_score im ON ts.aliquot_barcode = im.aliquot_barcode
WHERE total_tcr > 0
ORDER BY signature_name"

dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])

mycor <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub.dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	mycor[i] <- cor(sub.dat[,"enrichment_score"],sub.dat[,"total_tcr"],method="s")
}
names(mycor) <- cells

#               B.cells            CD4.mature          CD8.effector 
#             0.4693816             0.5126435             0.7331365 
# CD8.effector.NK.cells             Dendritic           Macrophages 
#             0.7015698             0.5736662             0.6263403 
#        Macrophages.M1        Macrophages.M2              NK.cells 
#             0.6849468             0.6133317             0.5167367 
#                 T.reg 
#             0.5328893 

plot_res <- dat[which(dat[,"signature_name"]=="CD8.effector"),]
plot_res[,"total_tcr"] <- log10(plot_res[,"total_tcr"] + 1)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/tcr_cd8_cor.pdf",width=2,height=2)
ggplot(plot_res,aes(x = total_tcr, y = enrichment_score)) +
geom_point(size=0.5) +
geom_smooth(method='lm',formula= y~x) +
labs(x = "TCR count (log)", y = "CD8+ T enrichment score", title="TCR abundance") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 
dev.off()

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
SELECT im.*,bs.total_bcr
FROM analysis.bcr_stats bs
LEFT JOIN analysis.davoli_immune_score im ON bs.aliquot_barcode = im.aliquot_barcode
--WHERE total_bcr > 0
ORDER BY signature_name"

dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])

mycor <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub.dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	mycor[i] <- cor(sub.dat[,"enrichment_score"],sub.dat[,"total_bcr"],method="s")
}
names(mycor) <- cells

#               B.cells            CD4.mature          CD8.effector 
#             0.4438014             0.4260259             0.5057637 
# CD8.effector.NK.cells             Dendritic           Macrophages 
#             0.5405564             0.4446273             0.4853530 
#        Macrophages.M1        Macrophages.M2              NK.cells 
#             0.5324960             0.4957067             0.4057343 
#                 T.reg 
#             0.4801683 

plot_res <- dat[which(dat[,"signature_name"]=="CD8.effector"),]
plot_res[,"total_tcr"] <- log10(plot_res[,"total_tcr"] + 1)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/tcr_cd8_cor.pdf",width=2,height=2)
ggplot(plot_res,aes(x = total_tcr, y = enrichment_score)) +
geom_point(size=0.5) +
geom_smooth(method='lm',formula= y~x) +
labs(x = "TCR count (log)", y = "CD8+ T enrichment score", title="TCR abundance") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 
dev.off()