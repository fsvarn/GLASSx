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