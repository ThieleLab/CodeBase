
library(RColorBrewer)
library(pheatmap)
library(readxl)
library(stringr)
setwd("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/mostImportantReactions_Plotted")

## plot reactions that were top hits in Random forests analyses

data = read.csv("Abundance_Adults_vs_infants_healthy.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

mycol <- c("white",brewer.pal(9,"Reds")[1:9])

png("Abundance_Adults_vs_infants_healthy.png", width = 10, height = 6, units = 'in', res = 300)
pheatmap(x,
         #col=mycol,
         scale="row",
         #cluster_rows=FALSE,
         cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         border=NA,
         show_legend=FALSE
         #show_rownames=FALSE,
         #show_colnames=FALSE
)
dev.off()