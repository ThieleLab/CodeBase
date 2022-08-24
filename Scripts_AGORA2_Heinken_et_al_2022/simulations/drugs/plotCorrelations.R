
library(RColorBrewer)
library(pheatmap)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd('..')
setwd('..')
setwd("Modeling_CRC/Plots")

############ Japanese diet ########################

#### Species correlations #####
taxonomy=read.csv("Taxonomy_JD_Species.csv",header=T,row.names=1,check.names = FALSE,stringsAsFactors=T)[c("Class","Phylum")]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(PhylumPalette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
ClassPalette = getPalette(length(levels(taxonomy$Class)))
ClassAnn.cols <- c(ClassPalette)
names(ClassAnn.cols) <- levels(taxonomy$Class)

ann_colors = list(
  Phylum=PhylumAnn.cols,
  Class=ClassAnn.cols)

data = read.csv("JD_Species.csv",header=T,row.names=1, check.names = FALSE)
mycol <- c(brewer.pal(9,"Blues")[5:1],"white",brewer.pal(9,"Reds")[1:9])
png("JD_SpeciesCorrelations.png", width = 10, height = 14, units = 'in', res = 300)
pheatmap(data,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
          annotation_row = taxonomy,
         annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()
