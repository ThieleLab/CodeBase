
#make sure all these packages are installed and loaded
library(NMF)
library(vegan)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(pheatmap)

setwd("C:/Users/Almut Heinken/National University of Ireland, Galway/Group_MSP - Documents/AGORA2/Modeling_CRC/Correlations")

############ Japanese diet ########################

#### Species correlations #####
taxonomy=read.csv("Taxonomy_JD_Species.csv",header=T,row.names=1)[c("Class","Phylum")]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(PhylumPalette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
ClassPalette = getPalette(length(levels(taxonomy$Class)))
ClassAnn.cols <- c(ClassPalette)
names(ClassAnn.cols) <- levels(taxonomy$Class)

ann_colors = list(
  # Metabolite_Subsystem=MetSubsystemAnn.cols,
  Phylum=PhylumAnn.cols,
  Class=ClassAnn.cols)

data = read.csv("JD_Species.csv",header=T,row.names=1, check.names = FALSE)
x <- t(as.matrix(data))
mycol <- c(brewer.pal(9,"Blues")[5:1],"white",brewer.pal(9,"Reds")[1:9])
png("JD_SpeciesCorrelations.png", width = 16, height = 8, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         #annotation_row=info,
         annotation_col=taxonomy,
         annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()

#### Genus correlations #####
taxonomy=read.csv("Taxonomy_JD_Genus.csv",header=T,row.names=1)[c("Class","Phylum")]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(PhylumPalette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
ClassPalette = getPalette(length(levels(taxonomy$Class)))
ClassAnn.cols <- c(ClassPalette)
names(ClassAnn.cols) <- levels(taxonomy$Class)

ann_colors = list(
  # Metabolite_Subsystem=MetSubsystemAnn.cols,
  Phylum=PhylumAnn.cols,
  Class=ClassAnn.cols)

data = read.csv("JD_Genus.csv",header=T,row.names=1, check.names = FALSE)
x <- t(as.matrix(data))
mycol <- c(brewer.pal(9,"Blues")[5:1],"white",brewer.pal(9,"Reds")[1:9])
png("JD_GenusCorrelations.png", width = 14, height = 8, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         #annotation_row=info,
         annotation_col=taxonomy,
         annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()

#### Family correlations #####
taxonomy=read.csv("Taxonomy_JD_Family.csv",header=T,row.names=1)[c("Class","Phylum")]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(PhylumPalette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
ClassPalette = getPalette(length(levels(taxonomy$Class)))
ClassAnn.cols <- c(ClassPalette)
names(ClassAnn.cols) <- levels(taxonomy$Class)

ann_colors = list(
  # Metabolite_Subsystem=MetSubsystemAnn.cols,
  Phylum=PhylumAnn.cols,
  Class=ClassAnn.cols)

data = read.csv("JD_Family.csv",header=T,row.names=1, check.names = FALSE)
x <- t(as.matrix(data))
mycol <- c(brewer.pal(9,"Blues")[5:1],"white",brewer.pal(9,"Reds")[1:9])
png("JD_FamilyCorrelations.png", width = 14, height = 8, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         #annotation_row=info,
         annotation_col=taxonomy,
         annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()

#### Order correlations #####
taxonomy=read.csv("Taxonomy_JD_Order.csv",header=T,row.names=1)[c("Class","Phylum")]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(PhylumPalette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
ClassPalette = getPalette(length(levels(taxonomy$Class)))
ClassAnn.cols <- c(ClassPalette)
names(ClassAnn.cols) <- levels(taxonomy$Class)

ann_colors = list(
  # Metabolite_Subsystem=MetSubsystemAnn.cols,
  Phylum=PhylumAnn.cols,
  Class=ClassAnn.cols)

data = read.csv("JD_Order.csv",header=T,row.names=1, check.names = FALSE)
x <- t(as.matrix(data))
mycol <- c(brewer.pal(9,"Blues")[5:1],"white",brewer.pal(9,"Reds")[1:9])
png("JD_OrderCorrelations.png", width = 14, height = 8, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         #annotation_row=info,
         annotation_col=taxonomy,
         annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()

#### Class correlations #####
taxonomy=read.csv("Taxonomy_JD_Class.csv",header=T,row.names=1)[c("Phylum")]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(PhylumPalette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)

ann_colors = list(
  # Metabolite_Subsystem=MetSubsystemAnn.cols,
  Phylum=PhylumAnn.cols)

data = read.csv("JD_Class.csv",header=T,row.names=1, check.names = FALSE)
x <- t(as.matrix(data))
mycol <- c(brewer.pal(9,"Blues")[5:1],"white",brewer.pal(9,"Reds")[1:9])
png("JD_ClassCorrelations.png", width = 14, height = 8, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         #annotation_row=info,
         annotation_col=taxonomy,
         annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()

#### Reaction correlations ####
data = read.csv("ReactionCorrelations.csv",header=T,row.names=1, check.names = FALSE)
x <- t(as.matrix(data))
mycol <- c(brewer.pal(9,"Blues")[6:1],"white",brewer.pal(9,"Reds")[1:9])
png("ReactionCorrelations.png", width = 10, height = 8, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         #annotation_row=info,
         #annotation_col=taxonomy,
         #annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()

