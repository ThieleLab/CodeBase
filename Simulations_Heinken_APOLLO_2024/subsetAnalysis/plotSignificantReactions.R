
library(RColorBrewer)
library(pheatmap)
library(readxl)
library(stringr)
setwd("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Heatmaps")

## plot different analyses of disease cases

mycol <- c("white",brewer.pal(4,"Reds")[1:4])

subsystems=read.delim("/Users/almut.heinken/Documents/Code/cobratoolbox/papers/2021_demeter/input/ReactionDatabase.txt",header=T,row.names=1)[c("Subsystem_general")]

#annotation <- annotation %>% replace_na(list(Disease.name = 'Healthy', Age.group = 'Adult'))

### Body sites ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Adults_body_sites_healthy/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Adults_body_sites_healthy/Adults_body_sites_healthy_samples.csv",header=T,row.names=1)[c("Body.site")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Body.site <- as.factor(annotation$Body.site)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Body.site)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Body.site)

subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Body.site=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("Adults_body_sites_healthy.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### Healthy adults vs. infants ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Adults_vs_infants_healthy/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Adults_vs_infants_healthy/Adults_vs_infants_healthy_samples.csv",header=T,row.names=1)[c("Age.group")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Age.group <- as.factor(annotation$Age.group)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Age.group)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Age.group)

subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Age.group=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("Adults_vs_infants_healthy.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### IBD ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/IBD_vs_healthy/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/IBD_vs_healthy/IBD_vs_healthy_samples.csv",header=T,row.names=1)[c("Disease.name")]

annotation$Disease.name <- str_replace_all(annotation$Disease.name,"Obesity","Healthy")

# define the colors
annotation$Disease.name <- as.factor(annotation$Disease.name)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Disease.name)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Disease.name)

subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Disease.name=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("IBD_vs_healthy.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### Infants premature vs. healthy ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infants_premature_vs_healthy/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infants_premature_vs_healthy/Infants_premature_vs_healthy_samples.csv",header=T,row.names=1)[c("Disease.name")]

# define the colors
annotation$Disease.name <- as.factor(annotation$Disease.name)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Disease.name)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Disease.name)

subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Disease.name=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("Infants_premature_vs_healthy.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()


### Infants undernourished vs. healthy ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infants_undernourished_vs_healthy/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infants_undernourished_vs_healthy/Infants_undernourished_vs_healthy_samples.csv",header=T,row.names=1)[c("Disease.name")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Disease.name <- as.factor(annotation$Disease.name)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Disease.name)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Disease.name)

subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Disease.name=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("Infants_undernourished_vs_healthy.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### Infection antibiotics vs no antibiotics ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infection_antibiotics_vs_no_antibiotics/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infection_antibiotics_vs_no_antibiotics/Infection_antibiotics_vs_no_antibiotics_samples.csv",header=T,row.names=1)[c("Antibiotics")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Antibiotics <- as.factor(annotation$Antibiotics)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Antibiotics)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Antibiotics)

subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Antibiotics=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("Infection_antibiotics_vs_no_antibiotics.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### Infection resistant vs susceptible ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infection_resistant_vs_susceptible/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infection_resistant_vs_susceptible/Infection_resistant_vs_susceptible_samples.csv",header=T,row.names=1)[c("Disease.severity","Stratification")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Stratification <- as.factor(annotation$Stratification)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Stratification)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Stratification)
annotation$Disease.severity <- as.factor(annotation$Disease.severity)
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
Palette = getPalette(length(levels(annotation$Disease.severity)))
Ann2.cols <- c(Palette)
names(Ann2.cols) <- levels(annotation$Disease.severity)
subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Stratification=Ann1.cols,
  Disease.severity=Ann2.cols,
  Subsystem_general=Ann3.cols
) 

png("Infection_resistant_vs_susceptible.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### Infection vs healthy ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infection_vs_healthy/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Infection_vs_healthy/Infection_vs_healthy_samples.csv",header=T,row.names=1)[c("Disease.name")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Disease.name <- as.factor(annotation$Disease.name)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Disease.name)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Disease.name)
subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Disease.name=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("Infection_vs_healthy.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()


### Obesity vs normalweight ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Obesity_vs_normalweight/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/Obesity_vs_normalweight/Obesity_vs_normalweight_samples.csv",header=T,row.names=1)[c("Disease.name")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Disease.name <- as.factor(annotation$Disease.name)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Disease.name)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Disease.name)
subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Disease.name=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("Obesity_vs_normalweight.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### PD_vs_healthy ###
data = read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/PD_vs_healthy/StatisticalAnalysis/SignificantFeatures_Reactions_presence.csv",header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)

annotation=read.csv("/Users/almut.heinken/National University of Ireland, Galway/Group_MSP - Dokumente/150k_Project/Analysis_microbiome_models/Subgroup_analysis/Subgroups/PD_vs_healthy/PD_vs_healthy_samples.csv",header=T,row.names=1)[c("Disease.name")]

subsystemsReduced <- subsystems
rowsToKeep <- intersect(rownames(subsystemsReduced), colnames(data))
subsystemsReduced = subsystemsReduced[rowsToKeep,, drop = FALSE]

# define the colors
annotation$Disease.name <- as.factor(annotation$Disease.name)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
Palette = getPalette(length(levels(annotation$Disease.name)))
Ann1.cols <- c(Palette)
names(Ann1.cols) <- levels(annotation$Disease.name)
subsystemsReduced$Subsystem_general <- as.factor(subsystemsReduced$Subsystem_general)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
Palette = getPalette(length(levels(subsystemsReduced$Subsystem_general)))
Ann3.cols <- c(Palette)
names(Ann3.cols) <- levels(subsystemsReduced$Subsystem_general)

ann_colors = list(
  Disease.name=Ann1.cols,
  Subsystem_general=Ann3.cols
) 

png("PD_vs_healthy.png", width = 10, height = 10, units = 'in', res = 300)
pheatmap(x,
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="ward.D",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize=10,
         fontsize_col=10,
         annotation_row=annotation,
         annotation_col=subsystemsReduced,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()
