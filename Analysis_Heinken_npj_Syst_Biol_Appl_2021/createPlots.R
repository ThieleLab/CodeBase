
library(NMF)
library(vegan)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(pheatmap)
library(bigmemory)
library(rstudioapi)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)

dir.create("R_plots")

resourcePath <- paste0(currPath, "/inputFiles/", collapse = NULL)
fluxPath <- paste0(currPath, "/NetFluxes/", collapse = NULL)
contrPath <- paste0(currPath, "/MicrobeContributions/", collapse = NULL)
corrPath <- paste0(currPath, "/Correlations/", collapse = NULL)
statPath <- paste0(currPath, "/Statistics/", collapse = NULL)
metProfPath <- paste0(currPath, "/MetaboliteProfiles/", collapse = NULL)
savePath <- paste0(currPath, "/R_plots/", collapse = NULL)
setwd(savePath)


#### create the plots for Figure 3 ####

#### Principal compoments analysis plot ###########
#### flux contributions ###########

# Stratification
metadata=read.csv(paste0(resourcePath,"metadata_IBD.csv", collapse = NULL))
data = read.csv(paste0(contrPath,"MicrobeContributions_Fluxes.csv", collapse = NULL),header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)
cover <- t(x)
pclass <- metadata

dis = as.character(pclass[,2])
names(dis) = pclass[,1]
dires <- capscale(cover~1,distance = "bray")
dires_sum <- summary(dires)

png("PCoA_fluxes.png", width = 14, height = 10, units = 'in', res = 300)
gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                   "Phenotype"=as.factor(dis[rownames(dires_sum$sites)]))
levels(gdat$Phenotype) = c("Healthy","IBD_dysbiotic", "IBD_nondysbiotic")
gdat[which(gdat$Phenotype=="IBD_dysbiotic"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="IBD_dysbiotic"),"PC1"])
gdat[which(gdat$Phenotype=="IBD_nondysbiotic"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="IBD_nondysbiotic"),"PC1"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
gdat[which(gdat$Phenotype=="IBD_dysbiotic"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="IBD_dysbiotic"),"PC2"])
gdat[which(gdat$Phenotype=="IBD_nondysbiotic"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="IBD_nondysbiotic"),"PC2"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
ggplot(data=gdat, aes(x=PC1, y=PC2, color=Phenotype)) +
  geom_hline(aes(yintercept=0),color="grey", size=1.5) +
  geom_vline(aes(xintercept=0),color="grey", size=1.5) +
  stat_ellipse(level = 0.95, size=2) +
  scale_colour_manual(values = c("royalblue2","firebrick2","yellow")) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=1.8) +
  geom_point(aes(x=PC1, y=PC2),size=5) +
  xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
  ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
  theme_bw(base_size = 20)
dev.off()


#### Correlations between net production fluxes and species correlations ###########

taxonomy=read.csv(paste0(corrPath,"Taxonomy_Correlations_Species.csv", collapse = NULL),header=T,row.names=1)[c("Genus","Phylum")]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
PhylumAnn.cols <- c(PhylumPalette)
names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
GenusPalette = getPalette(length(levels(taxonomy$Genus)))
GenusAnn.cols <- c(GenusPalette)
names(GenusAnn.cols) <- levels(taxonomy$Genus)

info=read.csv(paste0(resourcePath,"MetaboliteInformation.csv", collapse = NULL),header=T,row.names=1)[c("Metabolite_Subsystem")]
# Stratification
modtax <- info[order(info$Metabolite_Subsystem),]
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
SubsystemPalette = getPalette(length(levels(info$Metabolite_Subsystem)))

MetSubsystemAnn.cols <- c(SubsystemPalette)
names(MetSubsystemAnn.cols) <- levels(info$Metabolite_Subsystem)

ann_colors = list(
  Metabolite_Subsystem=MetSubsystemAnn.cols,
  Phylum=PhylumAnn.cols,
  Genus=GenusAnn.cols)
data = read.csv(paste0(corrPath,"Correlations_Species.csv", collapse = NULL),header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)
mycol <- c(brewer.pal(9,"Blues")[9:1],"white",brewer.pal(9,"Reds")[1:9])
png("NetProd_StrongSpeciesCorrelations.png", width = 12, height = 10, units = 'in', res = 300)
pheatmap(t(x),
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize_row=8,
         fontsize_col=8,
         annotation_row=info,
         annotation_col=taxonomy,
         annotation_colors=ann_colors,
         #show_rownames=FALSE,
         #show_colnames=FALSE,
         border=NA
)
dev.off()


#### reaction presence that differed between groups
# Stratification
metadata=read.csv(paste0(resourcePath,"metadata_IBD.csv", collapse = NULL),row.names=1)[c("Stratification")]

info=read.csv(paste0(resourcePath,"ReactionDatabase.csv", collapse = NULL),row.names=1)[c("Reaction_Subsystem")]
# Stratification
modtax <- info[order(info$Reaction_Subsystem),]
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
SubsystemPalette = getPalette(length(levels(info$Reaction_Subsystem)))
SubsystemAnn.cols <- c(SubsystemPalette)
names(SubsystemAnn.cols) <- levels(info$Reaction_Subsystem)

ann_colors = list(
  Stratification = c(Healthy = "royalblue2", IBD_nondysbiotic = "yellow", IBD_dysbiotic = "firebrick2"),
  Reaction_Subsystem=SubsystemAnn.cols
) 

data = read.delim(paste0(statPath,"ReactionPresence_SignificantFeatures.txt", collapse = NULL),header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)
mycol <- colorRampPalette(c("black", "red"))(100)
png("ReactionPresence.png", width = 10, height = 8, units = 'in', res = 300)
pheatmap(x,
         #scale="row",
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='binary',
         clustering_distance_cols='binary',
         fontsize=10,
         fontsize_col=8,
         annotation_col=metadata,
         annotation_row=info,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()

### plot subsystem abundances for supplementary material ###
##### reaction subsystem abundance ####
# Stratification
metadata=read.csv(paste0(resourcePath,"metadata_IBD.csv", collapse = NULL),row.names=1)[c("Stratification")]

data = read.delim(paste0(statPath,"SubsystemAbundance_SignificantFeatures.txt", collapse = NULL),header=T,row.names=1, check.names = FALSE)
x <- as.matrix(data)
#mycol <- c("white",brewer.pal(9,"Reds")[1:9])
mycol <- colorRampPalette(c("blue","white", "firebrick3"))(10000)
ann_colors = list(
  Stratification = c(Healthy = "royalblue2", IBD_nondysbiotic = "yellow", IBD_dysbiotic = "firebrick2")
) 
png("SubsystemsSummarized.png", width = 10, height = 8, units = 'in', res = 300)
pheatmap(x,
         scale="row",
         col=mycol,
         #cluster_rows=FALSE,
         #cluster_cols=FALSE,
         clustering_method="average",
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         fontsize_row=6,
         fontsize_col=8,
         annotation_col=metadata,
         annotation_colors=ann_colors,
         border=NA,
         show_legend=FALSE,
         #show_rownames=FALSE,
         show_colnames=FALSE
)
dev.off()


##### plot metabolite profiles for supplementary material ####

metadata=read.csv(paste0(resourcePath,"metadata_IBD.csv", collapse = NULL),header=T,row.names=1)[c("Stratification")]

mycol <- colorRampPalette(c("blue", "white", "red"))(10000)

setwd(metProfPath)
sheetnames=read.csv(paste0(metProfPath,"sheetnames_secretion_metabolites.csv", collapse = NULL),header=F,row.names=1);
######### read in all of the sheets with the data and generate the plots ###############
plotlist=list()
for (i in 1:nrow(sheetnames)){
  filename=as.character(sheetnames[i,1])
  plottitle=as.character(sheetnames[i,3])
  savename=as.character(sheetnames[i,4])
  data= read.csv(file=filename,header = T,row.names=1, check.names = FALSE)
  filename=as.character(sheetnames[i,2])
  taxonomy= read.csv(file=filename,header = T,row.names=2, check.names = FALSE)[c("Phylum","Class")]
  # Stratification
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))
  PhylumPalette = getPalette(length(levels(taxonomy$Phylum)))
  PhylumAnn.cols <- c(PhylumPalette)
  names(PhylumAnn.cols) <- levels(taxonomy$Phylum)
  getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
  ClassPalette = getPalette(length(levels(taxonomy$Class)))
  ClassAnn.cols <- c(ClassPalette)
  names(ClassAnn.cols) <- levels(taxonomy$Class)
  #
  ann_colors = list(
    Stratification = c(Healthy = "royalblue2", IBD_nondysbiotic = "yellow", IBD_dysbiotic = "firebrick2"),
    Phylum=PhylumAnn.cols,
    Class=ClassAnn.cols) 
  #
  x <- as.matrix(data)
  
  if (nrow(x) >100) {
    height_val = 20
  } else if (nrow(x) >50) {
    height_val = 16
  } else if (nrow(x) >20) {
    height_val = 14
  }  else {
    height_val = 8
  }
  
  plot <- pheatmap(x,
                   scale="row",
                   col=mycol,
                   #cluster_rows=FALSE,
                   cluster_cols=FALSE,
                   clustering_method="average",
                   clustering_distance_rows='euclidean',
                   clustering_distance_cols='euclidean',
                   fontsize=12,
                   annotation_col=metadata,
                   annotation_row=taxonomy,
                   annotation_colors=ann_colors,
                   border=NA,
                   show_legend=FALSE,
                   #show_rownames=FALSE,
                   show_colnames=FALSE,
                   main=plottitle
  )
  plotlist[[i]]=plot 
  #### save one plot by one ###########
  ggsave(plot = plotlist[[i]], file = paste(savename,".png",sep=""), width = 12, height = height_val, units = 'in', dpi = 300, path = metProfPath)
  
  dev.off()
}

