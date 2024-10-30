
## Plots the top stratifying reaction abundances between groups in microbiome 
## model scenarios as determined through random forests analysis.

library(RColorBrewer)
library(pheatmap)
library(bigmemory)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)
setwd("..")
setwd("..")

dir.create("Heatmaps")

## define the data sets to loop through

scenarios <- c("Adults_body_sites_healthy","Adults_healthy_by_country","Adults_vs_infants_healthy","IBD_vs_healthy","Infants_premature_vs_healthy","Infants_undernourished_vs_healthy","Infection_antibiotics_vs_no_antibiotics","Infection_resistant_vs_susceptible","Infection_vs_healthy","Obesity_vs_normalweight","PD_vs_healthy","T2D_vs_healthy")

datasets <- c("Reaction_abundance_")

for (i in 1:length(scenarios)) {
  infoPath = paste("data/analysis_MicrobiomeModels/Scenarios/",scenarios[i],"/",scenarios[i],"_samples.csv",sep="")
  if (i==1) {
  info=read.csv(infoPath,header=T,row.names=1, check.names =FALSE)[c("Body site")]
  }
  else if (i==2) {
    info=read.csv(infoPath,header=T,row.names=1, check.names =FALSE)[c("Country")]
  }
  else if (i==3) {
    info=read.csv(infoPath,header=T,row.names=1, check.names =FALSE)[c("Age group")]
  }
  else if (i==7) {
    info=read.csv(infoPath,header=T,row.names=1, check.names =FALSE)[c("Antibiotics")]
  }
  else if (i==8) {
    info=read.csv(infoPath,header=T,row.names=1, check.names =FALSE)[c("Stratification")]
  }
  else {
    info=read.csv(infoPath,header=T,row.names=1, check.names =FALSE)[c("Disease name")]
  }
  colnames(info) <- "Stratification"
  
  for (j in 1:length(datasets)) {
        
      dataFile = paste("results/microbiomes/Summary_for_figures/Feature_heatmaps/",datasets[j],scenarios[i],".csv",sep="")
      data=read.csv(dataFile,header=T,row.names=1, check.names = FALSE)
      x <- as.matrix(data)
      
      infoReduced <- info
      colsToKeep <- intersect(rownames(infoReduced), colnames(data))
      infoReduced = infoReduced[colsToKeep,, drop = FALSE]
      infoReduced$Stratification <- as.factor(infoReduced$Stratification)
      
      if (length(data[,1])>1) {
      
      ## add reaction annotations where applies
        subs = read.delim("input/ReactionDatabase.txt",header=T,row.names=1, check.names =FALSE)[c("Subsystem")]
        subsReduced <- subs
        rowsToKeep <- intersect(rownames(subsReduced), rownames(data))
        subsReduced = subsReduced[rowsToKeep,, drop = FALSE]
        subsReduced$Subsystem <- as.factor(subsReduced$Subsystem)
 
      mycol <- c("white",brewer.pal(9,"Oranges")[1:4],brewer.pal(9,"Reds")[5:9])
      #mycol <- c("white",brewer.pal(9,"Reds")[1:9])
      
      ### define colors
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      palette = getPalette(length(levels(infoReduced$Stratification)))
      infoReduced.cols <- c(palette)
      names(infoReduced.cols) <- levels(infoReduced$Stratification)
      ann_colors = list(
        Stratification=infoReduced.cols
      ) 
      
      getPalette = colorRampPalette(brewer.pal(8, "Accent"))
      palette = getPalette(length(levels(subsReduced$Subsystem)))
      subsReduced.cols <- c(palette)
      names(subsReduced.cols) <- levels(subsReduced$Subsystem)
      ann_colors = list(
        Stratification=infoReduced.cols,
        Subsystem=subsReduced.cols
      ) 
      imagePath =  paste("Heatmaps/TopFeatures_",datasets[j],scenarios[i],".png",sep="")
      if (nrow(x) > 51) {
        png(imagePath, width = 14, height = 16, units = 'in', res = 300)
        pheatmap(x,
                 scale="row",
                 col=mycol,
                 #cluster_rows=FALSE,
                 #cluster_cols=FALSE,
                 clustering_method="ward.D",
                 clustering_distance_rows='euclidean',
                 clustering_distance_cols='euclidean',
                 fontsize=12,
                 border=NA,
                 annotation_col=infoReduced,
                 annotation_row=subsReduced,
                 annotation_colors=ann_colors,
                 show_legend=TRUE,
                 #show_rownames=FALSE,
                 show_colnames=FALSE
        )
        dev.off()
      }
      else 
        if (nrow(x) > 31) {
      png(imagePath, width = 12, height = 12, units = 'in', res = 300)
      pheatmap(x,
               scale="row",
               col=mycol,
               #cluster_rows=FALSE,
               #cluster_cols=FALSE,
               clustering_method="ward.D",
               clustering_distance_rows='euclidean',
               clustering_distance_cols='euclidean',
               fontsize=12,
               border=NA,
               annotation_col=infoReduced,
               annotation_row=subsReduced,
               annotation_colors=ann_colors,
               show_legend=TRUE,
               #show_rownames=FALSE,
               show_colnames=FALSE
      )
      dev.off()
      }
      else
      if (nrow(x) > 3) {
      imagePath =  paste("Heatmaps/TopFeatures_",datasets[j],scenarios[i],".png",sep="")
      png(imagePath, width = 12, height = 8, units = 'in', res = 300)
      pheatmap(x,
                scale="row",
                col=mycol,
                #cluster_rows=FALSE,
                #cluster_cols=FALSE,
                clustering_method="ward.D",
                clustering_distance_rows='euclidean',
                clustering_distance_cols='euclidean',
                fontsize=12,
                border=NA,
                annotation_col=infoReduced,
                annotation_row=subsReduced,
                annotation_colors=ann_colors,
                show_legend=TRUE,
                #show_rownames=FALSE
                show_colnames=FALSE
      )
      dev.off()
      }
    }
  }
}
