
library(RColorBrewer)
library(pheatmap)
library(bigmemory)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)

## define the data sets to loop through

#versions <- c("150k","90k","150k_90k_combined")
versions <- c("150k_90k_combined")

infoFiles <- c("ExchangeAnnotations.txt","MetaboliteDatabase.txt","ReactionDatabase.txt")

datasets <- c("_uptake_secretion_results_","_internal_production_results_","_reaction_presence_results_")

#taxa <- c("Phylum","Class","Order","Family","Genus","Species")
taxa <- c("Phylum")

for (i in 1:length(versions)) {
  
  for (j in 1:length(datasets)) {
    infoPath = paste("input/",infoFiles[j],sep="")
    info=read.delim(infoPath,header=T,row.names=1, check.names =FALSE)[c("Subsystem")]
    
    for (k in 1:length(taxa)) {
      dataFile = paste("Model_properties_analysis/Summary_for_figures/Feature_heatmaps/",versions[i],datasets[j],taxa[k],".csv",sep="")
      data=read.csv(dataFile,header=T,row.names=1, check.names = FALSE)
      x <- as.matrix(data)
      
      infoReduced <- info
      rowsToKeep <- intersect(rownames(infoReduced), colnames(data))
      infoReduced = infoReduced[rowsToKeep,, drop = FALSE]
      infoReduced$Subsystem <- as.factor(infoReduced$Subsystem)
 
      if (j==1) {
        mycol <- c(brewer.pal(9,"Blues")[9:1],"white",brewer.pal(9,"Reds")[1:9])
      }
      else {
      mycol <- c("white",brewer.pal(9,"Reds")[1:9])
      }
      
      ### define colors
      getPalette = colorRampPalette(brewer.pal(8, "Set1"))
      palette = getPalette(length(levels(infoReduced$Subsystem)))
      infoReduced.cols <- c(palette)
      names(infoReduced.cols) <- levels(infoReduced$Subsystem)
      ann_colors = list(
       Subsystem=infoReduced.cols
      ) 
      
      imagePath =  paste("Heatmaps/Strains/TopFeatures_",versions[i],datasets[j],taxa[k],".png",sep="")
      
        png(imagePath, width = 11, height = 6, units = 'in', res = 300)
        pheatmap(x,
                 #scale="row",
                 col=mycol,
                 #cluster_rows=FALSE,
                 #cluster_cols=FALSE,
                 clustering_method="ward.D",
                 clustering_distance_rows='euclidean',
                 clustering_distance_cols='euclidean',
                 fontsize=10,
                 border=NA,
                 annotation_col=infoReduced,
                 annotation_colors=ann_colors,
                 show_legend=TRUE
                 #show_rownames=FALSE,
                 #show_colnames=FALSE
        )
        dev.off()
    }
  }
}
