

library(RColorBrewer)
library(pheatmap)
library(bigmemory)

currPath <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currPath)

## define the data sets to loop through

subsets <-
  c(
    "Adults_body_sites_healthy",
    "Adults_healthy_by_country",
    "Adults_vs_infants_healthy",
    "IBD_vs_healthy",
    "Infants_premature_vs_healthy",
    "Infants_undernourished_vs_healthy",
    "Infection_antibiotics_vs_no_antibiotics",
    "Infection_vs_healthy",
    "Obesity_vs_normalweight",
    "PD_vs_healthy",
    "T2D_vs_healthy"
  )

datasets <-
  c("Reaction_abundance_",
    "Reactions_presence_",
    "Subsystem_abundance_")

for (i in 1:length(subsets)) {
  infoPath = paste(
    "Analysis_microbiome_models/Subgroup_analysis/Subgroups/",
    subsets[i],
    "/",
    subsets[i],
    "_samples.csv",
    sep = ""
  )
  if (i == 1) {
    info = read.csv(
      infoPath,
      header = T,
      row.names = 1,
      check.names = FALSE
    )[c("Body site")]
  }
  else if (i == 2) {
    info = read.csv(
      infoPath,
      header = T,
      row.names = 1,
      check.names = FALSE
    )[c("Country")]
  }
  else if (i == 3) {
    info = read.csv(
      infoPath,
      header = T,
      row.names = 1,
      check.names = FALSE
    )[c("Age group")]
  }
  else if (i == 7) {
    info = read.csv(
      infoPath,
      header = T,
      row.names = 1,
      check.names = FALSE
    )[c("Antibiotics")]
  }
  else {
    info = read.csv(
      infoPath,
      header = T,
      row.names = 1,
      check.names = FALSE
    )[c("Disease name")]
  }
  colnames(info) <- "Stratification"
  
  for (j in 1:length(datasets)) {
  
    dataFile = paste(
      "Analysis_microbiome_models/Summary_for_figures/SignificantData/",
      datasets[j],
      subsets[i],
      ".csv",
      sep = ""
    )
    data = read.csv(
      dataFile,
      header = T,
      row.names = 1,
      check.names = FALSE
    )
    x <- as.matrix(data)
    
    infoReduced <- info
    colsToKeep <- intersect(rownames(infoReduced), colnames(data))
    infoReduced = infoReduced[colsToKeep, , drop = FALSE]
    colnames(infoReduced) <- "Stratification"
    infoReduced$Stratification <-
      as.factor(infoReduced$Stratification)
    
    if (length(data[, 1]) > 2) {
      ## add reaction annotations where applies
      subs = read.delim(
        "input/ReactionDatabase.txt",
        header = T,
        row.names = 1,
        check.names = FALSE
      )[c("Subsystem_general")]
      subsReduced <- subs
      rowsToKeep <-
        intersect(rownames(subsReduced), rownames(data))
      subsReduced = subsReduced[rowsToKeep, , drop = FALSE]
      subsReduced$Subsystem_general <- as.factor(subsReduced$Subsystem_general)
      
      if (j == 2) {
           mycol <- c("white",brewer.pal(9,"Reds")[1:9])
      }
      else {
        mycol <-
          c("white",
            brewer.pal(9, "Oranges")[1:4],
            brewer.pal(9, "Reds")[5:9])
      }

      ### define colors
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      palette = getPalette(length(levels(infoReduced$Stratification)))
      infoReduced.cols <- c(palette)
      names(infoReduced.cols) <- levels(infoReduced$Stratification)
      
      if (j == 1) {
        getPalette = colorRampPalette(brewer.pal(8, "Accent"))
        palette = getPalette(length(levels(subsReduced$Subsystem_general)))
        subsReduced.cols <- c(palette)
        names(subsReduced.cols) <- levels(subsReduced$Subsystem_general)
        ann_colors = list(Stratification = infoReduced.cols,
                          Subsystem_general = subsReduced.cols)
        imagePath =  paste(
          "Heatmaps/Microbiomes/SignificantFeatures_",
          datasets[j],
          subsets[i],
          ".png",
          sep = ""
        )

        png(
          imagePath,
          width = 14,
          height = 10,
          units = 'in',
          res = 300
        )
        pheatmap(
          x,
          scale = "row",
          col = mycol,
          #cluster_rows=FALSE,
          #cluster_cols=FALSE,
          clustering_method = "ward.D",
          clustering_distance_rows = 'euclidean',
          clustering_distance_cols = 'euclidean',
          fontsize = 12,
          border = NA,
          annotation_col = infoReduced,
          annotation_row = subsReduced,
          annotation_colors = ann_colors,
          show_legend = TRUE,
          show_rownames=FALSE,
          show_colnames = FALSE
        )
        dev.off()
      }
      
      if (j==2) {
        getPalette = colorRampPalette(brewer.pal(8, "Accent"))
        palette = getPalette(length(levels(subsReduced$Subsystem_general)))
        subsReduced.cols <- c(palette)
        names(subsReduced.cols) <- levels(subsReduced$Subsystem_general)
        ann_colors = list(Stratification = infoReduced.cols,
                          Subsystem_general = subsReduced.cols)
        
        imagePath =  paste(
          "Heatmaps/Microbiomes/SignificantFeatures_",
          datasets[j],
          subsets[i],
          ".png",
          sep = ""
        )
        png(
          imagePath,
          width = 10,
          height = 10,
          units = 'in',
          res = 300
        )
        pheatmap(
          x,
          #scale = "row",
          col = mycol,
          #cluster_rows=FALSE,
          #cluster_cols=FALSE,
          clustering_method = "ward.D",
          clustering_distance_rows = 'euclidean',
          clustering_distance_cols = 'euclidean',
          fontsize = 12,
          border = NA,
          annotation_col = infoReduced,
          annotation_row = subsReduced,
          annotation_colors = ann_colors,
          show_legend = TRUE,
          show_rownames=FALSE,
          show_colnames = FALSE
        )
        dev.off()
      }
      
      if (j==3) {
        ann_colors = list(Stratification = infoReduced.cols)
        
        imagePath =  paste(
          "Heatmaps/Microbiomes/SignificantFeatures_",
          datasets[j],
          subsets[i],
          ".png",
          sep = ""
        )
        
        if (length(data[, 1]) > 50) {
        png(
          imagePath,
          width = 12,
          height = 14,
          units = 'in',
          res = 300
        )
        pheatmap(
          x,
          scale = "row",
          col = mycol,
          #cluster_rows=FALSE,
          #cluster_cols=FALSE,
          clustering_method = "ward.D",
          clustering_distance_rows = 'euclidean',
          clustering_distance_cols = 'euclidean',
          fontsize = 10,
          border = NA,
          annotation_col = infoReduced,
          annotation_colors = ann_colors,
          show_legend = TRUE,
          #show_rownames=FALSE,
          show_colnames = FALSE
        )
        dev.off()
        }
        else if (length(data[, 1]) <= 50) {
          png(
            imagePath,
            width = 12,
            height = 8,
            units = 'in',
            res = 300
          )
          pheatmap(
            x,
            scale = "row",
            col = mycol,
            #cluster_rows=FALSE,
            #cluster_cols=FALSE,
            clustering_method = "ward.D",
            clustering_distance_rows = 'euclidean',
            clustering_distance_cols = 'euclidean',
            fontsize = 10,
            border = NA,
            annotation_col = infoReduced,
            annotation_colors = ann_colors,
            show_legend = TRUE,
            #show_rownames=FALSE,
            show_colnames = FALSE
          )
          dev.off()
        }
      }
    }
  }
}
