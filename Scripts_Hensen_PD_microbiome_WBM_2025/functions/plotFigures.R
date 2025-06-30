# Create upset plot from microbiome subsets

# - Define a helper function 'using' to load and install packages if needed
using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

setwd('C://Users//mspg//Documents')

# - Load required packages
using('tidyverse','readxl','UpSetR','here','grid','pheatmap','RColorBrewer','viridis','tidytext','ggalluvial')


################################################################################
# Upset plot for intersection of flux-assciated microbial species

# Load dataset
data <- read.csv(
  here('parkinson_recreated','outputs','microbeToFlux','fluxLimiterTable.csv'),
  check.names = FALSE)

# Start r graphics png
fileName <- here('parkinson_recreated','outputs','figures','microbeUpsetPlot.png')
png(file=fileName,width = 10, height = 5, units = "in", res = 300)
upset(data, 
      sets = rev(colnames(data)[2:7]),
      #nsets = 6,
      #order.by = "freq",
      point.size = 2.5, 
      line.size = 0.75, 
      mainbar.y.label = "Intersection size", 
      sets.x.label = "Microbial species",
      #nintersects  = 40,
      #group.by = "sets"
      #empty.intersections = "off",
      mb.ratio = c(0.7, 0.3),
      #sets.bar.color = "orange",
      #main.bar.color = 'blue',
      shade.color	= 'darkblue',
      decreasing = c(F, T),
      text.scale = 1.2,
      show.numbers = 'yes',
      keep.order = TRUE
      )
grid.text("Microbial bottlenecks to predicted fluxes in blood",
          x = 0.65, y=0.98, gp=gpar(fontsize=12))
dev.off()


################################################################################
# Produce bar plots for the mean flux sensitivity upon microbial abundances

# Load dataset
data <- read.csv(
  here('parkinson_recreated','outputs','microbeToFlux','microbeSensitivityTable.csv'),
  check.names = FALSE) %>%
  # Process data and make the negative values as positive (prediction reduction in flux)
  mutate(#Mean = Mean*-1, `2.5CI`= `2.5CI`*-1, `97.5CI` = `97.5CI`*-1,
         Metabolite = factor(Metabolite,
                             levels = c('Leucylleucine','L-leucine','Butyrate','Pantothenate','Nicotinic acid','Myristic acid')), 
         Taxa = str_replace(Taxa,'_',' '))
  

# Filter on the top n microbes per metabolite
n <- 20
dataPruned <- data %>% group_by(Metabolite) %>% 
  arrange(desc(Mean),.by_group = TRUE) %>%
  slice(1:n) %>% mutate(Taxa = factor(Taxa,levels = Taxa))


# Annotate microbial species with the direction of change in PD microbiomes
annotation <- read.csv(
  here('parkinson_recreated','outputs','wallen2022_speciesRes.csv')) %>%
  select(Species, Direction) %>% rename(Taxa = Species)

dataPruned1 <- dataPruned %>% 
  left_join(annotation) %>% 
  #group_by(Metabolite) %>% 
  arrange(Metabolite,Mean)# %>%
  #mutate(Taxa = factor(Taxa)) 

# Add fix to data
dataPruned1$Direction[dataPruned1$Taxa=='Anaerobutyricum hallii'] <- 'decreased'

# Change names for legend
dataPruned2 <- dataPruned1 %>%
  mutate(Direction = if_else(Direction=='decreased','Lower in PD',Direction),
         Direction = if_else(Direction=='increased','Higher in PD',Direction),
         Direction = if_else(Direction=='unchanged','Unchanged',Direction))

# TODO: Think making a heatmap of the flux contribution potentials. Such a visualisation might make it easier to compare 
# to flux contributions to the correlations.
# tableWide <- dataPruned2 %>% select(Metabolite,Taxa,Mean) %>% pivot_wider(names_from = Metabolite,values_from = Mean)


# Create bar plots

# Start r graphics png
fileName <- here('parkinson_recreated','outputs','figures','fluxSensitivityPlot.png')
#png(file=fileName,width = 8, height = 8, units = "in", res = 450)
png(file=fileName,width = 8, height = 9, units = "in", res = 450)
#png(file=fileName,width = 9, height = 6, units = "in", res = 450)
# horizontal
ggplot(dataPruned2) +
  geom_col( aes(x=Mean, y=reorder_within(Taxa,Mean,Metabolite),fill=Direction ),colour='black') +
  geom_errorbarh(aes(y = reorder_within(Taxa,Mean,Metabolite), xmax = `97.5CI`, xmin = `2.5CI`,),height = .3) +
  scale_y_reordered() +
  labs(y = '', 
       x = 'Predicted flux contributions in mmol/day/person',
       title = 'Microbial contributions towards predicted fluxes in blood',
       subtitle = 'Top 20 highest microbial contributors',
       fill = 'Change in relative abundance\nin PD microbiomes') +
  facet_wrap(vars(Metabolite),scales = "free",ncol=2) + 
  scale_fill_viridis(discrete = TRUE, option='viridis') +
  theme_bw()+
  theme(legend.position="bottom")
dev.off()


################################################################################
# Display heatmap of flux-microbe correlations

# Load dataset
data <- read.csv(
  here('parkinson_recreated','outputs','fluxRACorrTopMicrobes.csv'),
  check.names = FALSE,row.names = 1) %>%
  rownames_to_column(var = "Species") %>% 
  mutate(Species = str_replace(Species,'_',' ')) %>%
  replace(is.na(.), 0)

data <- data %>% filter(Species != 'Flux-associated taxa')

if (F) { # Do not apply now
  # Filter on microbes that correlate above the threshold
  threshold <- 0.15
  fun <- function(x,t) {abs(data[sapply(x, is.numeric)]) > t }
  data <- data[apply(fun(data,threshold), 1, any),]
}

# Load annotation data and merge with flux correlations

# Annotate microbial species with the direction of change in PD microbiomes
annotations <- read.csv(
  here('parkinson_recreated','outputs','wallen2022_speciesRes.csv')) %>%
  select(Phylum, Species, Direction)  %>%
  rename(`Shift in PD` = Direction)

mergedData <- left_join(data,annotations)

# Split correlation data and annotations
dataProcessed <- as.data.frame(mergedData[,colnames(data)])
rownames(dataProcessed) <- mergedData$Species
dataProcessed <- dataProcessed %>% select(-Species)

annotationProcessed <- mergedData[,colnames(annotations)]

# Create pheatmap
spec <- annotationProcessed %>% select(-Species) #%>% relocate(,'Phylum',.after = 'Genus')
rownames(spec) <- annotationProcessed$Species
#spec <- spec %>% select(-Order,-Family,-Genus)# %>% relocate(,'Phylum',.after = 'Family')

numCols <- 5
mycol <- c(brewer.pal(numCols,"Blues")[numCols:1],"white",brewer.pal(numCols,"Reds")[1:numCols])
breaks <- rev(c(0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4,-0.5,-0.6))+0.05

ann_colors = list(
  `Shift in PD` = c(`Higher in PD`="#440154", `Lower in PD`="#21908C",Unchanged="#FDE725"),
  Phylum = c(Actinobacteria="#fcfdbf", Bacteroidetes="#fe9f6d", Euryarchaeota="#de4968",
             Firmicutes="#8c2981",Proteobacteria="#3b0f70",Verrucomicrobia="#000004"))

# Create heatmap and save figure
filePath <- here('parkinson_recreated','outputs','figures','fluxMicrobeCorr.png')
#png(filePath, width = 5, height = 6, units = 'in', res = 450)
png(filePath, width = 6, height = 8, units = 'in', res = 450)
plt <- pheatmap(dataProcessed,
                #color = colorRampPalette(c("blue", "white", "red"))(10),
                color = mycol,
                breaks = breaks,
                display_numbers = T,
                cluster_rows = F, cluster_cols = F,
                clustering_method="average",
                clustering_distance_rows='euclidean',
                clustering_distance_cols='euclidean',         
                fontsize_row=7,
                fontsize_col=8,
                fontsize = 6,
                annotation_row = spec,
                annotation_colors = ann_colors,
                main = "Correlations of species relative abundances with\nmetabolite blood fluxes in mmol/day/person",
)
print(plt)
dev.off()


################################################################################
# Create chord diagram of top microbe-metabolite interactions

speciesToPlot <- c('Alistipes putredinis','Bacteroides dorei','Bacteroides ovatus','Bacteroides uniformis','Bacteroides vulgatus',
                   'Blautia wexlerae','Collinsella aerofaciens','Eubacterium rectale','Faecalibacterium prausnitzii',
                   'Phocaeicola plebeius','Prevotella copri','Roseburia intestinalis','Ruthenibacterium lactatiformans')

topClusterCombinations <- read.csv(
  here('parkinson_recreated','outputs','microbeToFlux','mostImportantMicrobialLimiters.csv'),check.names=F) %>%
  rename(Species = Row) %>%
  select(-Leucylleucine) %>%
  rename(`L-leucine / Leucylleucine` = `L-leucine`) %>%
  mutate(Species = gsub('_',' ',Species)) %>%
  filter(Species %in% speciesToPlot) %>%
  pivot_longer(-Species, names_to = 'Metabolite') %>%
  filter(value!=0) 

topClusterCombinations <- read.csv(
  here('parkinson_recreated','outputs','microbeToFlux','topClusterMicrobePresenceTable.csv'),
  check.names=F) %>% 
  #rename(Species = Row) %>%
  select(-Leucylleucine) %>%
  rename(`L-leucine / Leucylleucine` = `L-leucine`) %>%
  pivot_longer(-Species, names_to = 'Metabolite') %>%
  mutate(Species = gsub('_',' ',Species),
         Metabolite = gsub('_',' ',Metabolite)) %>%
  filter(value!=0) 

# Start r graphics png
fileName <- here('parkinson_recreated','outputs','figures','topMicrobeMets_alluvium.png')
png(file=fileName,width = 10, height = 6, units = "in", res = 450)

ggplot(data = topClusterCombinations,
       aes(axis1 = Species, axis2 = Metabolite, y=value)) +
  geom_alluvium(aes(fill = Metabolite)) +
  geom_stratum(fill='grey90') +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),size=4) +
  scale_x_discrete(limits = c("Species", "Metabolite"),expand = c(0, 0),
                   labels=c("Species" = "Gut microbial species", 
                            "Metabolite" = "Predicted metabolite\nin blood")) +
  scale_fill_viridis_d() +
  labs(title = 'Predicted key influencers of metabolic productions in blood',y='') +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()
