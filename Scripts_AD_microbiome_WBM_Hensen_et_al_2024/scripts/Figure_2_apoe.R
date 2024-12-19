# Figure 2: Creates multi-panel figures showing APOE-related metabolic analyses
#
# USAGE:
#    createAPOEFigures(fluxAllpruned, microbiome, metabolome)
#
# INPUTS:
#    fluxAllpruned          DataFrame containing metabolic flux data
#    microbiome       DataFrame containing microbial abundance data
#    metabolome       DataFrame containing serum metabolite concentrations
#
# OUTPUTS:
#    Multiple plot objects combined into a final figure
#    Saves two files: 'Fig_2.png' and 'Fig_2.svg' in the Figures directory
#
# .. Author: - Tim Hensen
#            - Annotated version created on 2024-11-18


updateBoxplotPvals <- function(object, input) {
# INPUTS:
# figure object: object
# Original dataset filtered on metabolite of interest

# Step 1: Calculate the correct p-values from a Dunns test

if ( any(variable.names(input)=='Flux') ) {
  dunnsRes <- dunn_test(input, Flux ~ apoePooled, detailed = TRUE)
}
  
if ( any(variable.names(input)=='Abundance') ) {
  dunnsRes <- dunn_test(input, Abundance ~ apoePooled, detailed = TRUE)
}
  
if ( any(variable.names(input)=='Concentration') ) {
  dunnsRes <- dunn_test(input, Concentration ~ apoePooled, detailed = TRUE)
}

dunnStats <- dunnsRes %>%
  select(group1, group2, p) %>% 
  arrange(c(1,3,2)) # Arrange for future merging

# Step 2: Obtain the p-value data from the figure
figureData <- ggplot_build(object)
objectData <- figureData$data[[2]]

# Find indices in the target (objectData) and source (dunnStats)

# Get target data
target <- objectData %>% select(annotation, group)

# Alter dunnStats so that it can be merged with the target
dunnStats1 <- dunnStats %>%
  mutate(group = paste(group1, group2, row_number(), sep = "-"))

# Merge with target
newData <- target %>%
  left_join(dunnStats1, by = "group") %>%
  mutate(annotation = ifelse(!is.na(p), p, annotation))

# Trim values so that the last two nonzero digits are preserved.

lastNzDigit <- floor(log10(newData$p)) # Find the first non-zero digit
roundingFactor <- abs(lastNzDigit)+1 # Add the second

# Trim values and add them to annotation
newData$annotation <- as.character(round(newData$p, digits = roundingFactor))

# Add new data to figureData
figureData$data[[2]]$annotation <- newData$annotation

# Transform figureData back to graph
objectNew <- ggplotify::as.ggplot(ggplot_gtable(figureData))
return(objectNew)
}


# Define statistical comparison groups for APOE genotypes
my_comparisons <- list(
  c('E2', 'E3'),  # Compare E2 vs E3
  c('E3', 'E4'),  # Compare E3 vs E4
  c("E2", "E4")   # Compare E2 vs E4
)

# Filter out E2/E4 genotype from analysis
fluxAllpruned <- fluxAll %>% filter(apoe != 24)

# Create custom theme for consistency across all plots
myTheme <- theme(
  legend.position = 'none',  # Remove legends as they're redundant
  axis.text = element_text(size=rel(1.1)),  # Increase axis text size
  axis.title.y = element_text(size=rel(1.2))  # Increase y-axis title size
)

## Panel A: Deoxycholate Flux Analysis
a <- fluxAllpruned %>% 
  filter(metabolite %in% 'Deoxycholate') %>% 
  # Test
  filter(!is.na(Flux)) %>%
  ggplot(aes(x=apoePooled, y=Flux, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +  # Semi-transparent outliers
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Flux in mmol/day/person',
    title = 'Deoxycholate fluxes in blood'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- fluxAllpruned %>% filter(metabolite == 'Deoxycholate')
a <- updateBoxplotPvals(a, input)


## Panel B: S-Adenosyl-L-methionine Flux Analysis
b <- fluxAllpruned %>% 
  filter(metabolite %in% 'S-Adenosyl-L-methionine') %>% 
  ggplot(aes(x=apoePooled, y=Flux, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Flux in mmol/day/person',
    title = 'S-adenosyl-L-methionine fluxes in blood'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- fluxAllpruned %>% filter(metabolite == 'S-Adenosyl-L-methionine')
b <- updateBoxplotPvals(b, input)


# Prune microbiome data
microbiomePruned <- microbiome  %>% filter(apoe != 24)

## Panel C: Eggerthella lenta Abundance Analysis
c <- microbiomePruned %>% 
  filter(Species %in% 'Eggerthella lenta') %>% 
  ggplot(aes(x=apoePooled, y=scale(log2(Abundance)), fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Normalised log2 relative abundance',
    title = expression(paste(italic('Eggerthella lenta'),' relative abundance'))
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- microbiomePruned %>% filter(Species == 'Eggerthella lenta')
c <- updateBoxplotPvals(c, input)


# Prune metabolome data
metabolomePruned <- metabolome %>% filter(apoe != 24)

## Panel D: Serum Deoxycholate Concentration Analysis
d <- metabolomePruned %>% 
  filter(metabolite %in% 'deoxycholate') %>% 
  ggplot(aes(x=apoePooled, y=Concentration, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Normalised concentration',
    title = 'Deoxycholate serum concentration'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- metabolomePruned %>% filter(metabolite == 'deoxycholate')
d <- updateBoxplotPvals(d, input)

## Additional Panels (E-I): Cholesterol-Related Metabolites Analysis
# Panel E: Serum Cholesterol
e <- metabolomePruned %>% 
  filter(metabolite %in% 'cholesterol') %>% 
  ggplot(aes(x=apoePooled, y=Concentration, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Normalised concentration',
    title = 'Cholesterol serum concentration'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- metabolomePruned %>% filter(metabolite == 'cholesterol')
e <- updateBoxplotPvals(e, input)

# Panel F: 4-cholesten-3-one
f <- metabolomePruned %>% 
  filter(metabolite %in% '4-cholesten-3-one') %>% 
  ggplot(aes(x=apoePooled, y=Concentration, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Normalised concentration',
    title = '4-cholesten-3-one serum concentration'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- metabolomePruned %>% filter(metabolite == '4-cholesten-3-one')
f <- updateBoxplotPvals(f, input)

# Panel G: 3beta-hydroxy-5-cholestenoate
g <- metabolomePruned %>% 
  filter(metabolite %in% '3beta-hydroxy-5-cholestenoate') %>% 
  ggplot(aes(x=apoePooled, y=Concentration, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Normalised concentration',
    title = '3-beta-hydroxy-5-cholestenoate\n serum concentration'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- metabolomePruned %>% filter(metabolite == '3beta-hydroxy-5-cholestenoate')
g <- updateBoxplotPvals(g, input)

# Panel H: 7-alpha-hydroxy-3-oxo-4-cholestenoate
h <- metabolome %>% 
  filter(metabolite %in% '7-alpha-hydroxy-3-oxo-4-cholestenoate (7-Hoca)') %>% 
  ggplot(aes(x=apoePooled, y=Concentration, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Normalised concentration',
    title = '7-alpha-hydroxy-3-oxo-4-cholestenoate\n serum concentration'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- metabolomePruned %>% filter(metabolite == '7-alpha-hydroxy-3-oxo-4-cholestenoate (7-Hoca)')
h <- updateBoxplotPvals(h, input)

# Panel I: Glycolithocholate
i <- metabolomePruned %>% 
  filter(metabolite %in% 'glycolithocholate') %>% 
  ggplot(aes(x=apoePooled, y=Concentration, fill=apoePooled)) + 
  geom_boxplot(outlier.alpha=0.2) +
  theme_bw() + 
  stat_compare_means(comparisons=my_comparisons, method='wilcox.test', vjust=0.15) +
  labs(
    x = '',
    y = 'Normalised concentration',
    title = 'Glycolithocholate serum\n concentrations'
  ) + 
  scale_fill_brewer(palette="Set2") + 
  myTheme

# Update p-values shown in a with those obtained from a dunns test
input <- metabolomePruned %>% filter(metabolite == 'glycolithocholate')
i <- updateBoxplotPvals(i, input)

# Define layout grid for combined plot
layout <- "
ABC
DEF
GF#
"

# Set figure dimensions
size <- 30 + 2  # 32 cm for both width and height

# Combine all plots and save
a + b + c + d + e + f + g + h + i +
  plot_annotation(tag_levels = 'A')  # Add panel labels A-I

# Save plots in both PNG and SVG formats
ggsave(here('Figures', 'Fig_2.png'), device='png', width=size, height=size, units='cm')
ggsave(here('Figures', 'Fig_2.svg'), device='svg', width=size, height=size, units='cm')