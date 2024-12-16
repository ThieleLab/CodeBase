# Figure 2: Creates multi-panel figures showing APOE-related metabolic analyses
#
# USAGE:
#    createAPOEFigures(fluxAll, microbiome, metabolome)
#
# INPUTS:
#    fluxAll          DataFrame containing metabolic flux data
#    microbiome       DataFrame containing microbial abundance data
#    metabolome       DataFrame containing serum metabolite concentrations
#
# OUTPUTS:
#    Multiple plot objects combined into a final figure
#    Saves two files: 'Fig_2.png' and 'Fig_2.svg' in the Figures directory
#
# .. Author: - Tim Hensen
#            - Annotated version created on 2024-11-18

# Define statistical comparison groups for APOE genotypes
my_comparisons <- list(
  c('E2', 'E3'),  # Compare E2 vs E3
  c('E3', 'E4'),  # Compare E3 vs E4
  c("E2", "E4")   # Compare E2 vs E4
)

# Create custom theme for consistency across all plots
myTheme <- theme(
  legend.position = 'none',  # Remove legends as they're redundant
  axis.text = element_text(size=rel(1.1)),  # Increase axis text size
  axis.title.y = element_text(size=rel(1.2))  # Increase y-axis title size
)

## Panel A: Deoxycholate Flux Analysis
a <- fluxAll %>% 
  filter(metabolite %in% 'Deoxycholate') %>% 
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

## Panel B: S-Adenosyl-L-methionine Flux Analysis
b <- fluxAll %>% 
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

## Panel C: Eggerthella lenta Abundance Analysis
c <- microbiome %>% 
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

## Panel D: Serum Deoxycholate Concentration Analysis
d <- metabolome %>% 
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

## Additional Panels (E-I): Cholesterol-Related Metabolites Analysis
# Panel E: Serum Cholesterol
e <- metabolome %>% 
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

# Panel F: 4-cholesten-3-one
f <- metabolome %>% 
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

# Panel G: 3beta-hydroxy-5-cholestenoate
g <- metabolome %>% 
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

# Panel I: Glycolithocholate
i <- metabolome %>% 
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