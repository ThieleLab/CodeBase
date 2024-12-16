# Figure 1: Multi-panel visualization of metabolic analysis
# This script creates a comprehensive figure with four panels:
# a) Volcano plot of age-related flux changes
# b) Bar plot showing microbial correlations with metabolites
# c) Forest plot of species abundance vs age
# d) Regression plot of serum arginine vs age
#
# INPUTS:
#    fluxreg          Structure containing flux regression results
#    metabolome       DataFrame containing metabolomics data
#    fluxMicrobeReg   DataFrame containing flux-microbe regression results
#    microbereg       Structure containing microbe regression results
#
# OUTPUTS:
#    plt             Combined plot object with all panels
#    Saves two files: 'Fig_1.png' and 'Fig_1.svg' in the Figures directory
#
# .. Author: - Tim Hensen
#            - Annotated version created on 2024-11-18


## Panel A: Volcano plot preparation
# Filter and prepare data for age-related flux changes
ageReg <- fluxreg$AD_risk %>% 
  filter(term == 'AgeERGO5') %>%      
  mutate(
    # Create significance categories based on p-values and FDR
    status = if_else(p.value < 0.05, 'p<0.05', 'p \u2265 0.05'),
    status = if_else(FDR < 0.05, 'FDR<0.05', status),
    status = factor(status, levels = c('FDR<0.05', 'p<0.05', 'p \u2265 0.05')),
    # Add labels for significant metabolites (p < 0.1)
    label = if_else(p.value < 0.1, metabolite, '')
  )

# Set seed for reproducible label positioning
set.seed(42)

# Create volcano plot with advanced formatting
a <- ageReg %>% 
  ggplot(aes(x = estimate, y = -log10(p.value), label = label)) +
  # Add points with significance-based coloring
  geom_point(size = 3, aes(colour = status)) +
  # Use a red-grey color scheme for significance levels
  scale_colour_manual(values = c('#E32528','#0072BD','#929292')) + 
  # Add non-overlapping labels
  geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.3,
    box.padding = 0.8,
    fontface = "plain",
    nudge_y = 0.15
  ) + 
  coord_cartesian(clip = "off") +
  # Add significance threshold line
  geom_hline(
    aes(yintercept = -log10(0.05)),
    col = "black",
    linetype = 'dashed',
    alpha = 0.4
  ) +
  # Add threshold label
  annotate(
    "text",
    x = max(ageReg$estimate) * 0.98,
    y = -log10(0.05) + 0.15,
    label = 'p=0.05'
  ) +
  # Customize plot appearance
  guides(color = guide_legend(
    ncol = length(unique(ageReg$status)),
    byrow = TRUE
  )) +
  labs(
    title = '\nPredicted flux value associations with age',
    x = 'Regression slope',
    y = '-log10(p-value)'
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

## Panel B: Bar plot preparation
# Filter and prepare data for microbial correlations
TOPX <- fluxMicrobeReg %>% 
  select(Species, `L-arginine`, `3-Dehydrocholate`) %>% 
  pivot_longer(cols = -Species, names_to = 'metabolite', values_to = 'R2') %>%
  arrange(R2) %>% 
  mutate(Species = fct_reorder(Species, R2)) %>% 
  filter(R2 > 0.1) %>%
  mutate(metabolite = gsub('3-Dehydrocholate', '3-dehydro-CA/CDCA', metabolite))

# Create bar plot
b <- TOPX %>%
  ggplot(aes(x = R2, y = fct_reorder(Species, R2), fill = metabolite)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  scale_fill_viridis(discrete = TRUE, option = "E", name = "Producer of:") +
  labs(
    title = 'Top microbial correlates of L-arginine \n and 3-dehydro-CA/CDCA fluxes',
    x = 'Explained variance (R2)',
    y = ''
  ) +
  xlim(c(0,1)) +
  theme_bw() +
  theme(legend.position = "right")

## Panel C: Forest plot preparation
# Create forest plot for selected species
c <- microbereg$AD_risk %>% 
  filter(term == 'AgeERGO5') %>% 
  rename(lower = `2.5 %`, upper = `97.5 %`, Species = metabolite) %>%
  filter(Species %in% TOPX$Species) %>%
  mutate(
    Species = factor(Species, levels = TOPX$Species),
    sig = if_else(p.value < 0.05, 'p<0.05', 'p \u2265 0.05')
  ) %>% 
  ggplot(aes(y = Species, x = estimate, fill = sig)) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_point(size = 3, shape = 21, fill = "#CDCECE") +
  geom_vline(xintercept = 0, color = "grey50", linetype = "solid", alpha = 0.5) +
  labs(
    title = 'Log2 relative species abundance \nassociated with age',
    x = 'Regression estimate',
    y = ''
  ) +
  theme_bw() +
  theme(legend.position = "right", legend.title = element_blank())

## Panel D: Regression plot preparation
# Filter arginine concentrations and create scatter plot
arginineConcentrations <- metabolome %>% filter(metabolite == 'arginine')
d <- arginineConcentrations %>% 
  ggplot(aes(x = AgeERGO5, y = Concentration)) +
  geom_point(color = "grey20", shape = 16, size = 3, alpha = 0.7) +
  stat_poly_line(color = "darkblue", size = 0.8, fill='yellow') +
  stat_poly_eq(
    aes(label = paste("p =", round(..p.value.., 3), ", R² =", round(..r.squared.., 4))),
    label.x = 1,
    label.y = 0.02,
    size = 4,
    color = "black",
    parse = FALSE
  ) +
  labs(
    title = 'Serum arginine concentration against age',
    y = 'Normalised serum\n concentration',
    x = 'Age in years'
  ) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust = -40))

# Define layout for combined plot
layout <- "
AAABB
AAACC
AAADD
"

# Combine all plots using patchwork
plt <- a + b + c + d +
  plot_annotation(tag_levels = 'A') +
  plot_layout(design = layout)

# Save plots in both PNG and SVG formats
# Define dimensions
h <- 23
w <- 35

# Save files
ggsave(here('Figures', 'Fig_1.png'), device = 'png', width = w, height = h, units = 'cm')
ggsave(plot = plt, here('Figures', 'Fig_1.svg'), device = 'svg', width = w, height = h, units = 'cm')