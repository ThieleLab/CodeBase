# Figure 4: Sex-specific cognition associations

# 4a: Stratified forest plot for sex-specific flux associations with cognition
# 4b: Stratified forest plot for Butyrivibrio crossotus and Bacteroides vulgatus

# Requirements:
# Flux regression results
# Microbe regression results

# 4a: stratified forest plot on sex-specific flux associations with cognition

# Figure 4: Creates a combined figure for sex-specific associations with cognition and
# their associated microbial correlates. 
#
# INPUTS:
#   fluxreg        - Data frame containing flux regression results
#   microbereg     - Data frame containing the microbe abundance regression results
#   metabolome     - Data frame containing metabolomic measurements
#
# OUTPUTS:
#   Fig_4.png      - Combined plots saved as PNG
#   Fig_4.svg      - Combined plots saved as SVG
#
# .. Author: - Tim Hensen
#            - Annotated version created on 2024-11-18


# Obtain metabolites to plot
mets <- fluxreg$Sex_moderation %>% filter(term == 'Flux:SexMale' & p.value<0.05) %>% select(metabolite)
mets <- mets$metabolite

# Obtain regression outcomes
fluxSexRes <- fluxreg$Sex_moderation %>% 
  filter(metabolite %in% mets & term=='Flux') %>%
  rename(lower=`2.5 %`,upper=`97.5 %`) %>%
  # Add classify metabolites by p-value
  mutate(Significance = if_else(p.value<0.05,'p<0.05','p \u2265 0.05'))

# Create forest plot
a <- fluxSexRes %>%  
  ggplot(aes(y = metabolite, x = estimate, fill = Significance, color = Cohort)) +
  
  # First add the horizontal lines
  geom_linerange(aes(xmin = lower, xmax = upper),linewidth=1.1,position=position_dodge(width = 0.5)) +
  
  # Then add the points
  geom_point(size = 3,shape=21,position=position_dodge(width = 0.5)) +
  
  # Now add a dashed line to indicate a beta of zero
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.4) +
  
  # Add x-axis limits
  xlim(c(min(fluxSexRes$lower),-1*min(fluxSexRes$lower))) +
  
  # Add titles
  labs(
    title = 'Sex-specific flux associations with the\nglobal cognition score',
    y='',
    x='Regression slope') +
  
  # Custom colors for Cohort and Significance
  scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +  # Adjust colors as needed
  scale_fill_manual(values = c("p<0.05" = "#E32528","p \u2265 0.05" = "#CDCECE")) +  # Adjust colors as needed
  
  
  # Add themes
  theme_bw() +
  theme(legend.position = "none") 
  #theme(legend.position = 'bottom',
  #      legend.direction = 'horizontal') 
a

# 4b: Explained variance of top species for sex-specific metabolites
mets[str_detect(mets,'ryptophan')] <- 'L-tryptophan'

TOPX <-fluxMicrobeReg %>% 
  select(Species,all_of(mets)) %>% 
  pivot_longer(cols=-Species,names_to = 'metabolite',values_to = 'R2') %>%
  arrange(R2) %>% 
  mutate(Species = fct_reorder(Species,R2),
         metabolite = if_else(metabolite=='L-tryptophan','L-tryptophan and\n downstream compounds',metabolite)) %>% 
  filter(R2>0.1)

# Calculate for each species the mean R2
specOrder <- TOPX %>% group_by(Species) %>% summarise(m=mean(R2)) %>% arrange(m)

b <- TOPX %>%
  # Factor the microbial species
  mutate(Species = factor(Species, levels = specOrder$Species)) %>%
  #slice_tail(n = 10) %>%
  
  # Create horizontal bar plot
  ggplot(aes(x = R2, y = Species,fill=metabolite)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  
  # Add colours
  scale_fill_viridis(discrete = TRUE, option = "E", name = "Producer of:") +
  #scale_fill_manual(values=c('#A5C8E8','#E3AEAF')) + 
  
  # Set limits on x-axis
  xlim(c(0,1)) +
  
  # Add titles
  labs(title = 'Top microbial correlates of metabolites with\nsex-specific associations with the global cognition score',
       x = 'Explained variance (R2)', 
       y = '') +

  theme_bw() + 
  theme(legend.position = "right")


# 4c: Stratified forest plot for main drivers of flux results

# Get the sex moderation results
microSexRes <- microbereg$Sex_moderation %>% 
  rename(lower=`2.5 %`,upper=`97.5 %`,Species=metabolite) %>%
  
  # Filter on metabolite associated with the top microbial species
  filter(Species %in% specOrder$Species & term=='Flux') %>%
  
  # Classify metabolites on significance
  mutate(Significance = if_else(p.value<0.05,'p<0.05','p \u2265 0.05'),
         Species = factor(Species, levels = specOrder$Species))

# Create forest plot
c <- microSexRes %>%  
  ggplot(aes(y = Species, x = estimate, fill = Significance, color = Cohort)) +

  # Add vertical line for beta = 0.
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.4) +
  #geom_vline(xintercept = 0, color = "blue", linetype = "dashed", cex = 1, alpha = 0.5) +
  
  # Add horizontal lines
  geom_linerange(aes(xmin = lower, xmax = upper),linewidth=1.1,position=position_dodge(width = 0.5)) +
  
  # Add points
  geom_point(size = 3,shape=21,position=position_dodge(width = 0.5)) +

  # Custom colors for Cohort and Significance
  scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +  # Adjust colors as needed
  scale_fill_manual(values = c("p<0.05" = "#E32528","p \u2265 0.05" = "#CDCECE")) +  # Adjust colors as needed
  
  # Add x limits
  xlim(c(min(microSexRes$lower),-1*min(microSexRes$lower))) +
  
  # Add titles
  labs(
    title = 'Sex-specific microbial associations\nwith the global cognition score',
    y='',
    x='Regression slope') +
  
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.text = element_text(size = rel(1.1)),  # Increase legend text size 
        legend.title = element_text(size = rel(1.1))) 
c


# Figure 4d: Sex-specific regressions for serum metabolite concentrations
# Obtain regression outcomes

# Manually select metabolites to visualise
mets <- c('creatine','arginine','tryptophan','kynurenine','butyrate/isobutyrate (4:0)')
metabolomeFiltered <- metabolome %>% filter(metabolite %in% mets)

# Perform regressions on metabolomics
file <-here('results',str_c('Metabolome_regressions_cognition_Sex_',Sys.Date(),'.xlsx'))
metabolomereg <-runRegressionAnalyses(metabolomeFiltered,'','Sex_moderation','','',file,'metabolome')


# Obtain regression results for metabolites of interest
metabSexRes <- metabolomereg$Sex_moderation %>% 
  filter(metabolite %in% mets & term=='Flux') %>%
  rename(lower=`2.5 %`,upper=`97.5 %`) %>%
  mutate(Significance = if_else(p.value<0.05,'p<0.05','p \u2265 0.05'))

# Make first letter uppercase
metabSexRes$metabolite <- str_to_title(metabSexRes$metabolite)

# Create forest plot for metabolites
d <- metabSexRes %>%  
  ggplot(aes(y = metabolite, x = estimate, fill = Significance, color = Cohort)) +

  # Add vertical line for beta = 0.
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.4) +
  
  # Add horizontal lines
  geom_linerange(aes(xmin = lower, xmax = upper),linewidth=1.1,position=position_dodge(width = 0.5)) +
  
  # Add regression points 
  geom_point(size = 3,shape=21,position=position_dodge(width = 0.5)) +
  
  # Custom colors for Cohort and Significance
  scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +  # Adjust colors as needed
  scale_fill_manual(values = c("p<0.05" = "#E32528","p \u2265 0.05" = "#CDCECE")) +  # Adjust colors as needed
  
  # Set x axis limits
  xlim(c(-1*max(metabSexRes$upper),max(metabSexRes$upper))) +
  
  # Add titles
  labs(
    title = 'Sex-specific serum metabolite associations\nwith the global cognition score',
    y='',
    x='Regression slope') +
  
  # Add themes
  theme_bw() +
  theme(legend.position = "none") 
d


# Define combined plot layout
layout <- "
AABB
AABB
CCBB
CCBB
CCBB
CCBB
CCDD
CCDD
"

# Combine plots and save
plt <- a+b+c+d + plot_layout(design=layout) + plot_annotation(tag_levels = 'A')
#plt <- a+b+c+d + plot_layout(design=layout, guides = 'collect') + plot_annotation(tag_levels = 'A') &
#  theme(legend.position = "right")

# Save figure
ggsave(here('Figures','Fig_4.png'), device = 'png',width = 30, height = 27, units = 'cm')
ggsave(plot = plt, here('Figures','Fig_4.svg'), device = 'svg',width = 35, height = 27, units = 'cm')



