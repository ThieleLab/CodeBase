# Figure 3: Creates volcano plots comparing metabolic fluxes and serum concentrations against cognition
#
# INPUTS:
#   fluxreg        - Data frame containing flux regression results
#   metabolome     - Data frame containing metabolomic measurements
#
# OUTPUTS:
#   Fig_3.png      - Combined volcano plots saved as PNG
#   Fig_3.svg      - Combined volcano plots saved as SVG
#
# .. Author: - Tim Hensen
#            - Annotated version created on 2024-11-18

# Figure 3: Combined volcano plots for flux and serum concentration analysis

## Figure 3a: Vulcanoplot of fluxes against cognition
cogReg <- fluxreg$AD_risk %>% 
  filter(term=='Flux') %>% 
  mutate(status = if_else(p.value<0.05, 'p<0.05','p \u2265 0.05'),
         status = if_else(FDR<0.05, 'FDR<0.05',status),
         status = factor(status,levels = c('FDR<0.05','p<0.05','p \u2265 0.05')),
         label = if_else(p.value<0.1,metabolite,'')) 

set.seed(42)

# Create vulcano plot for flux results
a <- 
  cogReg %>% 
  ggplot(aes(x=estimate,y=-log10(p.value),label = label)) +
  # Create scatter plot
  geom_point(size = 3,aes(colour = status)) +
  
  # Add labels to points
  geom_text_repel(max.overlaps = Inf) +
  
  # Add colours
  scale_colour_manual(values=c('#E32528','#0072BD','#929292')) +
  
  # Add line for p=0.05
  geom_hline(aes(yintercept = -log10(0.05)), col = "black", linetype = 'dashed', alpha = 0.4) +
  annotate("text", x = 0.5, y = -log10(0.05)+0.04, label = 'p=0.05') +

  # Add line for Beta = 0
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", alpha = 0.4) +
  # Add x-axis limits
  xlim(c(-0.5,0.5)) +
  ylim(c(0,3)) +  
  guides(color=guide_legend(ncol=length(unique(cogReg$status)), byrow=TRUE)) +
    
  # Add plot labels
  labs(
    title = '\nRegression outcomes of flux predictions against the global cognition score', 
    x= 'Regression slope', 
    y = '-log10(p-value)') +
    
  # Set theme
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.text=element_text(size=14))

# Print figure
a

## Figure 3B: Vulcanoplot of cholesterol and bile acid concentrations against cognition

# Filter on bile acids and cholesterol compounds
metabolomeFiltered <- metabolome %>% 
  filter(grepl('bile',SUB_PATHWAY,ignore.case = TRUE) |            
           grepl('cholesterol',metabolite,ignore.case = TRUE) | 
           grepl('4-cholesten-3-one',metabolite,ignore.case = TRUE) |
           grepl('7-alpha-hydroxy-3-oxo-4-cholestenoate',metabolite,ignore.case=TRUE)
  )

# Perform regressions on metabolomics
file <-here('results',str_c('Metabolome_regressions_Cholesterol_BA_age_cognition_',Sys.Date(),'.xlsx'))
metabolomereg <-runRegressionAnalyses(metabolomeFiltered,'AD_risk','','','',file,'metabolome')

# Find regression results and label them for visualisation
cogMetReg <- metabolomereg$AD_risk %>% 
  filter(term=='Flux') %>% 
  mutate(status = if_else(p.value<0.05, 'p<0.05','p \u2265 0.05'),
         status = if_else(FDR<0.05, 'FDR<0.05',status),
         status = factor(status,levels = c('FDR<0.05','p<0.05','p \u2265 0.05')),
         label = if_else(p.value<0.05,metabolite,''),
         label = gsub(' \\(1\\)|\\*','',label)) # Remove \ and 1 from metabolite names

# Create second vulcano plot
b <- cogMetReg %>% 
  ggplot(aes(x=estimate,y=-log10(p.value),label = label)) +
  
  # Create scatter plot
  geom_point(size = 3,aes(colour = status)) +
  
  # Add metabolite annotations
  geom_text_repel(max.overlaps = Inf) +
  
  # Add colours. Note that none of the results were FDR significant
  scale_colour_manual(values=c('#0072BD','#929292')) +
  
  # Add line for p=0.05 and annotate line
  geom_hline(aes(yintercept = -log10(0.05)), col = "black", linetype = 'dashed', alpha = 0.4) +
  annotate("text", x = 0.5, y = -log10(0.05)+0.04, label = 'p=0.05') +

  # Add vertical line for beta = 0
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", alpha = 0.4) +
  
  # Set x axis limits
  xlim(c(-0.5,0.5)) +
  ylim(c(0,3)) +

  # Set legend colours 
  guides(color=guide_legend(ncol=length(unique(cogMetReg$status)), byrow=TRUE)) +
  
  # Add plot titles
  labs(
    title = 'Regression outcomes of serum cholesterol and bile acid concentrations\nagainst the global cognition score',
       #subtitle = 'Cholesterol and bile acid metabolism',
       x= 'Regression slope', 
    y = '-log10(p-value)') +
  
  # Add plot theme
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.text=element_text(size=14))
b

# Combine plots and save
s <- 1.3  # Scaling factor for output dimensions
a+b +
  plot_annotation(tag_levels = 'A')
  ggsave(here('Figures','Fig_3.png'), device = 'png',width = 27*s, height = 15*s, units = 'cm')
  ggsave(here('Figures','Fig_3.svg'), device = 'svg',width = 27*s, height = 15*s, units = 'cm')
