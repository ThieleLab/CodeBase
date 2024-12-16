# Describe microbiome abundances and create figure

# Figure 1a: boxplots of phylum abundances 
# Figure 1b: Ordered boxplots of the most abundant species

# Figure 1b:
# 1: Calculate the mean abundance of each species
# 2: Sort by abundance
# 3: Filter on the top X species
# 4: Visualise 

# Save species abundances
speciesAbundance <- 
  microbiome %>% 
  mutate(Species = gsub('_',' ', Species)) %>% 
  group_by(Species) %>% 
  summarise(relAbundance = sum(Abundance,na.rm=TRUE)/sum(microbiome$Abundance)) %>% 
  arrange(desc(relAbundance))
write.table(speciesAbundance, here('results','speciesAbundance.txt'),row.names = FALSE)

# Order the species by the total amount of reads in the data and filter the top species that represent 80% of the reads
fig1b <-
  microbiome %>% 
  mutate(Species = gsub('_',' ', Species)) %>% 
  group_by(Species) %>% 
  summarise(relAbundance = sum(Abundance,na.rm=TRUE)/sum(microbiome$Abundance)) %>% 
  arrange(desc(relAbundance)) %>%
  filter(cumsum(relAbundance)<0.8) %>% 
  arrange(relAbundance) %>%
  mutate(Species = fct_inorder(Species)) %>% 
  ggplot(aes(y=Species,x=relAbundance)) + 
  geom_col(fill='cadetblue2',col='black') +
  labs(x='Total relative read abundance across all samples',y='',
       title = 'Most abundant species in the models',
       subtitle = 'Top species that together make up 80% of the reads') + 
  theme_classic()

# Figure 1a:
# 1: obtain all microbe annotations from Microbiome statistics-Relative abundance of lost sp
# 2: Match species in microbiome with phylum annotations
# 3: Create figure on the phylum abundances

# Load microbiome annotations
annotations <- read.csv(here('Data','microbiome_annotations.csv'))

# Obtain the number of taxa per classification before and after mapping
taxa <- annotations$microbe

# Obtain number of phyla
classifications <- list(Phylum=str_extract(taxa,'(?<=p_)(.*?)(?=_c)'),
                        Class=str_extract(taxa,'(?<=c_)(.*?)(?=_o)'),
                        Order=str_extract(taxa,'(?<=o_)(.*?)(?=_f)'),
                        Family=str_extract(taxa,'(?<=f_)(.*?)(?=_g)'),
                        Genus=str_extract(taxa,'(?<=g_)(.*?)(?=_s)'),
                        Species=str_extract(taxa,'(?=_g).*'))
classifications$Species <- gsub('_g_','',classifications$Species)
classifications$Species <- gsub('_s_',' ',classifications$Species)

# Create table and add phylum information to microbiome
taxaTable <- data.frame(classifications)
Phyla <- microbiome %>% left_join(taxaTable) %>% select(ID,Phylum,Abundance)

Phyla %>% 
  group_by(Phylum) %>% 
  summarise(Abundance=sum(Abundance,na.rm = TRUE)) %>% 
  mutate(relAbundance = Abundance/sum(Abundance)*100) %>% arrange(desc(relAbundance)) %>%
  write.table(., here('results','phylumAbundance.txt'),row.names = FALSE,sep = ',')


# Visuaise phylum read abundances
fig1a <- 
  Phyla %>% 
  group_by(Phylum) %>% 
  summarise(Abundance=sum(Abundance,na.rm = TRUE)) %>% 
  mutate(relAbundance = Abundance/sum(Abundance)) %>% arrange(desc(relAbundance)) %>% 
  filter(Phylum != 'NA') %>%
  mutate(Phylum = fct_inorder(Phylum)) %>%
  ggplot(aes(x=Phylum, y = relAbundance)) + 
  geom_col(fill='coral1',col='black') +
  labs(x='',y='Relative abundance',
       title = 'Phylum read abundances across all samples') + 
  theme_classic()

design <- "AAAABB"
fig1a + fig1b + plot_layout(design=design) + plot_annotation(tag_levels = 'A')


write.table(speciesAbundance, here('results','speciesAbundance.txt'),row.names = FALSE, sep = ',')
