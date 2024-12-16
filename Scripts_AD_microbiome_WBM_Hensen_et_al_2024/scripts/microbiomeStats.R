microbiomePhylumSummaryStats <- function(rawMicrobes, fluxAll){

## Step 1: Process the raw read counts

# Filter on flux samples
processedTotalReads <- rawMicrobes %>%
  
  # Convert wide format data to long format
  # Retains all columns except 'Taxon', pivoting others into 'ID' and 'value' columns.
  pivot_longer(cols = -c('Taxon'), names_to = 'ID') %>%
  
  # Convert back to wide format, using 'Taxon' as column names
  pivot_wider(names_from = 'Taxon') %>%
  
  # Clean up sample IDs by removing the prefix 'Sample'
  mutate(ID = gsub('Sample', '', ID)) %>%
  
  # Filter to retain only samples present in the fluxAll dataset
  filter(ID %in% fluxAll$ID) %>%
  
  # Convert back to long format for further transformations
  pivot_longer(cols = -c('ID'), names_to = 'Taxon') %>%
  
  # Replace zero values with NA to avoid including them in calculations
  mutate(value = if_else(value == 0, NA, value)) %>%
  
  # Extract species-level information from the 'Taxon' column
  group_by(ID) %>%
  mutate(
    Species = if_else(
      str_detect(Taxon, 's__'),                                      # Check if the 'Taxon' contains 's__'
      str_extract(Taxon, '(?<=g__).*') %>% gsub('; s__', ' ', .), 'NA')) %>% # Extract and clean species name
  ungroup() %>%
  
  
  
  # Identify whether each Taxon is in the mapped species list
  mutate(Species = if_else(Species == 'Ruminococcus champanellensis', 
                         'Ruminococcus chamellensis', Species),
         Species = if_else(Species == 'Prevotella__9 copri', 
                         'Prevotella_9 copri', Species),
         mapped = if_else(Species %in% mappedSpecies$mapped, 'mapped', 'unmapped')) %>%
  
  # Find the Phylum names
  mutate(Phylum = str_extract(Taxon, '(?<=p__)(.*?)(?=; c__)')) %>% select(-Taxon)


# Calculate the phylum level reads for the total and mapped reads

# First do it for all reads 
calculatePhylumReads <- function(x){
    x %>% group_by(ID, Phylum) %>% summarise(value = sum(value,na.rm = TRUE))
}

# Calculate the total phylum reads
phylumReadsTotal <- processedTotalReads %>% calculatePhylumReads()
# Calculate the mapped phylum reads
phylumReadsMapped <- processedTotalReads %>% filter(mapped=="mapped") %>% calculatePhylumReads()


## Step 2: Normalise the read counts 

normaliseReads <- function(x){
  x %>% group_by(ID) %>%
    mutate(
      # Scale each value to the total abundance for the sample
      value_norm = value / sum(value, na.rm = TRUE),
      # Remove extremely low values by replacing them with NA
      value_norm = if_else(value_norm < 1e-6, NA, value_norm))
}

# First normalise the species reads based on all reads
phylumReadsTotalnorm <- phylumReadsTotal %>% normaliseReads() %>% rename(all_reads = value_norm) %>% select(-value)
phylumReadsMappednorm <- phylumReadsMapped %>% normaliseReads() %>% rename(mapped_reads = value_norm) %>% select(-value)

# Combine data 
phylumReadsNorm <- phylumReadsTotalnorm %>% left_join(phylumReadsMappednorm)


## Step 3: Generate summary statistics

# Create function for formatting summary statistic data
meanSD <- function(x) {str_c(format(mean(x, na.rm = TRUE)*100,digits=4),' (',format(sd(x, na.rm = TRUE)*100,digits=4),')')}

raStatsPhyla <- phylumReadsNorm %>% 
  # Summarize abundance statistics for each species across all samples
  group_by(Phylum) %>%
  summarise(`R.A all reads: Mean (SD)` = meanSD(all_reads),
            `R.A mapped phyla: Mean (SD)` = meanSD(mapped_reads),
            mean_mapped = mean(mapped_reads,na.rm=TRUE)) %>%
  arrange(desc(mean_mapped)) %>% select(-mean_mapped)

return(raStatsPhyla)
}

