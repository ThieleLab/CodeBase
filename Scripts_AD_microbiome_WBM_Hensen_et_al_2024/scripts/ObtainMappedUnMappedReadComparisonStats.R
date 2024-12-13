## This script inputs the raw and mapped microbiome data and produces tables on 
# 1) A table with all taxonomies and mapping information 
# 2) A table with all species level reads.
# 3) A table with 
#     a) mean+sd of relative abundances for each species
#     b) mean+sd of relative abundances after filtering on species
#     c) mean+sd of relative abundances after filtering on mapped species
#
# 4) A table with all summed phylum level reads before mapping
# 5) A table with all summed MAPPED phylum level reads
# 6) A table with all summed UNMAPPED phylum level reads
# 7) A table with 
#     a) mean+sd of relative abundances for each phylum
#     c) mean+sd of relative abundances after filtering on mapped species

### Step 1: Load raw microbiome data and obtain taxonomy information ###

# Reads a taxonomy table and arranges it by the 'Taxon' column.
rawMicrobiomes <- read.csv(here('inputs', 'microbiome_abundances_and_taxonomy.csv')) %>% arrange(Taxon)

### Step 2: Process the microbiome data for analysis ###

# Filter on fluxAll IDs
rawMicrobiomesProcessed <-  rawMicrobiomes %>% 
    pivot_longer(cols = -c('Taxon'), names_to = 'ID') %>%
    pivot_wider(names_from = 'Taxon') %>% # Make long table
    mutate(ID = gsub('Sample','',ID)) %>% 
    filter(ID %in% fluxAll$ID) %>%
    pivot_longer(cols = -c('ID'), names_to = 'Taxon',values_to = 'All') %>%
    pivot_wider(names_from = ID,values_from = All)


### Step 2: Get taxa information ###
taxa <- rawMicrobiomesProcessed$Taxon

# Split taxonomic names based on ;
taxa_split <- strsplit(taxa, ";")

# Convert list to data frame
taxa_split <- do.call(rbind, lapply(taxa_split, function(x) {
  length(x) <- max(lengths(taxa_split))  # Ensure all vectors have the same length
  x
}))

# Set column names if known
colnames(taxa_split) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove all leading characters before "__"
taxa_processed <- taxa_split %>% as_tibble %>% mutate_all(function(x){sub(".*__", "",x)})

# Complete species names
taxa_processed <- taxa_processed %>%
  mutate(Species = if_else(Species != '', str_c(Genus, ' ', Species), Species),
         Species = if_else(Species == 'Ruminococcus champanellensis', 
                           'Ruminococcus chamellensis', Species))

# Add new information to 
rawMicrobiomesProcessedTaxa <- cbind(taxa_processed, rawMicrobiomesProcessed)# %>% select(-Taxon)


### Step 3: Load mapped species and identify the mapped species ###

# Load and process mapped species
mappedSpecies <- read.csv(here('inputs', 'MARS_species.csv')) %>%
  mutate(Species = gsub('pan', '', Species),         # Remove 'pan' prefix
         Species = gsub('_', ' ', Species)) %>%      # Replace underscores with spaces
  rename(mapped = Species) %>%                      # Rename column to 'mapped'
  select(mapped)                                    # Keep only the cleaned species column

# Add mapping info
rawMicrobiomesProcessedTaxaMappedInfo <- rawMicrobiomesProcessedTaxa %>% 
  mutate(AGORA2_mapped = if_else(Species %in% mappedSpecies$mapped, 'mapped', '')) %>%
  relocate(AGORA2_mapped, .after = Species)

# Sanity check. Are all mapped species included. The result should be 116
sum(rawMicrobiomesProcessedTaxaMappedInfo$AGORA2_mapped=='mapped')

# Rename rawMicrobiomesProcessedTaxaMappedInfo
microbiome <- rawMicrobiomesProcessedTaxaMappedInfo


### Step 4: Calculate the mapping coverages ###

getReads <- function(x) {summarise(x,across(`5713`:`6848`, sum)) %>% as.numeric()}

# Obtain read counts
readCounts <- list(
  # All reads
  microbiome %>% getReads(),
  # Species reads
  microbiome %>% filter(!is.na(Species)) %>% getReads(),
  # Mapped reads
  microbiome %>% filter(AGORA2_mapped == 'mapped') %>% getReads(),
  # Unmapped reads
  microbiome %>% filter(AGORA2_mapped != 'mapped') %>% getReads(),
  # Unmapped species reads
  microbiome %>% filter(AGORA2_mapped != 'mapped' & !is.na(Species)) %>% getReads()
)

# Function to compute mean and standard deviation as formatted text
meanSD <- function(x) {
  mean_x <- mean(x) * 100
  sd_x <- sd(x) * 100
  str_c(format(mean_x, digits = 4), " (", format(sd_x, digits = 4), ")")
}

# Calculate coverage statistics
coverage <- c(
  meanSD(readCounts[[2]] / readCounts[[1]]),  # Species coverage
  meanSD(readCounts[[3]] / readCounts[[1]]),  # Mapped species coverage
  meanSD(readCounts[[4]] / readCounts[[1]]),  # Unmapped coverage
  meanSD(readCounts[[5]] / readCounts[[1]])   # Unmapped species coverage
)

# Annotate coverage
coverageAnnotation <- c(
  "Mean (SD) read coverage of named species",
  "Mean (SD) read coverage of AGORA2 mapped species",
  "Mean (SD) read coverage of AGORA2 unmapped taxa",
  "Mean (SD) read coverage of AGORA2 unmapped species"
)

# Combine into a tibble
coverages <- tibble(
  Feature = coverageAnnotation,
  `Percentage of total reads covered` = coverage
)

### Step 5: Calculate the relative abundances ###

# Calculate relative abundances and filter on fluxAll IDs
calculateRelAb <- function(data){
  data %>% 
    pivot_longer(cols = -c('ID'), names_to = 'Taxon',values_to = 'All') %>%
    # Calculate the normalised reads and remove reads below 1e-6
    group_by(ID) %>% 
    mutate(All = All/sum(All,na.rm = TRUE), # Calculate the relative abundances
           All = if_else(All<1e-6,0,All)) %>% # Set all values below 1e-6 to zero
    ungroup() %>%
    pivot_wider(names_from = ID,values_from = All)
}

# Calculate relative abundances and filter on fluxAll IDs

# Set function for calculating the relative abundances
calculateRelAbun <- function(x){
  
  # Split annotation from data
  annotationVars <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Taxon","AGORA2_mapped")
  annotation <- x %>% select(annotationVars)
  data <- x %>% select(-annotationVars) 
  # Calculate column sums
  col_sums <- colSums(data)
  
  # Handle zero columns to avoid division by zero
  col_sums[col_sums == 0] <- NA  # Replace zero sums with NA
  
  # Normalize columns
  normalized_matrix <- sweep(data, 2, col_sums, "/")
  
  # Replace NA values with 0 for columns that had zero sum
  normalized_matrix[is.na(normalized_matrix)] <- 0
  
  # Add back annotation information
  relabun <- cbind(annotation,normalized_matrix)
  
  
  relabunLong <- relabun %>% pivot_longer(cols=`5713`:`6848`,names_to = 'ID',values_to = 'Abundance') %>%
  
  return(output)
}

# Calculate the relative abundances for the full data
relabun <- calculateRelAbun(microbiome)

# Calculate the relative abundances while filtering on the species
relabunSpecies <- microbiome %>% filter(!is.na(Species)) %>% calculateRelAbun()

# Calculate the relative abundances while filtering on the mapped species
relabunMapped <- microbiome %>% filter(AGORA2_mapped=='mapped') %>% calculateRelAbun()

# Next, calculate the mean+sd of reads before and aftere mapping
# Set up function
calculateTaxonAbundance <- function(x){
  x %>%  
    group_by(TaxonLevel) %>% 
    summarise(
      `Mean relative abundance` = mean(Abundance, na.rm = TRUE) * 100,  # Calculate mean abundance for each phylum, scaled to percentage
      `SD of relative abundance` = sd(Abundance, na.rm = TRUE) * 100  # Calculate standard deviation of abundance, scaled to percentage
    ) %>% 
    
    # Order by Mean Relative Abundance
    arrange(desc(`Mean relative abundance`))
}

# Calculate phylum level reads
getPhylumReads <- function(x){
  phylumReads <- x %>% 
    # Select Relevant Columns
    select(ID, Phylum, Abundance) %>% 
    mutate(Abundance = Abundance) %>%  # Redundant line; no transformation is applied here
    
    # Group by Phylum and ID, Then Aggregate Abundance
    group_by(Phylum, ID) %>% 
    summarise(
      Abundance = sum(Abundance, na.rm = TRUE))  # Sum abundance values for each 
}

# Get stats before mapping, while only considering species, and after mapping
speciesStatsBefore <- relabun %>% rename(TaxonLevel = Species) %>% calculateTaxonAbundance()
speciesStatsJustSpecies <- relabunSpecies %>% rename(TaxonLevel = Species) %>% calculateTaxonAbundance()
speciesStatsMapped <- relabunMapped %>% rename(TaxonLevel = Species) %>% calculateTaxonAbundance()

# Change table names for stats for species when considering all taxa
speciesStatsBefore <- rename(speciesStatsBefore,
                             `Mean relative abundance before mapping` =`Mean relative abundance`,
                             `SD of relative abundance before mapping` =`SD of relative abundance`)
# Change table names for stats for species when considering only species
speciesStatsJustSpecies <- rename(speciesStatsJustSpecies,
                             `Mean relative abundance before mapping filtered on species` =`Mean relative abundance`,
                             `SD of relative abundance before mapping filtered on species` =`SD of relative abundance`)
# Change table names for stats for mapped species
speciesStatsMapped <- rename(speciesStatsMapped,
                             `Mean relative abundance after mapping` =`Mean relative abundance`,
                             `SD of relative abundance after mapping` =`SD of relative abundance`)

# Combine tables
speciesStats <- merge(speciesStatsBefore,speciesStatsJustSpecies)
speciesStats <- merge(speciesStats,speciesStatsMapped)
speciesStats <- rename(speciesStats,Species=TaxonLevel)


# Now do the same for the phylum level reads

# Get phylum level reads before and after mapping
phylumReadsBefore <- relabun %>% getPhylumReads()
phylumReadsMapped <- relabunMapped %>% getPhylumReads()

# Calculate statistics
phylumStatsBefore <- phylumReadsBefore %>% rename(TaxonLevel = Phylum) %>% calculateTaxonAbundance()
phylumStatsAfter <- phylumReadsMapped %>% rename(TaxonLevel = Phylum) %>% calculateTaxonAbundance()

# Rename variables for merging
# Change table names for stats for species when considering all taxa
phylumStatsBefore <- rename(phylumStatsBefore,
                             `Mean relative abundance before mapping` =`Mean relative abundance`,
                             `SD of relative abundance before mapping` =`SD of relative abundance`)

# Change table names for stats for mapped species
phylumStatsAfter <- rename(phylumStatsAfter,
                             `Mean relative abundance after mapping` =`Mean relative abundance`,
                             `SD of relative abundance after mapping` =`SD of relative abundance`)

# Obtain combined phylum stats table
phylumStats <- merge(phylumStatsBefore,phylumStatsAfter)







# Get phylum level reads before and after mapping
relabunPhylum <- getPhylumReads(relabun)
relabunMappedPhylum <- getPhylumReads(relabunMapped)


# Get phylum level stats before and after mapping
phylumStatsBefore <-relabunPhylum %>% rename(TaxonLevel = Phylum) %>% calculateTaxonAbundance()
phylumStatsAfter <-relabunMappedPhylum %>% rename(TaxonLevel = Phylum) %>% calculateTaxonAbundance()

# Combine tables
renameTablesBeforeMapping <- function(x){x %>% rename(`Mean relative abundance before mapping` =`Mean relative abundance`,
                                         `SD of relative abundance before mapping` =`SD of relative abundance`)}
renameTablesAfterMapping <- function(x){x %>% rename(`Mean relative abundance after mapping` =`Mean relative abundance`,
                                                      `SD of relative abundance after mapping` =`SD of relative abundance`)}

# Rename tables
phylumStatsBefore <- renameTablesBeforeMapping(phylumStatsBefore)
phylumStatsAfter <- renameTablesAfterMapping(phylumStatsAfter)

# Obtain combined phylum stats table
phylumStats <- merge(phylumStatsBefore,phylumStatsAfter)


#     a) mean+sd of relative abundances for each species
#     b) mean+sd of relative abundances after filtering on species
#     c) mean+sd of relative abundances after filtering on mapped species

# Step 5a: Calculate the mean+sd for each species





