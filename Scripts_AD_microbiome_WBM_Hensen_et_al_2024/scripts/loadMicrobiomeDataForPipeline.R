loadMicrobiomeDataForPipeline <- function(fluxAll) {
#' This function loads, processes, and filters microbiome data for use in a pipeline, 
#' ensuring compatibility with other datasets such as fluxes and metadata. The 
#' microbiome data includes species abundance, taxonomic classifications (Phylum, Class, 
#' Order, Family, Genus, Species), and associated metadata. It removes species with low 
#' representation (less than 10% of samples) and returns a cleaned dataset.
#'
#' USAGE:
#'   loadMicrobiomeDataForPipeline(fluxAll)
#'
#' INPUTS:
#'   fluxAll   - A data frame containing flux data, used to filter samples based on their IDs.
#'               Only samples present in `fluxAll` are retained.
#'
#' OUTPUTS:
#'   microbiome - A data frame containing the cleaned microbiome data, including:
#'                - Species abundance data (with NA values for low abundance)
#'                - Taxonomic classifications (Phylum, Class, Order, Family, Genus, Species)
#'                - Metadata merged with the microbiome data.
#'
#' .. Author:
#'   - Tim Hensen, November 2024
#' 
# Load and process microbiome data for analysis

# Step 1: Load metadata and remove the unwanted column 'X'
metadata <- read.csv(here('inputs','metadata.csv')) %>% select(-X)

# Step 2: Load and process microbiome species abundance data
microbiome <- read.csv(here('inputs','MARS_species.csv')) %>% 
  mutate(Species = gsub('pan','',Species)) %>%  # Clean species names by removing 'pan'
  pivot_longer(cols = -c('Species'), names_to = 'ID') %>%  # Pivot from wide to long format
  pivot_wider(names_from = 'Species') %>%  # Pivot back to wide format
  mutate(ID = gsub('Sample','',ID)) %>%  # Clean sample IDs by removing 'Sample'
  filter(ID %in% fluxAll$ID) %>%  # Filter to keep only samples present in 'fluxAll'
  pivot_longer(cols = -ID, names_to = 'Species', values_to = 'Abundance') %>%  # Reformat data
  mutate(Abundance = if_else(Abundance < 1e-6, 0, Abundance),  # Replace low abundance values with 0
         Abundance = if_else(Abundance == 0, NA, Abundance),  # Convert zeros to NA
         ID = as.numeric(ID)) %>%  # Convert 'ID' to numeric type
  filter(!is.na(Abundance)) %>%  # Remove rows with NA abundance
  mutate(Species = gsub('_',' ', Species))  # Clean species names

#

# Step 3: Load microbiome annotations, which provide taxonomic information
annotations <- read.csv(here('inputs','microbiome_annotations.csv'))

# Step 4: Extract taxonomic classifications (Phylum, Class, Order, Family, Genus, Species)
taxa <- annotations$microbe  # Get list of microbiome taxa

# Extract different levels of taxonomy using regular expressions
classifications <- list(
  Phylum = str_extract(taxa, '(?<=p_)(.*?)(?=_c)'),  # Extract phylum from taxa names
  Class = str_extract(taxa, '(?<=c_)(.*?)(?=_o)'),  # Extract class from taxa names
  Order = str_extract(taxa, '(?<=o_)(.*?)(?=_f)'),  # Extract order from taxa names
  Family = str_extract(taxa, '(?<=f_)(.*?)(?=_g)'),  # Extract family from taxa names
  Genus = str_extract(taxa, '(?<=g_)(.*?)(?=_s)'),  # Extract genus from taxa names
  Species = str_extract(taxa, '(?=_g).*')  # Extract species from taxa names
)

# Clean species names by removing extra underscores
classifications$Species <- gsub('_g_','',classifications$Species)
classifications$Species <- gsub('_s_',' ',classifications$Species)

# Step 5: Create a taxonomy table and add phylum information to the microbiome data
taxaTable <- data.frame(classifications)

# Step 6: Merge the microbiome data with the taxonomy table
microbiome <- microbiome %>% left_join(taxaTable)

# Step 7: Filter out species present in fewer than 10% of samples
dropSpec <- microbiome %>%
  group_by(Species) %>%  # Group by species
  summarise(DROP = n() < (1065 / 10)) %>%  # Calculate if species is present in fewer than 10% of samples
  filter(DROP == TRUE) %>% select(Species)  # Filter species to be removed

# Step 8: Remove low-abundance species and merge with metadata
microbiome <- microbiome %>% filter(!Species %in% dropSpec$Species) %>% left_join(metadata)

# Return the processed microbiome data
return(microbiome)
}