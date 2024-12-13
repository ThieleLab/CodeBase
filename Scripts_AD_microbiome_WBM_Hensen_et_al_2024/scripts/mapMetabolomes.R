mapMetabolomes <- function(fluxAll) {
# Map metabolome data onto fluxes
#
# INPUT: 
#   fluxAll: A data frame containing flux data, assumed to have an 'ID' column

# OUTPUT:
#   A data frame containing mapped metabolome data, with associated metadata

# Load metabolome data and filter based on flux results
#
# .. Author:
#       - Tim Hensen       November, 2024
  
  
metabolome <- read.csv(
  here('inputs', 'KNN_imputed_metabolon_data_3rdMarch.csv'),
  check.names = FALSE
) 
colnames(metabolome)[1] <- 'ID'  # Rename the first column to 'ID' for consistency

# Load metabolome names and metadata
metnames <- read_excel(here('inputs', 'metabolon_named.xlsx')) %>% 
  rename(VMHID = VMHId)  # Rename column for consistency

# Extract relevant columns from metabolome names for mapping
namesMap <- metnames %>% 
  select(c('CHEM_ID', 'CHEMICAL_NAME', 'VMHID', 'SUPER_PATHWAY', 'SUB_PATHWAY'))

# Generate summaries of metabolomics pathways (Super and Sub pathways)
metSpecies <- list()
metSpecies[['metabolomicsSUPERPATHWAYS']] <- metnames %>% 
  group_by(SUPER_PATHWAY) %>% 
  summarise(
    Count = n(), # Count the number of occurances for each super pathway
    Fraction = percent(n() / nrow(metnames))  # Calculate percentage of total
  ) %>% 
  arrange(desc(Fraction))  # Order by Fraction descending

metSpecies[['metabolomicsSUBPATHWAYS']] <- metnames %>% 
  # Group the data by SUPER_PATHWAY and SUB_PATHWAY
  group_by(SUPER_PATHWAY, SUB_PATHWAY) %>% 
  summarise(
    Count = n(),  # Count the number of entries in each group
    Fraction = percent(n() / nrow(metnames))  # Calculate the fraction of the total rows for each group and format it as a percentage
  ) %>% 
  # Arrange the results in descending order of Fraction
  arrange(desc(Fraction))


# Save pathway abundance summaries to an Excel file
write_xlsx(
  metSpecies,
  here('results', 'metaboliteSpeciesAbundance.xlsx')
)

# Load a metabolite list, excluding the 'Pathway' column
metList <- read_excel(
  here('inputs', 'metaboliteList.xlsx')) %>% 
  select(-Pathway)

# Load metadata and remove the first column
metadata <- read_csv(here('inputs', 'metadata.csv')) %>% 
  select(-1)

# Process and map the metabolome data
mappedMetabolome <- metabolome %>% 
  filter(ID %in% fluxAll$ID) %>%  # Filter rows where ID matches those in fluxAll
  pivot_longer(cols = -c('ID'), names_to = 'CHEM_ID') %>%  # Reshape to long format
  pivot_wider(names_from = 'ID') %>%  # Reshape back to wide format
  mutate(CHEM_ID = as.numeric(CHEM_ID)) %>%  # Convert CHEM_ID to numeric
  left_join(namesMap) %>%  # Map metabolite names and pathways
  left_join(metList) %>%  # Join with metabolite list
  pivot_longer(
    cols = -c('VMHID', 'metabolite', 'CHEM_ID', 'CHEMICAL_NAME', 'SUPER_PATHWAY', 'SUB_PATHWAY'),
    names_to = 'ID',
    values_to = 'Concentration'
  ) %>%  # Reshape to long format
  select(-metabolite) %>% 
  rename(metabolite = CHEMICAL_NAME) %>%  # Rename for clarity
  mutate(ID = as.numeric(ID)) %>%  # Convert ID to numeric
  filter(!is.na(Concentration)) %>%  # Filter out rows with missing Concentration
  left_join(metadata)  # Add metadata

return(mappedMetabolome)  # Return the final mapped data
}