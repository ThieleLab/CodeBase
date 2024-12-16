prepareFluxes <- function(rescale, rounding) {
#' This function loads and preprocesses the flux data, including:
#' - Loading the base flux results
#' - Adding corrected bile acid fluxes
#' - Removing samples with mismatched gender
#' - Rounding fluxes to 6 decimal places
#' - Filtering on a list of metabolites of interest
#' - Calculating relative, fractional, or net microbiome influences
#' - Renaming VMH IDs to metabolite names
#
# INPUTS:
#   rescale   - String, specifies the type of flux rescaling. Options:
#                 'REL'  - Relative scaling (HM / GF).
#                 'FRAC' - Fractional scaling ((HM - GF) / HM).
#                 'NET'  - Absolute scaling (HM - GF).
#                 NULL   - No rescaling (default behavior).
#   rounding  - Integer, specifies the number of decimal places for rounding flux values.
#
# OUTPUTS:
#   A data frame containing:
#     - Sample IDs and sexes.
#     - Metabolic flux data mapped to metabolite names.
#
# .. Author:
#       - Tim Hensen       November, 2024

# Load metadata file with sample information
metadata <- read.csv(here('inputs', 'metadata.csv'))

# Load base flux data
fluxtable <- read.csv(here('inputs', 'flux_results.csv'), check.names = FALSE)

# Load corrected bile acid flux data
batable <- read.csv(here('inputs', 'Bile_acids_Corrected_table.csv'), check.names = FALSE)

# Merge bile acid data into the main flux table
fluxtable1 <- fluxtable %>% left_join(batable)

# Add specific bile acid fluxes from the new table into the original fluxtable
fluxtable$dchac <- fluxtable1$`dchac[bc]`
fluxtable$`3dhchol` <- fluxtable1$`3dhchol[bc]`
fluxtable$`3dhcdchol` <- fluxtable1$`3dhcdchol[bc]`
fluxtable$HC02194 <- fluxtable1$`HC02194[bc]`

# Process fluxtable: Filter out erroneous samples, round data, and recode sex
fluxtable <- fluxtable %>% 
  # Manually remove samples with methodological problems as found by dr. Frank J A van Rooij, Erasmus MC,
  # Department of epidemiology (see: MICROBIOME_RS_exclusions_for_Frank_20210122.csv)
  filter(!ID %in% c(5432, 7866, 7872, 7884, 7896, 7917, 6261)) %>% 
  mutate_at(vars(-c('ID', 'sex')), function(x) { round(x, digits = rounding) }) %>% # Round flux values
  mutate(sex = as.character(sex), sex = if_else(sex == 0, 'Male', 'Female')) %>% # Recode sex as 'Male' or 'Female'
  filter(ID %in% metadata$ID) # Retain samples present in metadata

# Handle rescaling based on the rescale parameter
if (rescale == 'REL') {
  gf <- read.csv(here('inputs', 'gf_fluxes.csv'), check.names = FALSE) %>%
    mutate_at(vars(-c('sex')), function(x) { round(x, digits = 6) }) %>% # Round germfree data
    mutate(sex = as.character(sex), sex = if_else(sex == 0, 'Male', 'Female')) %>%
    pivot_longer(cols = -c('sex'), names_to = 'VMHID', values_to = 'GF') # Reshape data
  
  # Calculate relative flux (HM / GF)
  fluxtable <- fluxtable %>%
    pivot_longer(cols = -c('sex', 'ID'), names_to = 'VMHID', values_to = 'HM') %>%
    right_join(gf) %>% mutate(M = HM / GF) %>% select(-HM, -GF) %>% 
    pivot_wider(names_from = VMHID, values_from = M) # Reshape back
}

if (rescale == 'FRAC') {
  gf <- read.csv(here('inputs', 'gf_fluxes.csv'), check.names = FALSE) %>%
    mutate_at(vars(-c('sex')), function(x) { round(x, digits = 6) }) %>%
    mutate(sex = as.character(sex), sex = if_else(sex == 0, 'Male', 'Female')) %>%
    pivot_longer(cols = -c('sex'), names_to = 'VMHID', values_to = 'GF')
  
  # Calculate fractional flux ((HM - GF) / HM)
  fluxtable <- fluxtable %>%
    pivot_longer(cols = -c('sex', 'ID'), names_to = 'VMHID', values_to = 'HM') %>%
    right_join(gf) %>% mutate(M = (HM - GF) / HM) %>% select(-HM, -GF) %>% 
    pivot_wider(names_from = VMHID, values_from = M)
}

if (rescale == 'NET') {
  gf <- read.csv(here('inputs', 'gf_fluxes.csv'), check.names = FALSE) %>%
    mutate_at(vars(-c('sex')), function(x) { round(x, digits = 6) }) %>%
    mutate(sex = as.character(sex), sex = if_else(sex == 0, 'Male', 'Female')) %>%
    pivot_longer(cols = -c('sex'), names_to = 'VMHID', values_to = 'GF')
  
  # Calculate net flux difference (HM - GF)
  fluxtable <- fluxtable %>%
    pivot_longer(cols = -c('sex', 'ID'), names_to = 'VMHID', values_to = 'HM') %>%
    right_join(gf) %>% mutate(M = HM - GF) %>% select(-HM, -GF) %>% 
    pivot_wider(names_from = VMHID, values_from = M)
}

# Load list of metabolites of interest
metList <- read_excel(here('inputs', 'metaboliteList.xlsx')) %>% 
  select(VMHID, metabolite)

# Replace VMHID identifiers with metabolite names and reshape table
fluxtable <- fluxtable %>%
  pivot_longer(cols = -c('ID', 'sex'), names_to = 'VMHID', values_to = 'Flux') %>%
  right_join(metList) %>% 
  select(-VMHID) %>% 
  pivot_wider(names_from = metabolite, values_from = Flux) %>% 
  select(ID, sex, sort(metList$metabolite)) %>% rename(Sex = sex)

# Return processed fluxtable
return(fluxtable)
}