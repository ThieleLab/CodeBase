regressFluxesAgainstMicrobes <- function(fluxesPruned,microbiome) {
#' regressFluxesAgainstMicrobes - Regresses microbe abundances against metabolic fluxes.
#'   This function performs linear regressions between microbe abundances and metabolic
#'   fluxes for each microbe and each flux. The R-squared values from these regressions 
#'   are stored in a result table, which is filtered to include only those correlations 
#'   between microbes and metabolites that the microbes are capable of producing.
#'   
#' USAGE:
#'   regressFluxesAgainstMicrobes(fluxesPruned, microbiome)
#'
#' INPUTS:
#'   fluxesPruned  - Data frame of pruned fluxes with samples as rows and fluxes as columns.
#'   microbiome    - Data frame of microbiome data with species abundance for each sample.
#'
#' OUTPUT:
#'   result_table1 - Data frame containing R-squared values from linear regressions
#'                  between each microbe and flux.
#'
#' .. Author:
#'   - Tim Hensen, November 2024

# Prepare the microbiome data: select species abundance and arrange by sample ID
microbes <- microbiome %>% 
  select(ID, Species, Abundance) %>%
  pivot_wider(names_from = Species, values_from = Abundance) %>%
  arrange(ID) %>%
  select(-ID)

# Prepare the flux data: arrange by sample ID and remove non-relevant columns
flux <- fluxesPruned %>%
  arrange(ID) %>%
  select(-ID, -Sex)

# Initialize an empty result table with rows for microbes and columns for fluxes
result_table <- matrix(data = NA, nrow = ncol(microbes), ncol = ncol(flux)) %>%
  as.data.frame()

rownames(result_table) <- names(microbes)
colnames(result_table) <- names(flux)

# Loop through each microbe (rows) and each flux (columns)
for (i in 1:nrow(result_table)) {
  for (j in 1:ncol(result_table)) {
    
    # Combine the flux data and microbe abundance data for regression
    DATA <- cbind(flux[, j], microbes[, i]) %>% as_tibble()
    names(DATA) <- make.names(names(DATA))
    
    # Create a formula for the linear regression
    formula <- paste0(names(DATA)[1], '~', names(DATA)[2])
    
    # Perform the linear regression and store the R-squared value
    fit <- lm(formula, data = DATA)
    result_table[i, j] <- summary(fit)$r.squared
  }
}

# Load metabolite-microbe links and process species names for consistency
metMicrobeLinks <- read_excel(here('inputs', 'metaboliteMicrobeLinks.xlsx')) %>%
  rename(Species = Row) %>%
  mutate(Species = gsub('pan', '', Species),
         Species = gsub('_', ' ', Species)) %>%
  pivot_longer(cols = -Species, names_to = 'VMHID') %>%
  left_join(
    read_excel(here('inputs', 'metaboliteList.xlsx')) %>%
      select(metabolite, VMHID)
  ) %>%
  filter(!is.na(metabolite)) %>%
  select(-VMHID) %>%
  pivot_wider(names_from = metabolite)

# Add known metabolite-microbe links not included in the original data
metMicrobeLinks <- metMicrobeLinks %>%
  mutate(`S-Adenosyl-L-methionine` = `L-methionine`,
         Creatine = `L-arginine`,
         Formaldehyde = rep(1, times = length(metMicrobeLinks$Formate)))

# Filter result_table by removing microbes that cannot produce a metabolite
result_table1 <- result_table %>%
  as_tibble(rownames = 'Species') %>%
  pivot_longer(cols = -Species, names_to = 'metabolite', values_to = 'Corr') %>%
  left_join(metMicrobeLinks %>% pivot_longer(cols = -Species, names_to = 'metabolite')) %>%
  mutate(Corr = if_else(value == 0, 0, Corr)) %>%
  select(-value) %>%
  pivot_wider(names_from = metabolite, values_from = Corr)

# Return the filtered result table with R-squared values for valid metabolite-microbe links
return(result_table1)
}