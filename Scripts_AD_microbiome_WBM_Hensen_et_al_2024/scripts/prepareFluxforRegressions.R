prepareFluxforRegressions <- function(fluxesPruned, SCALE) {
# Prepare flux data for regression analyses.
#
# INPUTS:
#   fluxesPruned: A data frame of pruned flux data. Expected to include 'ID' and 'Sex' columns.
#   SCALE: A logical flag indicating whether flux data should be log-transformed and scaled.
#
# OUTPUT:
#   fluxtable: A processed data frame combining flux data, metadata, and computed variables, 
#              optionally scaled for regression analysis.
#
# .. Author:
#       - Tim Hensen       November, 2024

# Step 1: Load Metadata
metadata <- read.csv(here('inputs', 'metadata.csv'))  # Load metadata containing demographic and genetic info.

# Step 2: Compute Z-scores and Principal Component Analysis (PCA)
# Identify valid rows based on non-missing values for specific metadata columns
id <- !is.na(metadata$g) & !is.na(metadata$AgeCollect) & !is.na(metadata$GRS_APOE)

# Subset and standardize specific risk-related variables using z-scores
adrisk <- metadata %>%
  select(LDST5, STR3T5, WFT5, WLTdel5, PPB_sum5, GRS_APOE) %>%  # Select columns related to AD risk
  filter(id) %>%  # Filter rows based on valid IDs
  mutate(across(everything(), scale))  # Standardize all selected columns (z-scores)

# Perform PCA on the standardized variables
pc <- prcomp(
  adrisk,  # PCA on Alzheimer's risk variables
  center = TRUE,  # Center the data
  scale. = TRUE   # Scale the data
)
# Add the first principal component (PC1) as a new variable 'AD_risk' in metadata
metadata$AD_risk[id] <- pc$x[, 1]

# Step 3: Prepare and Transform Metadata
metadata <- metadata %>%
  select(
    ID, Sex, AgeERGO5, AgeCollect, BMI, g, GRS_APOE, GRS, APOE_DIFF, apoe,
    apoePooled, apoe4, apoe2, education, AD_risk
  ) %>%  # Retain only relevant columns
  mutate(
    BMI = log10(BMI),  # Log-transform BMI
    AgeERGO5 = log10(AgeERGO5),  # Log-transform AgeERGO5
    AgeCollect = log10(AgeCollect)  # Log-transform AgeCollect
  )

# Step 4: Reshape and Merge Flux Data with Metadata
fluxtable <- fluxesPruned %>%
  pivot_longer(
    cols = -c('ID', 'Sex'),  # Reshape flux data from wide to long format, keeping 'ID' and 'Sex' intact
    names_to = 'metabolite',  # Column containing metabolite names
    values_to = 'Flux'  # Column containing flux values
  ) %>%
  left_join(metadata)  # Merge with metadata based on 'ID'

# Step 5: Log-Transform and Scale Flux Data (if SCALE is TRUE)
if (SCALE) {
  fluxtable <- fluxtable %>%
    group_by(metabolite) %>%  # Perform scaling group-wise for each metabolite
    mutate(
      Flux = log(Flux),  # Apply logarithmic transformation to flux values
      Flux = if_else(is.infinite(Flux), NaN, Flux),  # Handle cases where log transformation results in -Inf
      Flux = scale(Flux)  # Standardize flux values (z-scores)
    )
}

# Step 6: Return the Prepared Data Frame
return(fluxtable)  # Final output containing flux, metadata, and transformed variables
}
