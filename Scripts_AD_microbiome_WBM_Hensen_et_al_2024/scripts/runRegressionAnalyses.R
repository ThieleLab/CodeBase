runRegressionAnalyses <- function(fluxAll,AD_risk,Sex_moderation,e4_moderation,e2_moderation, file,type){
#' runRegressionAnalyses - Runs a series of regression analyses on metabolic fluxes or microbiome data.
#'   This function performs a series of regression analyses on the provided flux data, considering different 
#'   potential moderators (e.g., sex, APOE genotypes) and AD risk factors. It prepares the data, scales relevant 
#'   variables, and uses the performRegressions function to conduct the analyses. The results are organized into 
#'   separate groups based on the analysis type and saved to an Excel file.
#'   
#' USAGE:
#'   runRegressionAnalyses(fluxAll, AD_risk, Sex_moderation, e4_moderation, e2_moderation, file, type)
#'
#' INPUTS:
#'   fluxAll         - Data frame containing the flux data (e.g., metabolite concentrations or microbiome abundances).
#'   AD_risk         - String indicating if regressions against Alzheimer's disease (AD) risk factors should be performed ('AD_risk').
#'   Sex_moderation  - String indicating if sex interaction regressions should be performed ('Sex_moderation').
#'   e4_moderation   - String indicating if APOE4 moderation regressions should be performed ('e4_moderation').
#'   e2_moderation   - String indicating if APOE2 moderation regressions should be performed ('e2_moderation').
#'   file            - String specifying the file path to save the regression results as an Excel file.
#'   type            - String, either 'metabolome' or 'microbiome' to prepare the input data frame accordingly.
#'
#' OUTPUT:
#'   cogRegressions  - A list containing the regression results organized by the analysis type (e.g., AD risk, sex moderation).
#'
#' AUTHOR:
#'   - Tim Hensen, November 2024
#'
#' EXAMPLES:
#'   result <- runRegressionAnalyses(fluxAll, 'AD_risk', 'Sex_moderation', '', '', "results.xlsx", 'metabolome')


# Step 1: Prepare the input dataframe based on the type of data
if (type == 'metabolome') {
  fluxAll <- fluxAll %>% rename(Flux = Concentration)  # Rename 'Concentration' column to 'Flux' for metabolome data
} else if (type == 'microbiome') {
  fluxAll <- fluxAll %>% 
    mutate(Abundance = log2(Abundance)) %>%  # Log2 transform the relative abundances for microbiome data
    rename(Flux = Abundance, metabolite = Species)  # Rename 'Abundance' to 'Flux' and 'Species' to 'metabolite'
}

# Step 2: Scale the flux and age variables
fluxAll <- fluxAll %>% 
  group_by(metabolite) %>%  # Group by metabolite for scaling
  mutate(Flux = scale(Flux), AgeERGO5 = scale(AgeERGO5)) %>%  # Scale Flux and AgeERGO5 variables
  ungroup()

# Step 3: Source the regression function for performing regressions
source(here('scripts', 'performRegressions.R'))

# Step 4: Perform regression analyses if AD risk factor is specified
if (AD_risk == 'AD_risk') {
  # Perform regressions against AD risk factors
  ADreg <- list()
  ADreg[[1]] <- performRegressions(fluxAll, 'Flux ~  AgeERGO5 + Sex', 'AgeERGO5', '')  # Regression: Flux vs AgeERGO5 and Sex
  ADreg[[2]] <- performRegressions(fluxAll, 'g ~ Flux + AgeERGO5 + Sex + education + BMI', 'Flux', '')  # Regression: 'g' vs Flux, AgeERGO5, Sex, education, and BMI
  ADreg[[3]] <- performRegressions(fluxAll, 'Flux ~ GRS_APOE + Sex', 'GRS_APOE', '')  # Regression: Flux vs GRS_APOE and Sex
  AD_risk_regressions <- bind_rows(ADreg)  # Combine all AD risk regressions into a single dataframe
}

# Step 5: Perform regression analyses for sex moderation if specified
if (Sex_moderation == 'Sex_moderation') {
  # Perform regressions for sex interactions with fluxes and risk factors
  sexMod <- list()
  sexMod[[1]] <- performRegressions(fluxAll, 'Flux ~  AgeERGO5 + AgeERGO5*Sex', 'AgeERGO5:SexMale', '')  # Regression: Flux vs AgeERGO5 and interaction with Sex
  sexMod[[2]] <- performRegressions(fluxAll, 'g ~ Flux * Sex + AgeERGO5 + education + BMI', 'Flux:SexMale', '')  # Regression: 'g' vs Flux, Sex interaction, AgeERGO5, education, and BMI
  sexMod[[3]] <- performRegressions(fluxAll, 'Flux ~ GRS_APOE * Sex', 'GRS_APOE:SexMale', '')  # Regression: Flux vs GRS_APOE and interaction with Sex
  sexMod[[4]] <- performRegressions(fluxAll, 'Flux ~  AgeERGO5', 'AgeERGO5', 'Sex')  # Regression: Flux vs AgeERGO5, stratified by Sex
  sexMod[[5]] <- performRegressions(fluxAll, 'g ~ Flux + AgeERGO5 + education + BMI', 'Flux', 'Sex')  # Regression: 'g' vs Flux, stratified by Sex
  sexMod[[6]] <- performRegressions(fluxAll, 'Flux ~ GRS_APOE', 'GRS_APOE', 'Sex')  # Regression: Flux vs GRS_APOE, stratified by Sex
  AD_Sex_analyses <- bind_rows(sexMod)  # Combine all sex moderation analyses into a single dataframe
}

# Step 6: Organize regression results into a list
cogRegressions <- list()
if (AD_risk == 'AD_risk') {
  cogRegressions[['AD_risk']] <- AD_risk_regressions  # Add AD risk regressions to the result list
}
if (Sex_moderation == 'Sex_moderation') {
  cogRegressions[['Sex_moderation']] <- AD_Sex_analyses  # Add sex moderation regressions to the result list
}

# Step 7: Ungroup all results for further processing (if needed)
cogRegressions <- lapply(cogRegressions, ungroup)

# Step 8: Write the results to an Excel file
write_xlsx(cogRegressions, file)  # Save the results to an Excel file

# Return the regression results as a list
return(cogRegressions)
}