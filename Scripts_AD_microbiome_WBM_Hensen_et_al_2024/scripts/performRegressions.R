performRegressions <- function(flux, formula, Term, Filter) {
  # Generalized function for performing regression analyses on grouped data.
  # 
  # INPUTS:
  #   flux: A data frame containing the data to analyze.
  #   formula: A regression formula as a string (e.g., 'AgeERGO5 ~ {met} + apoe4 + BMI').
  #   Term: The term of interest in the regression model (e.g., the independent variable to focus on).
  #   Filter: A string indicating how to group the data (e.g., by 'Sex', 'apoe4', or 'apoe2').
  #
  # OUTPUT:
  # result: Data frame with regression results for each investigated compound
  # 
  # AUTHOR:
  #  - Tim Hensen, November 2024
  
  
  # Step 1: Define a helper function to perform linear regression.
  regress <- function(x, formula, Term) {
    fit <- lm(formula, data = x)  # Perform linear regression using the specified formula.
    ci95 <- confint(fit)  # Compute 95% confidence intervals for model coefficients.
    tfit <- tidy(fit)  # Summarize model results in a tidy format (from the broom package).
    mdl <- cbind(tfit, ci95) %>%  # Combine regression results with confidence intervals.
      as_tibble() %>% 
      filter(term == Term) %>%  # Filter to include only the term of interest.
      mutate(adj.R2 = summary(fit)$adj.r.squared)  # Add adjusted R-squared to the results.
    return(mdl)
  }
  
  # Step 2: Perform regression without filtering if Filter is an empty string.
  if (Filter == '') {
    result <- flux %>%
      group_by(metabolite) %>%  # Group data by each metabolite.
      do(regress(., formula, Term)) %>%  # Apply the regress function to each group.
      mutate(
        Formula = formula,  # Add the regression formula as a column.
        Cohort = rep('Complete', length(metabolite))  # Mark the cohort as 'Complete' for all results.
      ) %>%
      arrange(p.value) %>%  # Sort results by p-value.
      relocate(c('2.5 %', '97.5 %'), .after = estimate) %>%  # Move confidence intervals after the estimate column.
      relocate(Cohort, .after = term)  # Move cohort column after the term column.
    result$FDR <- p.adjust(result$p.value, method = 'fdr')  # Adjust p-values for false discovery rate (FDR).
    result <- result %>% relocate(FDR, .after = p.value)  # Relocate FDR column after p-value.
  }
  
  # Optional step: Perform regression on male and female sub cohorts.
  if (Filter == 'Sex') {
    result <- flux %>%
      group_by(metabolite, Sex) %>%  # Group data by metabolite and Sex.
      do(regress(., formula, Term)) %>%  # Apply the regress function to each group.
      mutate(Formula = formula) %>%  # Add the regression formula as a column.
      arrange(p.value) %>%  # Sort results by p-value.
      rename(Cohort = Sex) %>%  # Rename Sex column to Cohort for clarity.
      relocate(c('2.5 %', '97.5 %'), .after = estimate) %>%  # Move confidence intervals after the estimate column.
      relocate(Cohort, .after = term)  # Move cohort column after the term column.
    result$FDR <- p.adjust(result$p.value, method = 'fdr')  # Adjust p-values for FDR.
    result <- result %>% relocate(FDR, .after = p.value)  # Relocate FDR column after p-value.
  }
  
  # Optional step: Perform regression on apoe4 positive and negative sub cohorts.
  if (Filter == 'apoe4') {
    result <- flux %>%
      group_by(metabolite, apoe4) %>%  # Group data by metabolite and apoe4.
      do(regress(., formula, Term)) %>%  # Apply the regress function to each group.
      mutate(Formula = formula) %>%  # Add the regression formula as a column.
      arrange(p.value) %>%  # Sort results by p-value.
      rename(Cohort = apoe4) %>%  # Rename apoe4 column to Cohort for clarity.
      relocate(c('2.5 %', '97.5 %'), .after = estimate) %>%  # Move confidence intervals after the estimate column.
      relocate(Cohort, .after = term)  # Move cohort column after the term column.
    result$FDR <- p.adjust(result$p.value, method = 'fdr')  # Adjust p-values for FDR.
    result <- result %>% relocate(FDR, .after = p.value)  # Relocate FDR column after p-value.
  }
  
  # Optional step: Perform regression on apoe2 positive and negative sub cohorts.
  if (Filter == 'apoe2') {
    result <- flux %>%
      group_by(metabolite, apoe2) %>%  # Group data by metabolite and apoe2.
      do(regress(., formula, Term)) %>%  # Apply the regress function to each group.
      mutate(Formula = formula) %>%  # Add the regression formula as a column.
      arrange(p.value) %>%  # Sort results by p-value.
      rename(Cohort = apoe2) %>%  # Rename apoe2 column to Cohort for clarity.
      relocate(c('2.5 %', '97.5 %'), .after = estimate) %>%  # Move confidence intervals after the estimate column.
      relocate(Cohort, .after = term)  # Move cohort column after the term column.
    result$FDR <- p.adjust(result$p.value, method = 'fdr')  # Adjust p-values for FDR.
    result <- result %>% relocate(FDR, .after = p.value)  # Relocate FDR column after p-value.
  }
  
  return(result)  # Return the final result as a data frame.
}