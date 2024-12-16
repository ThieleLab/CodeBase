runAPOEanalyses <- function(input,type,file){
  ## runAPOEanalyses - Run comprehensive APOE genotype statistical analyses
  #
  # Performs statistical analyses of APOE genotype effects on metabolic or microbiome 
  # data, including Kruskal-Wallis tests, Dunn's post-hoc tests, and regression 
  # analyses while accounting for sex differences.
  #
  # USAGE:
  #
  #    apoeStats = runAPOEanalyses(input, type, file)
  #
  # INPUTS:
  #    input:         Data frame containing the analysis data with required columns:
  #                   * metabolites/species
  #                   * flux/concentration/abundance
  #                   * APOE genotype
  #                   * sex
  #    type:         Type of input data (char):
  #                   * 'metabolome' - metabolomic data
  #                   * 'microbiome' - microbiome data
  #    file:         File path for saving Excel results (char)
  #
  # OUTPUTS:
  #    apoeStats:    List containing statistical results:
  #                   * .kruskall_Gen - Kruskal-Wallis results for APOE genotypes
  #                   * .kruskall_Pooled - Kruskal-Wallis results for risk groups
  #                   * .Dunns_Gen - Dunn's test results for APOE genotypes
  #                   * .Dunns_Pooled - Dunn's test results for risk groups
  #                   * .AlleleComparisonReg - Regression results comparing variants
  #
  # .. Author:
  #       - Name: Tim Hensen
  #         Date: November 2024
  
  
  # Standardize column names based on input type
  if (type == 'metabolome') {
    # For metabolome data, rename Concentration to Flux
    input <- input %>% rename(Flux = Concentration)
  }
  
  if (type == 'microbiome') {
    # For microbiome data, rename Abundance to Flux and Species to metabolite
    input <- input %>% rename(Flux = Abundance, metabolite = Species)
  }
  
  # Filter out E2/E4 genotype from analysis
  input <- input %>% filter(apoe != 24)
  
  # Analysis 1: Kruskal-Wallis tests for six APOE genotypes
  GenotypeStats <-
    input %>% 
    group_by(metabolite) %>% 
    kruskal_test(Flux ~ apoe) %>%          # Non-parametric test for differences
    adjust_pvalue(method = "fdr") %>%       # Adjust for multiple comparisons
    add_significance() %>%                  # Add significance markers
    select(-c('.y.')) %>%                  # Remove unnecessary columns
    arrange(p)                             # Sort by p-value
  
  # Post-hoc analysis using Dunn's test for APOE genotypes
  DunnsStats <-
    input %>% 
    group_by(metabolite) %>% 
    dunn_test(Flux ~ apoe, detailed = TRUE) %>%             # Pairwise comparisons
    adjust_pvalue(method = "fdr") %>%       # Control false discovery rate
    add_significance() %>%
    select(-c('.y.')) %>% 
    arrange(p)
  
  # Analysis 2: Kruskal-Wallis tests for pooled risk groups (low/mid/high)
  PooledGenotypeStats <-
    input %>% 
    group_by(metabolite) %>% 
    kruskal_test(Flux ~ apoePooled) %>%     # Test pooled genotype groups
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    select(-c('.y.')) %>% 
    arrange(p)
  
  # Post-hoc Dunn's tests for pooled risk groups
  PooledDunnsStats <-
    input %>% 
    group_by(metabolite) %>%
    dunn_test(Flux ~ apoePooled, detailed = TRUE) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    select(-c('.y.')) %>% 
    arrange(p)
  
  # Analysis 3: Regression analyses accounting for sex as a covariate
  source(here('scripts','performRegressions.R'))
  
  # Compare APOE4 vs APOE2 with sex adjustment
  apoe2vsapoe4reg <- performRegressions(input, 'Flux ~ apoePooled + Sex', 'apoePooledE4', '')
  
  # Compare APOE3 vs APOE2 with sex adjustment
  apoe2vsapoe3reg <- performRegressions(input, 'Flux ~ apoePooled + Sex', 'apoePooledE3', '')
  
  # Analyze presence of E4 allele effect
  e4presence <- performRegressions(input, 'Flux ~ apoe4 + Sex', 'apoe4apoe4+', '')
  
  # Analyze presence of E2 allele effect
  e2presence <- performRegressions(input, 'Flux ~ apoe2 + Sex', 'apoe2apoe2+', '')
  
  # Combine all regression results
  apoeReg <- rbind(apoe2vsapoe4reg, apoe2vsapoe3reg, e4presence, e2presence)
  
  # Compile all analyses into a list
  apoeStats <- list(
    kruskall_Gen = GenotypeStats,
    kruskall_Pooled = PooledGenotypeStats, 
    Dunns_Gen = DunnsStats,
    Dunns_Pooled = PooledDunnsStats,
    AlleleComparisonReg = apoeReg
  )
  
  # Save results to Excel file
  write_xlsx(apoeStats, file)
  
  return(apoeStats)
}