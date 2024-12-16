# participantStatistics.R
#
# Description:
#   This script calculates and formats summary statistics for participant data, grouped by sex 
#   and across the entire cohort. It produces mean, standard deviation, and sample counts for 
#   various demographic and cognitive variables, as well as statistical comparisons between sexes. 
#   Outputs include an Excel file and a LaTeX table for reporting.
#
# USAGE:
#   Run this script after preparing the `fluxes2` dataset to ensure metadata alignment.
#
# INPUTS:
#   - Metadata file: 'Data/metadataTim.csv'
#   - Sample identifiers from the `fluxes2` dataset
#
# OUTPUTS:
#   - `participant_stats.xlsx`: Summary statistics in Excel format
#   - `participant_stats.txt`: Summary table formatted for LaTeX
#
# .. Author:
#       - Your Name       November, 2024

# Align metadata with samples in `fluxes2`
metadata <- read.csv(here('inputs', 'metadata.csv')) %>% 
  select(-X) %>%  # Remove redundant column
  filter(ID %in% fluxes2$ID)  # Keep rows where ID matches `fluxes2`

# Function to calculate mean and standard deviation as a formatted string
prepareNumbers <- function(x) {
  paste0(round(mean(x, na.rm = TRUE), 2), ' (', round(sd(x, na.rm = TRUE), 2), ')')
}

# Step 1: Calculate sex-specific summary statistics
meanSD <- metadata %>%
  group_by(Sex) %>%  # Group data by sex
  summarise(
    Samples = as.character(sum(!is.na(Sex))),  # Count non-NA samples
    Age = prepareNumbers(AgeERGO5),  # Mean (SD) of age
    BMI = prepareNumbers(BMI),  # Mean (SD) of BMI
    Education = paste0(table(education), collapse = "-"),  # Distribution of education levels
    MMSE = prepareNumbers(e5_13197),  # Mini-Mental State Examination
    LDST5 = prepareNumbers(LDST5),  # Letter-Digit Substitution Test
    STR3T5 = prepareNumbers(STR3T5),  # Stroop Test
    WFT5 = prepareNumbers(WFT5),  # Word Fluency Test
    WLTdel5 = prepareNumbers(WLTdel5),  # Delayed Word List Test
    PPB_sum5 = prepareNumbers(PPB_sum5),  # Purdue Pegboard Test
    Cognition = prepareNumbers(g),  # Global cognition score
    # Count APOE genotypes
    `APOE E2/E2` = as.character(sum(apoe == '22', na.rm = TRUE)),
    `APOE E2/E3` = as.character(sum(apoe == '23', na.rm = TRUE)),
    `APOE E2/E4` = as.character(sum(apoe == '24', na.rm = TRUE)),
    `APOE E3/E3` = as.character(sum(apoe == '33', na.rm = TRUE)),
    `APOE E3/E4` = as.character(sum(apoe == '34', na.rm = TRUE)),
    `APOE E4/E4` = as.character(sum(apoe == '44', na.rm = TRUE))
  ) %>% 
  pivot_longer(cols = -Sex, names_to = 'Variable') %>%  # Reshape data to long format
  pivot_wider(names_from = Sex)  # Reshape to wide format with sex as columns

# Step 2: Calculate total cohort statistics
getNonNA <- function(x) {
  as.character(sum(!is.na(x)))  # Count non-NA values
}
Samples <- metadata %>%
  summarise(
    Samples = getNonNA(Sex),
    Age = getNonNA(AgeERGO5),
    BMI = getNonNA(BMI),
    Education = getNonNA(education),
    MMSE = getNonNA(e5_13197),
    LDST5 = getNonNA(LDST5),
    STR3T5 = getNonNA(STR3T5),
    WFT5 = getNonNA(WFT5),
    WLTdel5 = getNonNA(WLTdel5),
    PPB_sum5 = getNonNA(PPB_sum5),
    Cognition = getNonNA(g),
    # Count APOE genotypes
    `APOE E2/E2` = as.character(sum(apoe == '22', na.rm = TRUE)),
    `APOE E2/E3` = as.character(sum(apoe == '23', na.rm = TRUE)),
    `APOE E2/E4` = as.character(sum(apoe == '24', na.rm = TRUE)),
    `APOE E3/E3` = as.character(sum(apoe == '33', na.rm = TRUE)),
    `APOE E3/E4` = as.character(sum(apoe == '34', na.rm = TRUE)),
    `APOE E4/E4` = as.character(sum(apoe == '44', na.rm = TRUE))
  ) %>% 
  pivot_longer(cols = everything(), names_to = 'Variable', values_to = 'N')  # Reshape to long format

# Step 3: Perform statistical tests for sex differences
pvals <- metadata %>%
  summarise(
    Samples = NA,  # Not applicable for p-values
    Age = scientific(t.test(AgeERGO5 ~ Sex, var.equal = TRUE)$p.value),
    BMI = scientific(t.test(BMI ~ Sex, var.equal = TRUE)$p.value),
    Education = scientific(fisher.test(education, Sex)$p.value),
    MMSE = scientific(t.test(e5_13197 ~ Sex, var.equal = TRUE)$p.value),
    LDST5 = scientific(t.test(LDST5 ~ Sex, var.equal = TRUE)$p.value),
    STR3T5 = scientific(t.test(STR3T5 ~ Sex, var.equal = TRUE)$p.value),
    WFT5 = scientific(t.test(WFT5 ~ Sex, var.equal = TRUE)$p.value),
    WLTdel5 = scientific(t.test(WLTdel5 ~ Sex, var.equal = TRUE)$p.value),
    PPB_sum5 = scientific(t.test(PPB_sum5 ~ Sex, var.equal = TRUE)$p.value),
    Cognition = scientific(t.test(g ~ Sex, var.equal = TRUE)$p.value),
    `APOE E2/E2` = NA,
    `APOE E2/E3` = NA,
    `APOE E2/E4` = NA,
    `APOE E3/E3` = NA,
    `APOE E3/E4` = NA,
    `APOE E4/E4` = NA
  ) %>% 
  pivot_longer(cols = everything(), names_to = 'Variable', values_to = '{p-valuedag}')

# Step 4: Merge statistics into a single table
participantStats <- merge(merge(Samples, meanSD, sort = FALSE), pvals, sort = FALSE) %>% 
  as_tibble() %>% 
  insertRows(., c(5, 13), new = NA)  # Add blank rows for formatting

# Step 5: Update variable names for clarity
participantStats$Variable[participantStats$Variable == 'Education'] <- 'Education**'
rows <- which(participantStats$Variable %in% c('MMSE', 'LDST5', 'STR3T5', 'WFT5', 'WLTdel5', 'PPB_sum5', 'Cognition'))
participantStats$Variable[rows] <- c(
  'Mini mental state examination score***',
  'Letter-digit substitution score',
  'Stroop colour and word test score', 
  'Word fluency score test',
  'Delayed word fluency test score',
  'Left hand Purdue pegboard score',
  'Global cognition score'
)
catvars <- c('Education**', 'APOE E2/E2', 'APOE E2/E3', 'APOE E2/E4', 'APOE E3/E3', 'APOE E3/E4', 'APOE E4/E4')
participantStats$Variable[which(participantStats$Variable %in% catvars)] <- paste(catvars, '(N)')

# Step 6: Save results to Excel and LaTeX formats
write.table(participantStats, here('results', 'participant_stats.xlsx'))
latex <- xtable(participantStats, align = c('l|', 'l', 'l', 'l', 'l', 'S[table-format=2.2e1]'),
                caption = '\\textbf{Participant characteristics}* 
  \\newline
  *Values reported as the mean (standard deviation) unless otherwise indicated...', 
                label = 'table:cohort'
)
print(latex, file = here('overleaf', 'participant_stats.txt'), include.rownames = FALSE)