groupFluxes <- function(fluxesPruned){
# groupFluxes identifies and groups metabolites with identical flux distributions based on near-perfect correlations.
#
# USAGE:
#   results <- groupFluxes(fluxesPruned)
#
# INPUTS:
#   fluxesPruned - A data frame of metabolite fluxes, with columns representing metabolites 
#                  and rows corresponding to samples. Columns "ID" and "Sex" are excluded 
#                  from correlation analysis.
#
# OUTPUTS:
#   results      - A list containing:
#                    * flux: A modified flux data frame where metabolites with identical 
#                            distributions are grouped, retaining only one representative per group.
#                    * groups: A named list where each key corresponds to a metabolite group name 
#                              and values are the original metabolite names in that group.
#
# .. Author:
#       - Tim Hensen       November, 2024


# Step 1: Prepare data for correlation analysis
f1 <- fluxesPruned %>% select(-ID, -Sex)  # Exclude non-metabolite columns
orinames <- names(f1)  # Preserve original metabolite names
names(f1) <- make.names(names(f1))  # Standardize column names for formula compatibility

# Initialize an empty matrix to store pairwise R-squared values
metcor <- matrix(data = NA, nrow = ncol(f1), ncol = ncol(f1)) %>% as.data.frame()
rownames(metcor) <- names(f1)  # Set row names
colnames(metcor) <- names(f1)  # Set column names

# Step 2: Compute pairwise R-squared values between metabolites
for (i in 1:nrow(metcor)) {  # Loop through rows
  for (j in 1:ncol(metcor)) {  # Loop through columns
    if (i != j) {  # Avoid self-comparison
      formula <- paste0(names(f1)[i], '~', names(f1)[j])  # Create regression formula
      fit <- lm(formula, data = f1)  # Perform linear regression
      metcor[i, j] <- summary(fit)$r.squared  # Extract R-squared value
    }
  }
}

# Restore original metabolite names
colnames(metcor) <- orinames
rownames(metcor) <- orinames

# Step 3: Identify strongly correlated metabolites (R-squared > 0.999)
your_network_dataframe <- metcor %>%
  as_tibble() %>% 
  mutate(Metabolite = colnames(metcor)) %>%
  pivot_longer(cols = -Metabolite, names_to = 'Target', values_to = 'Correlation') %>%
  filter(Correlation > 0.999) %>%  # Retain only near-perfect correlations
  rename(Source = Metabolite) %>%
  filter(Source != Target) %>%  # Exclude self-loops
  mutate(link = 1:length(Source))

# Step 4: Create a graph and identify subnetwork groups
graph <- graph_from_data_frame(your_network_dataframe, directed = FALSE)  # Create an undirected graph
subnetworks <- clusters(graph)  # Identify clusters (connected components)

# Map metabolites to subnetwork clusters
subnetwork_list <- subnetworks$membership
subnetwork_df <- data.frame(link = V(graph)$name, cluster = subnetwork_list)

# Organize metabolites into groups
fluxGroups <- list()
for (i in unique(subnetwork_list)) {
  subnetwork <- subnetwork_df[subnetwork_df$cluster == i, "link"]
  fluxGroups[[i]] <- toString(subnetwork)
}

# Step 5: Add descriptive group names (if available)
groupNames <- c(
  '3-dehydro-CA/CDCA', 
  '7-dehydro-CA/CDCA', 
  'Lysine & Aminoadipic acid', 
  'Tyrosine, Dopamine, and Norepinephrine',
  'Serine metabolism, GABA biosynthesis, and Polyamine biosynthesis', 
  'L-tryptophan and downstream compounds', 
  'Iso-CA/CDCA', 
  'L-threonine & L-alpha-aminobutyrate', 
  'Urso-CA/DCA'
)

fluxGroupNames <- list()
for (i in unique(subnetwork_list)) {
  subnetwork <- subnetwork_df[subnetwork_df$cluster == i, "link"]
  fluxGroupNames[[groupNames[i]]] <- subnetwork
}

# Step 6: Replace grouped metabolites with a single representative in the flux dataset
fluxes2 <- fluxesPruned
for (i in 1:length(fluxGroupNames)) {
  newName <- names(fluxGroupNames)[i]  # Use descriptive group name
  metabolites <- fluxGroupNames[[i]]  # Get metabolites in the group
  metAdd <- fluxes2[[metabolites[1]]]  # Use the flux values from the first metabolite
  fluxes2 <- fluxes2 %>% select(-all_of(metabolites))  # Remove all grouped metabolites
  fluxes2[[newName]] <- metAdd  # Add the representative metabolite
}

# Step 7: Return the results
results <- list(flux = fluxes2, groups = fluxGroupNames)
return(results)
}
