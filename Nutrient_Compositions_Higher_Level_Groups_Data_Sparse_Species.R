###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#   Nutrient Compositions for Higher-Level Groups and Data-Sparse Species     #
#                                                                             #
###############################################################################


# After we generated predicted nutrient compositions using the predicted ratios and simulations described in the Predicive Models script, we averaged our observed and predicted nutrient data at the species level to generate C, N, and P compositions. This was done to reduce the influence of over-represented species in our higher-level taxa calculations. 

# Then, for each taxonomic level above species, we calculated the mean and standard deviation of nutrient composition values based on each taxonomic level's corresponding species. 

# To account for compounding uncertainty among our observed and predicted values, we calculated a compounded standard deviation for each taxa.



#### NECESSARY PACKAGES ####

library(tidyverse)
library(data.table)
library(GGally)
library(nlme)
library(lme4)
library (MuMIn)
library(ggpubr)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(treemapify)




#### AVERAGING SPECIES-LEVEL DATA ####

##### Creating the Data Frame for Nutrient Composition Aggregation #####

# The following code generates a data frame that will be populated with the mean nutrient composition values for all species. Note: This process is specific to the species level. Due to the uneven distribution of nutrient values across species, we chose to first calculate averages at the species level. This approach ensures that when nutrient compositions are averaged at higher taxonomic levels, no single species disproportionately influences the overall average for any group.


# Creates data frame for averaging values.
NutrientContent_forAveraging <- Observed_and_Predicted_NutrientContent %>% select(1:17, 36:37, 40:41, 44:45)


# We had some data that was collected for species not listed in the Sea Around Us (SAU) dataset but that represented higher-level groups (e.g., families) that did have catch data in the SAU dataset. The following code is to make sure to extract that data.
NutrientContent_forAveraging <- NutrientContent_forAveraging %>% mutate(genus_species = paste(genus, species, sep = " "))


# The code below creates the species-level data frame for the mean nutrient compositions.
AggregatedNutrientContent_species <- NutrientContent_forAveraging %>% select (2:9)

# The code below combines the SAU Taxa Table with our data frame to check for correct taxonomic labels between our dataset and the SAU dataset.
AggregatedNutrientContent_species <- merge(AggregatedNutrientContent_species, SAU_TaxonTable[ , c(1, 4:6, 8, 10)], 
                                           by = "TaxonKey", all.x = T)


AggregatedNutrientContent_species <- AggregatedNutrientContent_species %>% 
  select("TaxonKey", "TaxonName", "CommonName", "phylum", "sub_phylum", "super_class",
         "class", "super_order", "order", "suborder_infraorder", "family", "genus", "species")


AggregatedNutrientContent_species <- AggregatedNutrientContent_species %>% mutate(genus_species = paste(genus, species, sep = " "))


# The code below breaks down the table so that any species that were accidentally repeated were eliminated.
AggregatedNutrientContent_species <- AggregatedNutrientContent_species[!duplicated(AggregatedNutrientContent_species$genus_species),]

# The code below changes the column names so that it matches how it appears in the SAU data frame.
colnames(AggregatedNutrientContent_species)[2] <- "scientific_name"
colnames(AggregatedNutrientContent_species)[3] <- "common_name"


# Creates the columns that will hold our mean composition values per species.
AggregatedNutrientContent_species$mean_WW_C <- NA
AggregatedNutrientContent_species$compound_WW_C_SD <- NA
AggregatedNutrientContent_species$mean_WW_N <- NA
AggregatedNutrientContent_species$compound_WW_N_SD <- NA
AggregatedNutrientContent_species$mean_WW_P <- NA
AggregatedNutrientContent_species$compound_WW_P_SD <- NA



##### For Loops for Averaging Nutrient Compositions - Species Level #####

# The for loops below average the nutrient composition values at the species level for all three nutrients. The equation used for the compound variance, and subsequent compound SD, can be found at https://www.emathzone.com/tutorials/basic-statistics/combined-variance.html


# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_species)){
  # This `ifelse` bracket check if, for row Z in the dataset, there is a C composition value available for the listed species. Otherwise, the for loop returns an 'NA'.
  {if(AggregatedNutrientContent_species$genus_species[Z] %in% 
      NutrientContent_forAveraging[!is.na(combined_WetWeightC)]$genus_species) {
    
    # The following code generates a data frame with the set of values that are available for the species listed in row Z.
    use_df <- NutrientContent_forAveraging[NutrientContent_forAveraging$genus_species == 
                                             AggregatedNutrientContent_species[Z]$genus_species]
    
    # Obtains the mean C composition value for the species in row Z.
    AggregatedNutrientContent_species$mean_WW_C[Z] <- mean(use_df$combined_WetWeightC, na.rm = T)
    
    # Obtains the compounded SD for the C composition values for the species in row Z.
    comp_SD <- sum((use_df$combined_WetWeight_C_SD^2 + (use_df$combined_WetWeightC - AggregatedNutrientContent_species$mean_WW_C[Z])^2), 
             na.rm = T)/nrow(use_df)
    
    AggregatedNutrientContent_species$compound_WW_C_SD[Z] <- sqrt(comp_SD)
    
  } else (NA)}
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_species)){
  # This `ifelse` bracket check if, for row Z in the dataset, there is a N composition value available for the listed species. Otherwise, the for loop returns an 'NA'.
  {if(AggregatedNutrientContent_species$genus_species[Z] %in% 
      NutrientContent_forAveraging[!is.na(combined_WetWeightN)]$genus_species) {
    
    # The following code generates a data frame with the set of values that are available for the species listed in row Z.
    use_df <- NutrientContent_forAveraging[NutrientContent_forAveraging$genus_species == 
                                             AggregatedNutrientContent_species[Z]$genus_species]
    
    # Obtains the mean N composition value for the species in row Z.
    AggregatedNutrientContent_species$mean_WW_N[Z] <- mean(use_df$combined_WetWeightN, na.rm = T)
    
    # Obtains the compounded SD for the N composition values for the species in row Z.
    comp_SD <- sum((use_df$combined_WetWeight_N_SD^2 + (use_df$combined_WetWeightN - AggregatedNutrientContent_species$mean_WW_N[Z])^2), 
             na.rm = T)/nrow(use_df)
    
    AggregatedNutrientContent_species$compound_WW_N_SD[Z] <- sqrt(comp_SD)
  }}
}



# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_species)){
  # This `ifelse` bracket check if, for row Z in the dataset, there is a P composition value available for the listed species. Otherwise, the for loop returns an 'NA'.
  {if(AggregatedNutrientContent_species$genus_species[Z] %in% 
      NutrientContent_forAveraging[!is.na(combined_WetWeightP)]$genus_species) {
    
    # The following code generates a data frame with the set of values that are available for the species listed in row Z.
    use_df <- NutrientContent_forAveraging[NutrientContent_forAveraging$genus_species == 
                                             AggregatedNutrientContent_species[Z]$genus_species]
    
    # Obtains the mean P composition value for the species in row Z.
    AggregatedNutrientContent_species$mean_WW_P[Z] <- mean(use_df$combined_WetWeightP, na.rm = T)
    
    # Obtains the compounded SD for the P composition values for the species in row Z.
    comp_SD <- sum((use_df$combined_WetWeight_P_SD^2 + (use_df$combined_WetWeightP - AggregatedNutrientContent_species$mean_WW_P[Z])^2), 
             na.rm = T)/nrow(use_df)
    
    AggregatedNutrientContent_species$compound_WW_P_SD[Z] <- sqrt(comp_SD)
  }}
}




#### AVERAGING HIGHER-LEVEL TAXA DATA ####


##### Creating the Data Frame for Nutrient Composition Aggregation #####

# The code below creates the data frame to fill out for all species and taxonomic groups that are included in the industrial SAU dataset.

# This line of code takes the species that were actually included in the fisheries dataset to the data frame we will create below. This does NOT include the species that we collected data for but for which there were no direct entries in the SAU dataset.
AggregatedNutrientContent_speciesonly <- AggregatedNutrientContent_species %>% filter(TaxonKey >= 600000)

# This code quickly subsets all the taxa in the SAU Taxon Table that is included among industrial catches. This code also works to retain only the taxonomic information we desired for each taxa.
AggregatedNutrientContent_AllLevels <- SAU_TaxonTable[ , c(1:13)]
AggregatedNutrientContent_AllLevels <- AggregatedNutrientContent_AllLevels[TaxonKey %in% SAU_Industrial_WithTaxaInfo$TaxonKey]

# This final bit adds the species that we have already averaged nutrient compositions for as well as add the columns to be filled out below.
AggregatedNutrientContent_AllLevels <- merge(AggregatedNutrientContent_AllLevels, 
                                             AggregatedNutrientContent_speciesonly[ , c(1, 15:20)], by = "TaxonKey", all.x = T)




##### Genus-Level #####

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the genera listed in the SAU taxa table (TaxonKey >= 500000 & <= 599999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 500000, 599999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a C composition value available for the listed genus. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$genus[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$genus) {
      
      # The following code generates a data frame with the set of species-level C values that are available for the genus listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$genus ==
                                                    AggregatedNutrientContent_AllLevels$genus[Z]]
      
      # Obtains the mean C composition value for the genus in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      # Obtains the compounded SD for the C composition values for the genus in row Z.
      comp_SD <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}




# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the genera listed in the SAU taxa table (TaxonKey >= 500000 & <= 599999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 500000, 599999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a N composition value available for the listed genus. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$genus[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$genus) {
      
      # The following code generates a data frame with the set of species-level N values that are available for the genus listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$genus == 
                                                    AggregatedNutrientContent_AllLevels$genus[Z]]
      
      # Obtains the mean N composition value for the genus in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      # Obtains the compounded SD for the N composition values for the genus in row Z.
      comp_SD <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}




# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the genera listed in the SAU taxa table (TaxonKey >= 500000 & <= 599999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 500000, 599999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a P composition value available for the listed genus. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$genus[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$genus) {
      
      # The following code generates a data frame with the set of species-level P values that are available for the genus listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$genus ==
                                                    AggregatedNutrientContent_AllLevels$genus[Z]]
      
      # Obtains the mean P composition value for the genus in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      # Obtains the compounded SD for the P composition values for the genus in row Z.
      comp_SD <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}





##### Family-Level #####

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the families listed in the SAU taxa table (TaxonKey >= 400000 & <= 499999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 400000, 499999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a C composition value available for the listed family Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$family[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$family) {
      
      # The following code generates a data frame with the set of species-level C values that are available for the family listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$family ==
                                                    AggregatedNutrientContent_AllLevels$family[Z]]
      
      # Obtains the mean C composition value for the family in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      # Obtains the compounded SD for the C composition values for the family in row Z.
      comp_SD <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the family listed in the SAU taxa table (TaxonKey >= 400000 & <= 499999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 400000, 499999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a N composition value available for the listed family Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$family[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$family) {
      
      # The following code generates a data frame with the set of species-level N values that are available for the family listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$family == 
                                                    AggregatedNutrientContent_AllLevels$family[Z]]
      
      # Obtains the mean N composition value for the family in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      # Obtains the compounded SD for the N composition values for the family in row Z.
      comp_SD <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the families listed in the SAU taxa table (TaxonKey >= 400000 & <= 499999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 400000, 499999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a P composition value available for the listed family Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$family[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$family) {
      
      # The following code generates a data frame with the set of species-level P values that are available for the family listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$family == 
                                                    AggregatedNutrientContent_AllLevels$family[Z]]
      
      # Obtains the mean P composition value for the family in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      # Obtains the compounded SD for the P composition values for the genus in row Z.
      v <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(v)
      
    }} else (NA)
}




##### Order-Level #####

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the orders listed in the SAU taxa table (TaxonKey >= 300000 & <= 399999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 300000, 399999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a C composition value available for the listed order. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$order[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$order) {
      
      # The following code generates a data frame with the set of species-level C values that are available for the order listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$order == 
                                                    AggregatedNutrientContent_AllLevels$order[Z]]
      
      # Obtains the mean C composition value for the order in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      # Obtains the compounded SD for the C composition values for the order in row Z.
      comp_SD <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the orders listed in the SAU taxa table (TaxonKey >= 300000 & <= 399999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 300000, 399999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a N composition value available for the listed order Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$order[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$order) {
      
      # The following code generates a data frame with the set of species-level N values that are available for the order listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$order == 
                                                    AggregatedNutrientContent_AllLevels$order[Z]]
      
      # Obtains the mean N composition value for the order in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      # Obtains the compounded SD for the N composition values for the order in row Z.
      comp_SD <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the orders listed in the SAU taxa table (TaxonKey >= 300000 & <= 399999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 300000, 399999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a P composition value available for the listed order. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$order[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$order) {
      
      # The following code generates a data frame with the set of species-level P values that are available for the order listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$order == 
                                                    AggregatedNutrientContent_AllLevels$order[Z]]
      
      # Obtains the mean P composition value for the order in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      # Obtains the compounded SD for the P composition values for the order in row Z.
      comp_SD <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}





##### Super-Class-, Class-, and Super-Order- Level #####

###### Class ######

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the class listed in the SAU taxa table (TaxonKey >= 200000 & <= 299999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 200000, 299999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a C composition value available for the listed class. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$class[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$class) {
      
      # The following code generates a data frame with the set of species-level C values that are available for the class listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$class == 
                                                    AggregatedNutrientContent_AllLevels$class[Z]]
      
      # Obtains the mean C composition value for the class in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      # Obtains the compounded SD for the C composition values for the class in row Z.
      comp_SD <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the class listed in the SAU taxa table (TaxonKey >= 200000 & <= 299999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 200000, 299999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a N composition value available for the listed class. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$class[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$class) {
      
      # The following code generates a data frame with the set of species-level N values that are available for the class listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$class == 
                                                    AggregatedNutrientContent_AllLevels$class[Z]]
      
      # Obtains the mean N composition value for the class in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      # Obtains the compounded SD for the N composition values for the class in row Z.
      comp_SD <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the class listed in the SAU taxa table (TaxonKey >= 200000 & <= 299999).
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 200000, 299999)) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a P composition value available for the listed class. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$class[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$class) {
      
      # The following code generates a data frame with the set of species-level P values that are available for the class listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$class == 
                                                    AggregatedNutrientContent_AllLevels$class[Z]]
      
      # Obtains the mean P composition value for the class in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      # Obtains the compounded SD for the P composition values for the genus in row Z.
      comp_SD <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}





###### Super-Order ######

# Some groups were listed as super-orders within the SAU taxa table so we had to specify which rows corresponded to those groups to average the mean nutrient composition values. The super-orders listed below are Pteriomorphia (TaxonKey = 290056) and Batoidea (TaxonKey = 300063).

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the super-orders listed in the SAU taxa table (specified below).
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 290056 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 300063) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a C composition value available for the listed super-order. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$super_order[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$super_order) {
      
      # The following code generates a data frame with the set of species-level C values that are available for the super-order listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$super_order == 
                                                    AggregatedNutrientContent_AllLevels$super_order[Z]]
      
      # Obtains the mean C composition value for the super-order in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      # Obtains the compounded SD for the C composition values for the super-order in row Z.
      comp_SD <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the super-orders listed in the SAU taxa table (specified below).
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 290056 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 300063) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a N composition value available for the listed super-order. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$super_order[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$super_order) {
      
      # The following code generates a data frame with the set of species-level N values that are available for the super-order listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$super_order == 
                                                    AggregatedNutrientContent_AllLevels$super_order[Z]]
      
      # Obtains the mean N composition value for the super-order in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      # Obtains the compounded SD for the N composition values for the super-order in row Z.
      comp_SD <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}



# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  # The first `ifelse` bracket filters out the super-order listed in the SAU taxa table (specified below).
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 290056 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 300063) {
    # The second `ifelse` bracket checks if, for row Z in the dataset, there is at least one species with a P composition value available for the listed super-order. Otherwise, the for loop returns an 'NA'.
    if (AggregatedNutrientContent_AllLevels$super_order[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$super_order) {
      
      # The following code generates a data frame with the set of species-level P values that are available for the super-order listed in row Z.
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$super_order == 
                                                    AggregatedNutrientContent_AllLevels$super_order[Z]]
      
      # Obtains the mean P composition value for the super-order in row Z.
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      # Obtains the compounded SD for the P composition values for the super-order in row Z.
      comp_SD <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(comp_SD)
      
    }} else (NA)
}




###### Super-Class ######

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 200538 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 210006) {
    if (AggregatedNutrientContent_AllLevels$super_class[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$super_class) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$super_class == 
                                                    AggregatedNutrientContent_AllLevels$super_class[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      v <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 200538 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 210006) {
    if (AggregatedNutrientContent_AllLevels$super_class[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$super_class) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$super_class == 
                                                    AggregatedNutrientContent_AllLevels$super_class[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      v <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(v)
      
    }} else (NA)
}





# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 200538 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 210006) {
    if (AggregatedNutrientContent_AllLevels$super_class[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$super_class) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$super_class == 
                                                    AggregatedNutrientContent_AllLevels$super_class[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      v <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



view(AggregatedNutrientContent_AllLevels[between(AggregatedNutrientContent_AllLevels$TaxonKey, 200000, 299999)])
hist(AggregatedNutrientContent_AllLevels[between(AggregatedNutrientContent_AllLevels$TaxonKey, 200000, 299999)]$mean_WW_C)
hist(AggregatedNutrientContent_AllLevels[between(AggregatedNutrientContent_AllLevels$TaxonKey, 200000, 299999)]$compound_WW_C_SD)
hist(AggregatedNutrientContent_AllLevels[between(AggregatedNutrientContent_AllLevels$TaxonKey, 200000, 299999)]$mean_WW_N)
hist(AggregatedNutrientContent_AllLevels[between(AggregatedNutrientContent_AllLevels$TaxonKey, 200000, 299999)]$compound_WW_N_SD)
hist(AggregatedNutrientContent_AllLevels[between(AggregatedNutrientContent_AllLevels$TaxonKey, 200000, 299999)]$mean_WW_P)
hist(AggregatedNutrientContent_AllLevels[between(AggregatedNutrientContent_AllLevels$TaxonKey, 200000, 299999)]$compound_WW_P_SD)





##### Phylum- and Sub-Phylum Level #####

###### Sub-phylum ######

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 100000, 199999)) {
    if (AggregatedNutrientContent_AllLevels$sub_phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$sub_phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$sub_phylum == 
                                                    AggregatedNutrientContent_AllLevels$sub_phylum[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      v <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 100000, 199999)) {
    if (AggregatedNutrientContent_AllLevels$sub_phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$sub_phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$sub_phylum == 
                                                    AggregatedNutrientContent_AllLevels$sub_phylum[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      v <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 100000, 199999)) {
    if (AggregatedNutrientContent_AllLevels$sub_phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$sub_phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$sub_phylum == 
                                                    AggregatedNutrientContent_AllLevels$sub_phylum[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      v <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



###### Phylum ######

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 100000, 199999)) {
    if (AggregatedNutrientContent_AllLevels$phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$phylum == 
                                                    AggregatedNutrientContent_AllLevels$phylum[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      v <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 100000, 199999)) {
    if (AggregatedNutrientContent_AllLevels$phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$phylum == 
                                                    AggregatedNutrientContent_AllLevels$phylum[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      v <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (between(AggregatedNutrientContent_AllLevels$TaxonKey[Z], 100000, 199999)) {
    if (AggregatedNutrientContent_AllLevels$phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$phylum == 
                                                    AggregatedNutrientContent_AllLevels$phylum[Z]]
      
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      v <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



###### Invertebrates ######

# Did not have data on any Echinodermata so the average for all invertebrates was used for that phylum.

# Carbon averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 100077 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 190277) {
    if (AggregatedNutrientContent_AllLevels$sub_phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_C)]$sub_phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$sub_phylum != "Vertebrata"]
      
      AggregatedNutrientContent_AllLevels$mean_WW_C[Z] <- mean(use_df$mean_WW_C, na.rm = T)
      
      v <- sum((use_df$compound_WW_C_SD^2 + (use_df$mean_WW_C - AggregatedNutrientContent_AllLevels$mean_WW_C[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_C_SD[Z] <- sqrt(v)
      
    }} else (NA)
}



# Nitrogen averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 100077 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 190277) {
    if (AggregatedNutrientContent_AllLevels$sub_phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_N)]$sub_phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$sub_phylum != "Vertebrata"]
      
      AggregatedNutrientContent_AllLevels$mean_WW_N[Z] <- mean(use_df$mean_WW_N, na.rm = T)
      
      v <- sum((use_df$compound_WW_N_SD^2 + (use_df$mean_WW_N - AggregatedNutrientContent_AllLevels$mean_WW_N[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_N_SD[Z] <- sqrt(v)
      
    }} else (NA)
}

# Phosphorus averaging
for (Z in 1:nrow(AggregatedNutrientContent_AllLevels)){
  if (AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 100077 | AggregatedNutrientContent_AllLevels$TaxonKey[Z] == 190277) {
    if (AggregatedNutrientContent_AllLevels$sub_phylum[Z] %in% AggregatedNutrientContent_species[!is.na(mean_WW_P)]$sub_phylum) {
      
      use_df <- AggregatedNutrientContent_species[AggregatedNutrientContent_species$sub_phylum != "Vertebrata"]
      
      AggregatedNutrientContent_AllLevels$mean_WW_P[Z] <- mean(use_df$mean_WW_P, na.rm = T)
      
      v <- sum((use_df$compound_WW_P_SD^2 + (use_df$mean_WW_P - AggregatedNutrientContent_AllLevels$mean_WW_P[Z])^2), 
               na.rm = T)/nrow(use_df)
      
      AggregatedNutrientContent_AllLevels$compound_WW_P_SD[Z] <- sqrt(v)
      
    }} else (NA)
}
