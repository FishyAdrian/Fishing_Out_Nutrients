###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#   Nutrient Compositions for Higher-Level Groups and Data-Sparse Species     #
#                                                                             #
###############################################################################


# After we generated predicted nutrient compositions using the predicted ratios and simulations described in the Predicive Models script, we averaged our observed and predicted nutrient data at the species level to generate C, N, and P compositions. To reduce the influence of over-represented species in our calculations, we first calculated the mean and standard deviations of nutrient composition values for each species represented in our dataset. Then, for each taxonomic level above species, we calculated the mean and standard deviation of nutrient composition values based on each taxonomic level's corresponding species. To account for compounding uncertainty among our observed and predicted values, we calculated the compounded standard deviation as described by the compound_SD equation below. NOTE: For observed values, we did not record SDs since these were not always reported. However, to account for uncertainty in observed nutrient values, we assigned SDs based on the in-species coefficients of variation (CV) for each nutrient observed in Czamanski et al. (2011) (https://link.springer.com/article/10.1007/s00227-011-1783-7). The coefficients of variation for each nutrient were as follows: CV of C = 0.0175, CV of N = 0.0097, and CV of P = 0.0528. This assured that each observed and predicted value had an associated SD to be used in calculating the compound standard deviation.



#### NECESSARY PACKAGES ####

library(tidyverse)
library(data.table)
library(ggplot2)
library(GGally)
library(nlme)
library (MuMIn)
library(ggpubr)
library(mvtnorm)




#### NUTRIENT COMPOSITIONS FOR HIGHER-LEVEL GROUPS AND DATA-SPARSE SPECIES ####

# After we generated predicted nutrient compositions using the predicted ratios and simulations described above, we averaged our observed and predicted nutrient data at the species level to generate mean C, N, and P compositions per species. To reduce the influence of over-represented species in our calculations, we first calculated the mean and standard deviations of nutrient composition values for each species represented in our initial dataset. Then, for each taxonomic level above species, we calculated the mean and standard deviation of nutrient composition values based on each taxonomic level's corresponding species. To account for compounding uncertainty among our observed and predicted values, we calculated a compounded standard deviation.

# In the following function, data is used to designate the data frame that has all of your observed and predicted values. Level designates the taxonomic level you are attempting to average. Nutrient and "nutrient_SD" designates the columns that hold your nutrient values and their corresponding SDs, respectively.

#### ASSIGNING NUTRIENT COMPOSITION VALUES FOR SPECIES AND GROUPS WITHOUT DATA ####

# Explanation here...

# Function to fill missing nutrient content based on taxonomic hierarchy
assign_nutrient_value <- function(data, taxonomic_levels, nutrient_col = NULL, sd_col = NULL) {
  
  for (i in 1:nrow(data)) {
    if (is.na(data[[nutrient_col]][i])) {
      for (level in taxonomic_levels) {
        taxon <- data[[level]][i]
        # Check if there is a non-NA nutrient value for the taxon in the same data frame
        matching_row <- data[!is.na(data[[nutrient_col]]) & data[[level]] == taxon, ]
        if (nrow(matching_row) > 0) {
          data[[nutrient_col]][i] <- matching_row[[nutrient_col]][1] # Assign the first matching value
          data[[sd_col]][i] <- matching_row[[sd_col]][1] # Assign the first matching SD value
          break
        }
      }
    }
  }
  
  return(data)
}


# You can assign the taxonomic levels you desire to draw averages from by assigning the levels to a vector like so:
taxonomic_levels <- c("genus", "family", "order", "super_order", "class", "super_class", "sub_phylum", "phylum")


