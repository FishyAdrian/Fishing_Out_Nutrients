###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#                              GitHub Test Code                               #
#                                                                             #
###############################################################################


#### NECESSARY PACKAGES ####
library(tidyverse)
library(data.table)
library(ggplot2)
library(GGally)
library(nlme)
library (MuMIn)
library(ggpubr)
library(mvtnorm)



#### Creating copy datasets for testing ####
test_PredRatios <- pred_Ratios %>% select(c(1:23))


#### Testing predictive model function - CHECKED ####

wet.NP_test_model <- ratio_model(Nutrient_Ratios, ratio = "N_P", nutrient = "WetWeight_Pwhole")

dry.NP_test_model <- ratio_model(Nutrient_Ratios, ratio = "N_P", nutrient = "DryWeight_Pwhole")

wet.NP_test_model <- ratio_model(Nutrient_Ratios, ratio = "N_P", nutrient = "WetWeight_N")


#### Testing nutrient ratio predictions - CHECKED ####

# Mean predictions remain the same, but the SDs in the function I've created for use here are drastically higher than the ones output by the 
pred_ratios_df <- predict_ratios(test_PredRatios, dry.NP_test_model, nutrient = "DryWeight_Pwhole")

pred_ratios_df <- predict_ratios(test_PredRatios, wet.NP_test_model, nutrient = "WetWeight_Pwhole")

