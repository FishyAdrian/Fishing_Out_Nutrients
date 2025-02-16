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



#### PREDICTIVE MODELS FOR NUTRIENT COMPOSITION - TESTED, All GOOD ####

test_PredRatios <- pred_Ratios %>% select(c(1:23))


# N:P from wet P
NP_from_wetP_test <- ratio_model(Nutrient_Ratios, ratio = "N_P", nutrient = "WetWeight_Pwhole")

# N:P from dry P
NP_from_dryP_test <- ratio_model(Nutrient_Ratios, ratio = "N_P", nutrient = "DryWeight_Pwhole")

# N:P from wet N
NP_from_wetN_test <- ratio_model(Nutrient_Ratios, ratio = "N_P", nutrient = "WetWeight_N")


#### NUTRIENT RATIO PREDICTIONS - TESTED, ALL GOOD ####

# N:P from Dry P
pred_ratios_df <- predict_ratios(test_PredRatios, NP_from_dryP_test, nutrient = "DryWeight_Pwhole")

hist(pred_ratios_df$pred_ratio - pred_Ratios$pred_log_NP_fromDryP)
hist(pred_ratios_df$pred_ratio_SD - pred_Ratios$predNP_fromDryP_SD)

# N:P from Wet P

pred_ratios_df <- predict_ratios(test_PredRatios, NP_from_wetP_test, nutrient = "WetWeight_Pwhole")

hist(pred_ratios_df$pred_ratio - pred_Ratios$pred_log_NP_fromWetP)
hist(pred_ratios_df$pred_ratio_SD - pred_Ratios$predNP_fromWetP_SD)





#### NUTRIENT COMPOSITION PREDICTIONS FROM PREDICTED RATIOS ####

# Nitrogen from N:P and Dry P

test_pred_N <- predict_N(pred_ratios_df, nutrient = "DryWeight_Pwhole", ratio_type = "N:P", pred_ratio = "pred_ratio", pred_ratio_SD = "pred_ratio_SD")

hist(test_pred_N$pred_nutrient - pred_Ratios_NutrientContent_fromStaticModel$derived_DryN)

# Nitrogen from N:P and Wet P

test_pred_N <- predict_N(pred_ratios_df, nutrient = "WetWeight_Pwhole", ratio_type = "N:P", pred_ratio = "pred_ratio", pred_ratio_SD = "pred_ratio_SD")

hist(test_pred_N$pred_nutrient - pred_Ratios_NutrientContent_fromStaticModel$derived_WetN)


