###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#                     Nutrient Extraction Molar Ratios                        #
#                                                                             #
###############################################################################

# This script was written after the first draft to calculate the molar ratios of the extractions. Therefore, that is why these code are seperate from the original nutrient extraction estimates script.

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


#### PER YEAR ####





#### PER MARINE REGION ####

# Firstly, the distributions per marine region for each time period were once again created.

# Carbon #

# Generates data frames to get C extraction distributions per time periods.
C_out_per_area <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))
C_out_per_area_1960s <- C_out_per_area
C_out_per_area_1990s <- C_out_per_area
C_out_per_area_2010s <- C_out_per_area


# C - All Years

for (A in 1:nrow(C_out_per_area)) {
  
  C_out_per_area[A, 2:101] <- colSums(C_out[which(C_out_per_area$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                    Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# C - 1960s

for (A in 1:nrow(C_out_per_area_1960s)) {
  
  C_out_per_area_1960s[A, 2:101] <- colSums(C_out[which(C_out_per_area_1960s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# C - 1990s

for (A in 1:nrow(C_out_per_area_1990s)) {
  
  C_out_per_area_1990s[A, 2:101] <- colSums(C_out[which(C_out_per_area_1990s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# C - 2010s

for (A in 1:nrow(C_out_per_area_2010s)) {
  
  C_out_per_area_2010s[A, 2:101] <- colSums(C_out[which(C_out_per_area_2010s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}



# Nitrogen #

# Re-generates data frames to get C extraction distributions per time periods.

N_out_per_area <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))
N_out_per_area_1960s <- N_out_per_area
N_out_per_area_1990s <- N_out_per_area
N_out_per_area_2010s <- N_out_per_area


# N - All Years

for (A in 1:nrow(N_out_per_area)) {
  
  N_out_per_area[A, 2:101] <- colSums(N_out[which(N_out_per_area$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                    Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# N - 1960s

for (A in 1:nrow(N_out_per_area_1960s)) {
  
  N_out_per_area_1960s[A, 2:101] <- colSums(N_out[which(N_out_per_area_1960s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# N - 1990s

for (A in 1:nrow(N_out_per_area_1990s)) {
  
  N_out_per_area_1990s[A, 2:101] <- colSums(N_out[which(N_out_per_area_1990s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# N - 2010s

for (A in 1:nrow(N_out_per_area_2010s)) {
  
  N_out_per_area_2010s[A, 2:101] <- colSums(N_out[which(N_out_per_area_2010s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}



# Phosphorus #

# Re-generates data frames to get N extraction distributions per time periods.# Phosphorus stochastic estimates per area
P_out_per_area <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))
P_out_per_area_1960s <- P_out_per_area
P_out_per_area_1990s <- P_out_per_area
P_out_per_area_2010s <- P_out_per_area


# P - All Years

for (A in 1:nrow(P_out_per_area)) {
  
  P_out_per_area[A, 2:101] <- colSums(P_out[which(P_out_per_area$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                    Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# P - 1960s

for (A in 1:nrow(P_out_per_area_1960s)) {
  
  P_out_per_area_1960s[A, 2:101] <- colSums(P_out[which(P_out_per_area_1960s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# P - 1990s

for (A in 1:nrow(P_out_per_area_1990s)) {
  
  P_out_per_area_1990s[A, 2:101] <- colSums(P_out[which(P_out_per_area_1990s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# P - 2010s

for (A in 1:nrow(P_out_per_area_2010s)) {
  
  P_out_per_area_2010s[A, 2:101] <- colSums(P_out[which(P_out_per_area_2010s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}





##### C:N Estimates #####

# Generates the data frame to hold the nutrient ratio estimates.
NutrientRatios_perArea <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name)) 

# Creates the columns to hold C:N estimates.
NutrientRatios_perArea$C.N_mean <- NA
NutrientRatios_perArea$C.N_SD <- NA
NutrientRatios_perArea$C.N_lowCI <- NA
NutrientRatios_perArea$C.N_highCI <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and N per marine region.
  C_output <- t(C_out_per_area %>% filter(area_name == A) %>% select(-1))
  
  N_output <- t(N_out_per_area %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and N distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, N_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (14/12))
  
  # Generates the mean C:N extraction molar ratios for 1960-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975, na.rm = TRUE))
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.N_mean[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.N_SD[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.N_lowCI[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.N_highCI[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}



# C:N estimates - 1960s

# Creates the columns to hold C:N estimates.
NutrientRatios_perArea$C.N_mean_1960_64 <- NA
NutrientRatios_perArea$C.N_SD_1960_64 <- NA
NutrientRatios_perArea$C.N_lowCI_1960_64 <- NA
NutrientRatios_perArea$C.N_highCI_1960_64 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and N per marine region.
  C_output <- t(C_out_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  N_output <- t(N_out_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and N distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, N_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (14/12))
  
  # Generates the mean C:N extraction molar ratios for 1960-1964. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.N_mean_1960_64[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.N_SD_1960_64[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.N_lowCI_1960_64[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.N_highCI_1960_64[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}



# C:N estimates - 1990s

# Creates the columns to hold C:N estimates.
NutrientRatios_perArea$C.N_mean_1993_97 <- NA
NutrientRatios_perArea$C.N_SD_1993_97 <- NA
NutrientRatios_perArea$C.N_lowCI_1993_97 <- NA
NutrientRatios_perArea$C.N_highCI_1993_97 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and N per marine region.
  C_output <- t(C_out_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  N_output <- t(N_out_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and N distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, N_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (14/12))
  
  # Generates the mean C:N extraction molar ratios for 1993-1997. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.N_mean_1993_97[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.N_SD_1993_97[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.N_lowCI_1993_97[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.N_highCI_1993_97[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}



# C:N estimates - 2010s

# Creates the columns to hold C:N estimates.
NutrientRatios_perArea$C.N_mean_2014_18 <- NA
NutrientRatios_perArea$C.N_SD_2014_18 <- NA
NutrientRatios_perArea$C.N_lowCI_2014_18 <- NA
NutrientRatios_perArea$C.N_highCI_2014_18 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and N per marine region.
  C_output <- t(C_out_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  N_output <- t(N_out_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and N distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, N_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (14/12))
  
  # Generates the mean C:N extraction molar ratios for 2014-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.N_mean_2014_18[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.N_SD_2014_18[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.N_lowCI_2014_18[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.N_highCI_2014_18[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}





##### C:P Estimates #####

# C:P estimates - All Years

NutrientRatios_perArea$C.P_mean <- NA
NutrientRatios_perArea$C.P_SD <- NA
NutrientRatios_perArea$C.P_lowCI <- NA
NutrientRatios_perArea$C.P_highCI <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  C_output <- t(C_out_per_area %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(C_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975, na.rm = TRUE))
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.P_mean[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.P_SD[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.P_lowCI[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.P_highCI[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}

# C:P estimates - 1960s

NutrientRatios_perArea$C.P_mean_1960_64 <- NA
NutrientRatios_perArea$C.P_SD_1960_64 <- NA
NutrientRatios_perArea$C.P_lowCI_1960_64 <- NA
NutrientRatios_perArea$C.P_highCI_1960_64 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  C_output <- t(C_out_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(C_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.P_mean_1960_64[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.P_SD_1960_64[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.P_lowCI_1960_64[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.P_highCI_1960_64[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}


# C:P estimates - 1990s

NutrientRatios_perArea$C.P_mean_1993_97 <- NA
NutrientRatios_perArea$C.P_SD_1993_97 <- NA
NutrientRatios_perArea$C.P_lowCI_1993_97 <- NA
NutrientRatios_perArea$C.P_highCI_1993_97 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  C_output <- t(C_out_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(C_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.P_mean_1993_97[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.P_SD_1993_97[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.P_lowCI_1993_97[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.P_highCI_1993_97[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}


# C:P estimates - 2010s

NutrientRatios_perArea$C.P_mean_2014_18 <- NA
NutrientRatios_perArea$C.P_SD_2014_18 <- NA
NutrientRatios_perArea$C.P_lowCI_2014_18 <- NA
NutrientRatios_perArea$C.P_highCI_2014_18 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  C_output <- t(C_out_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(C_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$C.P_mean_2014_18[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.P_SD_2014_18[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.P_lowCI_2014_18[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.P_highCI_2014_18[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}





##### N:P Estimates #####

# N:P estimates - All Years

NutrientRatios_perArea$N.P_mean <- NA
NutrientRatios_perArea$N.P_SD <- NA
NutrientRatios_perArea$N.P_lowCI <- NA
NutrientRatios_perArea$N.P_highCI <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  N_output <- t(N_out_per_area %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(N_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975, na.rm = TRUE))
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$N.P_mean[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$N.P_SD[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$N.P_lowCI[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$N.P_highCI[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}

# N:P estimates - 1960s

NutrientRatios_perArea$N.P_mean_1960_64 <- NA
NutrientRatios_perArea$N.P_SD_1960_64 <- NA
NutrientRatios_perArea$N.P_lowCI_1960_64 <- NA
NutrientRatios_perArea$N.P_highCI_1960_64 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  N_output <- t(N_out_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(N_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$N.P_mean_1960_64[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$N.P_SD_1960_64[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$N.P_lowCI_1960_64[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$N.P_highCI_1960_64[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}


# N:P estimates - 1990s

NutrientRatios_perArea$N.P_mean_1993_97 <- NA
NutrientRatios_perArea$N.P_SD_1993_97 <- NA
NutrientRatios_perArea$N.P_lowCI_1993_97 <- NA
NutrientRatios_perArea$N.P_highCI_1993_97 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  N_output <- t(N_out_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(N_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$N.P_mean_1993_97[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$N.P_SD_1993_97[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$N.P_lowCI_1993_97[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$N.P_highCI_1993_97[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}


# N:P estimates - 2010s

NutrientRatios_perArea$N.P_mean_2014_18 <- NA
NutrientRatios_perArea$N.P_SD_2014_18 <- NA
NutrientRatios_perArea$N.P_lowCI_2014_18 <- NA
NutrientRatios_perArea$N.P_highCI_2014_18 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  N_output <- t(N_out_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
  P_output <- t(P_out_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
  output <- as.data.frame(cbind(N_output, P_output))
  
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975), na.rm = TRUE)
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  NutrientRatios_perArea$N.P_mean_2014_18[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$N.P_SD_2014_18[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$N.P_lowCI_2014_18[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$N.P_highCI_2014_18[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}