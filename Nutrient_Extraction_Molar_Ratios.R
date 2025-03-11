###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#                     Nutrient Extraction Molar Ratios                        #
#                                                                             #
###############################################################################

# NOTE: This could should be used only after you have ran all the code from the 'Nutrient_Extraction_Estimates' so that the necessary data frames for the code below have been generated.

#### NECESSARY PACKAGES ####

library(data.table)  # For efficient data loading and writing
library(tidyverse)   # For data manipulation and visualization



#### PER TAXA ####

# For these, we simply calculated the mean molar ratios per taxa using the mean nutrient extractions for each taxa. This was done to have a sense of how C:N:P molar ratios were across the assessed taxa. These data were not used to create any figures.

# Creating data frame by taxa keeping trophic and functional group categories
NutrientExtraction.Ratios_perTaxa <- Fisheries_NutrientExtraction %>% 
  filter(year %in% c(1960:2018)) %>% 
  group_by(scientific_name, trophic_group, simp_functional_group, functional_group) %>% 
  summarize(tonnes = sum(tonnes, na.rm = TRUE),
            C_extracted = sum(C_extracted, na.rm = TRUE),
            N_extracted = sum(N_extracted, na.rm = TRUE),
            P_extracted = sum(P_extracted, na.rm = TRUE),
            mean_C = mean(mean_WW_C),
            mean_N = mean(mean_WW_N),
            mean_P = mean(mean_WW_P)) %>% 
  mutate(C.N = (C_extracted/N_extracted) * (14/12),
         C.P = (C_extracted/P_extracted) * (31/12),
         N.P = (N_extracted/P_extracted) * (31/14))


#### PER YEAR ####

##### C:N Estimates #####

# Generates the data frame to hold the C:N ratio estimates.
C.N_perYear <- data.frame(year = unique(Fisheries_NutrientExtraction$year)) %>%
  filter(year %in% 1960:2018) %>% arrange(year)

# Creates the columns to hold C:N estimates.
C.N_perYear$C.N_mean <- NA
C.N_perYear$C.N_SD <- NA
C.N_perYear$C.N_lowCI <- NA
C.N_perYear$C.N_highCI <- NA



for (A in c(1960:2018)) {
  
  # Transposes the distribution matrices of C and P per year
  C_output <- t(C_extractions_per_year %>% filter(year == A) %>% select(-1))
  N_output <- t(N_extractions_per_year %>% filter(year == A) %>% select(-1))
  
  # Combines the two transposed C and P distribution matrices per year
  output <- as.data.frame(cbind(C_output, N_output))
  
  # Calculates the distribution of molar ratios per year.
  output <- output %>% mutate(ratio = (V1/V2) * (14/12))
  
  # Generates the mean C:N extraction molar ratios per year for 1960-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  mean_output <- mean(output$ratio)
  sd_output <- sd(output$ratio)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975))
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  C.N_perYear$C.N_mean[C.N_perYear$year == A] <- mean_output
  C.N_perYear$C.N_SD[C.N_perYear$year == A] <- sd_output
  C.N_perYear$C.N_lowCI[C.N_perYear$year == A] <- lowCI_output
  C.N_perYear$C.N_highCI[C.N_perYear$year == A] <- highCI_output
  
}





##### C:P Estimates #####

# Generates the data frame to hold the C:P ratio estimates.
C.P_perYear <- data.frame(year = unique(Fisheries_NutrientExtraction$year)) %>%
  filter(year %in% 1960:2018) %>% arrange(year)

# Creates the columns to hold C:P estimates.
C.P_perYear$C.P_mean <- NA
C.P_perYear$C.P_SD <- NA
C.P_perYear$C.P_lowCI <- NA
C.P_perYear$C.P_highCI <- NA



for (A in c(1960:2018)) {
  
  # Transposes the distribution matrices of C and P per year.
  C_output <- t(C_extractions_per_year %>% filter(year == A) %>% select(-1))
  P_output <- t(P_extractions_per_year %>% filter(year == A) %>% select(-1))
  
  # Combines the two transposed C and P distribution matrices per year.
  output <- as.data.frame(cbind(C_output, P_output))
  
  # Calculates the distribution of molar ratios per year.
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  # Generates the mean C:P extraction molar ratios per year for 1960-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  mean_output <- mean(output$ratio)
  sd_output <- sd(output$ratio)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975))
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  C.P_perYear$C.P_mean[C.P_perYear$year == A] <- mean_output
  C.P_perYear$C.P_SD[C.P_perYear$year == A] <- sd_output
  C.P_perYear$C.P_lowCI[C.P_perYear$year == A] <- lowCI_output
  C.P_perYear$C.P_highCI[C.P_perYear$year == A] <- highCI_output
  
}





##### N:P Estimates #####

# Generates the data frame to hold the C:P ratio estimates.
N.P_perYear <- data.frame(year = unique(Fisheries_NutrientExtraction$year)) %>%
  filter(year %in% 1960:2018) %>% arrange(year)

# Creates the columns to hold N:P estimates.
N.P_perYear$N.P_mean <- NA
N.P_perYear$N.P_SD <- NA
N.P_perYear$N.P_lowCI <- NA
N.P_perYear$N.P_highCI <- NA



for (A in c(1960:2018)) {
  
  # Transposes the distribution matrices of N and P per year.
  N_output <- t(N_extractions_per_year %>% filter(year == A) %>% select(-1))
  P_output <- t(P_extractions_per_year %>% filter(year == A) %>% select(-1))
  
  # Combines the two transposed N and P distribution matrices per year.
  output <- as.data.frame(cbind(N_output, P_output))
  
  # Calculates the distribution of molar ratios per year.
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  # Generates the mean N:P extraction molar ratios per year for 1960-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  mean_output <- mean(output$ratio)
  sd_output <- sd(output$ratio)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975))
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  N.P_perYear$N.P_mean[N.P_perYear$year == A] <- mean_output
  N.P_perYear$N.P_SD[N.P_perYear$year == A] <- sd_output
  N.P_perYear$N.P_lowCI[N.P_perYear$year == A] <- lowCI_output
  N.P_perYear$N.P_highCI[N.P_perYear$year == A] <- highCI_output
  
}

# Merge all ratios into one data frame
NutrientRatios_perYear <- cbind(C.N_perYear, C.P_perYear, N.P_perYear)[ , -c(6,11)]




#### PER MARINE REGION ####

# Firstly, the distributions per marine region for each time period were once again created.

# Carbon #

# Generates data frames to get C extraction distributions per time periods.
C_extractions_per_area <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))
C_extractions_per_area_1960s <- C_extractions_per_area
C_extractions_per_area_1990s <- C_extractions_per_area
C_extractions_per_area_2010s <- C_extractions_per_area


# C - All Years

for (A in 1:nrow(C_extractions_per_area)) {
  
  C_extractions_per_area[A, 2:101] <- colSums(C_extractions[which(C_extractions_per_area$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                    Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# C - 1960s

for (A in 1:nrow(C_extractions_per_area_1960s)) {
  
  C_extractions_per_area_1960s[A, 2:101] <- colSums(C_extractions[which(C_extractions_per_area_1960s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# C - 1990s

for (A in 1:nrow(C_extractions_per_area_1990s)) {
  
  C_extractions_per_area_1990s[A, 2:101] <- colSums(C_extractions[which(C_extractions_per_area_1990s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# C - 2010s

for (A in 1:nrow(C_extractions_per_area_2010s)) {
  
  C_extractions_per_area_2010s[A, 2:101] <- colSums(C_extractions[which(C_extractions_per_area_2010s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}



# Nitrogen #

# Re-generates data frames to get C extraction distributions per time periods.

N_extractions_per_area <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))
N_extractions_per_area_1960s <- N_extractions_per_area
N_extractions_per_area_1990s <- N_extractions_per_area
N_extractions_per_area_2010s <- N_extractions_per_area


# N - All Years

for (A in 1:nrow(N_extractions_per_area)) {
  
  N_extractions_per_area[A, 2:101] <- colSums(N_extractions[which(N_extractions_per_area$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                    Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# N - 1960s

for (A in 1:nrow(N_extractions_per_area_1960s)) {
  
  N_extractions_per_area_1960s[A, 2:101] <- colSums(N_extractions[which(N_extractions_per_area_1960s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# N - 1990s

for (A in 1:nrow(N_extractions_per_area_1990s)) {
  
  N_extractions_per_area_1990s[A, 2:101] <- colSums(N_extractions[which(N_extractions_per_area_1990s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# N - 2010s

for (A in 1:nrow(N_extractions_per_area_2010s)) {
  
  N_extractions_per_area_2010s[A, 2:101] <- colSums(N_extractions[which(N_extractions_per_area_2010s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}



# Phosphorus #

# Re-generates data frames to get N extraction distributions per time periods.# Phosphorus stochastic estimates per area
P_extractions_per_area <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))
P_extractions_per_area_1960s <- P_extractions_per_area
P_extractions_per_area_1990s <- P_extractions_per_area
P_extractions_per_area_2010s <- P_extractions_per_area


# P - All Years

for (A in 1:nrow(P_extractions_per_area)) {
  
  P_extractions_per_area[A, 2:101] <- colSums(P_extractions[which(P_extractions_per_area$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                    Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# P - 1960s

for (A in 1:nrow(P_extractions_per_area_1960s)) {
  
  P_extractions_per_area_1960s[A, 2:101] <- colSums(P_extractions[which(P_extractions_per_area_1960s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# P - 1990s

for (A in 1:nrow(P_extractions_per_area_1990s)) {
  
  P_extractions_per_area_1990s[A, 2:101] <- colSums(P_extractions[which(P_extractions_per_area_1990s$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                          Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# P - 2010s

for (A in 1:nrow(P_extractions_per_area_2010s)) {
  
  P_extractions_per_area_2010s[A, 2:101] <- colSums(P_extractions[which(P_extractions_per_area_2010s$area_name[A] == Fisheries_NutrientExtraction$area_name &
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


# C:N Estimates - All Years (1960-2018)

for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and N per marine region.
  C_output <- t(C_extractions_per_area %>% filter(area_name == A) %>% select(-1))
  
  N_output <- t(N_extractions_per_area %>% filter(area_name == A) %>% select(-1))
  
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
  C_output <- t(C_extractions_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  N_output <- t(N_extractions_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
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
  C_output <- t(C_extractions_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  N_output <- t(N_extractions_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
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
  C_output <- t(C_extractions_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  N_output <- t(N_extractions_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
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

# All Years (1960-2018)

# Creates the columns to hold C:P estimates.
NutrientRatios_perArea$C.P_mean <- NA
NutrientRatios_perArea$C.P_SD <- NA
NutrientRatios_perArea$C.P_lowCI <- NA
NutrientRatios_perArea$C.P_highCI <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and P per marine region.
  C_output <- t(C_extractions_per_area %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and P distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  
  mean_output <- mean(output$ratio, na.rm = TRUE)
  sd_output <- sd(output$ratio, na.rm = TRUE)
  CI_output <- quantile(output$ratio, probs = c(0.025, 0.975, na.rm = TRUE))
  lowCI_output <- CI_output[1]
  highCI_output <- CI_output[2]
  
  # Generates the mean C:P extraction molar ratios for 1960-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
  NutrientRatios_perArea$C.P_mean[NutrientRatios_perArea$area_name == A] <- mean_output
  NutrientRatios_perArea$C.P_SD[NutrientRatios_perArea$area_name == A] <- sd_output
  NutrientRatios_perArea$C.P_lowCI[NutrientRatios_perArea$area_name == A] <- lowCI_output
  NutrientRatios_perArea$C.P_highCI[NutrientRatios_perArea$area_name == A] <- highCI_output
  
}



# C:P estimates - 1960s

# Creates the columns to hold C:P estimates.
NutrientRatios_perArea$C.P_mean_1960_64 <- NA
NutrientRatios_perArea$C.P_SD_1960_64 <- NA
NutrientRatios_perArea$C.P_lowCI_1960_64 <- NA
NutrientRatios_perArea$C.P_highCI_1960_64 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  # Transposes the distribution matrices of C and P per marine region.
  C_output <- t(C_extractions_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and P distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  # Generates the mean C:P extraction molar ratios for 1960-1964. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
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

# Creates the columns to hold C:P estimates.
NutrientRatios_perArea$C.P_mean_1993_97 <- NA
NutrientRatios_perArea$C.P_SD_1993_97 <- NA
NutrientRatios_perArea$C.P_lowCI_1993_97 <- NA
NutrientRatios_perArea$C.P_highCI_1993_97 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and P per marine region.
  C_output <- t(C_extractions_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and P distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  # Generates the mean C:P extraction molar ratios for 1993-1997. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
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

# Creates the columns to hold C:P estimates.
NutrientRatios_perArea$C.P_mean_2014_18 <- NA
NutrientRatios_perArea$C.P_SD_2014_18 <- NA
NutrientRatios_perArea$C.P_lowCI_2014_18 <- NA
NutrientRatios_perArea$C.P_highCI_2014_18 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of C and P per marine region.
  C_output <- t(C_extractions_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed C and P distribution matrices per marine region.
  output <- as.data.frame(cbind(C_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/12))
  
  # Generates the mean C:P extraction molar ratios for 2014-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
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

# N:P estimates - All Years (1960-2018)

# Creates the columns to hold N:P estimates.
NutrientRatios_perArea$N.P_mean <- NA
NutrientRatios_perArea$N.P_SD <- NA
NutrientRatios_perArea$N.P_lowCI <- NA
NutrientRatios_perArea$N.P_highCI <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of N and P per marine region.
  N_output <- t(N_extractions_per_area %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed N and P distribution matrices per marine region.
  output <- as.data.frame(cbind(N_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  # Generates the mean N:P extraction molar ratios for 1960-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
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

# Creates the columns to hold N:P estimates.
NutrientRatios_perArea$N.P_mean_1960_64 <- NA
NutrientRatios_perArea$N.P_SD_1960_64 <- NA
NutrientRatios_perArea$N.P_lowCI_1960_64 <- NA
NutrientRatios_perArea$N.P_highCI_1960_64 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of N and P per marine region.
  N_output <- t(N_extractions_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area_1960s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed N and P distribution matrices per marine region.
  output <- as.data.frame(cbind(N_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  # Generates the mean N:P extraction molar ratios for 1960-1964. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
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

# Creates the columns to hold N:P estimates.
NutrientRatios_perArea$N.P_mean_1993_97 <- NA
NutrientRatios_perArea$N.P_SD_1993_97 <- NA
NutrientRatios_perArea$N.P_lowCI_1993_97 <- NA
NutrientRatios_perArea$N.P_highCI_1993_97 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of N and P per marine region.
  N_output <- t(N_extractions_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area_1990s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed N and P distribution matrices per marine region.
  output <- as.data.frame(cbind(N_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  # Generates the mean N:P extraction molar ratios for 1993-1997. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
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

# Creates the columns to hold N:P estimates.
NutrientRatios_perArea$N.P_mean_2014_18 <- NA
NutrientRatios_perArea$N.P_SD_2014_18 <- NA
NutrientRatios_perArea$N.P_lowCI_2014_18 <- NA
NutrientRatios_perArea$N.P_highCI_2014_18 <- NA


for (A in unique(NutrientRatios_perArea$area_name)) {
  
  # Transposes the distribution matrices of N and P per marine region.
  N_output <- t(N_extractions_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  P_output <- t(P_extractions_per_area_2010s %>% filter(area_name == A) %>% select(-1))
  
  # Combines the two transposed N and P distribution matrices per marine region.
  output <- as.data.frame(cbind(N_output, P_output))
  
  # Calculates the distribution of molar ratios per marine region.
  output <- output %>% mutate(ratio = (V1/V2) * (31/14))
  
  # Generates the mean N:P extraction molar ratios for 2014-2018. Also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
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