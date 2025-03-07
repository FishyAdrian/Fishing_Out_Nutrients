###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#           Impact of Catch Uncertainties on Nutrient Extractions              #
#                                                                             #
###############################################################################


# This script was used to test the impact of catch uncertainties in our estimates. The Sea Around Us reports uncertainty scores in the catch data which provides approximate confidence intervals. The uncertainty scores were defined as follows: 1 (± 50% of the catch value), 2 (± 30%), 3 (± 20), and 4 (± 10). For catch values lacking an uncertainty score, we assigned the median score value of 3. Using the uncertainty scores, we calculated a standard deviation for each catch value. For example, if a catch value was 50 tonnes and had an uncertainty score of 2 (± 30%), then the resulting standard deviation would be 15 tonnes.

# Once again, we are not able to reproduce the SAU catch data in this repository but we reproduce the code below for reference.

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



#### SETTING UP DATA FRAME FOR UNCERTAINTIES ####

# We first created a copy of our fisheries dataset filtering for the years we considered in the study.
Fisheries_NutrientExtraction_ver2 <- Fisheries_NutrientExtraction %>% filter(year %in% 1960:2018)



# The following code makes the new column replacing all NA scores with the median uncertainty score that was found across the entire SAU database.
Fisheries_NutrientExtraction_ver2$uncertainty_landings <- replace_na(Fisheries_NutrientExtraction_ver2$uncertainty_score, 3)



# Then, the following code converts the scores to an uncertainty percentage according to the percentages assigned to each score on the SAU database (http://www.seaaroundus.org/catch-reconstruction-and-allocation-methods/).
Fisheries_NutrientExtraction_ver2$uncertainty_per <- ifelse(Fisheries_NutrientExtraction_ver2$uncertainty_landings == 1, 0.5, 
                                                            Fisheries_NutrientExtraction_ver2$uncertainty_landings)
Fisheries_NutrientExtraction_ver2$uncertainty_per <- ifelse(Fisheries_NutrientExtraction_ver2$uncertainty_per == 2.0,0.3,
                                                            Fisheries_NutrientExtraction_ver2$uncertainty_per)
Fisheries_NutrientExtraction_ver2$uncertainty_per <- ifelse(Fisheries_NutrientExtraction_ver2$uncertainty_per == 3.0,0.2,
                                                            Fisheries_NutrientExtraction_ver2$uncertainty_per)
Fisheries_NutrientExtraction_ver2$uncertainty_per <- ifelse(Fisheries_NutrientExtraction_ver2$uncertainty_per == 4.0,0.1,
                                                            Fisheries_NutrientExtraction_ver2$uncertainty_per)




#### FUNCTION TO PULL RANDOM FISHERIES VALUES ####

# The following function pulls a random fisheries value from a normal distribution determined by the reported catch data (acting as the mean for the distribution) and the SD associated with its uncertainty score.
norm.pull.fisheries <- function(x,y)
{
  rnorm(1,x,x*y)
}




# The following code runs the function on the fisheries landings for every row in the dataset.
Fisheries_NutrientExtraction_ver2$tonnes.w.error <- 
  mapply(norm.pull.fisheries, 
         Fisheries_NutrientExtraction_ver2$tonnes, Fisheries_NutrientExtraction_ver2$uncertainty_per)




#### NUTRIENT COMPOSITION VALUES ####

##### Making the Function to Pull a Random Nutrient Content Value #####

# Nutrient content values are produced from a beta distribution. This was to assure that we could have mean nutrient compositions and standard deviations that were kept positive and did not produce any negative nutrient compositions. This is how we got our nutrient composition values in the original nutrient extraction estimates.
norm.pull.nutrient <- function(x,y)
{
  rbeta(1, ((x/100)*(((x/100)*(1- (x/100)))/((y/100)^2)-1)),
        ((1-(x/100))*(((x/100)*(1- (x/100)))/((y/100)^2)-1))) * 100
}




#### NUTRIENT EXTRACTION ESTIMATES ####

##### Creating Distributions of Nutrient Extraction Estimates #####

# Creates vectors for the nutrient extraction estimates.
C_extraction_estimates.w.error <- vector()
N_extraction_estimates.w.error <- vector()
P_extraction_estimates.w.error <- vector()


# Create data frames to hold nutrient extraction estimate distributions per year
SAU_years <- c(1960:2018)
C_extraction_estimates_perYear.distribution <- as.data.frame(SAU_years)
N_extraction_estimates_perYear.distribution <- as.data.frame(SAU_years)
P_extraction_estimates_perYear.distribution <- as.data.frame(SAU_years)


# The following for loop creates a distribution of 100 nutrient extraction estimates, similar to how it was done in the original estimates. A distribution of 100 values was once again created for each row in the fisheries dataset and then these distributions were aggregated by year.

for (i in 1:100) {
  # The following code generates the random catch values from the corresponding distribution.
  Fisheries_NutrientExtraction_ver2$tonnes.w.error <- 
    mapply(norm.pull.fisheries, 
           Fisheries_NutrientExtraction_ver2$tonnes, Fisheries_NutrientExtraction_ver2$uncertainty_per)
  
  # The following code generates the random nutrient values from the corresponding distribution for each nutrient.
  ## Carbon
  Fisheries_NutrientExtraction_ver2$WW_C.w.error <- 
    mapply(norm.pull.nutrient, 
           Fisheries_NutrientExtraction_ver2$mean_WW_C, Fisheries_NutrientExtraction_ver2$compound_WW_C_SD)
  
  
  # Nitrogen
  Fisheries_NutrientExtraction_ver2$WW_N.w.error <- 
    mapply(norm.pull.nutrient, 
           Fisheries_NutrientExtraction_ver2$mean_WW_N, Fisheries_NutrientExtraction_ver2$compound_WW_N_SD)
  
  
  # Phosphorus
  Fisheries_NutrientExtraction_ver2$WW_P.w.error <- 
    mapply(norm.pull.nutrient, 
           Fisheries_NutrientExtraction_ver2$mean_WW_P, Fisheries_NutrientExtraction_ver2$compound_WW_P_SD)
  
  
  # Multiplies landings estimates by C, N, and P compositions to obtain extraction estimates per nutrient.
  Fisheries_NutrientExtraction_ver2 <- Fisheries_NutrientExtraction_ver2 %>% 
    mutate(C_extracted = tonnes.w.error * (WW_C.w.error/100),
           N_extracted = tonnes.w.error * (WW_N.w.error/100),
           P_extracted = tonnes.w.error * (WW_P.w.error/100))
  
  # The following code then adds up all the nutrient extraction estimates produced across the fisheries dataset for all 100 iterations. This produces a distribution of 100 total extraction values per nutrient.
  C_extraction_estimates.w.error[i] <- sum(Fisheries_NutrientExtraction_ver2$C_extracted)
  N_extraction_estimates.w.error[i] <- sum(Fisheries_NutrientExtraction_ver2$N_extracted)
  P_extraction_estimates.w.error[i] <- sum(Fisheries_NutrientExtraction_ver2$P_extracted)
  
  
  # Calculating nutrient extraction totals by year
  
  ## Carbon
  agreggated_C <- aggregate(Fisheries_NutrientExtraction_ver2$C_extracted, by = list(Fisheries_NutrientExtraction_ver2$year), FUN = sum)
  # Renames the first column to correspond with the distribution data frame.
  colnames(agreggated_C)[1] <- "SAU_years"
  
  # The C extractions for the iteration were then aggregated per year and added to the per year distribution data frame.
  C_extraction_estimates_perYear.distribution <- merge(C_extraction_estimates_perYear.distribution, agreggated_C, by = "SAU_years", sort = TRUE)
  
  ## Nitrogen
  agreggated_N <- aggregate(Fisheries_NutrientExtraction_ver2$N_extracted, by = list(Fisheries_NutrientExtraction_ver2$year), FUN = sum)
  # Renames the first column to correspond with the distribution data frame.
  colnames(agreggated_N)[1] <- "SAU_years"
  
  # The N extractions for the iteration were then aggregated per year and added to the per year distribution data frame.
  N_extraction_estimates_perYear.distribution <- merge(N_extraction_estimates_perYear.distribution, agreggated_N, by = "SAU_years", sort = TRUE)
  
  ## Phosphorus
  agreggated_P <- aggregate(Fisheries_NutrientExtraction_ver2$P_extracted, by = list(Fisheries_NutrientExtraction_ver2$year), FUN = sum)
  # Renames the first column to correspond with the distribution data frame.
  colnames(agreggated_P)[1] <- "SAU_years"
  
  # The P extractions for the iteration were then aggregated per year and added to the per year distribution data frame.
  P_extraction_estimates_perYear.distribution <- merge(P_extraction_estimates_perYear.distribution, agreggated_P, by = "SAU_years", sort = TRUE)
  
}



#### TOTAL EXTRACTION ESTIMATES ####

# Creates data frame to hold nutrient extraction estimates with fisheries uncertainties.
Total_NutrientExtraction <- data.frame(Nutrient, Extracted.mean, Extracted.SD, Extracted.median, Extracted.95Low, Extracted.95High)


# The code below obtains the mean and SD of total extractions for each nutrient for the period 1960-2018.
# Carbon
Total_NutrientExtraction[1,2] <- mean(C_extraction_estimates.w.error)
Total_NutrientExtraction[1,3] <- sd(C_extraction_estimates.w.error)
Total_NutrientExtraction[1,4] <- median(C_extraction_estimates.w.error)
Total_NutrientExtraction[1,5] <- quantile(C_extraction_estimates.w.error, probs = 0.025)
Total_NutrientExtraction[1,6] <- quantile(C_extraction_estimates.w.error, probs = 0.975)

# Nitrogen
Total_NutrientExtraction[2,2] <- mean(N_extraction_estimates.w.error)
Total_NutrientExtraction[2,3] <- sd(N_extraction_estimates.w.error)
Total_NutrientExtraction[2,4] <- median(N_extraction_estimates.w.error)
Total_NutrientExtraction[2,5] <- quantile(N_extraction_estimates.w.error, probs = 0.025)
Total_NutrientExtraction[2,6] <- quantile(N_extraction_estimates.w.error, probs = 0.975)

# Phosphorus
Total_NutrientExtraction[3,2] <- mean(P_extraction_estimates.w.error)
Total_NutrientExtraction[3,3] <- sd(P_extraction_estimates.w.error)
Total_NutrientExtraction[3,4] <- median(P_extraction_estimates.w.error)
Total_NutrientExtraction[3,5] <- quantile(P_extraction_estimates.w.error, probs = 0.025)
Total_NutrientExtraction[3,6] <- quantile(P_extraction_estimates.w.error, probs = 0.975)




#### EXTRACTIONS PER TIME PERIOD ####

##### Carbon #####
# Creates the data frame to hold the new carbon extractions estimates per time period.
Time.Period <- c("1960-1964", "1993-1997", "2014-2018", "Total")

Carbon_NutrientExtraction <- data.frame(Time.Period, Extracted.mean, Extracted.SD, Extracted.median, Extracted.95Low, Extracted.95High)

colnames(C_extraction_estimates_perYear.distribution)[2:101] <- c(1:100)



## 1960-1964
C_extraction_distribution_1960_64 <- C_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(1960:1964)) %>% 
  summarise(across(2:101, sum))

Carbon_NutrientExtraction[1,2] <- rowMeans(C_extraction_distribution_1960_64)
Carbon_NutrientExtraction[1,3] <- rowSds(as.matrix(C_extraction_distribution_1960_64))
Carbon_NutrientExtraction[1,4] <- rowMedians(as.matrix(C_extraction_distribution_1960_64))
Carbon_NutrientExtraction[1,5] <- rowQuantiles(as.matrix(C_extraction_distribution_1960_64), probs = 0.025)
Carbon_NutrientExtraction[1,6] <- rowQuantiles(as.matrix(C_extraction_distribution_1960_64), probs = 0.975)



## 1993-1997
C_extraction_distribution_1993_97 <- C_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(1993:1997)) %>% 
  summarise(across(2:101, sum))

Carbon_NutrientExtraction[2,2] <- rowMeans(C_extraction_distribution_1993_97)
Carbon_NutrientExtraction[2,3] <- rowSds(as.matrix(C_extraction_distribution_1993_97))
Carbon_NutrientExtraction[2,4] <- rowMedians(as.matrix(C_extraction_distribution_1993_97))
Carbon_NutrientExtraction[2,5] <- rowQuantiles(as.matrix(C_extraction_distribution_1993_97), probs = 0.025)
Carbon_NutrientExtraction[2,6] <- rowQuantiles(as.matrix(C_extraction_distribution_1993_97), probs = 0.975)



## 2014-2018
C_extraction_distribution_2014_18 <- C_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(2014:2018)) %>% 
  summarise(across(2:101, sum))

Carbon_NutrientExtraction[3,2] <- rowMeans(C_extraction_distribution_2014_18)
Carbon_NutrientExtraction[3,3] <- rowSds(as.matrix(C_extraction_distribution_2014_18))
Carbon_NutrientExtraction[3,4] <- rowMedians(as.matrix(C_extraction_distribution_2014_18))
Carbon_NutrientExtraction[3,5] <- rowQuantiles(as.matrix(C_extraction_distribution_2014_18), probs = 0.025)
Carbon_NutrientExtraction[3,6] <- rowQuantiles(as.matrix(C_extraction_distribution_2014_18), probs = 0.975)



## Total
Carbon_NutrientExtraction[4,2:6] <- Total_NutrientExtraction[1,2:6]



##### Nitrogen #####

# Creates the data frame to hold the new nitrogen extractions estimates per time period.
Time.Period <- c("1960-1964", "1993-1997", "2014-2018", "Total")

Nitrogen_NutrientExtraction <- data.frame(Time.Period, Extracted.mean, Extracted.SD, Extracted.median, Extracted.95Low, Extracted.95High)

colnames(N_extraction_estimates_perYear.distribution)[2:101] <- c(1:100)



## 1960-1964
N_extraction_distribution_1960_64 <- N_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(1960:1964)) %>% 
  summarise(across(2:101, sum))

Nitrogen_NutrientExtraction[1,2] <- rowMeans(N_extraction_distribution_1960_64)
Nitrogen_NutrientExtraction[1,3] <- rowSds(as.matrix(N_extraction_distribution_1960_64))
Nitrogen_NutrientExtraction[1,4] <- rowMedians(as.matrix(N_extraction_distribution_1960_64))
Nitrogen_NutrientExtraction[1,5] <- rowQuantiles(as.matrix(N_extraction_distribution_1960_64), probs = 0.025)
Nitrogen_NutrientExtraction[1,6] <- rowQuantiles(as.matrix(N_extraction_distribution_1960_64), probs = 0.975)



## 1993-1997
N_extraction_distribution_1993_97 <- N_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(1993:1997)) %>% 
  summarise(across(2:101, sum))

Nitrogen_NutrientExtraction[2,2] <- rowMeans(N_extraction_distribution_1993_97)
Nitrogen_NutrientExtraction[2,3] <- rowSds(as.matrix(N_extraction_distribution_1993_97))
Nitrogen_NutrientExtraction[2,4] <- rowMedians(as.matrix(N_extraction_distribution_1993_97))
Nitrogen_NutrientExtraction[2,5] <- rowQuantiles(as.matrix(N_extraction_distribution_1993_97), probs = 0.025)
Nitrogen_NutrientExtraction[2,6] <- rowQuantiles(as.matrix(N_extraction_distribution_1993_97), probs = 0.975)




## 2014-2018
N_extraction_distribution_2014_18 <- N_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(2014:2018)) %>% 
  summarise(across(2:101, sum))

Nitrogen_NutrientExtraction[3,2] <- rowMeans(N_extraction_distribution_2014_18)
Nitrogen_NutrientExtraction[3,3] <- rowSds(as.matrix(N_extraction_distribution_2014_18))
Nitrogen_NutrientExtraction[3,4] <- rowMedians(as.matrix(N_extraction_distribution_2014_18))
Nitrogen_NutrientExtraction[3,5] <- rowQuantiles(as.matrix(N_extraction_distribution_2014_18), probs = 0.025)
Nitrogen_NutrientExtraction[3,6] <- rowQuantiles(as.matrix(N_extraction_distribution_2014_18), probs = 0.975)



## Total
Nitrogen_NutrientExtraction[4,2:6] <- Total_NutrientExtraction[1,2:6]



##### Phosphorus #####

# Creates the data frame to hold the new phosphorus extractions estimates per time period.
Time.Period <- c("1960-1964", "1993-1997", "2014-2018", "Total")

Phosphorus_NutrientExtraction <- data.frame(Time.Period, Extracted.mean, Extracted.SD, Extracted.median, Extracted.95Low, Extracted.95High)

colnames(P_extraction_estimates_perYear.distribution)[2:101] <- c(1:100)



## 1960-1964
P_extraction_distribution_1960_64 <- P_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(1960:1964)) %>% 
  summarise(across(2:101, sum))

Phosphorus_NutrientExtraction[1,2] <- rowMeans(P_extraction_distribution_1960_64)
Phosphorus_NutrientExtraction[1,3] <- rowSds(as.matrix(P_extraction_distribution_1960_64))
Phosphorus_NutrientExtraction[1,4] <- rowMedians(as.matrix(P_extraction_distribution_1960_64))
Phosphorus_NutrientExtraction[1,5] <- rowQuantiles(as.matrix(P_extraction_distribution_1960_64), probs = 0.025)
Phosphorus_NutrientExtraction[1,6] <- rowQuantiles(as.matrix(P_extraction_distribution_1960_64), probs = 0.975)



## 1993-1997
P_extraction_distribution_1993_97 <- P_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(1993:1997)) %>% 
  summarise(across(2:101, sum))

Phosphorus_NutrientExtraction[2,2] <- rowMeans(P_extraction_distribution_1993_97)
Phosphorus_NutrientExtraction[2,3] <- rowSds(as.matrix(P_extraction_distribution_1993_97))
Phosphorus_NutrientExtraction[2,4] <- rowMedians(as.matrix(P_extraction_distribution_1993_97))
Phosphorus_NutrientExtraction[2,5] <- rowQuantiles(as.matrix(P_extraction_distribution_1993_97), probs = 0.025)
Phosphorus_NutrientExtraction[2,6] <- rowQuantiles(as.matrix(P_extraction_distribution_1993_97), probs = 0.975)



## 2014-2018
P_extraction_distribution_2014_18 <- P_extraction_estimates_perYear.distribution %>% filter(SAU_years %in% c(2014:2018)) %>% 
  summarise(across(2:101, sum))

Phosphorus_NutrientExtraction[3,2] <- rowMeans(P_extraction_distribution_2014_18)
Phosphorus_NutrientExtraction[3,3] <- rowSds(as.matrix(P_extraction_distribution_2014_18))
Phosphorus_NutrientExtraction[3,4] <- rowMedians(as.matrix(P_extraction_distribution_2014_18))
Phosphorus_NutrientExtraction[3,5] <- rowQuantiles(as.matrix(P_extraction_distribution_2014_18), probs = 0.025)
Phosphorus_NutrientExtraction[3,6] <- rowQuantiles(as.matrix(P_extraction_distribution_2014_18), probs = 0.975)



## Total
Phosphorus_NutrientExtraction[4,2:6] <- Total_NutrientExtraction[1,2:6]



