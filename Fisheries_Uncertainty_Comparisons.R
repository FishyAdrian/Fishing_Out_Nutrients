###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#           Impact of Catch Uncertainties on Nutrient Extractions              #
#                                                                             #
###############################################################################


# This script was used to test the impact of catch uncertainties in our estimates. The Sea Around Us reports uncertainty scores in the catch data which provides approximate confidence intervals. The uncertainty scores were defined as follows: 1 (± 50% of the catch value), 2 (± 30%), 3 (± 20), and 4 (± 10). For catch values lacking an uncertainty score, we assigned the median score value of 3. Using the uncertainty scores, we calculated a standard deviation for each catch value. For example, if a catch value was 50 tonnes and had an uncertainty score of 2 (± 30%), then the resulting standard deviation would be 15 tonnes.

# Once again, we are not able to provide the SAU catch data in this repository for this analysis but we reproduce the code below for reference. We have provided the generated distributions as an RData file titled "Extraction_Distributions_w_Fisheries_Uncerainties.RData" for the user to have as reference. The mean extraction estimates with fisheries uncertainties accounted for are reported in Supplementary Tables 8-10 in the manuscript.

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
library(matrixStats)



#### SETTING UP DATA FRAME FOR UNCERTAINTIES ####

# We first created a copy of our fisheries dataset filtering for the years we considered in the study (1960 - 2018).
Fisheries_NutrientExtraction_ver2 <- Fisheries_NutrientExtraction %>% filter(year %in% 1960:2018)



# The following code makes a new column to ensure all rows had an assigned uncertainty score. Any rows without uncertainty scores were applied a score the median uncertainty score of 3.
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

# The following function pulls a random fisheries value from a normal distribution determined by the tonnes landed (acting as the mean for the distribution) and the SD associated with its uncertainty score.
norm.pull.fisheries <- function(x,y)
{
  # The code below takes the landed amount (x) and multiplies it by its corresponding uncertainty percentage (y) to generate the appropriate standard deviation for the normal distribution.
  rnorm(1,x,x*y)
}




# The following code runs the norm.pull.fisheries function on the for every row in the fisheries dataset.
Fisheries_NutrientExtraction_ver2$tonnes.w.error <- 
  mapply(norm.pull.fisheries, 
         Fisheries_NutrientExtraction_ver2$tonnes, Fisheries_NutrientExtraction_ver2$uncertainty_per)




#### NUTRIENT COMPOSITION VALUES ####

##### Making the Function to Pull a Random Nutrient Content Value #####

# We recreated our methodology for creating distributions of nutrient composition values as a function. Nutrient composition values are produced from a beta distribution. This was to assure that we could have mean nutrient compositions and standard deviations that were kept positive and did not produce any negative nutrient compositions. The distribution uses the estimated mean nutrient composition as the mean and the estimated SD as the standard deviation for the distribution.
norm.pull.nutrient <- function(x,y)
{
  rbeta(1, ((x/100)*(((x/100)*(1- (x/100)))/((y/100)^2)-1)),
        ((1-(x/100))*(((x/100)*(1- (x/100)))/((y/100)^2)-1))) * 100
}




#### NUTRIENT EXTRACTION ESTIMATES ####

# For each nutrient, we created a vector to hold the extraction estimates.
C_extraction_estimates.w.error <- vector()
N_extraction_estimates.w.error <- vector()
P_extraction_estimates.w.error <- vector()


# Then, we created data frames to hold nutrient extraction estimate distributions for each year.
SAU_years <- c(1960:2018)
C_extraction_estimates_perYear.distribution <- as.data.frame(SAU_years)
N_extraction_estimates_perYear.distribution <- as.data.frame(SAU_years)
P_extraction_estimates_perYear.distribution <- as.data.frame(SAU_years)


# The following for loop creates a distribution of 100 nutrient extraction estimates for each row of the fisheries dataset. Each iteration of the for loop creates one value for each row of the dataset. Within each iteration, the total extractions per nutrient for each year were aggregated and summed up. This ran 100 times to generate a distribution of 100 extraction estimate values  

for (i in 1:100) {
  # The following code generates a random catch value for each row of the dataset from the distribution as described in the function above.
  Fisheries_NutrientExtraction_ver2$tonnes.w.error <- 
    mapply(norm.pull.fisheries, 
           Fisheries_NutrientExtraction_ver2$tonnes, Fisheries_NutrientExtraction_ver2$uncertainty_per)
  
  # The following code generates a random composition value for each nutrient and for each row in the dataset from the corresponding distributions as described in the function above.
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
  
  
  
  # Once random catch and nutrient composition values have been drawn, the following code multiplied the catch values by the C, N, and P compositions to obtain extraction estimates which also accounted for the fisheries uncertainties.
  Fisheries_NutrientExtraction_ver2 <- Fisheries_NutrientExtraction_ver2 %>% 
    mutate(C_extracted = tonnes.w.error * (WW_C.w.error/100),
           N_extracted = tonnes.w.error * (WW_N.w.error/100),
           P_extracted = tonnes.w.error * (WW_P.w.error/100))
  
  
  
  # The following code then takes all of the values calculated across the rows of the dataset to produce a total extraction amount for 1960-2018 for each nutrient. This process was repeated 100 times to create a distribution of 100 values per nutrient. The distributions used to assess the impact of fisheries uncertainties are provided in the RData file titled "Extraction_Distributions_w_Fisheries_Uncerainties.RData".
  C_extraction_estimates.w.error[i] <- sum(Fisheries_NutrientExtraction_ver2$C_extracted)
  N_extraction_estimates.w.error[i] <- sum(Fisheries_NutrientExtraction_ver2$N_extracted)
  P_extraction_estimates.w.error[i] <- sum(Fisheries_NutrientExtraction_ver2$P_extracted)
  
  
  
  
  # Calculating nutrient extraction totals by year
  
  # The following pieces of code repeated the process of adding up all extraction values across rows but did so for each year. Once again, this process was repeated 100 times for each nutrient such that each year had a distribution of 100 nutrient extraction values.
  
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

# Creates the columns that will be used in the data frames that will be made below.
Time.Period <- c("1960-1964", "1993-1997", "2014-2018", "Total")
Nutrient <- c("Carbon", "Nitrogen", "Phosphorus")
Extracted.mean <- NA
Extracted.median <- NA
Extracted.SD <- NA
Extracted.95Low <- NA
Extracted.95High <- NA

# Creates a data frame to hold total nutrient extraction estimates with fisheries uncertainties.
NutrientExtractions_w.fisheries.uncertainties <- data.frame(Nutrient, Extracted.mean, Extracted.SD, Extracted.median, Extracted.95Low, Extracted.95High)


# The code below obtains the mean and spread of total extractions for each nutrient for the period 1960-2018.
# Carbon
NutrientExtractions_w.fisheries.uncertainties[1,2] <- mean(C_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[1,3] <- sd(C_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[1,4] <- median(C_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[1,5] <- quantile(C_extraction_estimates.w.error, probs = 0.025)
NutrientExtractions_w.fisheries.uncertainties[1,6] <- quantile(C_extraction_estimates.w.error, probs = 0.975)

# Nitrogen
NutrientExtractions_w.fisheries.uncertainties[2,2] <- mean(N_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[2,3] <- sd(N_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[2,4] <- median(N_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[2,5] <- quantile(N_extraction_estimates.w.error, probs = 0.025)
NutrientExtractions_w.fisheries.uncertainties[2,6] <- quantile(N_extraction_estimates.w.error, probs = 0.975)

# Phosphorus
NutrientExtractions_w.fisheries.uncertainties[3,2] <- mean(P_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[3,3] <- sd(P_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[3,4] <- median(P_extraction_estimates.w.error)
NutrientExtractions_w.fisheries.uncertainties[3,5] <- quantile(P_extraction_estimates.w.error, probs = 0.025)
NutrientExtractions_w.fisheries.uncertainties[3,6] <- quantile(P_extraction_estimates.w.error, probs = 0.975)




#### EXTRACTIONS PER TIME PERIOD ####

##### Carbon #####

# Creates the columns that will be used in the data frames that will be made below.
Time.Period <- c("1960-1964", "1993-1997", "2014-2018", "Total")
Nutrient <- c("Carbon", "Nitrogen", "Phosphorus")
Extracted.mean <- NA
Extracted.median <- NA
Extracted.SD <- NA
Extracted.95Low <- NA
Extracted.95High <- NA

# Creates the data frame to hold total carbon extraction estimates with fisheries uncertainties for each time period.
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
Carbon_NutrientExtraction[4,2:6] <- NutrientExtractions_w.fisheries.uncertainties[1,2:6]



##### Nitrogen #####

# Creates the data frame to hold total nitrogen extraction estimates with fisheries uncertainties for each time period.
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
Nitrogen_NutrientExtraction[4,2:6] <- NutrientExtractions_w.fisheries.uncertainties[2,2:6]



##### Phosphorus #####

# Creates the data frame to hold total phosphorus extraction estimates with fisheries uncertainties for each time period.
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
Phosphorus_NutrientExtraction[4,2:6] <- NutrientExtractions_w.fisheries.uncertainties[3,2:6]



