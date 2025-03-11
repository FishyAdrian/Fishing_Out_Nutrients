###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#                       Nutrient Extraction Estimates                         #
#                                                                             #
###############################################################################


#### NECESSARY PACKAGES ####

library(data.table)  # For efficient data loading and writing
library(tidyverse)   # For data manipulation and visualization





#### LOADING NECESSARY DATA FRAMES ####

# Loads the Industrial Taxa Nutrient Content data frame. This data frame contains all of the estimated C, N, and P nutrient compositions and their associated SDs.
IndustrialTaxa_NutrientContent <- fread("https://raw.githubusercontent.com/FishyAdrian/Fishing_Out_Nutrients/refs/heads/main/IndustrialTaxa_NutrientContent.csv")


# To estimate nutrient extractions, we used data from the Sea Around Us (SAU) project (available at https://www.seaaroundus.org/data/) at both the Exclusive Economic Zone (EEZ) and High Seas levels. Specifically, we downloaded CSV files for each EEZ (n=283) and High Seas region (n=18), filtered for catches from industrial fisheries, and compiled the data into a dataset containing over 6.7 million records of landed catches, reported in tonnes. These data were stored in the data frame Fisheries_NutrientExtraction. We are not able to reproduce the SAU catch data in this repository. However, we include the Fisheries_NutrientExtraction_index CSV file which includes the identifying information for each row's taxa to allow our code to be reproduced with the extraction estimate distributions we produced in the analysis below. 


# Loads the Fisheries_NutrientExtraction data frame which holds the identifying information for each row of the Sea Around Us dataset.
Fisheries_NutrientExtraction <- fread("Fisheries_NutrientExtraction_index.csv")





#### GENERATING NUTRIENT COMPOSITION DISTRIBUTIONS ####

# For each row of the SAU dataset, named "Fisheries_NutrientExtraction", which comprised 6.7+ million rows of catch data, we generated 100 nutrient values drawn from a random normal distribution of each taxa's mean nutrient compositions. Then, each nutrient composition value generated in the distribution was multiplied by the landed amount reported to produce 100 nutrient extraciton estimates.


# For the first step, we wanted to make sure we could generate distributions of our nutrient compositions that were between 0 and 1 (i.e., no negative values). To accomplish this, we calculated alpha and beta values for our nutrient composition values in order to generate beta distributions below.


# Modifies the IndustrialTaxa_NutrientContent data frame to produce Beta values of our nutrient compositions.
IndustrialTaxa_NutrientContent <- IndustrialTaxa_NutrientContent %>% 
  mutate(alph_C = (mean_WW_C/100)*(((mean_WW_C/100)*(1- (mean_WW_C/100)))/((compound_WW_C_SD/100)^2)-1),
         bet_C = (1-(mean_WW_C/100))*(((mean_WW_C/100)*(1- (mean_WW_C/100)))/((compound_WW_C_SD/100)^2)-1),
         alph_N = (mean_WW_N/100)*(((mean_WW_N/100)*(1- (mean_WW_N/100)))/((compound_WW_N_SD/100)^2)-1),
         bet_N = (1-(mean_WW_N/100))*(((mean_WW_N/100)*(1- (mean_WW_N/100)))/((compound_WW_N_SD/100)^2)-1),
         alph_P = (mean_WW_P/100)*(((mean_WW_P/100)*(1- (mean_WW_P/100)))/((compound_WW_P_SD/100)^2)-1),
         bet_P = (1-(mean_WW_P/100))*(((mean_WW_P/100)*(1- (mean_WW_P/100)))/((compound_WW_P_SD/100)^2)-1))



# The following code created a vector of integers to identify which rows of the nutrient composition data frame had to be called when calculating nutrient extraction with the fisheries data frame.
taxaorder <- match(Fisheries_NutrientExtraction$TaxonKey, IndustrialTaxa_NutrientContent$TaxonKey)


# Using the alpha and beta values, we generated 100 nutrient composition values. Here, we generated a matrix where each row of the fisheries dataset had a distribution of 100 nutrient values based on the corresponding taxa for each row.

C_distributions <- matrix(
  rbeta(length(taxaorder) * 100, IndustrialTaxa_NutrientContent$alph_C[taxaorder], 
        IndustrialTaxa_NutrientContent$bet_C[taxaorder]), 
  nrow = length(taxaorder), ncol = 100
)

N_distributions <- matrix(
  rbeta(length(taxaorder) * 100, IndustrialTaxa_NutrientContent$alph_N[taxaorder], 
        IndustrialTaxa_NutrientContent$bet_N[taxaorder]), 
  nrow = length(taxaorder), ncol = 100
)

P_distributions <- matrix(
  rbeta(length(taxaorder) * 100, IndustrialTaxa_NutrientContent$alph_P[taxaorder], 
        IndustrialTaxa_NutrientContent$bet_P[taxaorder]), 
  nrow = length(taxaorder), ncol = 100
)


# Then, for each nutrient, we multiplied the matrix of generated composition values by the landed catch amounts (in tonnes) to produce a distribution of 100 extraction estimates for each row of the fisheries dataset.

C_extractions <- SAU_catch.data$tonnes * C_distributions
N_extractions <- SAU_catch.data$tonnes * N_distributions
P_extractions <- SAU_catch.data$tonnes * P_distributions





#### LOADING NUTRIENT EXTRACTION DISTRIBUTIONS ####

# Each row in the Fisheries_NutrientExtraction data frame corresponds to the rows in the nutrient extraction distribution matrices below. Therefore, the code provided can be used to recalculate the nutrient extraction estimates presented in the manuscript. 

# You can download the extraction distribution matrices from the following Fighsare repository: https://doi.org/10.6084/m9.figshare.28500593. It is recommended to load and work with one distribution matrix at a time given their large file size.

C_extractions <- fread("C_extraction_distributionmatrix.csv")

N_extractions <- fread("N_extraction_distributionmatrix.csv")

P_extractions <- fread("P_extraction_distributionmatrix.csv")





#### ESTIMATING MEAN NUTRIENT EXTRACTIONS ####

# The following code generates the mean nutrient extraction for each nutrient for each row in the Fisheries_NutrientExtraction data frame. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.

# Carbon
Fisheries_NutrientExtraction$C_extracted <- apply(C_extractions, 1, FUN = mean)
Fisheries_NutrientExtraction$C_extracted_SD <- apply(C_extractions, 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions, 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
Fisheries_NutrientExtraction$C_extracted_lowCI <- quantiles_CIs_C[1, ]
Fisheries_NutrientExtraction$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
Fisheries_NutrientExtraction$C_extracted_median <- quantiles_CIs_C[3, ]
Fisheries_NutrientExtraction$C_extracted_highIQR <- quantiles_CIs_C[4, ]
Fisheries_NutrientExtraction$C_extracted_highCI <- quantiles_CIs_C[5, ]


# Nitrogen
Fisheries_NutrientExtraction$N_extracted <- apply(N_extractions, 1, FUN = mean)
Fisheries_NutrientExtraction$N_extracted_SD <- apply(N_extractions, 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions, 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
Fisheries_NutrientExtraction$N_extracted_lowCI <- quantiles_CIs_N[1, ]
Fisheries_NutrientExtraction$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
Fisheries_NutrientExtraction$N_extracted_median <- quantiles_CIs_N[3, ]
Fisheries_NutrientExtraction$N_extracted_highIQR <- quantiles_CIs_N[4, ]
Fisheries_NutrientExtraction$N_extracted_highCI <- quantiles_CIs_N[5, ]


# Phosphorus
Fisheries_NutrientExtraction$P_extracted <- apply(P_extractions, 1, FUN = mean)
Fisheries_NutrientExtraction$P_extracted_SD <- apply(P_extractions, 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions, 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
Fisheries_NutrientExtraction$P_extracted_lowCI <- quantiles_CIs_P[1, ]
Fisheries_NutrientExtraction$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
Fisheries_NutrientExtraction$P_extracted_median <- quantiles_CIs_P[3, ]
Fisheries_NutrientExtraction$P_extracted_highIQR <- quantiles_CIs_P[4, ]
Fisheries_NutrientExtraction$P_extracted_highCI <- quantiles_CIs_P[5, ]





#### ESTIMATING NUTRIENT EXTRACTIONS ####

# For each group of extraction estimates (e.g., by year), we added up the distributions of 100 extraction estimates for all rows that matched the grouping category. For instance, to estimate annual nutrient extractions, we added up the distributions generated for all rows in the SAU database corresponding to a specific year. This process resulted in a single distribution of 100 estimates for each year. From this distribution, we calculated the mean to determine the estimated nutrient extraction. Additionally, the distribution was used to derive standard deviations (SDs) and 95% confidence intervals (CIs). This process was done for each year, time period, spatial region, trophic group, and functional group. This can be reproduced by using extraction distribution matrices provided in the Figshare repository (https://doi.org/10.6084/m9.figshare.28500593).


##### TOTAL NUTRIENT EXTRACTION #####

# Total
C_extractions_total <- colSums(C_extractions[which(Fisheries_NutrientExtraction$year %in% c(1960:2018))])
mean(C_extractions_total)

N_extractions_total <- colSums(N_extractions[which(Fisheries_NutrientExtraction$year %in% c(1960:2018))])
mean(N_extractions_total)

P_extractions_total <- colSums(P_extractions[which(Fisheries_NutrientExtraction$year %in% c(1960:2018))])
mean(P_extractions_total)

# The following estimates should result from the operations above.
# Carbon = 431,158,275 tonnes
# Nitrogen = 110,292,936 tonnes
# Phosphorus = 22,817,918 tonnes


# SDs
sd(C_extractions_total) # +/- 1,126,340
sd(N_extractions_total) # +/- 189,477.2
sd(P_extractions_total) # +/- 153,114.3


# 95% CIs
quantile(C_extractions_total, probs = c(0.025, 0.975)) # 428,865,830 - 433,477,727
quantile(N_extractions_total, probs = c(0.025, 0.975)) # 109,957,232 - 110,711,453
quantile(P_extractions_total, probs = c(0.025, 0.975)) # 22,513,904 - 23,110,214


# IQR
quantile(C_extractions_total, probs = c(0.25, 0.5, 0.75)) # 430,565,209; 431,271,544; 431,838,797
quantile(N_extractions_total, probs = c(0.25, 0.5, 0.75)) # 110,172,814; 110,286,014; 110,427,420 
quantile(P_extractions_total, probs = c(0.25, 0.5, 0.75)) # 22,715,600; 22,820,242; 22,913,843





##### PER YEAR #####

# Generates a data frame to summarize nutrient extraction estimates per year.
NutrientExtraction_perYear <- data.frame(year = c(1960:2018))



###### Carbon ######
# Creates data frame for holding C distributions per year.
C_extractions_per_year <- data.frame(year = c(1960:2018))

# For loop to add up C extraction distributions per year.
for (A in 1:nrow(C_extractions_per_year)) {
  
  C_extractions_per_year[A, 2:101] <- colSums(C_extractions[which(C_extractions_per_year$year[A] == Fisheries_NutrientExtraction$year)])
  
}


# The following code generates the mean C extraction for each year between 1960 and 2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perYear$C_extracted <- apply(C_extractions_per_year[ , 2:101], 1, FUN = mean)
NutrientExtraction_perYear$C_extracted_SD <- apply(C_extractions_per_year[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_year[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perYear$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perYear$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perYear$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perYear$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perYear$C_extracted_highCI <- quantiles_CIs_C[5, ]



###### Nitrogen ######
# Creates data frame for holding N distributions per year.
N_extractions_per_year <- data.frame(year = c(1960:2018))

# For loop to add up N extraction distributions per year.
for (A in 1:nrow(N_extractions_per_year)) {
  
  N_extractions_per_year[A, 2:101] <- colSums(N_extractions[which(N_extractions_per_year$year[A] == Fisheries_NutrientExtraction$year)])
  
}


# The following code generates the mean N extraction for each year between 1960 and 2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perYear$N_extracted <- apply(N_extractions_per_year[ , 2:101], 1, FUN = mean)
NutrientExtraction_perYear$N_extracted_SD <- apply(N_extractions_per_year[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_per_year[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perYear$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perYear$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perYear$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perYear$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perYear$N_extracted_highCI <- quantiles_CIs_N[5, ]



###### Phosphorus ######
# Creates data frame for holding P distributions per year.
P_extractions_per_year <- data.frame(year = c(1960:2018))

# For loop to add up P extraction distributions per year.
for (A in 1:nrow(P_extractions_per_year)) {
  
  P_extractions_per_year[A, 2:101] <- colSums(P_extractions[which(P_extractions_per_year$year[A] == Fisheries_NutrientExtraction$year)])
  
}


# The following code generates the mean P extraction for each year between 1960 and 2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perYear$P_extracted <- apply(P_extractions_per_year[ , 2:101], 1, FUN = mean)
NutrientExtraction_perYear$P_extracted_SD <- apply(P_extractions_per_year[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_per_year[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perYear$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perYear$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perYear$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perYear$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perYear$P_extracted_highCI <- quantiles_CIs_P[5, ]



##### PER TIME PERIOD #####

###### Carbon ######

# Creates data frame to hold time period C extraction distributions.
C_extractions_per_timeperiod <- data.frame(time_period = c("1960_64", "1993_97", "2014_18"))
colnames(C_extractions_per_timeperiod)[1] <- "time_period"

# Creates data frame to hold time period extraction estimates.
NutrientExtraction_perTimePeriod <- data.frame(time_period = c("1960_64", "1993_97", "2014_18"))
NutrientExtraction_perTimePeriod$C_extracted <- NA
NutrientExtraction_perTimePeriod$C_extracted_SD <- NA
NutrientExtraction_perTimePeriod$C_extracted_lowCI <- NA
NutrientExtraction_perTimePeriod$C_extracted_lowIQR <- NA
NutrientExtraction_perTimePeriod$C_extracted_median <- NA
NutrientExtraction_perTimePeriod$C_extracted_highIQR <- NA
NutrientExtraction_perTimePeriod$C_extracted_highCI <- NA



# 1960-1964 Period #

# For loop to add up C extraction distributions for 1960-1964.
for (A in 1:nrow(C_extractions_per_timeperiod)) {
  
  C_extractions_per_timeperiod[1, 2:101] <- colSums(C_extractions[which(Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean C extraction for the time period between 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[1, ]$C_extracted <- apply(C_extractions_per_timeperiod[1 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[1, ]$C_extracted_SD <- apply(C_extractions_per_timeperiod[1 , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_timeperiod[1 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[1, ]$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_highCI <- quantiles_CIs_C[5, ]



# 1993-1997 Period #

# For loop to add up C extraction distributions for 1993-1997.
for (A in 1:nrow(C_extractions_per_timeperiod)) {
  
  C_extractions_per_timeperiod[2, 2:101] <- colSums(C_extractions[which(Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean C extraction for the time period between 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[2, ]$C_extracted <- apply(C_extractions_per_timeperiod[2 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[2, ]$C_extracted_SD <- apply(C_extractions_per_timeperiod[2 , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_timeperiod[2 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[2, ]$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_highCI <- quantiles_CIs_C[5, ]



# 2014-2018 Period #

# For loop to add up C extraction distributions for 2014-2018.
for (A in 1:nrow(C_extractions_per_timeperiod)) {
  
  C_extractions_per_timeperiod[3, 2:101] <- colSums(C_extractions[which(Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean C extraction for the time period between 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[3, ]$C_extracted <- apply(C_extractions_per_timeperiod[3 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[3, ]$C_extracted_SD <- apply(C_extractions_per_timeperiod[3 , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_timeperiod[3 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[3, ]$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_highCI <- quantiles_CIs_C[5, ]




###### Nitrogen ######

# Creates data frame to hold time period distributions.
N_extractions_per_timeperiod <- data.frame(time_period = c("1960_64", "1993_97", "2014_18"))
colnames(N_extractions_per_timeperiod)[1] <- "time_period"

# Adds nitrogen columns to the time period estimates data frame.
NutrientExtraction_perTimePeriod$N_extracted <- NA
NutrientExtraction_perTimePeriod$N_extracted_SD <- NA
NutrientExtraction_perTimePeriod$N_extracted_lowCI <- NA
NutrientExtraction_perTimePeriod$N_extracted_lowIQR <- NA
NutrientExtraction_perTimePeriod$N_extracted_median <- NA
NutrientExtraction_perTimePeriod$N_extracted_highIQR <- NA
NutrientExtraction_perTimePeriod$N_extracted_highCI <- NA



# 1960-1964 Period #

# For loop to add up N extraction distributions for 1960-1964.
for (A in 1:nrow(N_extractions_per_timeperiod)) {
  
  N_extractions_per_timeperiod[1, 2:101] <- colSums(N_extractions[which(Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean N extraction for the time period between 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[1, ]$N_extracted <- apply(N_extractions_per_timeperiod[1 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[1, ]$N_extracted_SD <- apply(N_extractions_per_timeperiod[1 , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_per_timeperiod[1 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[1, ]$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_highCI <- quantiles_CIs_N[5, ]



# 1993-1997 Period #

# For loop to add up N extraction distributions for 1993-1997.
for (A in 1:nrow(N_extractions_per_timeperiod)) {
  
  N_extractions_per_timeperiod[2, 2:101] <- colSums(N_extractions[which(Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean N extraction for the time period between 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[2, ]$N_extracted <- apply(N_extractions_per_timeperiod[2 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[2, ]$N_extracted_SD <- apply(N_extractions_per_timeperiod[2 , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_per_timeperiod[2 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[2, ]$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_highCI <- quantiles_CIs_N[5, ]



# 2014-2018 Period #

# For loop to add up N extraction distributions for 2014-2018.
for (A in 1:nrow(N_extractions_per_timeperiod)) {
  
  N_extractions_per_timeperiod[3, 2:101] <- colSums(N_extractions[which(Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean N extraction for the time period between 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[3, ]$N_extracted <- apply(N_extractions_per_timeperiod[3 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[3, ]$N_extracted_SD <- apply(N_extractions_per_timeperiod[3 , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_per_timeperiod[3 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[3, ]$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_highCI <- quantiles_CIs_N[5, ]



###### Phosphorus ######

# Creates data frame to hold time period distributions.
P_extractions_per_timeperiod <- as.data.frame(c("1960_64", "1993_97", "2014_18"))
colnames(P_extractions_per_timeperiod)[1] <- "time_period"

# Adds phosphorus columns to the time period estimates data frame.
NutrientExtraction_perTimePeriod$P_extracted <- NA
NutrientExtraction_perTimePeriod$P_extracted_SD <- NA
NutrientExtraction_perTimePeriod$P_extracted_lowCI <- NA
NutrientExtraction_perTimePeriod$P_extracted_lowIQR <- NA
NutrientExtraction_perTimePeriod$P_extracted_median <- NA
NutrientExtraction_perTimePeriod$P_extracted_highIQR <- NA
NutrientExtraction_perTimePeriod$P_extracted_highCI <- NA

# 1960-1964 Period #

# For loop to add up P extraction distributions for 1960-1964.
for (A in 1:nrow(P_extractions_per_timeperiod)) {
  
  P_extractions_per_timeperiod[1, 2:101] <- colSums(P_extractions[which(Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean P extraction for the time period between 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[1, ]$P_extracted <- apply(P_extractions_per_timeperiod[1 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[1, ]$P_extracted_SD <- apply(P_extractions_per_timeperiod[1 , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_per_timeperiod[1 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[1, ]$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_highCI <- quantiles_CIs_P[5, ]



# 1993-1997 Period #

# For loop to add up P extraction distributions for 1993-1997.
for (A in 1:nrow(P_extractions_per_timeperiod)) {
  
  P_extractions_per_timeperiod[2, 2:101] <- colSums(P_extractions[which(Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean P extraction for the time period between 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[2, ]$P_extracted <- apply(P_extractions_per_timeperiod[2 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[2, ]$P_extracted_SD <- apply(P_extractions_per_timeperiod[2 , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_per_timeperiod[2 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[2, ]$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_highCI <- quantiles_CIs_P[5, ]



# 2014-2018 Period #

# For loop to add up P extraction distributions for 2014-2018.
for (A in 1:nrow(P_extractions_per_timeperiod)) {
  
  P_extractions_per_timeperiod[3, 2:101] <- colSums(P_extractions[which(Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean P extraction for the time period between 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTimePeriod[3, ]$P_extracted <- apply(P_extractions_per_timeperiod[3 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[3, ]$P_extracted_SD <- apply(P_extractions_per_timeperiod[3 , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_per_timeperiod[3 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[3, ]$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_highCI <- quantiles_CIs_P[5, ]





##### PER MARINE REGION #####

###### Carbon ######

# Creates data frame to hold marine region C extraction distributions.
C_extractions_perArea <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))

# Creates data frame to hold marine region extraction estimates.
NutrientExtraction_perArea <- data.frame(
  area_name = unique(Fisheries_NutrientExtraction$area_name),
  area_km2 = Fisheries_NutrientExtraction$area_km2[match(unique(Fisheries_NutrientExtraction$area_name), Fisheries_NutrientExtraction$area_name)]
)



# All Years (1960-2018) #

# For loop to add up C extraction distributions per marine region for 1960-2018.
for (A in 1:nrow(C_extractions_perArea)) {
  
  C_extractions_perArea[A, 2:101] <- colSums(C_extractions[which(C_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean C extraction per marine region for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$C_extracted <- apply(C_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$C_extracted_SD <- apply(C_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perArea$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perArea$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perArea$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perArea$C_extracted_highCI <- quantiles_CIs_C[5, ]



# 1960-1964 Period #

# For loop to add up C extraction distributions per marine region for 1960-1964.
for (A in 1:nrow(C_extractions_perArea)) {
  
  C_extractions_perArea[A, 2:101] <- colSums(C_extractions[which(C_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean C extraction per marine region for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$C_extracted_1960_64 <- apply(C_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$C_extracted_SD_1960_64 <- apply(C_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$C_extracted_lowCI_1960_64 <- quantiles_CIs_C[1, ]
NutrientExtraction_perArea$C_extracted_lowIQR_1960_64 <- quantiles_CIs_C[2, ]
NutrientExtraction_perArea$C_extracted_median_1960_64 <- quantiles_CIs_C[3, ]
NutrientExtraction_perArea$C_extracted_highIQR_1960_64 <- quantiles_CIs_C[4, ]
NutrientExtraction_perArea$C_extracted_highCI_1960_64 <- quantiles_CIs_C[5, ]



# 1993-1997 Period #

# For loop to add up C extraction distributions per marine region for 1993-1997.
for (A in 1:nrow(C_extractions_perArea)) {
  
  C_extractions_perArea[A, 2:101] <- colSums(C_extractions[which(C_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean C extraction per marine region for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$C_extracted_1993_97 <- apply(C_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$C_extracted_SD_1993_97 <- apply(C_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$C_extracted_lowCI_1993_97 <- quantiles_CIs_C[1, ]
NutrientExtraction_perArea$C_extracted_lowIQR_1993_97 <- quantiles_CIs_C[2, ]
NutrientExtraction_perArea$C_extracted_median_1993_97 <- quantiles_CIs_C[3, ]
NutrientExtraction_perArea$C_extracted_highIQR_1993_97 <- quantiles_CIs_C[4, ]
NutrientExtraction_perArea$C_extracted_highCI_1993_97 <- quantiles_CIs_C[5, ]




# 2014-2018 Period #

# For loop to add up C extraction distributions per marine region for 2014-2018.
for (A in 1:nrow(C_extractions_perArea)) {
  
  C_extractions_perArea[A, 2:101] <- colSums(C_extractions[which(C_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean C extraction per marine region for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$C_extracted_2014_18 <- apply(C_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$C_extracted_SD_2014_18 <- apply(C_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$C_extracted_lowCI_2014_18 <- quantiles_CIs_C[1, ]
NutrientExtraction_perArea$C_extracted_lowIQR_2014_18 <- quantiles_CIs_C[2, ]
NutrientExtraction_perArea$C_extracted_median_2014_18 <- quantiles_CIs_C[3, ]
NutrientExtraction_perArea$C_extracted_highIQR_2014_18 <- quantiles_CIs_C[4, ]
NutrientExtraction_perArea$C_extracted_highCI_2014_18 <- quantiles_CIs_C[5, ]




###### Nitrogen ######

# Creates data frame to hold marine region N extraction distributions.
N_extractions_perArea <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))


# All Years (1960-2018) #

# For loop to add up N extraction distributions per marine region for 1960-2018.
for (A in 1:nrow(N_extractions_perArea)) {
  
  N_extractions_perArea[A, 2:101] <- colSums(N_extractions[which(N_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                                   Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean N extraction per marine region for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$N_extracted <- apply(N_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$N_extracted_SD <- apply(N_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perArea$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perArea$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perArea$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perArea$N_extracted_highCI <- quantiles_CIs_N[5, ]



# 1960-1964 Period #

# For loop to add up N extraction distributions per marine region for 1960-1964.
for (A in 1:nrow(N_extractions_perArea)) {
  
  N_extractions_perArea[A, 2:101] <- colSums(N_extractions[which(N_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean N extraction per marine region for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$N_extracted_1960_64 <- apply(N_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$N_extracted_SD_1960_64 <- apply(N_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$N_extracted_lowCI_1960_64 <- quantiles_CIs_N[1, ]
NutrientExtraction_perArea$N_extracted_lowIQR_1960_64 <- quantiles_CIs_N[2, ]
NutrientExtraction_perArea$N_extracted_median_1960_64 <- quantiles_CIs_N[3, ]
NutrientExtraction_perArea$N_extracted_highIQR_1960_64 <- quantiles_CIs_N[4, ]
NutrientExtraction_perArea$N_extracted_highCI_1960_64 <- quantiles_CIs_N[5, ]



# 1993-1997 Period #

# For loop to add up N extraction distributions per marine region for 1993-1997.
for (A in 1:nrow(N_extractions_perArea)) {
  
  N_extractions_perArea[A, 2:101] <- colSums(N_extractions[which(N_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean N extraction per marine region for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$N_extracted_1993_97 <- apply(N_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$N_extracted_SD_1993_97 <- apply(N_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$N_extracted_lowCI_1993_97 <- quantiles_CIs_N[1, ]
NutrientExtraction_perArea$N_extracted_lowIQR_1993_97 <- quantiles_CIs_N[2, ]
NutrientExtraction_perArea$N_extracted_median_1993_97 <- quantiles_CIs_N[3, ]
NutrientExtraction_perArea$N_extracted_highIQR_1993_97 <- quantiles_CIs_N[4, ]
NutrientExtraction_perArea$N_extracted_highCI_1993_97 <- quantiles_CIs_N[5, ]



# 2014-2018 Period #

# For loop to add up N extraction distributions per marine region for 2014-2018.
for (A in 1:nrow(N_extractions_perArea)) {
  
  N_extractions_perArea[A, 2:101] <- colSums(N_extractions[which(N_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean N extraction per marine region for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$N_extracted_2014_18 <- apply(N_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$N_extracted_SD_2014_18 <- apply(N_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$N_extracted_lowCI_2014_18 <- quantiles_CIs_N[1, ]
NutrientExtraction_perArea$N_extracted_lowIQR_2014_18 <- quantiles_CIs_N[2, ]
NutrientExtraction_perArea$N_extracted_median_2014_18 <- quantiles_CIs_N[3, ]
NutrientExtraction_perArea$N_extracted_highIQR_2014_18 <- quantiles_CIs_N[4, ]
NutrientExtraction_perArea$N_extracted_highCI_2014_18 <- quantiles_CIs_N[5, ]




###### Phosphorus ######

# Creates data frame to hold marine region P extraction distributions.
P_extractions_perArea <- data.frame(area_name = unique(Fisheries_NutrientExtraction$area_name))


# All Years (1960-2018) #

# For loop to add up P extraction distributions per marine region for 1960-2018.
for (A in 1:nrow(P_extractions_perArea)) {
  
  P_extractions_perArea[A, 2:101] <- colSums(P_extractions[which(P_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name &
                                                                   Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean P extraction per marine region for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$P_extracted <- apply(P_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$P_extracted_SD <- apply(P_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perArea$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perArea$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perArea$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perArea$P_extracted_highCI <- quantiles_CIs_P[5, ]



# 1960-1964 Period #

# For loop to add up P extraction distributions per marine region for 1960-1964.
for (A in 1:nrow(P_extractions_perArea)) {
  
  P_extractions_perArea[A, 2:101] <- colSums(P_extractions[which(P_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean P extraction per marine region for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$P_extracted_1960_64 <- apply(P_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$P_extracted_SD_1960_64 <- apply(P_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$P_extracted_lowCI_1960_64 <- quantiles_CIs_P[1, ]
NutrientExtraction_perArea$P_extracted_lowIQR_1960_64 <- quantiles_CIs_P[2, ]
NutrientExtraction_perArea$P_extracted_median_1960_64 <- quantiles_CIs_P[3, ]
NutrientExtraction_perArea$P_extracted_highIQR_1960_64 <- quantiles_CIs_P[4, ]
NutrientExtraction_perArea$P_extracted_highCI_1960_64 <- quantiles_CIs_P[5, ]



# 1993-1997 Period #

# For loop to add up P extraction distributions per marine region for 1993-1997.
for (A in 1:nrow(P_extractions_perArea)) {
  
  P_extractions_perArea[A, 2:101] <- colSums(P_extractions[which(P_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean P extraction per marine region for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$P_extracted_1993_97 <- apply(P_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$P_extracted_SD_1993_97 <- apply(P_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$P_extracted_lowCI_1993_97 <- quantiles_CIs_P[1, ]
NutrientExtraction_perArea$P_extracted_lowIQR_1993_97 <- quantiles_CIs_P[2, ]
NutrientExtraction_perArea$P_extracted_median_1993_97 <- quantiles_CIs_P[3, ]
NutrientExtraction_perArea$P_extracted_highIQR_1993_97 <- quantiles_CIs_P[4, ]
NutrientExtraction_perArea$P_extracted_highCI_1993_97 <- quantiles_CIs_P[5, ]



# 2014-2018 Period #

# For loop to add up P extraction distributions per marine region for 2014-2018.
for (A in 1:nrow(P_extractions_perArea)) {
  
  P_extractions_perArea[A, 2:101] <- colSums(P_extractions[which(P_extractions_perArea$area_name[A] == Fisheries_NutrientExtraction$area_name & 
                                                                   Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean P extraction per marine region for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perArea$P_extracted_2014_18 <- apply(P_extractions_perArea[ , 2:101], 1, FUN = mean)
NutrientExtraction_perArea$P_extracted_SD_2014_18 <- apply(P_extractions_perArea[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perArea[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perArea$P_extracted_lowCI_2014_18 <- quantiles_CIs_P[1, ]
NutrientExtraction_perArea$P_extracted_lowIQR_2014_18 <- quantiles_CIs_P[2, ]
NutrientExtraction_perArea$P_extracted_median_2014_18 <- quantiles_CIs_P[3, ]
NutrientExtraction_perArea$P_extracted_highIQR_2014_18 <- quantiles_CIs_P[4, ]
NutrientExtraction_perArea$P_extracted_highCI_2014_18 <- quantiles_CIs_P[5, ]





###### Per Square Kilometer ######

# This piece is to calculate the total extractions per marine region on a per-square-kilometer basis. The estimates were done for the entire study period (1960-2018) and per time period (1960s, 1990s, 2010s).


# Carbon
NutrientExtraction_perArea <- NutrientExtraction_perArea %>% 
  mutate(C_extracted_perSqKm = C_extracted/area_km2,
         C_extracted_perSqKm_SD = C_extracted_SD/area_km2,
         C_extracted_perSqKm_1960_64 = C_extracted_1960_64/area_km2,
         C_extracted_1960_64_perSqKm_SD = C_extracted_SD_1960_64/area_km2,
         C_extracted_perSqKm_1993_97 = C_extracted_1993_97/area_km2,
         C_extracted_1993_97_perSqKm_SD = C_extracted_SD_1993_97/area_km2,
         C_extracted_perSqKm_2014_18 = C_extracted_2014_18/area_km2,
         C_extracted_2014_18_perSqKm_SD = C_extracted_SD_2014_18/area_km2)


# Nitrogen
NutrientExtraction_perArea <- NutrientExtraction_perArea %>% 
  mutate(N_extracted_perSqKm = N_extracted/area_km2,
         N_extracted_perSqKm_SD = N_extracted_SD/area_km2,
         N_extracted_perSqKm_1960_64 = N_extracted_1960_64/area_km2,
         N_extracted_1960_64_perSqKm_SD = N_extracted_SD_1960_64/area_km2,
         N_extracted_perSqKm_1993_97 = N_extracted_1993_97/area_km2,
         N_extracted_1993_97_perSqKm_SD = N_extracted_SD_1993_97/area_km2,
         N_extracted_perSqKm_2014_18 = N_extracted_2014_18/area_km2,
         N_extracted_2014_18_perSqKm_SD = N_extracted_SD_2014_18/area_km2)


# Phosphorus
NutrientExtraction_perArea <- NutrientExtraction_perArea %>% 
  mutate(P_extracted_perSqKm = P_extracted/area_km2,
         P_extracted_perSqKm_SD = P_extracted_SD/area_km2,
         P_extracted_perSqKm_1960_64 = P_extracted_1960_64/area_km2,
         P_extracted_1960_64_perSqKm_SD = P_extracted_SD_1960_64/area_km2,
         P_extracted_perSqKm_1993_97 = P_extracted_1993_97/area_km2,
         P_extracted_1993_97_perSqKm_SD = P_extracted_SD_1993_97/area_km2,
         P_extracted_perSqKm_2014_18 = P_extracted_2014_18/area_km2,
         P_extracted_2014_18_perSqKm_SD = P_extracted_SD_2014_18/area_km2)





###### Percent Change Between Periods ######

# This piece is to calculate the percent change in nutrient extractions between periods for each marine region.

NutrientExtraction_perArea <- NutrientExtraction_perArea %>% 
  mutate(C_extracted_perSqKm_60s_90s = 
           ((C_extracted_perSqKm_1993_97 - C_extracted_perSqKm_1960_64)/
              C_extracted_perSqKm_1960_64) *100,
         C_extracted_perSqKm_90s_10s = 
           ((C_extracted_perSqKm_2014_18 - C_extracted_perSqKm_1993_97)/
              C_extracted_perSqKm_1993_97) *100,
         N_extracted_perSqKm_60s_90s = 
           ((N_extracted_perSqKm_1993_97 - N_extracted_perSqKm_1960_64)/
              N_extracted_perSqKm_1960_64) *100,
         N_extracted_perSqKm_90s_10s = 
           ((N_extracted_perSqKm_2014_18 - N_extracted_perSqKm_1993_97)/
              N_extracted_perSqKm_1993_97) *100,
         P_extracted_perSqKm_60s_90s = 
           ((P_extracted_perSqKm_1993_97 - P_extracted_perSqKm_1960_64)/
              P_extracted_perSqKm_1960_64) *100,
         P_extracted_perSqKm_90s_10s = 
           ((P_extracted_perSqKm_2014_18 - P_extracted_perSqKm_1993_97)/
              P_extracted_perSqKm_1993_97) *100)





##### PER TROPHIC GROUP #####

###### Carbon ######

# Creates data frame to hold trophic group C extraction distributions.
C_extractions_perTG <- data.frame(trophic_group = unique(Fisheries_NutrientExtraction$trophic_group))

# Creates data frame to hold trophic group extraction estimates.
NutrientExtraction_perTG <- data.frame(trophic_group = unique(Fisheries_NutrientExtraction$trophic_group))



# All Years (1960-2018) #

# For loop to add up C extraction distributions per trophic group for 1960-2018.
for (A in 1:nrow(C_extractions_perTG)) {
  
  C_extractions_perTG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean C extraction per trophic group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$C_extracted <- apply(C_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$C_extracted_SD <- apply(C_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perTG$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perTG$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perTG$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perTG$C_extracted_highCI <- quantiles_CIs_C[5, ]



# 1960-1964 Period #

# For loop to add up C extraction distributions per trophic group for 1960-1964.
for (A in 1:nrow(C_extractions_perTG)) {
  
  C_extractions_perTG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean C extraction per trophic group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$C_extracted_1960_64 <- apply(C_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$C_extracted_SD_1960_64 <- apply(C_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$C_extracted_lowCI_1960_64 <- quantiles_CIs_C[1, ]
NutrientExtraction_perTG$C_extracted_lowIQR_1960_64 <- quantiles_CIs_C[2, ]
NutrientExtraction_perTG$C_extracted_median_1960_64 <- quantiles_CIs_C[3, ]
NutrientExtraction_perTG$C_extracted_highIQR_1960_64 <- quantiles_CIs_C[4, ]
NutrientExtraction_perTG$C_extracted_highCI_1960_64 <- quantiles_CIs_C[5, ]



# 1993-1997 Period #

# For loop to add up C extraction distributions per trophic group for 1993-1997.
for (A in 1:nrow(C_extractions_perTG)) {
  
  C_extractions_perTG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean C extraction per trophic group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$C_extracted_1993_97 <- apply(C_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$C_extracted_SD_1993_97 <- apply(C_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$C_extracted_lowCI_1993_97 <- quantiles_CIs_C[1, ]
NutrientExtraction_perTG$C_extracted_lowIQR_1993_97 <- quantiles_CIs_C[2, ]
NutrientExtraction_perTG$C_extracted_median_1993_97 <- quantiles_CIs_C[3, ]
NutrientExtraction_perTG$C_extracted_highIQR_1993_97 <- quantiles_CIs_C[4, ]
NutrientExtraction_perTG$C_extracted_highCI_1993_97 <- quantiles_CIs_C[5, ]




# 2014-2018 Period #

# For loop to add up C extraction distributions per trophic group for 2014-2018.
for (A in 1:nrow(C_extractions_perTG)) {
  
  C_extractions_perTG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean C extraction per trophic group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$C_extracted_2014_18 <- apply(C_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$C_extracted_SD_2014_18 <- apply(C_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$C_extracted_lowCI_2014_18 <- quantiles_CIs_C[1, ]
NutrientExtraction_perTG$C_extracted_lowIQR_2014_18 <- quantiles_CIs_C[2, ]
NutrientExtraction_perTG$C_extracted_median_2014_18 <- quantiles_CIs_C[3, ]
NutrientExtraction_perTG$C_extracted_highIQR_2014_18 <- quantiles_CIs_C[4, ]
NutrientExtraction_perTG$C_extracted_highCI_2014_18 <- quantiles_CIs_C[5, ]




###### Nitrogen ######

# Creates data frame to hold trophic group N extraction distributions.
N_extractions_perTG <- data.frame(trophic_group = unique(Fisheries_NutrientExtraction$trophic_group))


# All Years (1960-2018) #

# For loop to add up N extraction distributions per trophic group for 1960-2018.
for (A in 1:nrow(N_extractions_perTG)) {
  
  N_extractions_perTG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group &
                                                 Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean N extraction per trophic group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$N_extracted <- apply(N_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$N_extracted_SD <- apply(N_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perTG$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perTG$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perTG$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perTG$N_extracted_highCI <- quantiles_CIs_N[5, ]



# 1960-1964 Period #

# For loop to add up N extraction distributions per trophic group for 1960-1964.
for (A in 1:nrow(N_extractions_perTG)) {
  
  N_extractions_perTG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean N extraction per trophic group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$N_extracted_1960_64 <- apply(N_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$N_extracted_SD_1960_64 <- apply(N_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$N_extracted_lowCI_1960_64 <- quantiles_CIs_N[1, ]
NutrientExtraction_perTG$N_extracted_lowIQR_1960_64 <- quantiles_CIs_N[2, ]
NutrientExtraction_perTG$N_extracted_median_1960_64 <- quantiles_CIs_N[3, ]
NutrientExtraction_perTG$N_extracted_highIQR_1960_64 <- quantiles_CIs_N[4, ]
NutrientExtraction_perTG$N_extracted_highCI_1960_64 <- quantiles_CIs_N[5, ]



# 1993-1997 Period #

# For loop to add up N extraction distributions per trophic group for 1993-1997.
for (A in 1:nrow(N_extractions_perTG)) {
  
  N_extractions_perTG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean N extraction per trophic group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$N_extracted_1993_97 <- apply(N_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$N_extracted_SD_1993_97 <- apply(N_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$N_extracted_lowCI_1993_97 <- quantiles_CIs_N[1, ]
NutrientExtraction_perTG$N_extracted_lowIQR_1993_97 <- quantiles_CIs_N[2, ]
NutrientExtraction_perTG$N_extracted_median_1993_97 <- quantiles_CIs_N[3, ]
NutrientExtraction_perTG$N_extracted_highIQR_1993_97 <- quantiles_CIs_N[4, ]
NutrientExtraction_perTG$N_extracted_highCI_1993_97 <- quantiles_CIs_N[5, ]



# 2014-2018 Period #

# For loop to add up N extraction distributions per trophic group for 2014-2018.
for (A in 1:nrow(N_extractions_perTG)) {
  
  N_extractions_perTG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean N extraction per trophic group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$N_extracted_2014_18 <- apply(N_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$N_extracted_SD_2014_18 <- apply(N_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$N_extracted_lowCI_2014_18 <- quantiles_CIs_N[1, ]
NutrientExtraction_perTG$N_extracted_lowIQR_2014_18 <- quantiles_CIs_N[2, ]
NutrientExtraction_perTG$N_extracted_median_2014_18 <- quantiles_CIs_N[3, ]
NutrientExtraction_perTG$N_extracted_highIQR_2014_18 <- quantiles_CIs_N[4, ]
NutrientExtraction_perTG$N_extracted_highCI_2014_18 <- quantiles_CIs_N[5, ]




###### Phosphorus ######

# Creates data frame to hold trophic group P extraction distributions.
P_extractions_perTG <- data.frame(trophic_group = unique(Fisheries_NutrientExtraction$trophic_group))


# All Years (1960-2018) #

# For loop to add up P extraction distributions per trophic group for 1960-2018.
for (A in 1:nrow(P_extractions_perTG)) {
  
  P_extractions_perTG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group &
                                                 Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean P extraction per trophic group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$P_extracted <- apply(P_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$P_extracted_SD <- apply(P_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perTG$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perTG$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perTG$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perTG$P_extracted_highCI <- quantiles_CIs_P[5, ]



# 1960-1964 Period #

# For loop to add up P extraction distributions per trophic group for 1960-1964.
for (A in 1:nrow(P_extractions_perTG)) {
  
  P_extractions_perTG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean P extraction per trophic group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$P_extracted_1960_64 <- apply(P_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$P_extracted_SD_1960_64 <- apply(P_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$P_extracted_lowCI_1960_64 <- quantiles_CIs_P[1, ]
NutrientExtraction_perTG$P_extracted_lowIQR_1960_64 <- quantiles_CIs_P[2, ]
NutrientExtraction_perTG$P_extracted_median_1960_64 <- quantiles_CIs_P[3, ]
NutrientExtraction_perTG$P_extracted_highIQR_1960_64 <- quantiles_CIs_P[4, ]
NutrientExtraction_perTG$P_extracted_highCI_1960_64 <- quantiles_CIs_P[5, ]



# 1993-1997 Period #

# For loop to add up P extraction distributions per trophic group for 1993-1997.
for (A in 1:nrow(P_extractions_perTG)) {
  
  P_extractions_perTG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean P extraction per trophic group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$P_extracted_1993_97 <- apply(P_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$P_extracted_SD_1993_97 <- apply(P_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$P_extracted_lowCI_1993_97 <- quantiles_CIs_P[1, ]
NutrientExtraction_perTG$P_extracted_lowIQR_1993_97 <- quantiles_CIs_P[2, ]
NutrientExtraction_perTG$P_extracted_median_1993_97 <- quantiles_CIs_P[3, ]
NutrientExtraction_perTG$P_extracted_highIQR_1993_97 <- quantiles_CIs_P[4, ]
NutrientExtraction_perTG$P_extracted_highCI_1993_97 <- quantiles_CIs_P[5, ]



# 2014-2018 Period #

# For loop to add up P extraction distributions per trophic group for 2014-2018.
for (A in 1:nrow(P_extractions_perTG)) {
  
  P_extractions_perTG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perTG$trophic_group[A] == Fisheries_NutrientExtraction$trophic_group & 
                                                 Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean P extraction per trophic group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perTG$P_extracted_2014_18 <- apply(P_extractions_perTG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perTG$P_extracted_SD_2014_18 <- apply(P_extractions_perTG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perTG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTG$P_extracted_lowCI_2014_18 <- quantiles_CIs_P[1, ]
NutrientExtraction_perTG$P_extracted_lowIQR_2014_18 <- quantiles_CIs_P[2, ]
NutrientExtraction_perTG$P_extracted_median_2014_18 <- quantiles_CIs_P[3, ]
NutrientExtraction_perTG$P_extracted_highIQR_2014_18 <- quantiles_CIs_P[4, ]
NutrientExtraction_perTG$P_extracted_highCI_2014_18 <- quantiles_CIs_P[5, ]





##### PER SIMPLE FUNCTIONAL GROUP #####

# These are the functional groups that have been grouped together for simplicity as it was presented in the manuscript. The functional groups were classified in the fisheries data frame as follows:

Fisheries_NutrientExtraction$simp_functional_group <- "NA"
Fisheries_NutrientExtraction[functional_group_id %in% c(1:3)]$simp_functional_group <- "Pelagic"
Fisheries_NutrientExtraction[functional_group_id %in% c(4:6)]$simp_functional_group <- "Demersal"
Fisheries_NutrientExtraction[functional_group_id %in% c(7:9)]$simp_functional_group <- "Bathypelagic"
Fisheries_NutrientExtraction[functional_group_id %in% c(10:12)]$simp_functional_group <- "Bathydemersals"
Fisheries_NutrientExtraction[functional_group_id %in% c(13:15)]$simp_functional_group <- "Benthopelagic"
Fisheries_NutrientExtraction[functional_group_id %in% c(16:18)]$simp_functional_group <- "Reef Fish"
Fisheries_NutrientExtraction[functional_group_id %in% c(19:20)]$simp_functional_group <- "Sharks"
Fisheries_NutrientExtraction[functional_group_id %in% c(21:22)]$simp_functional_group <- "Rays"
Fisheries_NutrientExtraction[functional_group_id %in% c(23:24)]$simp_functional_group <- "Flatfish"
Fisheries_NutrientExtraction[functional_group_id == 25]$simp_functional_group <- "Cephalopods"
Fisheries_NutrientExtraction[functional_group_id == 26]$simp_functional_group <- "Shrimps"
Fisheries_NutrientExtraction[functional_group_id == 27]$simp_functional_group <- "Crabs/Lobsters"
Fisheries_NutrientExtraction[functional_group_id == 28]$simp_functional_group <- "Jellyfish"
Fisheries_NutrientExtraction[functional_group_id == 29]$simp_functional_group <- "Misc. Dem. Inverts"
Fisheries_NutrientExtraction[functional_group_id == 30]$simp_functional_group <- "Krill"



###### Carbon ######

# Creates data frame to hold simple functional group C extraction distributions.
C_extractions_perSFG <- data.frame(simp_functional_group = unique(Fisheries_NutrientExtraction$simp_functional_group))

# Creates data frame to hold simple functional group extraction estimates.
NutrientExtraction_perSFG <- data.frame(simp_functional_group = unique(Fisheries_NutrientExtraction$simp_functional_group))



# All Years (1960-2018) #

# For loop to add up C extraction distributions per simple functional group for 1960-2018.
for (A in 1:nrow(C_extractions_perSFG)) {
  
  C_extractions_perSFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean C extraction per simple functional group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$C_extracted <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$C_extracted_SD <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perSFG$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perSFG$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perSFG$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perSFG$C_extracted_highCI <- quantiles_CIs_C[5, ]



# 1960-1964 Period #

# For loop to add up C extraction distributions per simple functional group for 1960-1964.
for (A in 1:nrow(C_extractions_perSFG)) {
  
  C_extractions_perSFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean C extraction per simple functional group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$C_extracted_1960_64 <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$C_extracted_SD_1960_64 <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$C_extracted_lowCI_1960_64 <- quantiles_CIs_C[1, ]
NutrientExtraction_perSFG$C_extracted_lowIQR_1960_64 <- quantiles_CIs_C[2, ]
NutrientExtraction_perSFG$C_extracted_median_1960_64 <- quantiles_CIs_C[3, ]
NutrientExtraction_perSFG$C_extracted_highIQR_1960_64 <- quantiles_CIs_C[4, ]
NutrientExtraction_perSFG$C_extracted_highCI_1960_64 <- quantiles_CIs_C[5, ]



# 1993-1997 Period #

# For loop to add up C extraction distributions per simple functional group for 1993-1997.
for (A in 1:nrow(C_extractions_perSFG)) {
  
  C_extractions_perSFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean C extraction per simple functional group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$C_extracted_1993_97 <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$C_extracted_SD_1993_97 <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$C_extracted_lowCI_1993_97 <- quantiles_CIs_C[1, ]
NutrientExtraction_perSFG$C_extracted_lowIQR_1993_97 <- quantiles_CIs_C[2, ]
NutrientExtraction_perSFG$C_extracted_median_1993_97 <- quantiles_CIs_C[3, ]
NutrientExtraction_perSFG$C_extracted_highIQR_1993_97 <- quantiles_CIs_C[4, ]
NutrientExtraction_perSFG$C_extracted_highCI_1993_97 <- quantiles_CIs_C[5, ]




# 2014-2018 Period #

# For loop to add up C extraction distributions per simple functional group for 2014-2018.
for (A in 1:nrow(C_extractions_perSFG)) {
  
  C_extractions_perSFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean C extraction per simple functional group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$C_extracted_2014_18 <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$C_extracted_SD_2014_18 <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$C_extracted_lowCI_2014_18 <- quantiles_CIs_C[1, ]
NutrientExtraction_perSFG$C_extracted_lowIQR_2014_18 <- quantiles_CIs_C[2, ]
NutrientExtraction_perSFG$C_extracted_median_2014_18 <- quantiles_CIs_C[3, ]
NutrientExtraction_perSFG$C_extracted_highIQR_2014_18 <- quantiles_CIs_C[4, ]
NutrientExtraction_perSFG$C_extracted_highCI_2014_18 <- quantiles_CIs_C[5, ]




###### Nitrogen ######

# Creates data frame to hold simple functional group N extraction distributions.
N_extractions_perSFG <- data.frame(simp_functional_group = unique(Fisheries_NutrientExtraction$simp_functional_group))


# All Years (1960-2018) #

# For loop to add up N extraction distributions per simple functional group for 1960-2018.
for (A in 1:nrow(N_extractions_perSFG)) {
  
  N_extractions_perSFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group &
                                                                  Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean N extraction per simple functional group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$N_extracted <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$N_extracted_SD <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perSFG$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perSFG$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perSFG$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perSFG$N_extracted_highCI <- quantiles_CIs_N[5, ]



# 1960-1964 Period #

# For loop to add up N extraction distributions per simple functional group for 1960-1964.
for (A in 1:nrow(N_extractions_perSFG)) {
  
  N_extractions_perSFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean N extraction per simple functional group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$N_extracted_1960_64 <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$N_extracted_SD_1960_64 <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$N_extracted_lowCI_1960_64 <- quantiles_CIs_N[1, ]
NutrientExtraction_perSFG$N_extracted_lowIQR_1960_64 <- quantiles_CIs_N[2, ]
NutrientExtraction_perSFG$N_extracted_median_1960_64 <- quantiles_CIs_N[3, ]
NutrientExtraction_perSFG$N_extracted_highIQR_1960_64 <- quantiles_CIs_N[4, ]
NutrientExtraction_perSFG$N_extracted_highCI_1960_64 <- quantiles_CIs_N[5, ]



# 1993-1997 Period #

# For loop to add up N extraction distributions per simple functional group for 1993-1997.
for (A in 1:nrow(N_extractions_perSFG)) {
  
  N_extractions_perSFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean N extraction per simple functional group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$N_extracted_1993_97 <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$N_extracted_SD_1993_97 <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$N_extracted_lowCI_1993_97 <- quantiles_CIs_N[1, ]
NutrientExtraction_perSFG$N_extracted_lowIQR_1993_97 <- quantiles_CIs_N[2, ]
NutrientExtraction_perSFG$N_extracted_median_1993_97 <- quantiles_CIs_N[3, ]
NutrientExtraction_perSFG$N_extracted_highIQR_1993_97 <- quantiles_CIs_N[4, ]
NutrientExtraction_perSFG$N_extracted_highCI_1993_97 <- quantiles_CIs_N[5, ]



# 2014-2018 Period #

# For loop to add up N extraction distributions per simple functional group for 2014-2018.
for (A in 1:nrow(N_extractions_perSFG)) {
  
  N_extractions_perSFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean N extraction per simple functional group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$N_extracted_2014_18 <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$N_extracted_SD_2014_18 <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$N_extracted_lowCI_2014_18 <- quantiles_CIs_N[1, ]
NutrientExtraction_perSFG$N_extracted_lowIQR_2014_18 <- quantiles_CIs_N[2, ]
NutrientExtraction_perSFG$N_extracted_median_2014_18 <- quantiles_CIs_N[3, ]
NutrientExtraction_perSFG$N_extracted_highIQR_2014_18 <- quantiles_CIs_N[4, ]
NutrientExtraction_perSFG$N_extracted_highCI_2014_18 <- quantiles_CIs_N[5, ]




###### Phosphorus ######

# Creates data frame to hold simple functional group P extraction distributions.
P_extractions_perSFG <- data.frame(simp_functional_group = unique(Fisheries_NutrientExtraction$simp_functional_group))


# All Years (1960-2018) #

# For loop to add up P extraction distributions per simple functional group for 1960-2018.
for (A in 1:nrow(P_extractions_perSFG)) {
  
  P_extractions_perSFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group &
                                                                  Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean P extraction per simple functional group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$P_extracted <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$P_extracted_SD <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perSFG$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perSFG$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perSFG$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perSFG$P_extracted_highCI <- quantiles_CIs_P[5, ]



# 1960-1964 Period #

# For loop to add up P extraction distributions per simple functional group for 1960-1964.
for (A in 1:nrow(P_extractions_perSFG)) {
  
  P_extractions_perSFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean P extraction per simple functional group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$P_extracted_1960_64 <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$P_extracted_SD_1960_64 <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$P_extracted_lowCI_1960_64 <- quantiles_CIs_P[1, ]
NutrientExtraction_perSFG$P_extracted_lowIQR_1960_64 <- quantiles_CIs_P[2, ]
NutrientExtraction_perSFG$P_extracted_median_1960_64 <- quantiles_CIs_P[3, ]
NutrientExtraction_perSFG$P_extracted_highIQR_1960_64 <- quantiles_CIs_P[4, ]
NutrientExtraction_perSFG$P_extracted_highCI_1960_64 <- quantiles_CIs_P[5, ]



# 1993-1997 Period #

# For loop to add up P extraction distributions per simple functional group for 1993-1997.
for (A in 1:nrow(P_extractions_perSFG)) {
  
  P_extractions_perSFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean P extraction per simple functional group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$P_extracted_1993_97 <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$P_extracted_SD_1993_97 <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$P_extracted_lowCI_1993_97 <- quantiles_CIs_P[1, ]
NutrientExtraction_perSFG$P_extracted_lowIQR_1993_97 <- quantiles_CIs_P[2, ]
NutrientExtraction_perSFG$P_extracted_median_1993_97 <- quantiles_CIs_P[3, ]
NutrientExtraction_perSFG$P_extracted_highIQR_1993_97 <- quantiles_CIs_P[4, ]
NutrientExtraction_perSFG$P_extracted_highCI_1993_97 <- quantiles_CIs_P[5, ]



# 2014-2018 Period #

# For loop to add up P extraction distributions per simple functional group for 2014-2018.
for (A in 1:nrow(P_extractions_perSFG)) {
  
  P_extractions_perSFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perSFG$simp_functional_group[A] == Fisheries_NutrientExtraction$simp_functional_group & 
                                                                  Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean P extraction per simple functional group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perSFG$P_extracted_2014_18 <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perSFG$P_extracted_SD_2014_18 <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perSFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perSFG$P_extracted_lowCI_2014_18 <- quantiles_CIs_P[1, ]
NutrientExtraction_perSFG$P_extracted_lowIQR_2014_18 <- quantiles_CIs_P[2, ]
NutrientExtraction_perSFG$P_extracted_median_2014_18 <- quantiles_CIs_P[3, ]
NutrientExtraction_perSFG$P_extracted_highIQR_2014_18 <- quantiles_CIs_P[4, ]
NutrientExtraction_perSFG$P_extracted_highCI_2014_18 <- quantiles_CIs_P[5, ]





##### PER FUNCTIONAL GROUP #####

# These estimates are divided by all of the functional groups listed in the SAU dataset including those that are divided up by size groups. For example, there are seperate groups for small-, medium-, and large-sized pelagic species.

###### Carbon ######

# Creates data frame to hold functional group C extraction distributions.
C_extractions_perFG <- data.frame(functional_group_id = unique(Fisheries_NutrientExtraction$functional_group_id))

# Creates data frame to hold functional group extraction estimates.
NutrientExtraction_perFG <- data.frame(functional_group_id = unique(Fisheries_NutrientExtraction$functional_group_id), functional_group = unique(Fisheries_NutrientExtraction$functional_group))



# All Years (1960-2018) #

# For loop to add up C extraction distributions per functional group for 1960-2018.
for (A in 1:nrow(C_extractions_perFG)) {
  
  C_extractions_perFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean C extraction per functional group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$C_extracted <- apply(C_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$C_extracted_SD <- apply(C_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perFG$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perFG$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perFG$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perFG$C_extracted_highCI <- quantiles_CIs_C[5, ]



# 1960-1964 Period #

# For loop to add up C extraction distributions per functional group for 1960-1964.
for (A in 1:nrow(C_extractions_perFG)) {
  
  C_extractions_perFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean C extraction per functional group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$C_extracted_1960_64 <- apply(C_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$C_extracted_SD_1960_64 <- apply(C_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$C_extracted_lowCI_1960_64 <- quantiles_CIs_C[1, ]
NutrientExtraction_perFG$C_extracted_lowIQR_1960_64 <- quantiles_CIs_C[2, ]
NutrientExtraction_perFG$C_extracted_median_1960_64 <- quantiles_CIs_C[3, ]
NutrientExtraction_perFG$C_extracted_highIQR_1960_64 <- quantiles_CIs_C[4, ]
NutrientExtraction_perFG$C_extracted_highCI_1960_64 <- quantiles_CIs_C[5, ]



# 1993-1997 Period #

# For loop to add up C extraction distributions per functional group for 1993-1997.
for (A in 1:nrow(C_extractions_perFG)) {
  
  C_extractions_perFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean C extraction per functional group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$C_extracted_1993_97 <- apply(C_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$C_extracted_SD_1993_97 <- apply(C_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$C_extracted_lowCI_1993_97 <- quantiles_CIs_C[1, ]
NutrientExtraction_perFG$C_extracted_lowIQR_1993_97 <- quantiles_CIs_C[2, ]
NutrientExtraction_perFG$C_extracted_median_1993_97 <- quantiles_CIs_C[3, ]
NutrientExtraction_perFG$C_extracted_highIQR_1993_97 <- quantiles_CIs_C[4, ]
NutrientExtraction_perFG$C_extracted_highCI_1993_97 <- quantiles_CIs_C[5, ]




# 2014-2018 Period #

# For loop to add up C extraction distributions per functional group for 2014-2018.
for (A in 1:nrow(C_extractions_perFG)) {
  
  C_extractions_perFG[A, 2:101] <- colSums(C_extractions[which(C_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean C extraction per functional group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$C_extracted_2014_18 <- apply(C_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$C_extracted_SD_2014_18 <- apply(C_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$C_extracted_lowCI_2014_18 <- quantiles_CIs_C[1, ]
NutrientExtraction_perFG$C_extracted_lowIQR_2014_18 <- quantiles_CIs_C[2, ]
NutrientExtraction_perFG$C_extracted_median_2014_18 <- quantiles_CIs_C[3, ]
NutrientExtraction_perFG$C_extracted_highIQR_2014_18 <- quantiles_CIs_C[4, ]
NutrientExtraction_perFG$C_extracted_highCI_2014_18 <- quantiles_CIs_C[5, ]




###### Nitrogen ######

# Creates data frame to hold functional group N extraction distributions.
N_extractions_perFG <- data.frame(functional_group_id = unique(Fisheries_NutrientExtraction$functional_group_id))


# All Years (1960-2018) #

# For loop to add up N extraction distributions per functional group for 1960-2018.
for (A in 1:nrow(N_extractions_perFG)) {
  
  N_extractions_perFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id &
                                                                 Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean N extraction per functional group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$N_extracted <- apply(N_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$N_extracted_SD <- apply(N_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perFG$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perFG$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perFG$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perFG$N_extracted_highCI <- quantiles_CIs_N[5, ]



# 1960-1964 Period #

# For loop to add up N extraction distributions per functional group for 1960-1964.
for (A in 1:nrow(N_extractions_perFG)) {
  
  N_extractions_perFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean N extraction per functional group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$N_extracted_1960_64 <- apply(N_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$N_extracted_SD_1960_64 <- apply(N_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$N_extracted_lowCI_1960_64 <- quantiles_CIs_N[1, ]
NutrientExtraction_perFG$N_extracted_lowIQR_1960_64 <- quantiles_CIs_N[2, ]
NutrientExtraction_perFG$N_extracted_median_1960_64 <- quantiles_CIs_N[3, ]
NutrientExtraction_perFG$N_extracted_highIQR_1960_64 <- quantiles_CIs_N[4, ]
NutrientExtraction_perFG$N_extracted_highCI_1960_64 <- quantiles_CIs_N[5, ]



# 1993-1997 Period #

# For loop to add up N extraction distributions per functional group for 1993-1997.
for (A in 1:nrow(N_extractions_perFG)) {
  
  N_extractions_perFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean N extraction per functional group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$N_extracted_1993_97 <- apply(N_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$N_extracted_SD_1993_97 <- apply(N_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$N_extracted_lowCI_1993_97 <- quantiles_CIs_N[1, ]
NutrientExtraction_perFG$N_extracted_lowIQR_1993_97 <- quantiles_CIs_N[2, ]
NutrientExtraction_perFG$N_extracted_median_1993_97 <- quantiles_CIs_N[3, ]
NutrientExtraction_perFG$N_extracted_highIQR_1993_97 <- quantiles_CIs_N[4, ]
NutrientExtraction_perFG$N_extracted_highCI_1993_97 <- quantiles_CIs_N[5, ]



# 2014-2018 Period #

# For loop to add up N extraction distributions per functional group for 2014-2018.
for (A in 1:nrow(N_extractions_perFG)) {
  
  N_extractions_perFG[A, 2:101] <- colSums(N_extractions[which(N_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean N extraction per functional group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$N_extracted_2014_18 <- apply(N_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$N_extracted_SD_2014_18 <- apply(N_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$N_extracted_lowCI_2014_18 <- quantiles_CIs_N[1, ]
NutrientExtraction_perFG$N_extracted_lowIQR_2014_18 <- quantiles_CIs_N[2, ]
NutrientExtraction_perFG$N_extracted_median_2014_18 <- quantiles_CIs_N[3, ]
NutrientExtraction_perFG$N_extracted_highIQR_2014_18 <- quantiles_CIs_N[4, ]
NutrientExtraction_perFG$N_extracted_highCI_2014_18 <- quantiles_CIs_N[5, ]




###### Phosphorus ######

# Creates data frame to hold functional group P extraction distributions.
P_extractions_perFG <- data.frame(functional_group_id = unique(Fisheries_NutrientExtraction$functional_group_id))


# All Years (1960-2018) #

# For loop to add up P extraction distributions per functional group for 1960-2018.
for (A in 1:nrow(P_extractions_perFG)) {
  
  P_extractions_perFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id &
                                                                 Fisheries_NutrientExtraction$year %in% c(1960:2018))])
  
}

# The following code generates the mean P extraction per functional group for 1960-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$P_extracted <- apply(P_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$P_extracted_SD <- apply(P_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perFG$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perFG$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perFG$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perFG$P_extracted_highCI <- quantiles_CIs_P[5, ]



# 1960-1964 Period #

# For loop to add up P extraction distributions per functional group for 1960-1964.
for (A in 1:nrow(P_extractions_perFG)) {
  
  P_extractions_perFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(1960:1964))])
  
}

# The following code generates the mean P extraction per functional group for 1960-1964. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$P_extracted_1960_64 <- apply(P_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$P_extracted_SD_1960_64 <- apply(P_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$P_extracted_lowCI_1960_64 <- quantiles_CIs_P[1, ]
NutrientExtraction_perFG$P_extracted_lowIQR_1960_64 <- quantiles_CIs_P[2, ]
NutrientExtraction_perFG$P_extracted_median_1960_64 <- quantiles_CIs_P[3, ]
NutrientExtraction_perFG$P_extracted_highIQR_1960_64 <- quantiles_CIs_P[4, ]
NutrientExtraction_perFG$P_extracted_highCI_1960_64 <- quantiles_CIs_P[5, ]



# 1993-1997 Period #

# For loop to add up P extraction distributions per functional group for 1993-1997.
for (A in 1:nrow(P_extractions_perFG)) {
  
  P_extractions_perFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

# The following code generates the mean P extraction per functional group for 1993-1997. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$P_extracted_1993_97 <- apply(P_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$P_extracted_SD_1993_97 <- apply(P_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$P_extracted_lowCI_1993_97 <- quantiles_CIs_P[1, ]
NutrientExtraction_perFG$P_extracted_lowIQR_1993_97 <- quantiles_CIs_P[2, ]
NutrientExtraction_perFG$P_extracted_median_1993_97 <- quantiles_CIs_P[3, ]
NutrientExtraction_perFG$P_extracted_highIQR_1993_97 <- quantiles_CIs_P[4, ]
NutrientExtraction_perFG$P_extracted_highCI_1993_97 <- quantiles_CIs_P[5, ]



# 2014-2018 Period #

# For loop to add up P extraction distributions per functional group for 2014-2018.
for (A in 1:nrow(P_extractions_perFG)) {
  
  P_extractions_perFG[A, 2:101] <- colSums(P_extractions[which(P_extractions_perFG$functional_group_id[A] == Fisheries_NutrientExtraction$functional_group_id & 
                                                                 Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

# The following code generates the mean P extraction per functional group for 2014-2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perFG$P_extracted_2014_18 <- apply(P_extractions_perFG[ , 2:101], 1, FUN = mean)
NutrientExtraction_perFG$P_extracted_SD_2014_18 <- apply(P_extractions_perFG[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_extractions_perFG[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perFG$P_extracted_lowCI_2014_18 <- quantiles_CIs_P[1, ]
NutrientExtraction_perFG$P_extracted_lowIQR_2014_18 <- quantiles_CIs_P[2, ]
NutrientExtraction_perFG$P_extracted_median_2014_18 <- quantiles_CIs_P[3, ]
NutrientExtraction_perFG$P_extracted_highIQR_2014_18 <- quantiles_CIs_P[4, ]
NutrientExtraction_perFG$P_extracted_highCI_2014_18 <- quantiles_CIs_P[5, ]



