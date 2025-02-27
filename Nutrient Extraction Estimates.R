###############################################################################
#                                                                             #
#                           Fishing Out Nutrients                             #
#                       Nutrient Extraction Estimates                         #
#                                                                             #
###############################################################################


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


# Then, we multiplied the matrix of the generated composition values values by the landed amount (tonnes) to generate a distribution of 100 estimates of nutrient extractions. We used the Sea Around Us data (https://www.seaaroundus.org/data/) at the EEZ and High Seas level meaning we download the CSV data files for each EEZ (n=283) and High Seas region (n=18) and compiled all of the data into one data frame (Fisheries_NutrientExtraction) that comprised 6.7+ million records of landed catches (reported in tonnes in the SAU data).

C_extractions <- SAU_catch.data$tonnes * C_distributions
N_extractions <- SAU_catch.data$tonnes * N_distributions
P_extractions <- SAU_catch.data$tonnes * P_distributions

# These extraction distribution matrices have been saved as CSVs which can be found in the following Figshare repository: 10.6084/m9.figshare.28500593.


# Finally, we took C_extracted to be the mean of the 100 extraction estimates produced above. We also obtained the SD, IQR, and 95% CIs. The 'Fisheries_NutrientExtraction' represents our finalized data frame which had all the estimates of C, N, and P extraction per row of the fisheries record (n=)





#### LOADING NUTRIENT EXTRACTION DISTRIBUTIONS ####

# It is recommended to load and work with each distribution matrix one at a time given their size. It is also recommended to use the data.table package to load it.

C_extractions <- fread("C_extraction_distributionmatrix.csv")

N_extractions <- fread("N_extraction_distributionmatrix.csv")

P_extractions <- fread("P_extraction_distributionmatrix.csv")





#### ESTIMATING MEAN NUTRIENT EXTRACTIONS ####

# The following code generates the mean nutrient extraction for each nutrient for each row in the compiled catch dataset. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.

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

##### PER YEAR #####

# Carbon
C_extractions_per_year <- data.frame(year = c(1960:2018))

# For loop to add up carbon extraction estimates and generate a distribution of 100 carbon extraction estimates per year.
for (A in 1:nrow(C_extractions_per_year)) {
  
  C_extractions_per_year[A, 2:101] <- colSums(C_extractions[which(C_extractions_per_year$year[A] == Fisheries_NutrientExtraction$year)])
  
}


# Generates a data frame to summarize nutrient extraction estimates per year.
NutrientExtraction_perYear <- data.frame(year = c(1960:2018))


# The following code generates the mean nutrient extraction for each nutrient for each year between 1960 and 2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perYear$C_extracted <- apply(C_extractions_per_year[ , 2:101], 1, FUN = mean)
NutrientExtraction_perYear$C_extracted_SD <- apply(C_extractions_per_year[ , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_year[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perYear$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perYear$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perYear$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perYear$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perYear$C_extracted_highCI <- quantiles_CIs_C[5, ]



# Nitrogen
N_extractions_per_year <- data.frame(year = c(1960:2018))

# For loop to add up nitrogen extraction estimates and generate a distribution of 100 nitrogen extraction estimates per year.
for (A in 1:nrow(N_extractions_per_year)) {
  
  N_extractions_per_year[A, 2:101] <- colSums(N_extractions[which(N_extractions_per_year$year[A] == Fisheries_NutrientExtraction$year)])
  
}


# The following code generates the mean nutrient extraction for each nutrient for each year between 1960 and 2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perYear$N_extracted <- apply(N_extractions_per_year[ , 2:101], 1, FUN = mean)
NutrientExtraction_perYear$N_extracted_SD <- apply(N_extractions_per_year[ , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_extractions_per_year[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perYear$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perYear$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perYear$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perYear$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perYear$N_extracted_highCI <- quantiles_CIs_N[5, ]



# Phosphorus
P_extractions_per_year <- data.frame(year = c(1960:2018))

# For loop to add up phosphorus extraction estimates and generate a distribution of 100 phosphorus extraction estimates per year.
for (A in 1:nrow(P_out_per_year)) {
  
  P_out_per_year[A, 2:101] <- colSums(P_extractions[which(P_out_per_year$year[A] == Fisheries_NutrientExtraction$year)])
  
}


# The following code generates the mean nutrient extraction for each nutrient for each year between 1960 and 2018. It also generates corresponding SDs, the upper and lower bounds of the IQR, and the upper and lower bounds of the 95% confidence intervals.
NutrientExtraction_perYear$P_extracted <- apply(P_out_per_year[ , 2:101], 1, FUN = mean)
NutrientExtraction_perYear$P_extracted_SD <- apply(P_out_per_year[ , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_out_per_year[ , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perYear$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perYear$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perYear$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perYear$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perYear$P_extracted_highCI <- quantiles_CIs_P[5, ]
