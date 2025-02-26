#### ESTIMATING NUTRIENT EXTRACTIONS ####

# To account for uncertainties around producing nutrient extraction estimates, we generated 100 nutrient values drawn from a random normal distribution of each taxa's mean nutrient compositions for each row of the SAU dataset. The SAU data comprised 6.7+ million rows of data where each row of data had a distribution of 100 nutrient composition values for the row's corresponding taxa.

# We wanted to make sure we could generate distributions that were between 0 and 100. Since some of our values are very close to 0 (especially in the case of phosphorus), we wanted to ensure that no negative values were accidentally generated during our random draws. To accomplish this, we transformed our nutrient compositions into alpha and beta values so we could generate beta distributions. The transformations we made are recreated below.

taxa_nutrient.data <- taxa_nutrient.data %>% 
  mutate(alph_C = (mean_C/100)*(((mean_C/100)*(1- (mean_C/100)))/((SD_C/100)^2)-1),
         bet_C = (1-(mean_C/100))*(((mean_C/100)*(1- (mean_C/100)))/((SD_C/100)^2)-1),
         alph_N = (mean_N/100)*(((mean_N/100)*(1- (mean_N/100)))/((SD_N/100)^2)-1),
         bet_N = (1-(mean_N/100))*(((mean_N/100)*(1- (mean_N/100)))/((SD_N/100)^2)-1),
         alph_P = (mean_P/100)*(((mean_P/100)*(1- (mean_P/100)))/((SD_P/100)^2)-1),
         bet_P = (1-(mean_P/100))*(((mean_P/100)*(1- (mean_P/100)))/((SD_P/100)^2)-1))


# The following code created a vector of integers to identify which rows of the nutrient composition data frame had to be called when calculating nutrient extraction with the fisheries data frame. TaxonKey is a unique ID code for each taxa. These ID codes were provided by the SAU and used to match taxa between our nutrient composition data frame and the catch data frame.
taxaorder <- match(taxa_nutrient.data$TaxonKey, SAU_catch.data$TaxonKey)


# Then, using the newly transformed values, we generated 100 nutrient composition values taken from beta distributions. Here, we generated a matrix where each row of the fisheries dataset had a distribution of 100 nutrient values based on the corresponding taxa for each row. Each of the 100 generate values were then multiplied by the landed amount reported in the fisheries dataset. The code used to generate these matrices is reproduced below.

C_distributions <- matrix(
  rbeta(length(taxaorder) * 100, taxa_nutrient.data$alph_C[taxaorder], 
        taxa_nutrient.data$bet_C[taxaorder]), 
  nrow = length(taxaorder), ncol = 100
)

N_distributions <- matrix(
  rbeta(length(taxaorder) * 100, taxa_nutrient.data$alph_N[taxaorder], 
        taxa_nutrient.data$bet_N[taxaorder]), 
  nrow = length(taxaorder), ncol = 100
)

P_distributions <- matrix(
  rbeta(length(taxaorder) * 100, taxa_nutrient.data$alph_P[taxaorder], 
        taxa_nutrient.data$bet_P[taxaorder]), 
  nrow = length(taxaorder), ncol = 100
)

# Then, we multiplied the matrix of the generated values by the landed amount (tonnes) to generate a distribution of 100 estimates of nutrient extractions. We used the Sea Around Us data (https://www.seaaroundus.org/data/) at the EEZ and High Seas level meaning we download the CSV data files for each EEZ (n=283) and High Seas region (n=18) and compiled all of the data into one data frame that comprised 6.7+ million records of landed catches (reported in tonnes in the SAU data).

C_extractions <- SAU_catch.data$tonnes * C_distributions
N_extractions <- SAU_catch.data$tonnes * N_distributions
P_extractions <- SAU_catch.data$tonnes * P_distributions


# Finally, we took C_extracted to be the mean of the 100 extraction estimates produced above. We also obtained the SD, IQR, and 95% CIs. The 'Fisheries_NutrientExtraction' represents our finalized data frame which had all the estimates of C, N, and P extraction per row of the fisheries record (n=)

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



#### ESTIMATING NUTRIENT EXTRACTIONS BY CATEGORIES ####

##### YEAR #####
Nutrient_extractions_perYear <- data.frame(year = 1950:2018)

C_extractions_perYear <- C_extractions %>%
  as.data.frame() %>%
  bind_cols(year = Fisheries_NutrientExtraction$year) %>%
  group_by(year) %>%
  summarise(across(starts_with("V"), sum)) # the 'starts_with("V)' piece is to summarize across the columns that have the simulated values from the matrix.

N_extractions_perYear <- N_extractions %>%
  as.data.frame() %>%
  bind_cols(year = Fisheries_NutrientExtraction$year) %>%
  group_by(year) %>%
  summarise(across(starts_with("V"), sum))

P_extractions_perYear <- P_extractions %>%
  as.data.frame() %>%
  bind_cols(year = Fisheries_NutrientExtraction$year) %>%
  group_by(year) %>%
  summarise(across(starts_with("V"), sum)) 

Nutrient_extractions_perYear <- C_extractions_perYear %>%
  rowwise() %>%
  mutate(
    C_extracted = mean(c_across(starts_with("V")), na.rm = TRUE),
    C_extracted_SD = sd(c_across(starts_with("V")), na.rm = TRUE),
    quantiles_CIs_C = list(quantile(c_across(starts_with("V")), probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))
  ) %>%
  unnest_wider(quantiles_CIs_C, names_sep = "_") %>%
  rename(
    C_extracted_lowCI = quantiles_CIs_C_1,
    C_extracted_lowIQR = quantiles_CIs_C_2,
    C_extracted_median = quantiles_CIs_C_3,
    C_extracted_highIQR = quantiles_CIs_C_4,
    C_extracted_highCI = quantiles_CIs_C_5
  ) %>%
  select(year, everything())

N_extractions_perYear <- N_extractions_perYear %>%
  rowwise() %>%
  mutate(
    N_extracted = mean(c_across(starts_with("V")), na.rm = TRUE),
    N_extracted_SD = sd(c_across(starts_with("V")), na.rm = TRUE),
    quantiles_CIs_N = list(quantile(c_across(starts_with("V")), probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))
  ) %>%
  unnest_wider(quantiles_CIs_N, names_sep = "_") %>%
  rename(
    N_extracted_lowCI = quantiles_CIs_N_1,
    N_extracted_lowIQR = quantiles_CIs_N_2,
    N_extracted_median = quantiles_CIs_N_3,
    N_extracted_highIQR = quantiles_CIs_N_4,
    N_extracted_highCI = quantiles_CIs_N_5
  )

P_extractions_perYear <- P_extractions_perYear %>%
  rowwise() %>%
  mutate(
    P_extracted = mean(c_across(starts_with("V")), na.rm = TRUE),
    P_extracted_SD = sd(c_across(starts_with("V")), na.rm = TRUE),
    quantiles_CIs_P = list(quantile(c_across(starts_with("V")), probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))
  ) %>%
  unnest_wider(quantiles_CIs_P, names_sep = "_") %>%
  rename(
    P_extracted_lowCI = quantiles_CIs_P_1,
    P_extracted_lowIQR = quantiles_CIs_P_2,
    P_extracted_median = quantiles_CIs_P_3,
    P_extracted_highIQR = quantiles_CIs_P_4,
    P_extracted_highCI = quantiles_CIs_P_5
  )

NutrientExtraction_perYear <- Nutrient_extractions_perYear %>%
  full_join(N_extractions_perYear, by = "year") %>%
  full_join(P_extractions_perYear, by = "year") %>%
  select(year, starts_with("C_"), starts_with("N_"), starts_with("P_"))





##### TIME PERIODS #####

C_extractions_per_timeperiod <- data.frame(
  time_period = c("1950_54", "1993_97", "2014_18"),
  stringsAsFactors = FALSE
)

N_extractions_per_timeperiod <- data.frame(
  time_period = c("1950_54", "1993_97", "2014_18"),
  stringsAsFactors = FALSE
)

P_extractions_per_timeperiod <- data.frame(
  time_period = c("1950_54", "1993_97", "2014_18"),
  stringsAsFactors = FALSE
)



###### 1950-1954 ######
# Carbon
for (A in 1:nrow(C_extractions_per_timeperiod)) {
  
  C_extractions_per_timeperiod[1, 2:101] <- colSums(C_extractions[which(Fisheries_NutrientExtraction$year %in% c(1950:1954))])
  
}

NutrientExtraction_perTimePeriod <- as.data.frame(c("1950_54", "1993_97", "2014_18"))
colnames(NutrientExtraction_perTimePeriod)[1] <- "time_period"
NutrientExtraction_perTimePeriod$C_extracted <- NA
NutrientExtraction_perTimePeriod$C_extracted_SD <- NA
NutrientExtraction_perTimePeriod$C_extracted_lowCI <- NA
NutrientExtraction_perTimePeriod$C_extracted_lowIQR <- NA
NutrientExtraction_perTimePeriod$C_extracted_median <- NA
NutrientExtraction_perTimePeriod$C_extracted_highIQR <- NA
NutrientExtraction_perTimePeriod$C_extracted_highCI <- NA


NutrientExtraction_perTimePeriod[1, ]$C_extracted <- apply(C_extractions_per_timeperiod[1 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[1, ]$C_extracted_SD <- apply(C_extractions_per_timeperiod[1 , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_timeperiod[1 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[1, ]$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perTimePeriod[1, ]$C_extracted_highCI <- quantiles_CIs_C[5, ]

# cleanup
remove(quantiles_CIs_C)

# Nitrogen
for (A in 1:nrow(N_out_per_timeperiod)) {
  
  N_out_per_timeperiod[1, 2:101] <- colSums(N_out[which(Fisheries_NutrientExtraction$year %in% c(1950:1954))])
  
}

NutrientExtraction_perTimePeriod$N_extracted <- NA
NutrientExtraction_perTimePeriod$N_extracted_SD <- NA
NutrientExtraction_perTimePeriod$N_extracted_lowCI <- NA
NutrientExtraction_perTimePeriod$N_extracted_lowIQR <- NA
NutrientExtraction_perTimePeriod$N_extracted_median <- NA
NutrientExtraction_perTimePeriod$N_extracted_highIQR <- NA
NutrientExtraction_perTimePeriod$N_extracted_highCI <- NA


NutrientExtraction_perTimePeriod[1, ]$N_extracted <- apply(N_out_per_timeperiod[1 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[1, ]$N_extracted_SD <- apply(N_out_per_timeperiod[1 , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_out_per_timeperiod[1 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[1, ]$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perTimePeriod[1, ]$N_extracted_highCI <- quantiles_CIs_N[5, ]

remove(quantiles_CIs_N)


# Phosphorus
for (A in 1:nrow(P_out_per_timeperiod)) {
  
  P_out_per_timeperiod[1, 2:101] <- colSums(P_out[which(Fisheries_NutrientExtraction$year %in% c(1950:1954))])
  
}

NutrientExtraction_perTimePeriod$P_extracted <- NA
NutrientExtraction_perTimePeriod$P_extracted_SD <- NA
NutrientExtraction_perTimePeriod$P_extracted_lowCI <- NA
NutrientExtraction_perTimePeriod$P_extracted_lowIQR <- NA
NutrientExtraction_perTimePeriod$P_extracted_median <- NA
NutrientExtraction_perTimePeriod$P_extracted_highIQR <- NA
NutrientExtraction_perTimePeriod$P_extracted_highCI <- NA


NutrientExtraction_perTimePeriod[1, ]$P_extracted <- apply(P_out_per_timeperiod[1 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[1, ]$P_extracted_SD <- apply(P_out_per_timeperiod[1 , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_out_per_timeperiod[1 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[1, ]$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perTimePeriod[1, ]$P_extracted_highCI <- quantiles_CIs_P[5, ]

remove(quantiles_CIs_P)




###### 1993-1997 ######

# Carbon
for (A in 1:nrow(C_extractions_per_timeperiod)) {
  
  C_extractions_per_timeperiod[2, 2:101] <- colSums(C_extractions[which(Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

NutrientExtraction_perTimePeriod[2, ]$C_extracted <- apply(C_extractions_per_timeperiod[2 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[2, ]$C_extracted_SD <- apply(C_extractions_per_timeperiod[2 , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_timeperiod[2 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[2, ]$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perTimePeriod[2, ]$C_extracted_highCI <- quantiles_CIs_C[5, ]

remove(quantiles_CIs_C)

# Nitrogen
for (A in 1:nrow(N_out_per_timeperiod)) {
  
  N_out_per_timeperiod[2, 2:101] <- colSums(N_out[which(Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

NutrientExtraction_perTimePeriod[2, ]$N_extracted <- apply(N_out_per_timeperiod[2 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[2, ]$N_extracted_SD <- apply(N_out_per_timeperiod[2 , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_out_per_timeperiod[2 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[2, ]$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perTimePeriod[2, ]$N_extracted_highCI <- quantiles_CIs_N[5, ]

remove(quantiles_CIs_N)

# Phosphorus
for (A in 1:nrow(P_out_per_timeperiod)) {
  
  P_out_per_timeperiod[2, 2:101] <- colSums(P_out[which(Fisheries_NutrientExtraction$year %in% c(1993:1997))])
  
}

NutrientExtraction_perTimePeriod[2, ]$P_extracted <- apply(P_out_per_timeperiod[2 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[2, ]$P_extracted_SD <- apply(P_out_per_timeperiod[2 , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_out_per_timeperiod[2 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[2, ]$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perTimePeriod[2, ]$P_extracted_highCI <- quantiles_CIs_P[5, ]

remove(quantiles_CIs_P)



###### 2014-2018 ######

# Carbon
for (A in 1:nrow(C_extractions_per_timeperiod)) {
  
  C_extractions_per_timeperiod[3, 2:101] <- colSums(C_extractions[which(Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

NutrientExtraction_perTimePeriod[3, ]$C_extracted <- apply(C_extractions_per_timeperiod[3 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[3, ]$C_extracted_SD <- apply(C_extractions_per_timeperiod[3 , 2:101], 1, FUN = sd)
quantiles_CIs_C <- apply(C_extractions_per_timeperiod[3 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[3, ]$C_extracted_lowCI <- quantiles_CIs_C[1, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_lowIQR <- quantiles_CIs_C[2, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_median <- quantiles_CIs_C[3, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_highIQR <- quantiles_CIs_C[4, ]
NutrientExtraction_perTimePeriod[3, ]$C_extracted_highCI <- quantiles_CIs_C[5, ]

remove(quantiles_CIs_C)


# Nitrogen
for (A in 1:nrow(N_out_per_timeperiod)) {
  
  N_out_per_timeperiod[3, 2:101] <- colSums(N_out[which(Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

NutrientExtraction_perTimePeriod[3, ]$N_extracted <- apply(N_out_per_timeperiod[3 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[3, ]$N_extracted_SD <- apply(N_out_per_timeperiod[3 , 2:101], 1, FUN = sd)
quantiles_CIs_N <- apply(N_out_per_timeperiod[3 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[3, ]$N_extracted_lowCI <- quantiles_CIs_N[1, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_lowIQR <- quantiles_CIs_N[2, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_median <- quantiles_CIs_N[3, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_highIQR <- quantiles_CIs_N[4, ]
NutrientExtraction_perTimePeriod[3, ]$N_extracted_highCI <- quantiles_CIs_N[5, ]

remove(quantiles_CIs_N)


# Phosphorus
for (A in 1:nrow(P_out_per_timeperiod)) {
  
  P_out_per_timeperiod[3, 2:101] <- colSums(P_out[which(Fisheries_NutrientExtraction$year %in% c(2014:2018))])
  
}

NutrientExtraction_perTimePeriod[3, ]$P_extracted <- apply(P_out_per_timeperiod[3 , 2:101], 1, FUN = mean)
NutrientExtraction_perTimePeriod[3, ]$P_extracted_SD <- apply(P_out_per_timeperiod[3 , 2:101], 1, FUN = sd)
quantiles_CIs_P <- apply(P_out_per_timeperiod[3 , 2:101], 1, FUN = quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
NutrientExtraction_perTimePeriod[3, ]$P_extracted_lowCI <- quantiles_CIs_P[1, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_lowIQR <- quantiles_CIs_P[2, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_median <- quantiles_CIs_P[3, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_highIQR <- quantiles_CIs_P[4, ]
NutrientExtraction_perTimePeriod[3, ]$P_extracted_highCI <- quantiles_CIs_P[5, ]

remove(quantiles_CIs_P)




##### FUNCTIONAL GROUPS #####

