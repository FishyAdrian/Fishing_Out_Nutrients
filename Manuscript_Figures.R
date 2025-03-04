###############################################################################
#                                                                             #
#                            Manuscript Figures                               #
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
library(patchwork)
library(xtable)
library(cowplot)
library(grid)
library(gridExtra)




#### FIGURE 1 ####

# Defines color palette for the three nutrients (C, N, and P).
cols <- c("#56B4E9", "#009E73", "#D55E00")



##### Extraction Plots #####

# NOTE: We cannot reproduce the Sea Around Us landings data so the landings data used in Figure 1a is not reproduced in the repository, but the code is provided below.


## Panel a - Landings per year

# Creates data frame for annual landings for industrial fisheries.
annualLandings <- Fisheries_NutrientExtraction %>% group_by(year) %>% 
  summarize(annual_landings = sum(tonnes))

annualLandings_plot <- ggplot(annualLandings, aes(x = year, y = annual_landings)) +
  geom_line(lwd = 0.8) +
  scale_x_continuous(name = NULL, limits = c(1960,2018),
                     breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(0, 17e+06, 34e+06, 51e+06, 68e+06, 85e+06),
                     labels = c("0", "17", "34", "51", "68", "85"), limits = c(0, 87.937e+06)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "a  Landings") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_line(),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))



# Panel b - Carbon extractions per year
annualCE_plot <- ggplot(NutrientExtraction_perYear, aes(x = year, y = C_extracted)) +
  geom_ribbon(aes(ymin = C_extracted_lowCI, ymax = C_extracted_highCI), alpha = 0.6, fill = "#56B4E9") +
  geom_line(lwd = 0.8) +
  scale_x_continuous(name = NULL, breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(0, 2e+06, 4e+06, 6e+06, 8e+06, 10e+06), 
                     labels = c("0", "2", "4", "6", "8", "10"), limits = c(0, 10e+06)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "b  Carbon") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_line(),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))



# Panel c - Nitrogen extractions per year
annualNE_plot <- ggplot(NutrientExtraction_perYear, aes(x = year, y = N_extracted)) +
  geom_ribbon(aes(ymin = N_extracted_lowCI, ymax = N_extracted_highCI), alpha = 0.6, fill = "#009E73") +
  geom_line(lwd = 0.8) +
  scale_x_continuous(name = NULL, breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(0, 0.5e+06, 1e+06, 1.5e+06, 2e+06, 2.5e+06), 
                     labels = c("0", "0.5", "1", "1.5", "2", "2.5"), limits = c(0, 2.53e+06)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "c  Nitrogen") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_line(),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))



# Panel d - Phosphorus extractions per year
annualPE_plot <- ggplot(NutrientExtraction_perYear, aes(x = year, y = P_extracted)) +
  geom_ribbon(aes(ymin = P_extracted_lowCI, ymax = P_extracted_highCI), alpha = 0.6, fill = "#D55E00") +
  geom_line(lwd = 0.8) +
  scale_x_continuous(name = NULL, breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(0, 1.2e+05, 2.4e+05, 3.6e+05, 4.8e+05, 6e+05), 
                     labels = c("0", "0.12", "0.24", "0.36", "0.48", "0.6"), limits = c(0, 6.25e+05)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "d  Phosphorus") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_line(),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))





##### Nutrient Ratio Plots #####

# Panel e - C:N extraction per year
annualC.N_plot <- ggplot(NutrientRatios_perYear, aes(x = year, y = C.N_mean)) +
  geom_ribbon(aes(ymin = C.N_lowCI, ymax = C.N_highCI), alpha = 0.6, fill = "#2ca9af") +
  geom_line(lwd = 0.8) +
  geom_line(aes(x = c(1960:2018), y = 5.208333), linetype = "longdash", color = "black", lwd = 0.6) +
  scale_x_continuous(name = NULL, breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(seq(4, 5.25, by = 0.25)),
                     limits = c(3.95, 5.3)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "e  C:N") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_line(),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))


# Panel f - C:P extraction per year
annualC.P_plot <- ggplot(NutrientRatios_perYear, aes(x = year, y = C.P_mean)) +
  geom_ribbon(aes(ymin = C.P_lowCI, ymax = C.P_highCI), alpha = 0.6, fill = "#978975") +
  geom_line(lwd = 0.8) +
  geom_line(aes(x = c(1960:2018), y = 53.81944), linetype = "longdash", color = "black", lwd = 0.6) +
  scale_x_continuous(name = NULL, breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(seq(35, 60, by = 5)),
                     limits = c(33, 62)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "f  C:P") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))


# Panel g - N:P extraction per year
annualN.P_plot <- ggplot(NutrientRatios_perYear, aes(x = year, y = N.P_mean)) +
  geom_ribbon(aes(ymin = N.P_lowCI, ymax = N.P_highCI), alpha = 0.6, fill = "#6b7e3a") +
  geom_line(lwd = 0.8) +
  geom_line(aes(x = c(1960:2018), y = 10.33333), linetype = "longdash", color = "black", lwd = 0.6) +
  scale_x_continuous(name = NULL, breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(seq(7.5, 13.5, by = 1.5)),
                     limits = c(7, 13.5)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "g  N:P") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))





##### Exporting Figure 1 #####
jpeg("Figure1_manuscript.jpeg", width = 18, height = 19.38, units = 'cm', res = 600)

# Sets up blank
blank <- ggplot() + geom_blank() + theme_void()

# Sets up left-side panels
Figure1.w.ratios_1 <- ggarrange(annualLandings_plot,
                                annualNE_plot,
                                annualC.N_plot,
                                annualN.P_plot,
                                ncol = 1, nrow = 4,
                                align = c("v"),
                                heights = c(1, 1, 1, 1.16))
Figure1.w.ratios_1 <- annotate_figure(Figure1.w.ratios_1,
                                      left = text_grob("Molar ratio                                                      Million tonnes", rot = 90, size = 12))

# Sets up right-side panels
Figure1.w.ratios_2 <- ggarrange(annualCE_plot,
                                annualPE_plot,
                                annualC.P_plot,
                                blank,
                                ncol = 1, nrow = 4,
                                align = c("v"),
                                heights = c(1, 1, 1.16, 1))

# Combines two sides
Figure1.w.ratios <- ggarrange(Figure1.w.ratios_1, Figure1.w.ratios_2,
                              ncol = 2,
                              align = "v")
# Labels the x-axis
annotate_figure(Figure1.w.ratios,
                bottom = text_grob("Year", size = 12))

dev.off()




#### FIGURE 2 ####

# Maps for Figure 2 were made in ArcGIS Pro 3.1.2 using the data included in the data frames listed below. Data frames were exported as CSVs and imported into ArcGIS Pro to create the figures.

# For the right-side panels, the following data frame was used:
NutrientExtraction_perArea

# For the right-hand panels which demonstrate the percent difference between our estimates and those produced using the uniform values, the following data frame was used:
GloAv_NutrientExtraction_perArea





#### FIGURE 3 ####

# Maps for Figure 3 were made in ArcGIS Pro 3.1.2 using the data included in the data frames listed below. Data frames were exported as CSVs and imported into ArcGIS Pro to create the figures.

# All panels for this figure all came from data included in the following data frame:
NutrientRatios_perArea





#### FIGURE 4 ####

# Maps for Figure 4 were made in ArcGIS Pro 3.1.2 using the data included in the data frame listed below. Data frame was exported as a CSV and imported into ArcGIS Pro to create the figures.

# All panels for this figure all came from data included in the following data frame:
NutrientExtraction_perArea






#### FIGURE 5 ####

# Sets the color legend for the different periods in the following order: 1960-1964, 1993-1997, 2014-2018, and 1960-2018.
cols_years <- c("#332288", "#DDCC77", "#AA4499", "#117733")


# Transforms the NutrientExtraction_perTG data frame from wide to long formats. The first data frame contains the estimates and the second contains the standard deviations.
NutrientExtraction_perTG_long1 <- as.data.table(NutrientExtraction_perTG %>% select(c(1, seq(2, 79, by = 7))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "metric_tons"))
NutrientExtraction_perTG_long2 <- as.data.table(NutrientExtraction_perTG %>% select (c(1, seq(3, 80, by = 7))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "metric_tons_SD"))

# Merges the two long data frames to create the data frame needed for the plot.
NE_perTG_forPlot <- cbind(NutrientExtraction_perTG_long1, NutrientExtraction_perTG_long2)
NE_perTG_forPlot <- NE_perTG_forPlot %>% select(-c(4:5))


# Creates the column to designate the nutrient in each row.
NE_perTG_forPlot$nutrient <- "NA"
NE_perTG_forPlot[years %in% c("C_extracted", "C_extracted_1960_64", "C_extracted_1993_97", "C_extracted_2014_18")]$nutrient <- "Carbon"
NE_perTG_forPlot[years %in% c("N_extracted", "N_extracted_1960_64", "N_extracted_1993_97", "N_extracted_2014_18")]$nutrient <- "Nitrogen"
NE_perTG_forPlot[years %in% c("P_extracted", "P_extracted_1960_64", "P_extracted_1993_97", "P_extracted_2014_18")]$nutrient <- "Phosphorous"

# Changes the values to correspond to the correct time periods.
NE_perTG_forPlot[years %in% c("C_extracted", "N_extracted", "P_extracted")]$years <- "All Years"
NE_perTG_forPlot[years %in% c("C_extracted_1960_64", "N_extracted_1960_64", "P_extracted_1960_64")]$years <- "1960-64"
NE_perTG_forPlot[years %in% c("C_extracted_1993_97", "N_extracted_1993_97", "P_extracted_1993_97")]$years <- "1993-97"
NE_perTG_forPlot[years %in% c("C_extracted_2014_18", "N_extracted_2014_18", "P_extracted_2014_18")]$years <- "2014-18"


##### Panels a-b #####

# Panel a was created as a diagram in PowerPoint.

# Panel b was created as a simple pie chart using Excel with the data from data frame NutrientExtraction_perTG.


##### Extraction Plots #####

## Panel c - Carbon extractions
CE_perTG_plot <- NE_perTG_forPlot %>% filter(years != "All Years" & nutrient == "Carbon") %>% 
  ggplot(aes(x = trophic_group, y=metric_tons, 
             fill=factor(years, levels= rev(c("2014-18", "1993-97", "1960-64", "All Years"))))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators"),
                   labels = NULL) +
  scale_y_continuous(name = "Million tonnes", expand = c(0.02, 0.02),
                     limits = c(0, 28000000),
                     breaks = c(0, 5.5e+6, 11e+6, 16.5e+6, 22e+6, 27.5e+6),
                     labels = c("0", "5.5", "11.0", "16.5", "22.0", "27.5")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "c") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title.y = element_text(size=12),
        legend.position = c(.225, .7),
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

## Panel d - Nitrogen extractions
NE_perTG_plot <- NE_perTG_forPlot %>% filter(years != "All Years" & nutrient == "Nitrogen") %>% ggplot(aes(x = trophic_group, y=metric_tons,
                                                                                                           fill=factor(years, levels= rev(c("2014-18", "1993-97", "1960-64", "All Years"))))) +
  geom_bar(color = "black",stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators"),
                   labels = NULL) +
  scale_y_continuous(name = "Million tonnes", expand = c(0.02, 0.02),
                     limits = c(0, 7030826),
                     breaks = c(0, 1.4e+6, 2.8e+6, 4.2e+6, 5.6e+6, 7e+6),
                     labels = c("0", "1.4", "2.8", "4.2", "5.6", "7.0")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "d") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")

## Panel e - Phosphorus extractions
PE_perTG_plot <- NE_perTG_forPlot %>% filter(years != "All Years" & nutrient == "Phosphorous") %>% ggplot(aes(x = trophic_group, y=metric_tons, 
                                                                                                              fill=factor(years, levels= rev(c("2014-18", "1993-97", "1960-64", "All Years"))))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  scale_y_continuous(name = "Million tonnes", expand = c(0.02, 0.02),
                     limits = c(0, 1553681),
                     breaks = c(0, 0.3e+6, 0.6e+6, 0.9e+6, 1.2e+6, 1.5e+6),
                     labels = c("0", "0.3", "0.6", "0.9", "1.2", "1.5")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "e") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")

##### Comparison Plot #####

test1 <- as.data.table(GloAv_NutrientExtraction_perTG %>% select(c(1, seq(102, 190, by = 8))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "per_change"))
test2 <- as.data.table(GloAv_NutrientExtraction_perTG %>% select(c(1, seq(103, 191, by = 8))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "per_change_SD"))

GloAv_NE_perTG_forPlot <- cbind(test1, test2)
GloAv_NE_perTG_forPlot <- GloAv_NE_perTG_forPlot %>% select(-c(4:5))


GloAv_NE_perTG_forPlot$nutrient <- "NA"
GloAv_NE_perTG_forPlot[years %in% c("C_dif_per", "C_dif_1960_64_per", "C_dif_1993_97_per", "C_dif_2014_18_per")]$nutrient <- "Carbon"
GloAv_NE_perTG_forPlot[years %in% c("N_dif_per", "N_dif_1960_64_per", "N_dif_1993_97_per", "N_dif_2014_18_per")]$nutrient <- "Nitrogen"
GloAv_NE_perTG_forPlot[years %in% c("P_dif_per", "P_dif_1960_64_per", "P_dif_1993_97_per", "P_dif_2014_18_per")]$nutrient <- "Phosphorus"

GloAv_NE_perTG_forPlot[years %in% c("C_dif_per", "N_dif_per", "P_dif_per")]$years <- "All Years"
GloAv_NE_perTG_forPlot[years %in% c("C_dif_1960_64_per", "N_dif_1960_64_per", "P_dif_1960_64_per")]$years <- "1960-64"
GloAv_NE_perTG_forPlot[years %in% c("C_dif_1993_97_per", "N_dif_1993_97_per", "P_dif_1993_97_per")]$years <- "1993-97"
GloAv_NE_perTG_forPlot[years %in% c("C_dif_2014_18_per", "N_dif_2014_18_per", "P_dif_2014_18_per")]$years <- "2014-18"


## Panel f
perTG_per_dif_plot <- GloAv_NE_perTG_forPlot %>% filter(years == "All Years") %>% 
  ggplot(aes(x = trophic_group, y=per_change, 
             fill=factor(nutrient, levels= c("Carbon", "Nitrogen", "Phosphorus")))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  scale_y_continuous(name = "Percent difference", expand = c(0.02, 0.02),
                     limits = c(-10, 30),
                     breaks = c(-10, 0, 10, 20, 30)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(title = "Nutrient")) +
  labs(title = "f") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=60, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        legend.key.size = unit(4, 'mm'),
        axis.title =element_text(size=12),
        legend.position = c(.75, .85),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))





##### Exporting Figure 4 #####

jpeg("Figure4_manuscript.jpeg", width = 6.5, height = 5.75, units = 'in', res = 600)
Figure4.4panel <- ggarrange(CE_perTG_plot, NE_perTG_plot, PE_perTG_plot, perTG_per_dif_plot,
                            heights = c(0.55, 1))
annotate_figure(Figure4.4panel,
                bottom = text_grob("Trophic group", size = 12))
dev.off()







#### FIGURE 6 ####

# Maps for Figure 6 were made in ArcGIS Pro 3.1.2. The code below outlines how the proportions were calculated. Once calculated, the data frames below were exported as CSVs and imported into ArcGIS Pro to create the panels.

##### Trophic Group Proportions #####

# Proportion of Carbon Extraction through Mesopredators per Area, All Years
NutrientExtraction_perArea.mesopredators <- Fisheries_NutrientExtraction %>% 
  filter(trophic_group == "Mesopredators" & year %in% c(1960:2018)) %>% 
  group_by(area_name) %>% 
  summarise(C_extracted_mesopredators = sum(C_extracted))

NutrientExtraction_perArea.TG.prop <- NutrientExtraction_perArea[ , c(1,2)]
NutrientExtraction_perArea.TG.prop <- merge(NutrientExtraction_perArea.TG.prop, NutrientExtraction_perArea.mesopredators, 
                                            by = "area_name", 
                                            all.x = TRUE)
NutrientExtraction_perArea.TG.prop <- NutrientExtraction_perArea.TG.prop %>% 
  mutate(CE_meso_prop = C_extracted_mesopredators/C_extracted)



# Proportion of Carbon Extraction through High-Level Predators (trophic level > 4) per Area, All Years
NutrientExtraction_perArea.highpred <- Fisheries_NutrientExtraction %>% 
  filter((trophic_group == "High-level predators" | trophic_group == "Top predators") & year %in% c(1960:2018)) %>% 
  group_by(area_name) %>% 
  summarise(C_extracted_highpred = sum(C_extracted))

NutrientExtraction_perArea.TG.prop <- merge(NutrientExtraction_perArea.TG.prop, NutrientExtraction_perArea.highpred,
                                            by = "area_name",
                                            all.x = TRUE)

NutrientExtraction_perArea.TG.prop <- NutrientExtraction_perArea.TG.prop %>% 
  mutate(CE_highpred_prop = C_extracted_highpred/C_extracted)



# Proportion of Carbon Extraction through Low-level consumers (trophic level = 2-2.8) per Area, All Years
NutrientExtraction_perArea.low.level.con <- Fisheries_NutrientExtraction %>% 
  filter(trophic_group == "Low-level consumers" & year %in% c(1960:2018)) %>% 
  group_by(area_name) %>% 
  summarise(C_extracted_low.level.con = sum(C_extracted))

NutrientExtraction_perArea.TG.prop <- merge(NutrientExtraction_perArea.TG.prop, NutrientExtraction_perArea.low.level.con, 
                                            by = "area_name",
                                            all.x = TRUE)
NutrientExtraction_perArea.TG.prop <- NutrientExtraction_perArea.TG.prop %>% 
  mutate(CE_low.level_prop = C_extracted_low.level.con/C_extracted)



# Verify if proportions add up

NutrientExtraction_perArea.TG.prop %>% mutate(Total_prop = CE_meso_prop + CE_highpred_prop + CE_low.level_prop) %>% filter(Total_prop > 1)



##### Functional Group Proportions #####

# Proportion of Carbon Extraction through Pelagics per Area, All Years
NutrientExtraction_perArea.pelagics <- Fisheries_NutrientExtraction %>% 
  filter(simp_functional_group == "Pelagic" & year %in% c(1960:2018)) %>% 
  group_by(area_name) %>% 
  summarise(C_extracted_pelagic = sum(C_extracted))

NutrientExtraction_perArea.FG.prop <- NutrientExtraction_perArea[ , c(1,2)]
NutrientExtraction_perArea.FG.prop <- merge(NutrientExtraction_perArea.FG.prop, NutrientExtraction_perArea.pelagics, 
                                            by = "area_name", 
                                            all.x = TRUE)
NutrientExtraction_perArea.FG.prop <- NutrientExtraction_perArea.FG.prop %>% 
  mutate(CE_pel_prop = C_extracted_pelagic/C_extracted)



# Proportion of Carbon Extraction through Demersals per Area, All Years
NutrientExtraction_perArea.demersals <- Fisheries_NutrientExtraction %>% 
  filter(simp_functional_group == "Demersal" & year %in% c(1960:2018)) %>% 
  group_by(area_name) %>% 
  summarise(C_extracted_demersal = sum(C_extracted))

NutrientExtraction_perArea.FG.prop <- merge(NutrientExtraction_perArea.FG.prop, NutrientExtraction_perArea.demersals, 
                                            by = "area_name",
                                            all.x = TRUE)
NutrientExtraction_perArea.FG.prop <- NutrientExtraction_perArea.FG.prop %>% 
  mutate(CE_dem_prop = C_extracted_demersal/C_extracted)



# Proportion of Carbon Extraction through Benthopelagics per Area, All Years
NutrientExtraction_perArea.benthopelagics <- Fisheries_NutrientExtraction %>% 
  filter(simp_functional_group == "Benthopelagic" & year %in% c(1960:2018)) %>% 
  group_by(area_name) %>% 
  summarise(C_extracted_benthopel = sum(C_extracted))

NutrientExtraction_perArea.FG.prop <- merge(NutrientExtraction_perArea.FG.prop, NutrientExtraction_perArea.benthopelagics, 
                                            by = "area_name",
                                            all.x = TRUE)
NutrientExtraction_perArea.FG.prop <- NutrientExtraction_perArea.FG.prop %>% 
  mutate(CE_benthopel_prop = C_extracted_benthopel/C_extracted)



# Proportion of Carbon Extraction through all other functional groups (Shrimps, Jellyfish, Sharks, Cephalopods, Rays, Flatfish, Crabs/Lobsters, Bathydemersal, Reef Fish, Bathypelagic, Krill) per Area, All Years
NutrientExtraction_perArea.miscfun <- Fisheries_NutrientExtraction %>% 
  filter(!(simp_functional_group %in% c("Pelagic", "Demersal", "Benthopelagic")) & year %in% c(1960:2018)) %>% 
  group_by(area_name) %>% 
  summarise(C_extracted_miscfun = sum(C_extracted))

NutrientExtraction_perArea.FG.prop <- merge(NutrientExtraction_perArea.FG.prop, NutrientExtraction_perArea.miscfun, 
                                            by = "area_name",
                                            all.x = TRUE)
NutrientExtraction_perArea.FG.prop <- NutrientExtraction_perArea.FG.prop %>% 
  mutate(CE_miscfun_prop = C_extracted_miscfun/C_extracted)



# Verify if proportions add up

NutrientExtraction_perArea.FG.prop %>% 
  mutate(Total_prop = CE_pel_prop + CE_dem_prop + CE_benthopel_prop + CE_miscfun_prop) %>% 
  filter(Total_prop > 1)







# Saving the proportion tables

NutrientExtraction_perArea_group.props <- merge(NutrientExtraction_perArea.TG.prop, 
                                                NutrientExtraction_perArea.FG.prop[ , c(1,3:10)], 
                                                by = "area_name")

fwrite(NutrientExtraction_perArea_group.props, "G:/My Drive/Utah State/Thesis/Thesis_Manuscripts/C, N, P and Ecology Manuscript/Tables/R Output Tables/NutrientExtraction_perArea_groupprops_ver2.csv")
fwrite(NutrientExtraction_perArea_group.props, "NutrientExtraction_perArea_groupprops_ver2.csv")









#### FIGURE 7 ####

# The following code is a second attempt which featured a nutrient (and landings) per plot with each time period featured (but not All Years). This also did not feature a secondary y-axis with the percentages.

test1 <- as.data.table(NutrientExtraction_perSFG %>% select(c(1, seq(2, 79, by = 7))) %>% pivot_longer(!simp_functional_group, names_to = "years", values_to = "metric_tons"))
test2 <- as.data.table(NutrientExtraction_perSFG %>% select (c(1, seq(3, 80, by = 7))) %>% pivot_longer(!simp_functional_group, names_to = "years", values_to = "metric_tons_SD"))

NE_perSFG_forPlot <- cbind(test1, test2)
NE_perSFG_forPlot <- NE_perSFG_forPlot %>% select(-c(4:5))
view(NE_perSFG_forPlot)


NE_perSFG_forPlot$nutrient <- "NA"
NE_perSFG_forPlot[years %in% c("C_extracted", "C_extracted_1960_64", "C_extracted_1993_97", "C_extracted_2014_18")]$nutrient <- "Carbon"
NE_perSFG_forPlot[years %in% c("N_extracted", "N_extracted_1960_64", "N_extracted_1993_97", "N_extracted_2014_18")]$nutrient <- "Nitrogen"
NE_perSFG_forPlot[years %in% c("P_extracted", "P_extracted_1960_64", "P_extracted_1993_97", "P_extracted_2014_18")]$nutrient <- "Phosphorous"


NE_perSFG_forPlot[years %in% c("C_extracted", "N_extracted", "P_extracted")]$years <- "All Years"
NE_perSFG_forPlot[years %in% c("C_extracted_1960_64", "N_extracted_1960_64", "P_extracted_1960_64")]$years <- "1960-64"
NE_perSFG_forPlot[years %in% c("C_extracted_1993_97", "N_extracted_1993_97", "P_extracted_1993_97")]$years <- "1993-97"
NE_perSFG_forPlot[years %in% c("C_extracted_2014_18", "N_extracted_2014_18", "P_extracted_2014_18")]$years <- "2014-18"
view(NE_perSFG_forPlot)



new_labs <- c("Benthopelagic" = "Benthopelagic", "Cephalopods" = "Cephalopods","Demersal" = "Demersal", 
              "Misc. Dem. Inverts" = "Misc. dem. inverts", "Pelagic" = "Pelagic", "Reef Fish" = "Reef fish", "Shrimps" = "Shrimps")





##### Extraction Plots #####
CE_perSFG_plot2 <- NE_perSFG_forPlot %>% 
  filter(simp_functional_group %in% c("Benthopelagic", "Cephalopods" ,"Demersal", 
                                      "Misc. Dem. Inverts", "Pelagic", "Reef Fish", "Shrimps"),
         years != "All Years",
         nutrient == "Carbon") %>% 
  ggplot(aes(x = simp_functional_group, y=metric_tons, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, limits = c(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                             "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  scale_y_continuous(name = "Million tonnes", expand = c(0.02, 0.02),
                     limits = c(0, 25759737),
                     breaks = c(0, 5e+6, 10e+6, 15e+6, 20e+6, 25e+6),
                     labels = c("0", "5", "10", "15", "20", "25")) +
  scale_fill_manual(values = cols_years) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "a", subtitle = "Carbon") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title.y = element_text(size=12),
        legend.position = c(.83, .65),
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))



NE_perSFG_plot2 <- NE_perSFG_forPlot %>% 
  filter(simp_functional_group %in% c("Benthopelagic", "Cephalopods" ,"Demersal", 
                                      "Misc. Dem. Inverts", "Pelagic", "Reef Fish", "Shrimps"),
         years != "All Years",
         nutrient == "Nitrogen") %>% 
  ggplot(aes(x = simp_functional_group, y=metric_tons, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, limits = c(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                             "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  scale_y_continuous(name = "Million tonnes", expand = c(0.02, 0.02),
                     limits = c(0, 6537389),
                     breaks = c(0, 1.3e+6, 2.6e+6, 3.9e+6, 5.2e+6, 6.5e+6),
                     labels = c("0", "1.3", "2.6", "3.9", "5.2", "6.5")) +
  scale_fill_manual(values = cols_years) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "b", subtitle = "Nitrogen") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")



PE_perSFG_plot2 <- NE_perSFG_forPlot %>% 
  filter(simp_functional_group %in% c("Benthopelagic", "Cephalopods" ,"Demersal", 
                                      "Misc. Dem. Inverts", "Pelagic", "Reef Fish", "Shrimps"),
         years != "All Years",
         nutrient == "Phosphorous") %>% 
  ggplot(aes(x = simp_functional_group, y=metric_tons, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18"))), ) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, limits = c(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                             "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  scale_y_continuous(name = "Million tonnes", expand = c(0.02, 0.02),
                     limits = c(0, 1400000),
                     breaks = c(0, 0.25e+6, 0.5e+6, 0.75e+6, 1.0e+6, 1.25e+6),
                     labels = c("0", "0.25", "0.50", "0.75", "1.00", "1.25")) +
  scale_fill_manual(values = cols_years) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "c", subtitle = "Phosphorous") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")




##### Comparison Plot #####
test1 <- as.data.table(GloAv_NutrientExtraction_perSFG %>% select(c(1, seq(102, 190, by = 8))) %>% pivot_longer(!simp_functional_group, names_to = "years", values_to = "per_change"))
test2 <- as.data.table(GloAv_NutrientExtraction_perSFG %>% select(c(1, seq(103, 191, by = 8))) %>% pivot_longer(!simp_functional_group, names_to = "years", values_to = "per_change_SD"))

GloAv_NE_perSFG_forPlot <- cbind(test1, test2)
GloAv_NE_perSFG_forPlot <- GloAv_NE_perSFG_forPlot %>% select(-c(4:5))


GloAv_NE_perSFG_forPlot$nutrient <- "NA"
GloAv_NE_perSFG_forPlot[years %in% c("C_dif_per", "C_dif_1960_64_per", "C_dif_1993_97_per", "C_dif_2014_18_per")]$nutrient <- "Carbon"
GloAv_NE_perSFG_forPlot[years %in% c("N_dif_per", "N_dif_1960_64_per", "N_dif_1993_97_per", "N_dif_2014_18_per")]$nutrient <- "Nitrogen"
GloAv_NE_perSFG_forPlot[years %in% c("P_dif_per", "P_dif_1960_64_per", "P_dif_1993_97_per", "P_dif_2014_18_per")]$nutrient <- "Phosphorous"

GloAv_NE_perSFG_forPlot[years %in% c("C_dif_per", "N_dif_per", "P_dif_per")]$years <- "All Years"
GloAv_NE_perSFG_forPlot[years %in% c("C_dif_1960_64_per", "N_dif_1960_64_per", "P_dif_1960_64_per")]$years <- "1960-64"
GloAv_NE_perSFG_forPlot[years %in% c("C_dif_1993_97_per", "N_dif_1993_97_per", "P_dif_1993_97_per")]$years <- "1993-97"
GloAv_NE_perSFG_forPlot[years %in% c("C_dif_2014_18_per", "N_dif_2014_18_per", "P_dif_2014_18_per")]$years <- "2014-18"



# The following plots are showing the percentage of extraction composed by each trophic level where each individual plot is a nutrient and each group is a time period.

plot_SFG <- c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
              "Shrimps", "Misc. Dem. Inverts", "Reef Fish")

# Comparison plot
perSFG_per_dif_plot <- GloAv_NE_perSFG_forPlot %>% filter(years == "All Years") %>% 
  ggplot(aes(x = simp_functional_group, y=per_change, 
             fill=factor(nutrient, levels= c("Carbon", "Nitrogen", "Phosphorous")))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, limits = plot_SFG,
                   labels = new_labs) +
  scale_y_continuous(name = "Percent difference", expand = c(0.02, 0.02),
                     limits = c(-20, 200),
                     breaks = c(0, 50, 100, 150, 200)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(title = "Nutrient")) +
  labs(title = "d") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        legend.key.size = unit(4, 'mm'),
        axis.title =element_text(size=12),
        legend.position = c(.25, .7),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))



# Export figure
setwd("G:/My Drive/Utah State/Thesis/Thesis_Manuscripts/C, N, P and Ecology Manuscript/Figures")
jpeg("Figure7_ver2_manuscript.jpeg", width = 18, height = 18, units = 'cm', res = 600)
Figure7 <- ggarrange(CE_perSFG_plot2, NE_perSFG_plot2, PE_perSFG_plot2, perSFG_per_dif_plot)
annotate_figure(Figure7,
                bottom = text_grob("Functional group", size = 12))
dev.off()












#### SUPPLEMENTARY FIGURE 1 ####

annualCE_dif_plot <- ggplot(GloAv_NutrientExtraction_perYear, aes(x = year, y = C_dif)) +
  geom_ribbon(aes(ymin = C_dif_lowCI, ymax = C_dif_highCI), alpha = 0.6, fill = "#56B4E9") +
  geom_line(lwd = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  scale_x_continuous(name = NULL, breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1950", "1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Thousand tonnes", breaks = c(0, 0.4e+06, 0.8e+06, 1.2e+06, 1.6e+06, 2e+06), 
                     labels = c("0", "400", "800", "1,200", "1,600", "2,000"), limits = c(0, 2e+06)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "a", subtitle = "Carbon") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))


annualNE_dif_plot <- ggplot(GloAv_NutrientExtraction_perYear, aes(x = year, y = N_dif)) +
  geom_ribbon(aes(ymin = N_dif_lowCI, ymax = N_dif_highCI), alpha = 0.6, fill = "#009E73") +
  geom_line(lwd = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  scale_x_continuous(name = NULL, breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1950", "1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Thousand tonnes", breaks = c(-135e+03, -90e+03, -45e+03, 0, 45e+03, 90e+03, 135e+03),
                     labels = c("-135", "-90", "-45", "0", "45", "90", "135"), limits = c(-110e+3, 135e+03)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "b", subtitle = "Nitrogen") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))


annualPE_dif_plot <- ggplot(GloAv_NutrientExtraction_perYear, aes(x = year, y = P_dif)) +
  geom_ribbon(aes(ymin = P_dif_lowCI, ymax = P_dif_highCI), alpha = 0.6, fill = "#D55E00") +
  geom_line(lwd = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  scale_x_continuous(name = NULL, breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1950", "1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = "Thousand tonnes", breaks = c(-150e+03, -100e+03, -50e+03, 0, 50e+03, 100e+03, 1500e+03), 
                     labels = c("-150", "-100", "-50", "0", "50", "100", "150"), limits = c(-150e+3, 80e+03)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "c", subtitle = "Phosphorus") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))



setwd("G:/My Drive/Utah State/Thesis/Thesis_Manuscripts/C, N, P and Ecology Manuscript/Figures")
jpeg("SupFig1_manuscript.jpeg", width = 18, height = 13.71, units = 'cm', res = 600)
Figure1_comp <- ggarrange(annualCE_dif_plot, annualNE_dif_plot, annualPE_dif_plot)
annotate_figure(Figure1_comp, bottom = text_grob("Year", size = 12))
dev.off()






#### SUPPLEMENTARY FIGURES 2-4 ####

# Maps for Supplementary Figures 2-4 were made in ArcGIS Pro 3.1.2 using the data included in the 'NutrientRatios_perArea' data frame listed below. This data frame was exported as a CSV and imported into ArcGIS Pro to create the figures.

