###############################################################################
#                                                                             #
#                            Manuscript Figures                               #
#                                                                             #
###############################################################################



#### Necessary Libraries ####

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





#### Exporting Figure 1 ####
jpeg("Figure1_ver2_manuscript.jpeg", width = 18, height = 19.38, units = 'cm', res = 600)

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














#### Figure 4 - 6 panel - With Comparisons ####
cols_years <- c("#332288", "#DDCC77", "#AA4499", "#117733")


test1 <- as.data.table(NutrientExtraction_perTG %>% select(c(1, seq(2, 79, by = 7))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "metric_tons"))
test2 <- as.data.table(NutrientExtraction_perTG %>% select (c(1, seq(3, 80, by = 7))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "metric_tons_SD"))

NE_perTG_forPlot <- cbind(test1, test2)
NE_perTG_forPlot <- NE_perTG_forPlot %>% select(-c(4:5))


NE_perTG_forPlot$nutrient <- "NA"
NE_perTG_forPlot[years %in% c("C_extracted", "C_extracted_1960_64", "C_extracted_1993_97", "C_extracted_2014_18")]$nutrient <- "Carbon"
NE_perTG_forPlot[years %in% c("N_extracted", "N_extracted_1960_64", "N_extracted_1993_97", "N_extracted_2014_18")]$nutrient <- "Nitrogen"
NE_perTG_forPlot[years %in% c("P_extracted", "P_extracted_1960_64", "P_extracted_1993_97", "P_extracted_2014_18")]$nutrient <- "Phosphorous"


NE_perTG_forPlot[years %in% c("C_extracted", "N_extracted", "P_extracted")]$years <- "All Years"
NE_perTG_forPlot[years %in% c("C_extracted_1960_64", "N_extracted_1960_64", "P_extracted_1960_64")]$years <- "1960-64"
NE_perTG_forPlot[years %in% c("C_extracted_1993_97", "N_extracted_1993_97", "P_extracted_1993_97")]$years <- "1993-97"
NE_perTG_forPlot[years %in% c("C_extracted_2014_18", "N_extracted_2014_18", "P_extracted_2014_18")]$years <- "2014-18"

# The following plots are showing the percentage of extraction composed by each trophic level where each individual plot is a nutrient and each group is a time period.

CE_perTG_plot <- NE_perTG_forPlot %>% filter(years != "All Years" & nutrient == "Carbon") %>% 
  ggplot(aes(x = trophic_group, y=metric_tons, 
             fill=factor(years, levels=c("2014-18", "1993-97", "1960-64", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 28000000),
                     breaks = c(0, 5.5e+6, 11e+6, 16.5e+6, 22e+6, 27.5e+6),
                     labels = c("0", "5.5", "11.0", "16.5", "22.0", "27.5")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "c") +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5, size = 11, color = "black"), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), axis.title =element_text(size=12),
        legend.position = "none")

NE_perTG_plot <- NE_perTG_forPlot %>% filter(years != "All Years" & nutrient == "Nitrogen") %>% ggplot(aes(x = trophic_group, y=metric_tons,
                                                                                                           fill=factor(years, levels=c("2014-18", "1993-97", "1960-64", "All Years")))) +
  geom_bar(color = "black",stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 7030826),
                     breaks = c(0, 1.4e+6, 2.8e+6, 4.2e+6, 5.6e+6, 7e+6),
                     labels = c("0", "1.4", "2.8", "4.2", "5.6", "7.0")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "d") +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5, size = 11, color = "black"), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), axis.title =element_text(size=12),
        legend.position = "none")

PE_perTG_plot <- NE_perTG_forPlot %>% filter(years != "All Years" & nutrient == "Phosphorous") %>% ggplot(aes(x = trophic_group, y=metric_tons, 
                                                                                                              fill=factor(years, levels=c("2014-18", "1993-97", "1960-64", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 1553681),
                     breaks = c(0, 0.3e+6, 0.6e+6, 0.9e+6, 1.2e+6, 1.5e+6),
                     labels = c("0", "0.3", "0.6", "0.9", "1.2", "1.5")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "e") +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5, size = 11, color = "black"), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), axis.title =element_text(size=12),
        legend.position = "none")

##### Figure 4 - Comparison Plot #####

CE_perTG_per_dif_plot <- GloAv_NE_perTG_forPlot %>% filter(nutrient == "Carbon") %>% 
  ggplot(aes(x = trophic_group, y=per_change, 
             fill=factor(years, levels=c("All Years", "2014-18", "1993-97", "1960-64")))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators"),
                   labels = NULL) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(-20, 40),
                     breaks = c(-20, -10, 0, 10, 20, 30, 40)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = NULL) +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5, size = 11, color = "black"), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), axis.title =element_text(size=12), 
        legend.position = "none")

NE_perTG_per_dif_plot <- GloAv_NE_perTG_forPlot %>% filter(nutrient == "Nitrogen") %>% 
  ggplot(aes(x = trophic_group, y=per_change, 
             fill=factor(years, levels=c("All Years", "2014-18", "1993-97", "1960-64")))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators"),
                   labels = NULL) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(-20, 40),
                     breaks = c(-20, -10, 0, 10, 20, 30, 40)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = NULL) +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5, size = 11, color = "black"), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), axis.title =element_text(size=12), 
        legend.position = c(.75, .7),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

PE_perTG_per_dif_plot <- GloAv_NE_perTG_forPlot %>% filter(nutrient == "Phosphorous") %>% 
  ggplot(aes(x = trophic_group, y=per_change, 
             fill=factor(years, levels=c("All Years", "2014-18", "1993-97", "1960-64")))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators"),
                   labels = NULL) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(-20, 40),
                     breaks = c(-20, -10, 0, 10, 20, 30, 40)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = NULL) +
  theme_classic() +
  theme(axis.text.x=element_text(hjust=0.5, size = 11, color = "black"), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), axis.title =element_text(size=12), 
        legend.position = "none")









#### Figure 4 - 4 panel - with Comparisons ####
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

##### Figure 4 - Comparison Plot #####

perTG_per_dif_plot <- GloAv_NE_perTG_forPlot %>% filter(years == "All Years") %>% 
  ggplot(aes(x = trophic_group, y=per_change, 
             fill=factor(nutrient, levels= c("Carbon", "Nitrogen", "Phosphorous")))) +
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












#### TG - Nutrient Ratios ####


test1 <- as.data.table(NutrientRatios_perTG %>% select(c(1, seq(2, 46, by = 4))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "ratios"))
test2 <- as.data.table(NutrientRatios_perTG %>% select (c(1, seq(3, 47, by = 4))) %>% pivot_longer(!trophic_group, names_to = "years", values_to = "ratios_SD"))

Ratios_perTG_forPlot <- cbind(test1, test2)
Ratios_perTG_forPlot <- Ratios_perTG_forPlot %>% select(-c(4:5))



Ratios_perTG_forPlot <- Ratios_perTG_forPlot %>% 
  mutate(ratio_type = case_when(
    str_detect(years, "^C.N_mean") ~ "C:N",
    str_detect(years, "^C.P_mean") ~ "C:P",
    str_detect(years, "^N.P_mean") ~ "N:P",
    TRUE ~ "NA"
  ),
  years = case_when(
    str_detect(years, "1960_64") ~ "1960-64",
    str_detect(years, "1993_97") ~ "1993-97",
    str_detect(years, "2014_18") ~ "2014-18",
    TRUE ~ "All Years"
  ))


# C:N Plot
C.N_perTG_plot <- Ratios_perTG_forPlot %>% filter(ratio_type == "C:N") %>% 
  ggplot(aes(x = trophic_group, y=ratios, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = ratios - ratios_SD, 
                    ymax = ratios + ratios_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  # scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
  #                    limits = c(0, 25759737),
  #                    breaks = c(0, 5e+6, 10e+6, 15e+6, 20e+6, 25e+6),
  #                    labels = c("0", "5", "10", "15", "20", "25")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "a C:N") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")


C.P_perTG_plot <- Ratios_perTG_forPlot %>% filter(ratio_type == "C:P") %>% 
  ggplot(aes(x = trophic_group, y=ratios, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = ratios - ratios_SD, 
                    ymax = ratios + ratios_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 70)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "b C:P") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")



N.P_perTG_plot <- Ratios_perTG_forPlot %>% filter(ratio_type == "N:P") %>% 
  ggplot(aes(x = trophic_group, y=ratios, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = ratios - ratios_SD, 
                    ymax = ratios + ratios_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, 
                   limits = c("Low-level consumers", "Mesopredators", "High-level predators")) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 15)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "c N:P") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")


ggarrange(C.N_perTG_plot, C.P_perTG_plot, N.P_perTG_plot)






##### Exporting Figure 4 #####
# Original version
setwd("G:/My Drive/Utah State/Thesis/Thesis_Manuscripts/C, N, P and Ecology Manuscript/Figures")
jpeg("Figure4_ver2_manuscript.jpeg", width = 124, height = 130, units = 'mm', res = 600)
Figure4 <- ggarrange(CE_perTG_plot, NE_perTG_plot, PE_perTG_plot, ncol = 1, common.legend = TRUE, legend = "right")
annotate_figure(Figure4,
                bottom = text_grob("Million tonnes", size = 12),
                left = text_grob("Trophic group", rot = 90, size = 12))
dev.off()


# 6-panel version
jpeg("Figure4_ver2_w.Comparisons_6.panel.jpeg", width = 6.5, height = 5.75, units = 'in', res = 600)
plot_grid <- ggarrange(
  CE_perTG_plot, CE_perTG_per_dif_plot,
  NE_perTG_plot, NE_perTG_per_dif_plot,
  PE_perTG_plot, PE_perTG_per_dif_plot,
  ncol = 2, nrow = 3,
  widths = c(1.75, 1),
  align = "h")

# Create the "Million tonnes" and "Percent difference" labels at the bottom
bottom_labels <- arrangeGrob(
  textGrob("Million tonnes", gp = gpar(fontsize = 12), hjust = -0.15),
  textGrob("Percent difference", gp = gpar(fontsize = 12), hjust = 0.4),
  ncol = 2,
  widths = c(2, 1)  # Match the plot column widths
)

# Create the "Trophic group" label on the left side
left_label <- textGrob("Trophic group", rot = 90, gp = gpar(fontsize = 12))

# Combine everything into the final layout with left, plot grid, and bottom labels
final_plot <- grid.arrange(
  arrangeGrob(left_label, plot_grid, ncol = 2, widths = c(0.5, 10)),
  bottom_labels,
  ncol = 1,
  heights = c(10, 0.5)  # Adjust heights for spacing
)
dev.off()


# 4-panel version
jpeg("Figure4_ver2_w.Comparisons_4.panel.jpeg", width = 6.5, height = 5.75, units = 'in', res = 600)
Figure4.4panel <- ggarrange(CE_perTG_plot, NE_perTG_plot, PE_perTG_plot, perTG_per_dif_plot,
                            heights = c(0.55, 1))
annotate_figure(Figure4.4panel,
                bottom = text_grob("Trophic group", size = 12))
dev.off()







#### Figure 5 - ver2 ####

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













#### Figure 6 - ver2 ####

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
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
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
        axis.title =element_text(size=12),
        legend.position = c(.83, .75),
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
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
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
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
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

# Export figure
setwd("G:/My Drive/Utah State/Thesis/Thesis_Manuscripts/C, N, P and Ecology Manuscript/Figures")
jpeg("Figure6_ver2_Manuscript.jpeg", width = 6.5, height = 6.5, units = 'in', res = 600)
Figure6 <- ggarrange(CE_perSFG_plot2, NE_perSFG_plot2, PE_perSFG_plot2)
annotate_figure(Figure6,
                bottom = text_grob("Functional group", size = 12),
                left = text_grob("Million tonnes", rot = 90, size = 12))
dev.off()



#### Figure 6 - 6 panel - with Comparisons ####

CE_perSFG_plot3 <- NE_perSFG_forPlot %>% 
  filter(simp_functional_group %in% c("Benthopelagic", "Cephalopods" ,"Demersal", 
                                      "Misc. Dem. Inverts", "Pelagic", "Reef Fish", "Shrimps"),
         years != "All Years",
         nutrient == "Carbon") %>% 
  ggplot(aes(x = simp_functional_group, y=metric_tons, 
             fill=factor(years, levels = rev(c("1960-64", "1993-97", "2014-18"))))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, limits = rev(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                               "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 25759737),
                     breaks = c(0, 5e+6, 10e+6, 15e+6, 20e+6, 25e+6),
                     labels = c("0", "5", "10", "15", "20", "25")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "a") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")



NE_perSFG_plot3 <- NE_perSFG_forPlot %>% 
  filter(simp_functional_group %in% c("Benthopelagic", "Cephalopods" ,"Demersal", 
                                      "Misc. Dem. Inverts", "Pelagic", "Reef Fish", "Shrimps"),
         years != "All Years",
         nutrient == "Nitrogen") %>% 
  ggplot(aes(x = simp_functional_group, y=metric_tons, 
             fill=factor(years, levels= rev(c("1960-64", "1993-97", "2014-18"))))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, limits = rev(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                               "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 6537389),
                     breaks = c(0, 1.3e+6, 2.6e+6, 3.9e+6, 5.2e+6, 6.5e+6),
                     labels = c("0", "1.3", "2.6", "3.9", "5.2", "6.5")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "c") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")



PE_perSFG_plot3 <- NE_perSFG_forPlot %>% 
  filter(simp_functional_group %in% c("Benthopelagic", "Cephalopods" ,"Demersal", 
                                      "Misc. Dem. Inverts", "Pelagic", "Reef Fish", "Shrimps"),
         years != "All Years",
         nutrient == "Phosphorous") %>% 
  ggplot(aes(x = simp_functional_group, y=metric_tons, 
             fill=factor(years, levels= rev(c("1960-64", "1993-97", "2014-18"))))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = metric_tons - metric_tons_SD, 
                    ymax = metric_tons + metric_tons_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, limits = rev(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                               "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(0, 1400000),
                     breaks = c(0, 0.25e+6, 0.5e+6, 0.75e+6, 1.0e+6, 1.25e+6),
                     labels = c("0", "0.25", "0.50", "0.75", "1.00", "1.25")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "e") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")


##### Figure 6 - Comparison Plot #####
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

CE_perSFG_per_dif_plot <- GloAv_NE_perSFG_forPlot %>% filter(nutrient == "Carbon") %>% 
  ggplot(aes(x = simp_functional_group, y=per_change, 
             fill=factor(years, levels= rev(c("1960-64", "1993-97", "2014-18", "All Years"))))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, limits = rev(plot_SFG),
                   labels = NULL) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(-20, 260),
                     breaks = c(0, 50, 100, 150, 200, 250)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "b") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = c(.7, .65),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

NE_perSFG_per_dif_plot <- GloAv_NE_perSFG_forPlot %>% filter(nutrient == "Nitrogen") %>% 
  ggplot(aes(x = simp_functional_group, y=per_change, 
             fill=factor(years, levels= rev(c("1960-64", "1993-97", "2014-18", "All Years"))))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, limits = rev(plot_SFG),
                   labels = NULL) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(-20, 260),
                     breaks = c(0, 50, 100, 150, 200, 250)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "d") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")

PE_perSFG_per_dif_plot <- GloAv_NE_perSFG_forPlot %>% filter(nutrient == "Phosphorous") %>% 
  ggplot(aes(x = simp_functional_group, y=per_change, 
             fill=factor(years, levels= rev(c("1960-64", "1993-97", "2014-18", "All Years"))))) +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = per_change - per_change_SD, 
                    ymax = per_change + per_change_SD), 
                width = 0.5, position=position_dodge(.9)) +
  coord_flip() +
  scale_x_discrete(name = NULL, limits = rev(plot_SFG),
                   labels = NULL) +
  scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
                     limits = c(-20, 260),
                     breaks = c(0, 50, 100, 150, 200, 250)) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "f") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")





# Export figure
setwd("G:/My Drive/Utah State/Thesis/Thesis_Manuscripts/C, N, P and Ecology Manuscript/Figures")
jpeg("Figure6_ver2_w.Comparisons_6.panel.jpeg", width = 6.5, height = 7.5, units = 'in', res = 600)
plot_grid <- grid.arrange(
  CE_perSFG_plot3, CE_perSFG_per_dif_plot,
  NE_perSFG_plot3, NE_perSFG_per_dif_plot,
  PE_perSFG_plot3, PE_perSFG_per_dif_plot,
  ncol = 2, nrow = 3,
  widths = c(2, 1)  # Set column widths for plots
)

# Create the "Million tonnes" and "Percent change" labels at the bottom
bottom_labels <- arrangeGrob(
  textGrob("Million tonnes", gp = gpar(fontsize = 12), hjust = -0.15),
  textGrob("Percent change", gp = gpar(fontsize = 12), hjust = 0.4),
  ncol = 2,
  widths = c(2, 1)  # Match the plot column widths
)

# Create the "Functional group" label on the left side
left_label <- textGrob("Functional group", rot = 90, gp = gpar(fontsize = 12))

# Combine everything into the final layout with left, plot grid, and bottom labels
final_plot <- grid.arrange(
  arrangeGrob(left_label, plot_grid, ncol = 2, widths = c(0.5, 10)),
  bottom_labels,
  ncol = 1,
  heights = c(10, 0.5)  # Adjust heights for spacing
)
dev.off()


#### Figure 7 ####

# Extraction Plots
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
Figure6 <- ggarrange(CE_perSFG_plot2, NE_perSFG_plot2, PE_perSFG_plot2, perSFG_per_dif_plot)
annotate_figure(Figure6,
                bottom = text_grob("Functional group", size = 12))
dev.off()













##### Nutrient Ratios - SFG #####

test1 <- as.data.table(NutrientRatios_perSFG %>% select(c(1, seq(2, 46, by = 4))) %>% pivot_longer(!simp_functional_group, names_to = "years", values_to = "ratios"))
test2 <- as.data.table(NutrientRatios_perSFG %>% select (c(1, seq(3, 47, by = 4))) %>% pivot_longer(!simp_functional_group, names_to = "years", values_to = "ratios_SD"))

Ratios_perSFG_forPlot <- cbind(test1, test2)
Ratios_perSFG_forPlot <- Ratios_perSFG_forPlot %>% select(-c(4:5))



Ratios_perSFG_forPlot <- Ratios_perSFG_forPlot %>% 
  mutate(ratio_type = case_when(
    str_detect(years, "^C.N_mean") ~ "C:N",
    str_detect(years, "^C.P_mean") ~ "C:P",
    str_detect(years, "^N.P_mean") ~ "N:P",
    TRUE ~ "NA"
  ),
  years = case_when(
    str_detect(years, "1960_64") ~ "1960-64",
    str_detect(years, "1993_97") ~ "1993-97",
    str_detect(years, "2014_18") ~ "2014-18",
    TRUE ~ "All Years"
  ))


# C:N Plot
C.N_perSFG_plot <- Ratios_perSFG_forPlot %>% filter(ratio_type == "C:N") %>% 
  ggplot(aes(x = simp_functional_group, y=ratios, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = ratios - ratios_SD, 
                    ymax = ratios + ratios_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, limits = c(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                             "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  # scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
  #                    limits = c(0, 25759737),
  #                    breaks = c(0, 5e+6, 10e+6, 15e+6, 20e+6, 25e+6),
  #                    labels = c("0", "5", "10", "15", "20", "25")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "a C:N") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")


C.P_perSFG_plot <- Ratios_perSFG_forPlot %>% filter(ratio_type == "C:P") %>% 
  ggplot(aes(x = simp_functional_group, y=ratios, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = ratios - ratios_SD, 
                    ymax = ratios + ratios_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, limits = c(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                             "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  # scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
  #                    limits = c(0, 25759737),
  #                    breaks = c(0, 5e+6, 10e+6, 15e+6, 20e+6, 25e+6),
  #                    labels = c("0", "5", "10", "15", "20", "25")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "b C:P") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")


N.P_perSFG_plot <- Ratios_perSFG_forPlot %>% filter(ratio_type == "N:P") %>% 
  ggplot(aes(x = simp_functional_group, y=ratios, 
             fill=factor(years, levels=c("1960-64", "1993-97", "2014-18", "All Years")))) +
  geom_bar(color = "black", stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin = ratios - ratios_SD, 
                    ymax = ratios + ratios_SD), 
                width = 0.5, position=position_dodge(.9)) +
  scale_x_discrete(name = NULL, limits = c(c("Pelagic", "Demersal", "Benthopelagic", "Cephalopods", 
                                             "Shrimps", "Misc. Dem. Inverts", "Reef Fish")),
                   labels = new_labs) +
  # scale_y_continuous(name = NULL, expand = c(0.02, 0.02),
  #                    limits = c(0, 25759737),
  #                    breaks = c(0, 5e+6, 10e+6, 15e+6, 20e+6, 25e+6),
  #                    labels = c("0", "5", "10", "15", "20", "25")) +
  scale_fill_manual(values = cols_years,
                    breaks = c("1960-64", "1993-97", "2014-18", "All Years")) +
  guides(fill = guide_legend(title = "Time period")) +
  labs(title = "c N:P") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 11), 
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12),
        legend.position = "none")

ggarrange(C.N_perSFG_plot, C.P_perSFG_plot, N.P_perSFG_plot)



#### Nutrient Ratios - Area ####

# QUANTILES FOR MAPS

quantile(NutrientRatios_perArea$C.N_mean, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1))
# 3.363591 4.242760 4.319635 4.557575 4.651687 4.701315 4.792739 4.984906 5.621544 7.039742 
quantile(NutrientRatios_perArea$C.N_mean_1960_64, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE)
quantile(NutrientRatios_perArea$C.N_mean_1993_97, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE) # Max was highest here
quantile(NutrientRatios_perArea$C.N_mean_2014_18, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE) # Lowest value was from here

quantile(NutrientRatios_perArea$C.P_mean, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1))
# 31.61255  39.51315  41.05745  44.75972  47.92143  51.85496  56.32065  68.17303  71.67475 124.00140

quantile(NutrientRatios_perArea$C.P_mean_1960_64, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE)
quantile(NutrientRatios_perArea$C.P_mean_1993_97, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE) # Highest value here
quantile(NutrientRatios_perArea$C.P_mean_2014_18, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE) # lowest value here

quantile(NutrientRatios_perArea$N.P_mean, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1))
# 5.521145  8.700568  9.068076  9.618509 10.217537 10.747252 12.193586 14.303808 14.918744 26.414102  
quantile(NutrientRatios_perArea$N.P_mean_1960_64, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE)
quantile(NutrientRatios_perArea$N.P_mean_1993_97, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE) # lowest and highest value here
quantile(NutrientRatios_perArea$N.P_mean_2014_18, probs = c(0, 0.05, 0.1, 0.3, 0.45, 0.55, 0.7, 0.9, 0.95, 1),na.rm = TRUE)

# BOXPLOTS
Ratios_perArea_forPlot <- as.data.table(NutrientRatios_perArea %>% select(c(1, seq(2, 46, by = 4))) %>% pivot_longer(!area_name, names_to = "years", values_to = "ratios"))

Ratios_perArea_forPlot <- Ratios_perArea_forPlot %>% 
  mutate(ratio_type = case_when(
    str_detect(years, "^C.N_mean") ~ "C:N",
    str_detect(years, "^C.P_mean") ~ "C:P",
    str_detect(years, "^N.P_mean") ~ "N:P",
    TRUE ~ "NA"
  ),
  years = case_when(
    str_detect(years, "1960_64") ~ "1960-64",
    str_detect(years, "1993_97") ~ "1993-97",
    str_detect(years, "2014_18") ~ "2014-18",
    TRUE ~ "All Years"
  ))



ggplot(Ratios_perArea_forPlot, aes(x = years, y = ratios)) +
  geom_boxplot() +
  facet_wrap(~ ratio_type, scales = "free_y") +
  labs(x = "Time period", y = "Nutrient ratio", title = "Mean nutrient ratios per time period") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))






#### Supplementary Figure 1 - Annual extraction comparisons ####

annualCE_dif_plot <- ggplot(GloAv_NutrientExtraction_perYear, aes(x = year, y = C_dif)) +
  geom_line(lwd = 0.8, color = "#56B4E9") +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  scale_x_continuous(name = NULL, breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1950", "1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(0, 0.4e+06, 0.8e+06, 1.2e+06, 1.6e+06, 2e+06), 
                     labels = c("0", "400", "800", "1,200", "1,600", "2,000"), limits = c(0, 1.7e+06)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "a", subtitle = "Carbon") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))


annualNE_dif_plot <- ggplot(GloAv_NutrientExtraction_perYear, aes(x = year, y = N_dif)) +
  geom_line(lwd = 0.8, color = "#009E73") +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  scale_x_continuous(name = NULL, breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1950", "1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(-11e+03, 0, 11e+03, 22e+03, 33e+03, 44e+03), 
                     labels = c("-11", "0", "11", "22", "33", "44"), limits = c(-11e+3, 45e+03)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "b", subtitle = "Nitrogen") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))




annualPE_dif_plot <- ggplot(GloAv_NutrientExtraction_perYear, aes(x = year, y = P_dif)) +
  geom_line(lwd = 0.8, color = "#D55E00") +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, color = "black") +
  scale_x_continuous(name = NULL, breaks = c(1950, 1960, 1970, 1980, 1990, 2000, 2010, 2018),
                     labels =  c("1950", "1960", "1970", "1980", "1990", "2000", "2010", "2018"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = NULL, breaks = c(0, 5e+03, 10e+03, 15e+03, 20e+03, 25e+03), 
                     labels = c("0", "5", "10", "15", "20", "25"), limits = c(-1e+3, 28e+03)) +
  guides(color = guide_legend(title = "")) +
  labs(title = "c", subtitle = "Phosphorous") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 11), 
        plot.title = element_text(size = 13, face = "bold"),
        axis.text.y=element_text(size=11), 
        axis.title =element_text(size=12))


jpeg("SupplementaryFigure1_AnnualComparisons_manuscript.jpeg", width = 18, height = 10.83, units = 'cm', res = 600)
Figure1_comp <- ggarrange(annualCE_dif_plot, annualNE_dif_plot, annualPE_dif_plot)
annotate_figure(Figure1_comp,
                bottom = text_grob("Year", size = 12),
                left = text_grob("Thousand tonnes", rot = 90, size = 12))
dev.off()




#### Supplementary Figures  and  - vers2 ####

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









