## Overview
This repository contains the code and data for *'Fisheries disrupt marine nutrient cycles through biomass extraction'* by González Ortiz et al. (in review). 
The analysis was conducted in **R**, and the necessary R packages for replicating the analyses are outlined at the top of each script.

---

## Code Structure  
The scripts are organized to follow the method subsections of the manuscript:  

### **Predictive Models for Nutrient Composition**  
- Script: `Predictive_Model_Nutrient_Composition.R`  

### **Nutrient Compositions for Higher-Level Groups and Data-Sparse Species**  
- Script: `Nutrient_Compositions_Higher_Level_Groups_Data_Sparse_Species.R`  

### **Estimating Nutrient Extraction**  
- Scripts:  
  - `Nutrient_Extraction_Estimates.R`  
  - `Nutrient_Extraction_Uniform_Value_Comparisons.R`  
  - `Nutrient_Extraction_Molar_Ratios.R`  

### **Uncertainties and Limitations**  
- Script: `Fisheries_Uncertainty_Comparisons.R`  

### **Reproducing Figures**  
- Script: `Manuscript_Figures.R`  

---

## Data Availability  
Due to copyright restrictions, the raw nutrient and fisheries data cannot be reproduced in this repository. However, the following resources are provided:  

### **Nutrient Data Sources**  
- A list of sources from which the raw nutrient data was compiled is available as **Supplementary Data 1** in the manuscript and as a spreadsheet (`List_of_Sources.xlsx`) in the `data` folder of this repository.  

### **Fisheries Data**  
- The fisheries data was obtained from the **Sea Around Us database** (Pauly et al., 2020). Catch data was downloaded by Exclusive Economic Zone (EEZ) and High Seas region and compiled into a unified dataset.  
- Only catch data from industrial fisheries (both reported and reconstructed) was considered. Artisanal, subsistence, and recreational fisheries were not considered. Bycatch was also excluded.  

### **Extraction Estimates**  
- Total catches (tonnes year⁻¹) for each species and taxonomic group caught in an EEZ or high seas region between 1960 and 2018 resulted in a dataset of over 6.7 million entries.  
- For each entry, 100 random nutrient composition values were sampled from a normal distribution based on species/taxonomic group means and standard deviations. These values were multiplied by the total landed amount to generate 100 nutrient extraction estimates per entry.  
- Mean extraction values for carbon (C), nitrogen (N), and phosphorus (P) were calculated, along with standard deviations and 95% confidence intervals.  
- The full distribution matrices of extraction estimates for each nutrient are too large for GitHub but are available on **Figshare**: [https://doi.org/10.6084/m9.figshare.28500593](https://doi.org/10.6084/m9.figshare.28500593).  
- Extraction estimates were grouped by year, time period, marine region, trophic group, and functional group. The final estimates per category can be reproduced by executing the nutrient extraction scripts and are also available in the `data` folder.  

---

## Citation  
If you use this code or data, please cite the original manuscript once published.

---

## Contact  
Please contact **Adrian Gonzalez Ortiz** ([adgon@umich.edu](mailto:adgon@umich.edu)) for questions regarding the code or data. 
