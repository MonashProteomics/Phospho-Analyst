# import libraries
library(tidyverse)
library(plyr)
library(rstudioapi)
library(DEP)
library(SummarizedExperiment)
# # Select the folder
# selected_folder <- selectDirectory(
#   caption = "Select Directory",
#   label = "Select",
#   path = getActiveProject()
# )
# Set the path
setwd("/Users/hailey/Documents/UNITS/Thesis/MaxQuant/LFQ-Phospho")
getwd()

# 1. Load the Phospho(STY)site.txt file in R
list.files()
Phospho_data <- read.delim("Phospho (STY)Sites.txt")

# 2. Load experimental design file
exp_design <- read.delim("exp_design_p16_0022.txt")

# 3. Remove Reverse sequences
# 4. Remove potential contaminants
# drop the two columns
Phospho_data <- Phospho_data %>% filter(Reverse != '+') %>% filter(Potential.contaminant != '+') %>%
  select(- 'Reverse',-'Potential.contaminant')

# 5. Expand Site table
# get all intensity columns
intensity <- grep("^Intensity.+|Intensity", colnames(Phospho_data)) 
# get the required intensity columns
intensity_cols <- grep("^Intensity.+___\\d", colnames(Phospho_data))
intensity_names <- colnames(Phospho_data[,intensity_cols])
intensity_names
# get the intensity columns need to be dropped
drop_cols <- setdiff(intensity, intensity_cols)
# drop columns
Phospho_data_new <- subset(Phospho_data, select = -drop_cols)

# expand site table
Phospho_data_ex <- Phospho_data_new %>% tidyr::pivot_longer(cols = contains(intensity_names),
                                                        names_to = c('.value','Multiplicity'),
                                                        names_pattern = '(.*)___(.)',
                                                        values_drop_na = TRUE)
# check the size of data
dim(Phospho_data_new)
dim(Phospho_data_ex)

# 6. Log2 transform the data
intensity_cols_new <- grep("^Intensity.", colnames(Phospho_data_ex))
Phospho_data_ex[,intensity_cols_new] <- log2(as.matrix(Phospho_data_ex[,intensity_cols_new]))
# change infinite value to 0
Phospho_data_ex[mapply(is.infinite, Phospho_data_ex)] <- 0

# 7. Filter based on the localization probability ( remove rows with <0.75 localization probability
Phospho_data_ex <- Phospho_data_ex %>% filter(Localization.prob >= 0.75)
dim(Phospho_data_ex)

# 8. Create new column for peptide sequence without scores and localisation probability
# use regex to extract peptide sequence from column "Phospho..STY..Probabilities"
peptide.sequence <- Phospho_data_ex$Phospho..STY..Probabilities %>% gsub("[^[A-Z]+","",.)
Phospho_data_pre <- dplyr::mutate(Phospho_data_ex,peptide.sequence, .after = "Phospho..STY..Score.diffs")

# # Write the data to a csv file
# write.csv(Phospho_data_pre, paste0(selected_folder,"/",'test_phospho_sites.csv'),
#           row.names = FALSE)

# 9. Convert the data into SummarisedExperiment object.
Phospho_data_pre <- Phospho_data_pre %>% 
  mutate(name = paste(Proteins,Positions.within.proteins,Multiplicity, sep = '_'))
Phospho_data_pre <- Phospho_data_pre %>% 
  mutate(ID = paste(id,Multiplicity, sep = '_'))

Phospho_data_unique_names <- make_unique(Phospho_data_pre, 'name','ID', delim = ";")

intensity_ints <- grep("^Intensity.", colnames(Phospho_data_unique_names))
Phospho_data_se <- make_se(Phospho_data_unique_names, intensity_ints, exp_design)

# plot_frequency(Phospho_data_se)
