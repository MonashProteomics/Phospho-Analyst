# import libraries
library(tidyverse)
library(plyr)
library(rstudioapi)
library(DEP)
library(SummarizedExperiment)
library(msImpute)
library(reticulate)
library(limma)

# # Select the folder
# selected_folder <- selectDirectory(
#   caption = "Select Directory",
#   label = "Select",
#   path = getActiveProject()
# )
# Set the path
setwd("/Users/hailey/Documents/UNITS/Thesis/MaxQuant/LFQ-Phospho")
getwd()

# Load the Phospho(STY)site.txt and experimental design file in R
Phospho_data <- read.delim("Phospho (STY)Sites.txt")
exp_design <- read.delim("exp_design_p16_0022.txt")

# Apply the pre-processing function
result_list <-pre_precessing(Phospho_data, exp_design)

# get required se for visualization
Phospho_data_filter <- result_list$data_filter
Phospho_data_norm <- result_list$data_norm
Phospho_diff_all <- result_list$diff_all
Phospho_diff_all_1 <- result_list$diff_all_1

# QC plots
# Plot a barplot of the number of identified proteins per samples
plot_numbers(Phospho_data_norm)
# Plot a barplot of the protein identification overlap between samples
plot_coverage(Phospho_data_norm)

# normalization plot before and after
plot_normalization(Phospho_data_filter,
                   Phospho_data_norm)
# check missing value
plot_missval(Phospho_data_filter)
# with and without missing values density
plot_detect(Phospho_data_norm)

# befor and after imputation
plot_imputation(Phospho_data_norm,Phospho_diff_all)

# histogram
# plot_p_hist(dep)
plot_p_hist(Phospho_diff_all)
plot_p_hist(Phospho_diff_all_1) 

# correlation
# plot_cor(dep, significant = FALSE)
plot_cor(Phospho_diff_all, significant = FALSE)
plot_cor(Phospho_diff_all_1, significant = FALSE)

# coefficient of variation (CV)
plot_cvs(Phospho_diff_all)

# 14. PCA plot
pca_plot <- DEP::plot_pca(Phospho_diff_all, point_size = 4, indicate = "condition") 
pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                    size = 4,
                                    box.padding = unit(0.1, 'lines'),
                                    point.padding = unit(0.1, 'lines'),
                                    segment.size = 0.5)