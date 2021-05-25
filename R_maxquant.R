# import libraries
library(tidyverse)
library(plyr)
library(rstudioapi)
# Select the folder
selected_folder <- selectDirectory(
  caption = "Select Directory",
  label = "Select",
  path = getActiveProject()
)
# Set the path
setwd(selected_folder)
getwd()

# 1. Load the Phospho(STY)site.txt file in R
list.files()
Phospho_data <- read.delim("Phospho (STY)Sites.txt")

# 2. Load experimental design file
exp_design <- read.delim("exp_design_p16_0022.txt")

# 3. Remove Reverse sequences
# 4. Remove potential contaminants
re_Phospho_data <- Phospho_data %>% filter(Reverse != '+') %>% filter(Potential.contaminant != '+')

# save col_names for using
col_names <- colnames(re_Phospho_data)
col_names
# 5. Expand Site table
intensity_cols <- grep("^Intensity.+___\\d", col_names)
intensity_names <- colnames(re_Phospho_data[,intensity_cols])
intensity_names

# change zero value of intensities into NA
for (f in intensity_names){
  re_Phospho_data[f][re_Phospho_data[f] == 0] <-NA
}

# test
new_df <- data.frame()
expand_function <- function(a_row){
  rows <- matrix(rep(a_row,each=length(intensity_names)),nrow=length(intensity_names),
                 dimnames = list(intensity_names,col_names))
  for (i in 1:length(intensity_names)){
    for (j in 1:length(col_names)){
      if ( colnames(rows)[j] %in% intensity_names &
           row.names(rows)[i] != colnames(rows)[j]){
        rows[row.names(rows)[i],colnames(rows)[j]] = NA
      }
    }
  }
  new_df <- rbind(new_df, rows)
  return(new_df)
}

# Apply the function to insert value to the new dataFrame
result <- apply(re_Phospho_data,1,expand_function)

# unlist the result
result_df1 <-plyr::ldply (result, data.frame)

dim(result_df1)
# remove rows if selected columns(intensity_cols) are all NA
result_df2 <- result_df1[apply(result_df1[,intensity_names],1,function(x) any(!is.na(x))),]
dim(result_df2)

# change intensity back to numeric
to_number <- function(x){
  return(x %>% as.numeric())
}
for (f in intensity_names){
  result_df2[f] <- lapply(result_df2[f],to_number)
}

# 6. Log2 transform the data
log_Phospho_data <- log2(as.matrix(result_df2[,intensity_cols]))

# # check intensity_cols
# typeof(result_df2[1,'Intensity.C1___1'])
# intensity_cols
# length(intensity_cols)
# re_Phospho_data[1,intensity_cols]
# 7. Filter based on the localization probability ( remove rows with <0.75 localization probability
result_df3 <- result_df2 %>% filter(Localization.prob >= 0.75)
dim(result_df3)

# Write the data to a csv file
write.csv(result_df3, paste0(selected_folder,"/",'test_phospho_sites.csv'),
          row.names = FALSE)


# # Install bionic packages
# library(BiocManager)
# BiocManager::install('DEP',force = TRUE)
# library(DEP)
