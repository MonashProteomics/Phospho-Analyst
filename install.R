# specify the required packages
packages = c("tidyverse", 
             "testthat",
             "shiny",
             "shinydashboard",
             "shinyjs",
             "shinyalert",
             "DT",
             "ggrepel",
             "httr",
             "rjson",
             "svglite",
             "fresh",
             "ggplot2",
             "ggExtra",
             "ggpubr",
             "shinyWidgets",
             "stats",
             "matrixStats",
             "tinytex",
             "shinyBS",
             "KSEAapp",
             "rhandsontable",
             "shinycssloaders",
             "shiny.info",
             "ggVennDiagram",
             "scales",
             "UpSetR",
             "SummarizedExperiment",
             "DEP",
             "ComplexHeatmap",
             "limma")


bio_packages = c("SummarizedExperiment",
                 "DEP",
                 "ComplexHeatmap",
                 "limma")

# use this function to check if each package is on the local machine
# if any are not, the missing package(s) will be installed
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    if (x %in% bio_packages){
      install.packages("BiocManager")
      library('BiocManager')
      BiocManager::install(x, dependencies = TRUE)
    } else {
      install.packages(x, dependencies = TRUE)
    }
  }
})
