# syntax=docker/dockerfile:1
FROM rocker/shiny-verse:4.2.1

RUN apt-get update && apt-get install -yq \
    libhdf5-dev libnetcdf-dev build-essential libgd-dev libbz2-dev libudunits2-dev libproj-dev libgdal-dev \
    texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra

RUN Rscript -e 'install.packages(c("devtools", "tidyverse", "ggrepel", "httr", "rjson", "mvtnorm", "tmvtnorm", \
"imputeLCMD", "plotly", "DT", "BiocManager","testthat", "RColorBrewer", "shiny","shinyalert","shinydashboard", \
"shinyjs", "svglite", "rhandsontable", "shinyBS", "shinyWidgets", "ggVennDiagram", "shinycssloaders","shiny.info", \
"fresh", "ggplot2", "ggExtra", "ggpubr", "reshape2", "generics", "stats", "matrixStats", "tinytex", "KSEAapp"), dependencies=TRUE)'

#FROM bioconductor/bioconductor_docker:RELEASE_3_15
RUN Rscript -e 'BiocManager::install(pkgs=c("DEP", "SummarizedExperiment", "limma", "ComplexHeatmap","pcaMethods","impute"), ask=F, dependencies=TRUE)'

COPY ./ /srv/shiny-server/phospho-analyst
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
RUN chmod -R +r /srv/shiny-server/phospho-analyst
