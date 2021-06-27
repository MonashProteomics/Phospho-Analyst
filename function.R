# import libraries
library(tidyverse)
library(plyr)
library(rstudioapi)
library(DEP)
library(SummarizedExperiment)
library(msImpute)
library(reticulate)
library(limma)

#### ===== Pre-processing data ===== #####

pre_precessing <- function(maxquant_data,expdesign){
  # 1. Load the Phospho(STY)site.txt file in R
  maxquant_output<- maxquant_data
  
  # 2. Load experimental design file
  exp_design <- expdesign
  
  # 3. Remove Reverse sequences
  # 4. Remove potential contaminants
  # drop the two columns
  maxquant_output <- maxquant_output %>% filter(Reverse != '+') %>% filter(Potential.contaminant != '+') %>%
    select(- 'Reverse',-'Potential.contaminant')
  
  # 5. Expand Site table
  # get all intensity columns
  intensity <- grep("^Intensity.+|Intensity", colnames(maxquant_output)) 
  # get the required intensity columns
  intensity_cols <- grep("^Intensity.+___\\d", colnames(maxquant_output))
  intensity_names <- colnames(maxquant_output[,intensity_cols])
  intensity_names
  # get the intensity columns need to be dropped
  drop_cols <- setdiff(intensity, intensity_cols)
  # drop columns
  data_new <- subset(maxquant_output, select = -drop_cols)
  
  # expand site table
  data_ex <- data_new %>% tidyr::pivot_longer(cols = contains(intensity_names),
                                                              names_to = c('.value','Multiplicity'),
                                                              names_pattern = '(.*)___(.)',
                                                              values_drop_na = TRUE)
  # check the size of data
  dim(data_new)
  dim(data_ex)
  
  # 6. Log2 transform the data
  # intensity_cols_new <- grep("^Intensity.", colnames(data_ex))
  
  # data_ex[,intensity_cols_new] <- log2(as.matrix(data_ex[,intensity_cols_new]))
  
  # # change infinite value to 0
  # data_ex[mapply(is.infinite, data_ex)] <- 0
  
  # 7. Filter based on the localization probability ( remove rows with <0.75 localization probability
  data_ex <- data_ex %>% filter(Localization.prob >= 0.75)
  dim(data_ex)
  
  # 8. Create new column for peptide sequence without scores and localisation probability
  # use regex to extract peptide sequence from column "Phospho..STY..Probabilities"
  peptide.sequence <- data_ex$Phospho..STY..Probabilities %>% gsub("[^[A-Z]+","",.)
  data_pre <- dplyr::mutate(data_ex,peptide.sequence, .after = "Phospho..STY..Score.diffs")
  
  # 9. Convert the data into SummarisedExperiment object.
  data_pre <- data_pre %>% 
    mutate(name = paste(Proteins,Positions.within.proteins,Multiplicity, sep = '_'))
  data_pre <- data_pre %>% 
    mutate(ID = paste(id,Multiplicity, sep = '_'))
  
  data_unique_names <- make_unique(data_pre, 'name','ID', delim = ";")
  
  intensity_ints <- grep("^Intensity.", colnames(data_unique_names))
  data_se <- make_se(data_unique_names, intensity_ints, exp_design)
  
  # 11. Remove rows containing many missing values
  ## If thr=0, it means all the replicates in each sample should have valid values, no missing values is allowed
  ## thr=1 means at least 2 out 3 replicates should have valid values in each sample 
  # Phospho_data_filter<-DEP::filter_missval(Phospho_data_se, thr = 1)
  data_filter<-DEP::filter_proteins(data_se,type = "fraction", min = 0.8)
  
  # 12. Normalise the data ( Median normalisation)
  data_norm<-DEP::normalize_vsn(data_filter)
  
  # meanSdPlot(Phospho_data_filter)
  # meanSdPlot(Phospho_data_norm)
  
  # 13. Replace missing values (imputation)
  assay_df <- assay(data_norm)
  rowData <- rowData(data_norm)
  colData <- colData(data_norm)
  
  # plotCV2(assay_df, trend = FALSE)
  
  # selectFeatures(assay_df, method="ebm", group=c("Control","Mutant")) 
  
  xcomplete <- msImpute(assay_df, method="v2-mnar", group=c("Control","Mutant")) 
  
  # put the imputation result back to se object
  data_impute <- SummarizedExperiment(assays=list(xcomplete),colData=colData,rowData=rowData)
  
  # Differential enrichment test, type of contrasts will be tested
  diff_all <-test_limma(data_impute,type='all')
#  diff_all <-test_diff(data_impute,type='all')
  dep <- add_rejections(diff_all) # Mark significant peptide sequences
  
  # result list
  result_list <- list('data_se' = data_se, 'data_filter' = data_filter,'data_norm' = data_norm, 
                      'data_impute' = data_impute,'diff_all' = diff_all, 'dep' = dep)
  return(result_list)
}

#### ===== limma BH FDR ===== #####

test_limma <- function(se, type = c("control", "all", "manual"),
                       control = NULL, test = NULL,
                       design_formula = formula(~ 0 + condition),
                       paired = FALSE) {
  #require("dplyr", "tidyr", "purrr")
  
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type),
                          class(design_formula) == "formula")
  if (paired == FALSE){
    design_formula <- design_formula
  }else{
    design_formula<-formula(~ 0 + condition + replicate)
  }
  
  
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  
  col_data <- colData(se)
  raw <- assay(se)
  
  if(any(!c("name", "ID") %in% colnames(rowData(se)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if(any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)), "'")
  }
  
  if(!is.null(control)) {
    # Show error if control input is not valid
    assertthat::assert_that(is.character(control),
                            length(control) == 1)
    if(!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"), "'",
           call. = FALSE)
    }
  }
  
  # variables in formula
  variables <- terms.formula(design_formula) %>%
    attr(., "variables") %>%
    as.character() %>%
    .[-1]
  
  # Throw error if variables are not col_data columns
  if(any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if(variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }
  
  # Obtain variable factors
  for(var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  
  # Make an appropriate design matrix
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  
  # Generate contrasts to be tested
  # Either make all possible combinations ("all"),
  # only the contrasts versus the control sample ("control") or
  # use manual contrasts
  conditions <- as.character(unique(condition))
  if(type == "all") {
    # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")
    
    if(!is.null(control)) {
      # Make sure that contrast containing
      # the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if(length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>%
          gsub(paste(control, "- ", sep = " "), "", .) %>%
          paste(" - ", control, sep = "")
      }
    }
    
  }
  if(type == "control") {
    # Throw error if no control argument is present
    if(is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")
    
    # Make contrasts
    cntrst <- paste(conditions[!conditions %in% control],
                    control,
                    sep = " - ")
  }
  if(type == "manual") {
    # Throw error if no test argument is present
    if(is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))
    
    if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"),
           "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
           "'.", call. = FALSE)
    }
    
    cntrst <- gsub("_vs_", " - ", test)
    
  }
  # Print tested contrasts
  message("Tested contrasts: ",
          paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))
  
  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  
  if(any(is.na(raw))) {
    for(i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
    }
  }
  
  eB_fit <- eBayes(contrast_fit)
  
  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", adjust.method="BH", coef = comp,
                    number = Inf, confint = TRUE)
    # res <- res[!is.na(res$t),]
    #fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    # res$qval <- res$adj.P.Value
    #res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }
  
  #limma_res<- topTable(eB_fit, sort.by = 'B', adjust.method="BH", coef = cntrst, number = Inf, confint = T )
  # limma_res$comparison <- rep(cntrst, dim(limma_res)[1])
  #limma_res <- rownames_to_column(limma_res)
  # Retrieve the differential expression test results
  limma_res <- purrr::map_df(cntrst, retrieve_fun)
  
  # Select the logFC, CI and qval variables
  table <- limma_res %>%
    dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, adj.P.Val, comparison) %>%
    dplyr::mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
    tidyr::gather(variable, value, -c(rowname,comparison)) %>%
    dplyr::mutate(variable = dplyr::recode(variable, logFC = "diff", P.Value = "p.val", adj.P.Val = "p.adj")) %>%
    tidyr::unite(temp, comparison, variable) %>%
    tidyr::spread(temp, value)
  rowData(se) <- merge(rowData(se), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE)
  return(se)
  #return(table)
}

coef_variation<-function(x){
  coef=sd(x)/mean(x)
}

#### Plot CVs

plot_cvs<-function(se) {
  
  ## backtransform data
  untransformed_intensity<- 2^(assay(se))
  exp_design<-colData(se)
  
  ### merge untransformed to exp design and calculate cvs
  
  cvs_group<- untransformed_intensity %>% data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather("ID", "Intensity", -rowname) %>%
    dplyr::left_join(.,data.frame(exp_design), by="ID") %>%
    dplyr::group_by(rowname,condition) %>%
    dplyr::summarise(cvs=coef_variation(Intensity)) %>%
    dplyr::group_by(condition)%>%
    dplyr::mutate(condition_median=median(cvs))
  
  p1 <-  ggplot(cvs_group, aes(cvs, color=condition, fill=condition)) +
    geom_histogram(alpha=.5, bins= 20, show.legend = FALSE) +
    facet_wrap(~condition) +
    geom_vline(aes(xintercept=condition_median, group=condition),color='grey40',
               linetype="dashed") +
    labs(title= 'Sample Coefficient of Variation', x="Coefficient of Variation", y="Count") +
    theme_DEP2() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold")) 
  
  p1 +geom_text(aes(x=max(cvs_group$cvs)-0.6,
                    y=max(ggplot_build(p1)$data[[1]]$ymax*1.1), 
                    label=paste0("Median =",round(condition_median,2)*100,"%",by="")),
                show.legend = FALSE, size=4)
  
}
