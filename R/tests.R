### Test if column names are proper in experiment design file

exp_design_test<-function(exp_design){
  col_names<-colnames(exp_design)
  ## 
  if(!"label" %in% col_names){
    stop(safeError("The column 'label'(case sensitive) is not found in the Experimental Design File"))
  }
  
  else if (!"condition" %in% col_names){
    stop(safeError("The column 'condition' (case sensitive) is not found in the Experimental Design File"))
  }
  
  else if (!"replicate" %in% col_names){
    stop(safeError("The column 'replicate' (case sensitive) is not found in the Experimental Design File"))
  }
  
}

### Test if column names are proper in maxquant Phospho(STY)sites file
phospho_input_test<-function(maxquant_input){
  col_names<-colnames(maxquant_input)
  ## 
  if(!"Gene.names" %in% col_names){
    # stop(safeError("The column 'Gene names' is not found in the MaxQuant Phospho(STY)sites File"))
    "The column 'Gene names' is not found in the MaxQuant Phospho(STY)sites File"
  }
  
  else if (any(grepl("^Intensity.+___\\d", col_names))==FALSE){
    stop(safeError("Columns starting with 'Intensity' are not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  # else if (!"Protein.IDs" %in% col_names){
  #   stop(safeError("The column 'Protein IDs' is not found in the MaxQuant Phospho(STY)sites File"))
  # }
  
  else if (!"Reverse" %in% col_names){
    stop(safeError("The column 'Reverse' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Potential.contaminant" %in% col_names){
    stop(safeError("The column 'Potential contaminant' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Localization.prob" %in% col_names){
    stop(safeError("The column 'Localization prob' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Phospho..STY..Probabilities" %in% col_names){
    stop(safeError("The column 'Phospho (STY) Probabilities' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Positions.within.proteins" %in% col_names){
    stop(safeError("The column 'Positions within proteins' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Amino.acid" %in% col_names){
    stop(safeError("The column 'Amino acid' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Position" %in% col_names){
    stop(safeError("The column 'Position' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"id" %in% col_names){
    stop(safeError("The column 'id' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Protein" %in% col_names){
    stop(safeError("The column 'Protein' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
  else if (!"Protein.names" %in% col_names){
    stop(safeError("The column 'Protein names' is not found in the MaxQuant Phospho(STY)sites File"))
  }
  
}

### Test if column names are proper in maxquant ProteinGroups file
maxquant_input_test<-function(maxquant_input){
  col_names<-colnames(maxquant_input)
  ## 
  if(!"Gene.names" %in% col_names){
    stop(safeError("The column 'Gene names' is not found in the MaxQuant proteinGroups File"))
  }
  
  else if (any(grepl("LFQ", col_names))==FALSE){
    stop(safeError("Columns starting with 'LFQ' are not found in the MaxQuant proteinGroups File"))
  }
  
  else if (!"Protein.IDs" %in% col_names){
    stop(safeError("The column 'Protein IDs' is not found in the MaxQuant proteinGroups File"))
  }
  
  else if (!"Reverse" %in% col_names){
    stop(safeError("The column 'Reverse' is not found in the MaxQuant proteinGroups File"))
  }
  
  else if (!"Potential.contaminant" %in% col_names){
    stop(safeError("The column 'Potential contaminant' is not found in the MaxQuant proteinGroups File"))
  }
  
  else if (!"Only.identified.by.site" %in% col_names){
    stop(safeError("The column 'Only identified by site' is not found in the MaxQuant proteinGroups File"))
  }
  
  else if (!"Razor...unique.peptides" %in% col_names){
    stop(safeError("The column 'Razor + unique peptides' is not found in the MaxQuant proteinGroups File"))
  }
  
  else if (!"Protein.names" %in% col_names){
    stop(safeError("The column 'Protein names' is not found in the MaxQuant proteinGroups File"))
  }
  
}


### Test if experimental design names and LFQ column names match
test_match_lfq_column_design_phospho<-function(unique_data, lfq_columns, exp_design){
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(unique_data),
                          is.integer(lfq_columns),
                          is.data.frame(exp_design))
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(unique_data))) {
    stop(safeError("'Gene name', 'Position', and/or 'id' columns are not present in
          Phospho(STY)sites input file"
    ))
  }
  
  if(any(!c("label", "condition", "replicate") %in% colnames(exp_design))) {
    stop(safeError("'label', 'condition' and/or 'replicate' columns
         are not present in the experimental design"))
  }
  
  if(any(!apply(unique_data[, lfq_columns], 2, is.numeric))) {
    stop(safeError("specified 'columns' should be numeric
         Run make_se_parse() with the appropriate columns as argument"))
  }
  
  raw <- unique_data[, lfq_columns]
  
  expdesign <- mutate(exp_design, condition = make.names(condition)) %>%
    tidyr::unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  
  if(any(is.na(matched))) {
    stop(safeError("The labels/'run names' in the experimental design DID NOT match
         with intensity column names in maxquants Phospho(STY)sites file
         Run Phospho-Analyst with correct labels in the experimental design"))
  }
}

test_match_lfq_column_design<-function(unique_data, lfq_columns, exp_design){
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(unique_data),
                          is.integer(lfq_columns),
                          is.data.frame(exp_design))
  
  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(unique_data))) {
    stop(safeError("'Gene name' and/or 'Protein ID' columns are not present in
          protein groups input file"
        ))
  }
  
  if(any(!c("label", "condition", "replicate") %in% colnames(exp_design))) {
    stop(safeError("'label', 'condition' and/or 'replicate' columns
         are not present in the experimental design"))
  }
  
  if(any(!apply(unique_data[, lfq_columns], 2, is.numeric))) {
    stop(safeError("specified 'columns' should be numeric
         Run make_se_parse() with the appropriate columns as argument"))
  }
  
  raw <- unique_data[, lfq_columns]
  
  expdesign <- mutate(exp_design, condition = make.names(condition)) %>%
    tidyr::unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  
  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  
  if(any(is.na(matched))) {
    stop(safeError("The labels/'run names' in the experimental design DID NOT match
         with lfq column names in maxquants proteinGroups file
         Run Phospho-Analyst with correct labels in the experimental design"))
  }
}



enrichment_output_test<-function(dep, database){
  significant <- SummarizedExperiment::rowData(dep) %>%
    as.data.frame() %>%
    dplyr::select(name, significant) %>%
    dplyr::filter(significant) %>%
    dplyr::mutate(name = gsub("[.|_].*", "", name))
  test_enrichment_output<-enrichr_mod(significant$name, databases = database)
  if(nrow(test_enrichment_output[[1]])==0)
    stop(safeError("Enrichment analysis failed. 
                   Please check if the gene names are in Entrenz Gene Symbol format. 
                   (eg. ASM24, MYO6)"))
}

null_enrichment_test<-function(gsea_result,alpha=0.05){
  gsea_df<-gsea_result %>% group_by(contrast, var) %>% dplyr::filter(Adjusted.P.value <= alpha)
  if(nrow(gsea_df)==0){
    stop(safeError("No enriched term found at FDR cutoff 0.05. 
                   Enrichment plot could not be displayed. 
                   However, the results (non-significant hits) can still be accessed 
                   through 'Download table' tab."))
  }
}

ids_test<-function(filtered_data){
  if("Evidence.IDs" %in% colnames(filtered_data)){
    filtered_data$`Evidence.IDs`<-stringr::str_trunc(as.character(filtered_data$`Evidence.IDs`), 25000)
  }
  if("MS.MS.IDs" %in% colnames(filtered_data)){
    filtered_data$`MS.MS.IDs`<-stringr::str_trunc(as.character(filtered_data$`MS.MS.IDs`), 25000)
  }
  
  return(filtered_data)
  
}