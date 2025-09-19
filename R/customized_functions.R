# missing value heatmap with percentage
plot_missval_new <- function(se, full_dataset = FALSE) {
  # Show error if input is not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  
  se_assay <- assay(se)
  # # Show error if there are no missing values
  # if(!any(is.na(se_assay))) {
  #   stop("No missing values in '", deparse(substitute(se)), "'",
  #        call. = FALSE)
  # }
  
  # Make assay data binary (1 = valid value, 0 = missing value)
  df <- se_assay %>% data.frame(.)
  colnames(df) <- paste0(colnames(df), " (", round((colMeans(is.na(df)))*100,2), "%)")
  missval <- df[apply(df, 1, function(x) any(is.na(x))), ]
  missval <- ifelse(is.na(missval), 0, 1)
  
  
  
  # Plot binary heatmap()
  # Show message if there are no missing values
  if(!any(is.na(se_assay))) {
    p <- ggplot() +
      annotate("text", x = 1,  y = 1,
               size = 7,
               label = "No missing values in dataset") + 
      theme_void()
  } else {
    if(full_dataset){ # full dataset
      df[is.na(df)] = 0
      colnames(df) <- colnames(se_assay) # remove percentage of missing values in full dataset names
      col_range <- seq(max(df),0, -(max(df)/7))
      used_color = circlize::colorRamp2(col_range,c("#800f26", "#bd1a26","#e3211c","#fc4e2a","#fd8d3c","#fed976","#ffeda0","#ffffff"))
      scale_list <- seq(0, plyr::round_any(max(df),5,ceiling),10)
      if(ncol(df) == 1) {
        col_clust = FALSE
      } else {
        col_clust = TRUE
      }
      if(nrow(df) == 1) {
        row_clust = FALSE
      } else {
        row_clust = TRUE
      }
      Color_heatmap <- Heatmap(df,
                               col = used_color,
                               cluster_rows = col_clust,
                               cluster_columns = row_clust,
                               show_row_names = FALSE,
                               show_column_names = TRUE,
                               column_names_side = "top",
                               column_names_gp = gpar(fontsize = 12),
                               name = "Log2 transformed",
                               heatmap_legend_param = list(color_bar = "continuous",
                                                           legend_direction = "horizontal",
                                                           legend_width = unit(5, "cm"),
                                                           title_position = "lefttop",
                                                           at = scale_list,
                                                           labels = scale_list),
                               top_annotation = get_annotation(se,"condition")
      )
      p <- draw(Color_heatmap, heatmap_legend_side = "top")
    } else { # subset dataset
      ht2 = Heatmap(missval,
                    col = c("white", "black"),
                    column_names_side = "top",
                    show_row_names = FALSE,
                    show_column_names = TRUE,
                    name = "Missing values pattern",
                    column_names_gp = gpar(fontsize = 14),
                    heatmap_legend_param = list(at = c(0, 1),
                                                labels = c("Missing value", "Valid value")))
      p <- draw(ht2, heatmap_legend_side = "top")
    }
  }
  return(p)
}

# convert protein id to 
convert_uniprot <- function(df, Accession){
  df$UniProt <- df[[Accession]] %>% gsub(";.*", "", .) 
  uniprot_df <- mapUniProt("UniProtKB_AC-ID", "UniProtKB", query = df$UniProt)
  uniprot_df <- uniprot_df  %>% dplyr::select("From","Gene.Names") %>% group_by(From) %>% summarize(Gene.Names = paste(Gene.Names, collapse = ";"))
  uniprot_df$Gene.Names <- uniprot_df$Gene.Names  %>% gsub("[[:space:]].*","",.)
  df <- df %>% left_join(uniprot_df, by =  c("UniProt" = "From"))
  df <- df %>% dplyr::relocate(Gene.Names, .after = Accession)
  colnames(df)[colnames(df) == "Gene.Names"] = "Gene.names"
  return(df)
}
