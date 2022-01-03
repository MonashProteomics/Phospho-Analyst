#Define server logic to read selected file ----
server <- function(input, output,session){
  options(shiny.maxRequestSize=100*1024^2)## Set maximum upload size to 100MB
  #  Show elements on clicking Start analysis button
  observeEvent(input$analyze ,{ 
    if(input$analyze==0){
      return()
    }
    shinyjs::hide("quickstart_info")
    shinyjs::show("panels")
    shinyjs::show("qc_tab")
  })
  
  # Hide other pages if only upload one of the phosphosite or protein data file
  observeEvent(input$analyze ,{ 
    if (is.null(input$file1)){
      hideTab(inputId = "panel_list", target = "Phosphosite")
      hideTab(inputId = "panel_list", target = "Comparison")
      hideTab(inputId = "panel_list", target = "Phosphosite(corrected)")
    }
    
    if (is.null(input$file2)){
      hideTab(inputId = "panel_list", target = "ProteinGroup")
      hideTab(inputId = "panel_list", target = "Comparison")
      hideTab(inputId = "panel_list", target = "Phosphosite(corrected)")
    }
    
  })
  
  
  observeEvent(input$analyze ,{
    if(input$analyze==0 ){
      return()
    }
    shinyalert("In Progress!", "Data analysis has started, wait until table and plots
                appear on the screen", type="info",
               closeOnClickOutside = TRUE,
               closeOnEsc = TRUE,
               timer = 25000) # timer in miliseconds (10 sec)

    # if(phospho_exp_data$Data$condition != protein_exp_data$Data_pr$condition){
    #   shinyalert("Warning!","The conditions in phosphosite and protein experimental design table are different,
    #              please check it.", type="warning",
    #              closeOnClickOutside = TRUE,
    #              closeOnEsc = TRUE)
    #   # httpuv:::service()
    #   # continue <<- !isTRUE(session$input$stopThis)
    # }
    # else(
    #   shinyalert("In Progress!", "Data analysis has started, wait until table and plots
    #             appear on the screen", type="info",
    #              closeOnClickOutside = TRUE,
    #              closeOnEsc = TRUE,
    #              timer = 25000) # timer in miliseconds (10 sec)
    # )

  })
  
  observe({
    if (input$tabs_selected=="demo"){
      shinyalert("Demo results loading!...", "Wait until table and plots
                appear on the screen", type="info",
                 closeOnClickOutside = TRUE,
                 closeOnEsc = TRUE,
                 timer = 6000)
    }
  })
  
  #### Phosphosite page logic ========== #############
  ####======= Render Functions
  output$volcano_cntrst <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("volcano_cntrst",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  ##comparisons
  output$contrast <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$contrast_1 <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_1",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$downloadTable <- renderUI({
    if(!is.null(dep())){
      selectizeInput("dataset",
                     "Choose a dataset to save" ,
                     c("Results","Original_matrix",
                       "Imputed_matrix",
                       "Full_dataset"))
    }
  })
  
  output$downloadButton <- renderUI({
    if(!is.null(dep())){
      downloadButton('downloadData', 'Save')
    }
  })
  
  output$downloadZip <- renderUI({
    if(!is.null(dep())){
      downloadButton('downloadZip1', 'Download result plots')
    }
  })
  output$downloadreport <- renderUI({
    if(!is.null(dep())){
      downloadButton('downloadReport', 'Download Report')
    }
  })
  
  output$downloadPlots <- renderUI({
    if(!is.null(dep())){
      downloadButton('downloadPlots1', 'Download Plots')
    }
  })
  
  ## make reactive elements of the input data
  phospho_data_input<-reactive({NULL})
  protein_data_input<-reactive({NULL})
  exp_design_input<-reactive({NULL})
  
  phospho_data_input<-eventReactive(input$analyze,{
    inFile<-input$file1
    if(is.null(inFile))
      return(NULL)
    temp_data<-read.delim(inFile$datapath,
                          header = TRUE,
                          fill= TRUE, # to fill any missing data
                          sep = "\t"
    )
    #validate(maxquant_input_test(temp_data))
    return(temp_data)
  })
  
  protein_data_input<-eventReactive(input$analyze,{
    inFile<-input$file2
    if(is.null(inFile))
      return(NULL)
    temp_data1<-read.delim(inFile$datapath,
                           header = TRUE,
                           fill= TRUE, # to fill any missing data
                           sep = "\t"
    )
    #validate(maxquant_input_test(temp_data))
    return(temp_data1)
  })
  
  ####======= interactive exp_design (phosphosite)=======####
  phospho_exp_data <- reactiveValues(Data = NULL)
  
  observe({
    req(input$file1)
    # read file
    tryCatch(
      {
        Data <- read.delim(input$file1$datapath,
                           header = TRUE,
                           fill= TRUE, # to fill any missing data
                           sep = "\t")
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    # set value to created exp_design table
    phospho_exp_data$Data <- get_exp_design(Data)
  })
  
  # overview changes in exp_design table
  output$print_phospho <- renderPrint({
    phospho_exp_data$Data
  })
  
  # editable exp_design table
  output$exp_phospho<-DT::renderDataTable(phospho_exp_data$Data,
                                          selection = 'none',
                                          # rownames = FALSE,
                                          edit = list(target = "all", 
                                                      disable = list(columns = c(0,1)),
                                                      numeric = list(columns = 3)),
                                          class = "display"
  )
  
  # proxy_exp <- DT::dataTableProxy("exp_phospho")
  
  observeEvent(input$exp_phospho_cell_edit, {
    phospho_exp_data$Data <<- editData(phospho_exp_data$Data, input$exp_phospho_cell_edit, 'exp_phospho')
  })
  
  # download edited exp_design table
  output$download_exp <- downloadHandler("phosphosite_exp_design.csv",
                                     content = function(file){
                                       write.csv(phospho_exp_data$Data, file, row.names = F)
                                     },
                                     contentType = "text/csv")
  
  # phosphosite exp_desgin file
  exp_design_input<-eventReactive(input$analyze,{
    inFile<-input$file3
    if (is.null(inFile) ||input$save_exp>0){
      temp_df <- phospho_exp_data$Data
      
    }
    else{
      temp_df<-read.delim(inFile$datapath,
                          header = TRUE,
                          sep="\t",
                          stringsAsFactors = FALSE)
      exp_design_test(temp_df)
      temp_df$label<-as.character(temp_df$label)
      temp_df$condition<-trimws(temp_df$condition, which = "left")
    }
    return(temp_df)
  })
  
  ####======= interactive exp_design (protein)=======####
  protein_exp_data <- reactiveValues(Data_pr = NULL)
  
  observe({
    req(input$file2)
    
    # read file
    tryCatch(
      {
        Data_pr <- read.delim(input$file2$datapath,
                         header = TRUE,
                         fill= TRUE, # to fill any missing data
                         sep = "\t")
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    # set value to created exp_design table
    protein_exp_data$Data_pr <- get_exp_design_pr(Data_pr)
  })
  
  # # overview changes in exp_design table
  # output$print_protein <- renderPrint({
  #   protein_exp_data$Data_pr
  # })
  
  # editable exp_design table
  output$exp_protein<-DT::renderDataTable(
    protein_exp_data$Data_pr,selection = 'none', 
    edit = list(target = "all", 
                disable = list(columns = c(0,1)),
                numeric = list(columns = 3)),
    class = "display"
    )
  
  proxy_pr <- DT::dataTableProxy("exp_protein")
  
  observeEvent(input$exp_protein_cell_edit, {
    protein_exp_data$Data_pr <<- editData(protein_exp_data$Data_pr, input$exp_protein_cell_edit, 'exp_protein')
  })
  
  # download edited exp_design table
  output$download_exp_pr <- downloadHandler("protein_exp_design.csv", 
                                     content = function(file){
                                       write.csv(protein_exp_data$Data_pr, file, row.names = F)
                                     },
                                     contentType = "text/csv")
  
  # proteinGroup exp_desgin file
  exp_design_input_1<-eventReactive(input$analyze,{
    inFile<-input$file4
    if (is.null(inFile)||input$save_exp>0){
      temp_df <- protein_exp_data$Data_pr
      
    }
    else{
      temp_df<-read.delim(inFile$datapath,
                          header = TRUE,
                          sep="\t",
                          stringsAsFactors = FALSE)
      exp_design_test(temp_df)
      temp_df$label<-as.character(temp_df$label)
      temp_df$condition<-trimws(temp_df$condition, which = "left")
    }
    return(temp_df)
  })
  
  ### Reactive components
  # Data cleaning
  cleaned_data<- reactive({
    ## check which dataset
    if(!is.null (phospho_data_input() )){
      phospho_data <- reactive({phospho_data_input()})
    }
    
    if(grepl('+',phospho_data()$Reverse)){
      filtered_data<-dplyr::filter(phospho_data(),Reverse!="+")
    }
    else{filtered_data<-phospho_data()}
    if(grepl('+',filtered_data$Potential.contaminant)){
      filtered_data<-dplyr::filter(filtered_data,Potential.contaminant!="+")
    }
    
    filtered_data<-ids_test(filtered_data)
    
    ## Expand Site table
    # get all intensity columns
    intensity <- grep("^Intensity.+|Intensity", colnames(filtered_data)) 
    # get the required intensity columns
    intensity_cols <- grep("^Intensity.+___\\d", colnames(filtered_data))
    intensity_names <- colnames( filtered_data[,intensity_cols])
    
    # ensure all intensity columns are numeric type
    filtered_data[,intensity_cols] <- sapply(filtered_data[,intensity_cols],as.numeric)
    
    # get the intensity columns need to be dropped
    drop_cols <- setdiff(intensity, intensity_cols)
    # drop columns
    data_new <- subset(filtered_data, select = -drop_cols)
    
    # expand site table
    data_ex <- data_new %>% tidyr::pivot_longer(cols = contains(intensity_names),
                                                names_to = c('.value','Multiplicity'),
                                                names_pattern = '(.*)___(.)',
                                                values_drop_na = TRUE)
    
    data_ex <- data_ex %>% dplyr::filter(Localization.prob >= 0.75)
    peptide.sequence <- data_ex$Phospho..STY..Probabilities %>% gsub("[^[A-Z]+","",.)
    data_pre <- dplyr::mutate(data_ex,peptide.sequence, .after = "Phospho..STY..Score.diffs")
    data_pre$Residue.Both <- map2(data_pre$Positions.within.proteins, data_pre$Amino.acid,create_Residue.Both_func)
    data_pre$Gene.names <- data_pre$Gene.names %>% toupper()
    
    # create unique name and ID columns for the data
    data_pre <- data_pre %>%
      mutate(name = paste(Gene.names,Positions.within.proteins, Multiplicity,sep = '_'))
    # data_pre <- data_pre %>% 
    #   mutate(name = paste(Protein,Positions.within.proteins, Multiplicity,sep = '_'))
    
    data_pre <- data_pre %>% 
      mutate(ID = paste(id,Multiplicity, sep = '_'))
    return(data_pre)
  })
  
  
  ## Convert the data into SummarisedExperiment object,removing multipple missing value rows
  processed_data<- reactive({
    if(!is.null (exp_design_input() )){
      exp_design<-reactive({exp_design_input()})
    }
    message(exp_design())
    data_pre <- cleaned_data()
    data_unique_names <- DEP::make_unique(data_pre, 'name','ID', delim = ";")
    
    intensity_ints <- grep("^Intensity.", colnames(data_unique_names))
    
    #test_match_lfq_column_design(data_unique,intensity_ints, exp_design())
    data_se <- DEP::make_se(data_unique_names, intensity_ints, exp_design())
    
    as.numeric(rowData(data_se)$Score.diff) # add test
    
    ## Check for matching columns in maxquant and experiment design file
    # Check number of replicates
    if(max(exp_design()$replicate)<3){
      threshold<-0
    } else  if(max(exp_design()$replicate)==3){
      threshold<-1
    } else if(max(exp_design()$replicate)<6 ){
      threshold<-2
    } else if (max(exp_design()$replicate)>=6){
      threshold<-trunc(max(exp_design()$replicate)/2)
    }
    
    # removing multipple missing value rows
    filter_missval(data_se,thr = threshold)
  })
  
  unimputed_table<-reactive({
    temp<-assay(processed_data())
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"original_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) 
    #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  imputed_data<-reactive({
    DEP::impute(processed_data(),input$imputation)
  })
  
  normalised_data<-reactive({
    normalize_vsn(imputed_data())
  })
  
  imputed_table<-reactive({
    temp<-assay(imputed_data())
    #tibble::rownames_to_column(temp,var = "ProteinID")
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"imputed_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  diff_all<-reactive({
    test_diff(normalised_data(),type = 'all')
  })
  
  dep<-reactive({
    if(input$fdr_correction=="BH"){
      diff_all<-test_limma(normalised_data(),type='all', paired = input$paired)
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    else{
      diff_all<-test_diff(normalised_data(),type='all')
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    
    if(length(unique(exp_design_input()$condition)) <= 2){
      return(diff_all_rej)
      
    }
    else if(length(unique(exp_design_input()$condition)) >= 3){
      anova_dep <- diff_all
      # get assay data
      intensity <- assay(anova_dep)
      exp_design <- exp_design_input()
      exp_design_rename<-exp_design
      exp_design_rename$label<-paste(exp_design_rename$condition, exp_design_rename$replicate, sep = "_")
      
      # reshape intensity columns
      data_reshape<-reshape2::melt(intensity,value.name = "intensity", variable.name = "label")
      colnames(data_reshape)<-c("uid", "label", "intensity")
      
      # Join the table
      data_experiment<-left_join(data_reshape, exp_design_rename, by="label")
      
      # apply anova function
      anova<-data_experiment %>%
        group_by(`uid`) %>%
        do(anova_function(.)) %>% dplyr::select(p.value) %>%
        ungroup()
      
      colnames(anova)<-c("name", "anova_p.val")
      
      # add anova p.value to row data
      rowData(anova_dep) <- merge(rowData(anova_dep), anova, by = 'name', sort = FALSE)
      anova_dep_rej <- DEP::add_rejections(anova_dep,alpha = input$p, lfc= input$lfc)
      
      # calculate adjusted anova p.value to data
      anova_adj <-anova
      anova_adj$anova_p.adj <- p.adjust(anova_adj$anova_p.val,method = input$fdr_correction)
      anova_adj <- anova_adj %>% select(-anova_p.val)
      
      # add adjusted anova p.value to row data
      rowData(anova_dep_rej) <- merge(rowData(anova_dep_rej), anova_adj, by = 'name', sort = FALSE)
      return(anova_dep_rej)
    }
    
  })

  comparisons<-reactive ({
    temp<-capture.output(DEP::test_diff(normalised_data(),type='all'),type = "message")
    gsub(".*: ","",temp)
    ## Split conditions into character vector
    unlist(strsplit(temp,","))
    ## Remove leading and trailing spaces
    trimws(temp)
  })
  
  
  ## Results plot inputs
  
  ## PCA Plot
  pca_input<-eventReactive(input$analyze ,{ 
    if(input$analyze==0 ){
      return()
    }
    if (num_total()<=500){
      if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
        pca_plot<-DEP::plot_pca(dep(), n=num_total(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
      else{
        pca_plot<-DEP::plot_pca(dep(), n=num_total(), point_size = 4, indicate = "condition") 
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
    }
    else{
      if(length(levels(as.factor(colData(dep())$replicate))) <= 6){
        pca_plot<-DEP::plot_pca(dep(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
      else{
        #pca_label<-SummarizedExperiment::colData(dep())$replicate
        pca_plot<-DEP::plot_pca(dep(), point_size = 4, indicate = "condition")
        #pca_plot<-pca_plot + geom_point()
        pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                                      size = 4,
                                                      box.padding = unit(0.1, 'lines'),
                                                      point.padding = unit(0.1, 'lines'),
                                                      segment.size = 0.5)
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        
        #        pca_plot<-DEP::plot_pca(dep(), point_size = 4, indicate = "condition")
        #         pca_plot + ggrepel::geom_text_repel(aes(label=SummarizedExperiment::colData(dep())$replicate),
        #                                        size = 5,
        #                                           box.padding = unit(0.1, 'lines'),
        #                                          point.padding = unit(0.1, 'lines'),
        #                                         segment.size = 0.5)
        return(pca_plot)
      }
    }
    
  })
  
  ### Heatmap Differentially expressed proteins
  heatmap_input<-eventReactive(input$analyze ,{ 
    if(input$analyze==0 ){
      return()
    }
    get_cluster_heatmap(dep(),
                        type="centered",kmeans = TRUE,
                        k=input$k_number, col_limit = 6,
                        indicate = "condition"
    )
  })
  
  ### Volcano Plot
  volcano_input <- reactive({
    if(!is.null(input$volcano_cntrst)) {
      if(length(unique(exp_design_input()$condition)) <= 2) {
        plot_volcano_new(dep(),
                         input$volcano_cntrst,
                         FALSE,
                         input$check_names,
                         input$p_adj)
      } else {
        plot_volcano_new(dep(),
                         input$volcano_cntrst,
                         input$check_anova,
                         input$check_names,
                         input$p_adj)
      }
      
    }
  })
  
  volcano_df<- reactive({
    if(!is.null(input$volcano_cntrst)) {
      get_volcano_df(dep(),
                     input$volcano_cntrst) 
      
    }
  })
  
  volcano_input_selected<-reactive({
    if(!is.null(input$volcano_cntrst)){
      
      if (!is.null(input$contents_rows_selected)){
        proteins_selected<-data_result()[c(input$contents_rows_selected),]## get all rows selected
      }
      else if(!is.null(input$protein_brush)){
        proteins_selected<-data_result()[data_result()[["Phosphosite"]] %in% protein_name_brush(), ] 
      }
      print(proteins_selected)
      ## convert contrast to x and padj to y
      diff_proteins <- grep(paste(input$volcano_cntrst, "_log2", sep = ""),
                            colnames(proteins_selected))
      if(input$p_adj=="FALSE"){
        padj_proteins <- grep(paste(input$volcano_cntrst, "_p.val", sep = ""),
                              colnames(proteins_selected))
      }
      else{
        padj_proteins <- grep(paste(input$volcano_cntrst, "_p.adj", sep = ""),
                              colnames(proteins_selected))
      }
      
      df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                               y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                               name = proteins_selected$Phosphosite)
      # print(df_protein)
      if(length(unique(exp_design_input()$condition)) <= 2){
        p<-plot_volcano_new(dep(),
                            input$volcano_cntrst,
                            FALSE,
                            input$check_names,
                            input$p_adj)
        
      } else {
        p<-plot_volcano_new(dep(),
                            input$volcano_cntrst,
                            input$check_anova,
                            input$check_names,
                            input$p_adj)
      }
      p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
        ggrepel::geom_text_repel(data = df_protein,
                                 aes(x, y, label = name),
                                 size = 4,
                                 box.padding = unit(0.1, 'lines'),
                                 point.padding = unit(0.1, 'lines'),
                                 segment.size = 0.5)## use the dataframe to plot points
      
    }
  })
  
  protein_input<-reactive({ 
    protein_selected  <- data_result()[input$contents_rows_selected,1]
    protein_selected <-as.character(protein_selected)
    if(length(levels(as.factor(colData(dep())$replicate))) <= 8){
      plot_protein(dep(), protein_selected, as.character(input$type))
    }
    else{
      protein_plot<-plot_protein(dep(), protein_selected, as.character(input$type))
      protein_plot + scale_color_brewer(palette = "Paired")
    }
    
  })
  
  
  ## QC Inputs
  norm_input <- reactive({
    plot_normalization(processed_data(),
                       normalised_data())
  })
  
  missval_input <- reactive({
    plot_missval(processed_data())
  })
  
  detect_input <- reactive({
    plot_detect(processed_data())
  })
  
  imputation_input <- reactive({
    plot_imputation(processed_data(),
                    diff_all())
  })
  
  p_hist_input <- reactive({
    plot_p_hist(dep())
  })
  
  numbers_input <- reactive({
    plot_numbers(processed_data()) +
      labs(title= "Phosphosites per sample", y = "Number of phosphosites")
  })
  
  coverage_input <- reactive({
    plot_coverage(processed_data())+
      labs(title= "Phosphosites per sample", y = "Number of phosphosites")
  })
  
  correlation_input<-reactive({
    plot_cor(dep(),significant = FALSE)
  })
  
  cvs_input<-reactive({
    plot_cvs(dep())
  })
  
  num_total<-reactive({
    dep() %>%
      nrow()
  }) 
  
  ## Enrichment inputs
  go_input<-eventReactive(input$go_analysis,{
    withProgress(message = 'Gene ontology enrichment is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    
    if(!is.null(input$contrast)){
      enrichment_output_test(dep(), input$go_database)
      go_results<- test_gsea_mod_phospho(dep(), databases = input$go_database, contrasts = TRUE)
      null_enrichment_test(go_results)
      if (input$go_database == "KEGG" | input$go_database == "Reactome"){
        plot_go<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast,
                                      databases=input$pathway_database, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
          xlab(NULL)
      }
      else{
        plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast,
                                  databases = input$go_database, nrow = 2, term_size = 8) + aes(stringr::str_wrap(Term, 60)) +
          xlab(NULL)
      }
      
      go_list<-list("go_result"=go_results, "plot_go"=plot_go)
      return(go_list)
    }
  })
  # 
  # pathway_input<-eventReactive(input$pathway_analysis,{
  #   progress_indicator("Pathway Analysis is running....")
  #   enrichment_output_test(dep(), input$pathway_database)
  #   pathway_results<- test_gsea_mod_phospho(dep(), databases=input$pathway_database, contrasts = TRUE)
  #   null_enrichment_test(pathway_results)
  #   plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast,
  #                                 databases=input$pathway_database, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
  #     xlab(NULL)
  #   pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
  #   return(pathway_list)
  # })
  
  KSEA_input<-eventReactive(input$KSEA_analysis,{
    progress_indicator("Kinase-Substrate Analysis is running....")
    
    result_df <- get_results_phospho(dep(),FALSE)
    print(input$contrast_1)  #test
    col_selected <- c('Protein','Gene.names','peptide.sequence', 'Residue.Both',
                      paste(input$contrast_1, "_p.val", sep = ""),
                      paste(input$contrast_1, "_log2 fold change", sep = ""))
    print(col_selected) #test
    
    # select required columns and rename them
    column_names <- c('Protein','Gene','Peptide','Residue.Both','p','FC')
    PX <- result_df %>% dplyr::select (col_selected)
    names(PX) <- column_names
    KSData <- KSEAapp::KSData 
    
    # Generate a summary bar plot using the KSEA.Barplot() function
    plot_KSEA <- KSEAapp::KSEA.Barplot(KSData, PX, NetworKIN=TRUE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.01, export=FALSE)
    return(plot_KSEA)
  })
  
  
  #### Interactive UI
  output$significantBox <- renderInfoBox({
    num_total <- dep() %>%
      nrow()
    num_signif <- dep() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
    frac <- num_signif / num_total
    
    info_box <- infoBox("Significant phosphosites",
                          paste0(num_signif,
                                 " out of ",
                                 num_total),
                          paste0(signif(frac * 100, digits = 3),
                                 "% of phosphosites differentially expressed across all conditions"),
                          icon = icon("stats", lib = "glyphicon"),
                          color = "olive",
                          # fill = TRUE,
                          width = 4)
    
    return(info_box)
  })
  
  ##### Get results dataframe from Summarizedexperiment object
  data_result<-reactive({
    if(length(unique(exp_design_input()$condition)) <= 2){
      get_results_phospho(dep(),FALSE) %>% dplyr::select (-Residue.Both,-peptide.sequence,-Protein)
    } else {
      get_results_phospho(dep(),TRUE) %>% dplyr::select (-Residue.Both,-peptide.sequence,-Protein)
    }
  })
  
  
  #### Data table
  output$contents <- DT::renderDataTable({
    withProgress(message = 'Result table calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    df<- data_result()
    return(df)
  },
  options = list(scrollX = TRUE,
                 autoWidth=TRUE,
                 columnDefs= list(list(width = '450px', targets = c(-1))))
  )
  
  ## Deselect all rows button
  proxy <- dataTableProxy("contents")
  
  observeEvent(input$clear,{
    proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$original,{
    output$contents <- DT::renderDataTable({
      df<- data_result()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '450px', targets = c(-1))))
    )
  })
  
  # Name brush function
  protein_name_brush<- reactive({
    #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
                               xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  protein_name_click<- reactive({
    protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
    #xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  
  brush <- NULL
  makeReactiveBinding("brush")
  
  observeEvent(input$protein_brush,{
    output$contents <- DT::renderDataTable({
      
      df<- data_result()[data_result()[["Phosphosite"]] %in% protein_name_brush(), ]
      return(df)
    },
    options = list(scrollX= TRUE)
    )
    
    proteins_selected<-data_result()[data_result()[["Phosphosite"]] %in% protein_name_brush(), ] ## get all rows selected
    ## convert contrast to x and padj to y
    diff_proteins <- grep(paste(input$volcano_cntrst, "_log2", sep = ""),
                          colnames(proteins_selected))
    if(input$p_adj=="FALSE"){
      padj_proteins <- grep(paste(input$volcano_cntrst, "_p.val", sep = ""),
                            colnames(proteins_selected))
    }
    else{
      padj_proteins <- grep(paste(input$volcano_cntrst, "_p.adj", sep = ""),
                            colnames(proteins_selected))
    }
    df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                             y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                             name = proteins_selected$Phosphosite)
    #print(df_protein)
    
    if(length(unique(exp_design_input()$condition)) <= 2){
      p<-plot_volcano_new(dep(),
                          input$volcano_cntrst,
                          FALSE,
                          input$check_names,
                          input$p_adj)
      
    } else {
      p<-plot_volcano_new(dep(),
                          input$volcano_cntrst,
                          input$check_anova,
                          input$check_names,
                          input$p_adj)
    }
    
    p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
      ggrepel::geom_text_repel(data = df_protein,
                               aes(x, y, label = name),
                               size = 4,
                               box.padding = unit(0.1, 'lines'),
                               point.padding = unit(0.1, 'lines'),
                               segment.size = 0.5)
    
    output$volcano <- renderPlot({
      withProgress(message = 'Volcano Plot calculations are in progress',
                   detail = 'Please wait for a while', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      p
    })
    return(p)
  })
  
  observeEvent(input$resetPlot,{
    session$resetBrush("protein_brush")
    brush <<- NULL
    
    output$contents <- DT::renderDataTable({
      df<- data_result()
      return(df)
    },
    options = list(scrollX= TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '450px', targets = c(-1))))
    )
    
    output$volcano <- renderPlot({
      volcano_input()
    })
  })
  
  ## Render Result Plots
  output$pca_plot<-renderPlot({
    pca_input()
  })
  
  output$heatmap<-renderPlot({
    withProgress(message = 'Heatmap rendering is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    heatmap_input()
  })
  
  output$volcano <- renderPlot({
    withProgress(message = 'Volcano Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_rows_selected)){
      volcano_input()
    }
    else if(!is.null(input$volcano_cntrst)){
      volcano_input_selected()
    }# else close
  })
  
  output$protein_plot<-renderPlot({
    if(!is.null(input$contents_rows_selected)){
      protein_input()
    }
  })
  
  
  ### QC Outputs
  output$sample_corr <-renderPlot({
    correlation_input()
  })
  
  output$sample_cvs <- renderPlot({
    cvs_input()
  })
  
  output$norm <- renderPlot({
    norm_input()
  })
  
  output$missval <- renderPlot({
    missval_input()
  })
  
  output$detect <- renderPlot({
    detect_input()
  })
  
  output$imputation <- renderPlot({
    imputation_input()
  })
  
  # output$p_hist <- renderPlot({
  #   p_hist_input()
  # })
  
  output$numbers <- renderPlot({
    numbers_input()
  })
  
  output$coverage <- renderPlot({
    coverage_input()
  })
  
  ## Enrichment Outputs
  output$go_enrichment<-renderPlot({
    go_input()$plot_go
  })
  
  output$KSEA_enrichment<-renderPlot({
    KSEA_input()
  })
  
  ##### Download Functions
  datasetInput <- reactive({
    switch(input$dataset,
           "Results" = get_results_phospho(dep()),
           "Original_matrix"= unimputed_table(),
           # "significant_proteins" = get_results(dep()) %>%
           #   filter(significant) %>%
           #   select(-significant),
           "Imputed_matrix" = imputed_table(),
           "Full_dataset" = get_df_wide(dep()))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$dataset, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(datasetInput(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ### === Cluster Download ==== ####
  
  individual_cluster <- reactive({
    cluster_number <- input$cluster_number
    cluster_all <- heatmap_input()
    data_result()[cluster_all[[cluster_number]],]
  })
  
  # output$text1 <- renderPrint({
  #   paste(individual_cluster())
  # })
  
  output$downloadCluster <- downloadHandler(
    filename = function() { paste("Cluster_info_",input$cluster_number, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(individual_cluster(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  output$downloadVolcano <- downloadHandler(
    filename = function() {
      paste0("Volcano_", input$volcano_cntrst, ".pdf")
    },
    content = function(file) {
      pdf(file)
      if(is.null(input$protein_brush)){
        print(volcano_input())
        dev.off()
      }
      else{
        observeEvent(input$protein_brush,{
          print(p)
        })
        print(volcano_input_selected())
        dev.off()
      }
    }
  )
  
  
  ## Protein plot download
  output$downloadProtein <- downloadHandler(
    filename = function() {
      paste0(input$type,".pdf")
    },
    content = function(file) {
      pdf(file)
      print(protein_input())
      dev.off()
    }
  )
  
  ###### ==== DOWNLOAD GO TABLE ==== ####
  output$downloadGO <- downloadHandler(
    filename = function() { paste("GO_enrichment_",input$go_database, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(go_input()$go_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
  output$downloadPA <- downloadHandler(
    filename = function() { paste("Pathway_enrichment_",input$pathway_database, ".csv", sep = "") }, 
    ## use = instead of <-
    content = function(file) {
      write.table(pathway_input()$pa_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  output$download_hm_svg<-downloadHandler(
    filename = function() { "heatmap.svg" }, 
    ## use = instead of <-
    content = function(file) {
      heatmap_plot<-DEP::plot_heatmap(dep(),"centered", k=6, indicate = "condition")
      svg(file)
      print(heatmap_plot)
      
      
      
      
      dev.off()
    }
  )
  
  #####===== Download Report (phosphosite)=====#####
  output$downloadReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "Phospho-Analyst(phosphosite)report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "Phosphosite_report.Rmd")
      file.copy("Phosphosite_report.Rmd", tempReport, overwrite = TRUE)
      
      sig_proteins<-dep() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()
      
      tested_contrasts<- gsub("_significant", "", 
                              colnames(SummarizedExperiment::rowData(dep()))[grep("_significant", 
                                                                                  colnames(SummarizedExperiment::rowData(dep())))])
      pg_width<- ncol(normalised_data()) / 2.5
      # Set up parameters to pass to Rmd document
      params <- list(data = processed_data,
                     alpha = input$p,
                     lfc = input$lfc,
                     num_signif= sig_proteins,
                     pg_width = pg_width,
                     tested_contrasts= tested_contrasts,
                     numbers_input= numbers_input,
                     detect_input = detect_input,
                     imputation_input = imputation_input,
                     missval_input = missval_input,
                     p_hist_input = p_hist_input,
                     pca_input = pca_input,
                     coverage_input= coverage_input,
                     correlation_input =correlation_input,
                     heatmap_input = heatmap_input,
                     cvs_input= cvs_input,
                     dep = dep
      )
      
      # Knit the document, passing in the `params` list
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  ###### ==== DOWNLOAD QC plots svg ==== ####
  
  output$download_pca_svg<-downloadHandler(
    filename = function() { "PCA_plot.svg" }, 
    content = function(file) {
      svg(file)
      print(pca_input())
      dev.off()
    }
  )
  
  output$download_corr_svg<-downloadHandler(
    filename = function() { "Correlation_plot.svg" }, 
    content = function(file) {
      svg(file)
      print(correlation_input())
      dev.off()
    }
  )
  
  output$download_cvs_svg<-downloadHandler(
    filename = function() { "Sample_CV.svg" }, 
    content = function(file) {
      svg(file)
      print(cvs_input())
      dev.off()
    }
  )
  
  output$download_num_svg<-downloadHandler(
    filename = function() { "Proteins_plot.svg" }, 
    content = function(file) {
      svg(file)
      print(numbers_input())
      dev.off()
    }
  )
  
  output$download_cov_svg<-downloadHandler(
    filename = function() { "Coverage_plot.svg" }, 
    content = function(file) {
      svg(file)
      print(coverage_input())
      dev.off()
    }
  )
  
  output$download_norm_svg<-downloadHandler(
    filename = function() { "Normalization_plot.svg" }, 
    content = function(file) {
      svg(file)
      print(norm_input())
      dev.off()
    }
  )
  
  output$download_missval_svg<-downloadHandler(
    filename = function() { "Missing_value_heatmap.svg" }, 
    content = function(file) {
      svg(file)
      print(missval_input())
      dev.off()
    }
  )
  
  output$download_imp_svg<-downloadHandler(
    filename = function() { "Imputation_plot.svg" }, 
    content = function(file) {
      svg(file)
      print(imputation_input())
      dev.off()
    }
  )
  
  
  #### Protein page logic ========== #############
  
  ####======= Render Functions
  
  output$volcano_cntrst_pr <- renderUI({
    if (!is.null(comparisons_pr())) {
      df <- SummarizedExperiment::rowData(dep_pr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("volcano_cntrst_pr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  ##comparisons
  output$contrast_pr <- renderUI({
    if (!is.null(comparisons_pr())) {
      df <- SummarizedExperiment::rowData(dep_pr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_pr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$contrast_1_pr <- renderUI({
    if (!is.null(comparisons_pr())) {
      df <- SummarizedExperiment::rowData(dep_pr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_1_pr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$downloadTable_pr <- renderUI({
    if(!is.null(dep_pr())){
      selectizeInput("dataset_pr",
                     "Download data table" ,
                     c("Results",
                       "Full dataset"))
    }
  })
  
  output$downloadButton_pr <- renderUI({
    if(!is.null(dep_pr())){
      downloadButton('downloadData_pr', 'Save')
    }
  })
  
  
  output$downloadreport_pr <- renderUI({
    if(!is.null(dep_pr())){
      downloadButton('downloadReport_pr', 'Download Report')
    }
  })
  
  output$downloadPlots_pr <- renderUI({
    if(!is.null(dep_pr())){
      downloadButton('downloadPlots1_pr', 'Download Plots')
    }
  })
  
  comparisons_pr<-reactive({
    temp<-capture.output(DEP::test_diff(imputed_data_pr(),type='all'),type = "message")
    gsub(".*: ","",temp)
    ## Split conditions into character vector
    unlist(strsplit(temp,","))
    ## Remove leading and trailing spaces
    trimws(temp)
  })
  
  
  ## Results plot inputs
  
  ## PCA Plot
  
  pca_label_pr<-reactive({
    pca_lable<-levels(as.factor(colData(dep_pr())$replicate))
    print(pca_label)
  })
  
  pca_input_pr<-eventReactive(input$analyze ,{
    if(input$analyze==0 ){
      return()
    }
    if (num_total_pr()<=500){
      if(length(levels(as.factor(colData(dep_pr())$replicate))) <= 6){
        pca_plot<- DEP::plot_pca(dep_pr(), n=num_total_pr(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
      else{
        pca_plot<-DEP::plot_pca(dep_pr(), n=num_total_pr(), point_size = 4, indicate = "condition")
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
    }
    else{
      if(length(levels(as.factor(colData(dep_pr())$replicate))) <= 6){
        pca_plot<-DEP::plot_pca(dep_pr(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }else{
        pca_label<-SummarizedExperiment::colData(dep_pr())$replicate
        pca_plot<-DEP::plot_pca(dep_pr(), point_size = 4, indicate = "condition")
        #pca_plot<-pca_plot + geom_point()
        pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                                      size = 4,
                                                      box.padding = unit(0.1, 'lines'),
                                                      point.padding = unit(0.1, 'lines'),
                                                      segment.size = 0.5)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
    }
    
  })
  
  ### Heatmap Differentially expressed proteins
  heatmap_input_pr<-eventReactive(input$analyze ,{
    if(input$analyze==0 ){
      return()
    }
    get_cluster_heatmap(dep_pr(),
                        type="centered",kmeans = TRUE,
                        k=6, col_limit = 6,
                        indicate = "condition"
    )
  })
  
  ### Volcano Plot
  volcano_input_pr <- reactive({
    if(!is.null(input$volcano_cntrst_pr)) {
      if(length(unique(exp_design_input_1()$condition)) <= 2){
        plot_volcano_new(dep_pr(),
                            input$volcano_cntrst_pr,
                            FALSE,
                            input$check_names_pr,
                            input$p_adj_pr)
        
      } else {
        plot_volcano_new(dep_pr(),
                         input$volcano_cntrst_pr,
                         input$check_anova_pr,
                         input$check_names_pr,
                         input$p_adj_pr)
      }
      # plot_volcano_new(dep_pr(),
      #                  input$volcano_cntrst_pr,
      #                  input$check_anova_pr,
      #                  input$check_names_pr,
      #                  input$p_adj_pr)
      # 
    }
  })
  
  volcano_df_pr<- reactive({
    if(!is.null(input$volcano_cntrst_pr)) {
      get_volcano_df(dep_pr(),
                     input$volcano_cntrst_pr)
      
    }
  })
  
  
  volcano_input_selected_pr<-reactive({
    if(!is.null(input$volcano_cntrst_pr)){
      if (!is.null(input$contents_nr_rows_selected)){
        proteins_selected<-data_result_pr()[c(input$contents_nr_rows_selected),]## get all rows selected
      }
      else if(!is.null(input$protein_brush_pr)){
        proteins_selected<-data_result_pr()[data_result_pr()[["Gene Name"]] %in% protein_name_brush_pr(), ] 
      }
      ## convert contrast to x and padj to y
      diff_proteins <- grep(paste(input$volcano_cntrst_pr, "_log2", sep = ""),
                            colnames(proteins_selected))
      if(input$p_adj_pr=="FALSE"){
        padj_proteins <- grep(paste(input$volcano_cntrst_pr, "_p.val", sep = ""),
                              colnames(proteins_selected))
      }
      else{
        padj_proteins <- grep(paste(input$volcano_cntrst_pr, "_p.adj", sep = ""),
                              colnames(proteins_selected))
      }
      
      df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                               y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                               name = proteins_selected$`Gene Name`)
      #print(df_protein)
      if(length(unique(exp_design_input_1()$condition)) <= 2){
        p<-plot_volcano_new(dep_pr(),
                         input$volcano_cntrst_pr,
                         FALSE,
                         input$check_names_pr,
                         input$p_adj_pr)
        
      } else {
        p<-plot_volcano_new(dep_pr(),
                         input$volcano_cntrst_pr,
                         input$check_anova_pr,
                         input$check_names_pr,
                         input$p_adj_pr)
      }
      # 
      # p<-plot_volcano_new(dep_pr(),
      #                     input$volcano_cntrst_pr,
      #                     input$check_anova_pr,
      #                     input$check_names_pr,
      #                     input$p_adj_pr)
      
      p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
        ggrepel::geom_text_repel(data = df_protein,
                                 aes(x, y, label = name),
                                 size = 4,
                                 box.padding = unit(0.1, 'lines'),
                                 point.padding = unit(0.1, 'lines'),
                                 segment.size = 0.5)## use the dataframe to plot points
      
    }
  })
  
  protein_input_pr<-reactive({ 
    
    protein_selected  <- data_result_pr()[input$contents_pr_rows_selected,1]
    
    #protein<-row_selected$name
    if(length(levels(as.factor(colData(dep_pr())$replicate))) <= 8){
      plot_protein(dep_pr(), protein_selected, input$type_pr)
    }
    else{
      protein_plot<-plot_protein(dep_pr(), protein_selected, input$type_pr)
      protein_plot + scale_color_brewer(palette = "Paired")
    }
    
  })
  
  ## Get processed data
  cleaned_data_pr<- reactive({
    ## check which dataset
    if(!is.null (protein_data_input() )){
      protein_data <- reactive({protein_data_input()})
    }
    
    if(grepl('+',protein_data()$Reverse)){
      filtered_data<-dplyr::filter(protein_data(),Reverse!="+")
    }
    else{filtered_data<-protein_data()}
    if(grepl('+',filtered_data$Potential.contaminant)){
      filtered_data<-dplyr::filter(filtered_data,Potential.contaminant!="+")
    }
    if(grepl('+',filtered_data$Only.identified.by.site)){
      filtered_data<-dplyr::filter(filtered_data,Only.identified.by.site!="+") 
    }
    if(input$single_peptide==TRUE){
      filtered_data <-filtered_data
    }
    else{filtered_data<-dplyr::filter(filtered_data,as.numeric(Razor...unique.peptides)>=2)}
    
    filtered_data<-ids_test(filtered_data)
    filtered_data$Gene.names <- filtered_data$Gene.names %>% toupper()
    return(filtered_data)
  })
  
  
  processed_data_pr<-reactive({
    if(!is.null (exp_design_input_1() )){
      exp_design<-reactive({exp_design_input_1()})
    }
    message(exp_design())
    
    data_pre <- cleaned_data_pr()
    data_unique<- DEP::make_unique(data_pre,"Gene.names","Protein.IDs",delim=";")
    lfq_columns<-grep("LFQ.", colnames(data_unique))
    
    ## Check for matching columns in maxquant and experiment design file
    test_match_lfq_column_design(data_unique,lfq_columns, exp_design())
    data_se<-DEP:::make_se(data_unique,lfq_columns,exp_design())
    
    # Check number of replicates
    if(max(exp_design()$replicate)<3){
      threshold<-0
    } else  if(max(exp_design()$replicate)==3){
      threshold<-1
    } else if(max(exp_design()$replicate)<6 ){
      threshold<-2
    } else if (max(exp_design()$replicate)>=6){
      threshold<-trunc(max(exp_design()$replicate)/2)
    }
    
    
    filter_missval(data_se,thr = threshold)
  })
  
  unimputed_table_pr<-reactive({
    temp<-assay(processed_data_pr())
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"original_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) 
    #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  normalised_data_pr<-reactive({
    print(assay(processed_data_pr()))
    DEP::normalize_vsn(processed_data_pr())
  })
  
  
  imputed_data_pr<-reactive({
    DEP::impute(processed_data_pr(),input$imputation)
  })
  
  imputed_table_pr<-reactive({
    temp<-assay(imputed_data_pr())
    #tibble::rownames_to_column(temp,var = "ProteinID")
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"imputed_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  
  diff_all_pr<-reactive({
    test_diff(imputed_data_pr(),type = 'all')
  })
  
  dep_pr<-reactive({
    if(input$fdr_correction=="BH"){
      diff_all<-test_limma(imputed_data_pr(),type='all', paired = input$paired)
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    else{
      diff_all<-test_diff(imputed_data_pr(),type='all')
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    
    if(length(unique(exp_design_input_1()$condition)) <= 2){
      return(diff_all_rej)
      
    }
    else if(length(unique(exp_design_input_1()$condition)) >= 3){
      anova_dep <- diff_all
      # get assay data
      intensity <- assay(anova_dep)
      exp_design <- exp_design_input_1()
      exp_design_rename<-exp_design
      exp_design_rename$label<-paste(exp_design_rename$condition, exp_design_rename$replicate, sep = "_")
      
      # reshape intensity columns
      data_reshape<-reshape2::melt(intensity,value.name = "intensity", variable.name = "label")
      colnames(data_reshape)<-c("uid", "label", "intensity")
      
      # Join the table
      data_experiment<-left_join(data_reshape, exp_design_rename, by="label")
      
      # apply anova function
      anova<-data_experiment %>%
        group_by(`uid`) %>%
        do(anova_function(.)) %>% dplyr::select(p.value) %>%
        ungroup()
      
      colnames(anova)<-c("name", "anova_p.val")
      
      # add anova p.value to row data
      rowData(anova_dep) <- merge(rowData(anova_dep), anova, by = 'name', sort = FALSE)
      anova_dep_rej <- DEP::add_rejections(anova_dep,alpha = input$p, lfc= input$lfc)
      
      # calculate adjusted anova p.value to data
      anova_adj <-anova
      anova_adj$anova_p.adj <- p.adjust(anova_adj$anova_p.val,method = input$fdr_correction)
      anova_adj <- anova_adj %>% select(-anova_p.val)

      # add adjusted anova p.value to row data
      rowData(anova_dep_rej) <- merge(rowData(anova_dep_rej), anova_adj, by = 'name', sort = FALSE)
      return(anova_dep_rej)
    }
  })
  
  ## QC Inputs
  norm_input_pr <- reactive({
    plot_normalization(processed_data_pr(),
                       normalised_data_pr())
  })
  
  missval_input_pr <- reactive({
    plot_missval(processed_data_pr())
  })
  
  detect_input_pr <- reactive({
    plot_detect(processed_data_pr())
  })
  
  imputation_input_pr <- reactive({
    plot_imputation(processed_data_pr(),
                    diff_all_pr())
  })
  
  p_hist_input_pr <- reactive({
    plot_p_hist(dep_pr())
  })
  
  numbers_input_pr <- reactive({
    plot_numbers(normalised_data_pr())
  })
  
  coverage_input_pr <- reactive({
    plot_coverage(normalised_data_pr())
  })
  
  correlation_input_pr<-reactive({
    plot_cor(dep_pr())
  })
  
  cvs_input_pr<-reactive({
    plot_cvs(dep_pr())
  })
  
  num_total_pr<-reactive({
    dep_pr() %>%
      nrow()
  }) 
  
  ## Enrichment inputs
  
  go_input_pr<-eventReactive(input$go_analysis_pr,{
    withProgress(message = 'Gene ontology enrichment is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    
    if(!is.null(input$contrast_pr)){
      enrichment_output_test(dep_pr(), input$go_database_pr)
      go_results<- test_gsea_mod(dep_pr(), databases = input$go_database_pr, contrasts = TRUE)
      plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast_pr,
                                databases = input$go_database_pr, nrow = 2, term_size = 8) + 
        aes(stringr::str_wrap(Term, 60)) +
        xlab(NULL)
      go_list<-list("go_result"=go_results, "plot_go"=plot_go)
      return(go_list)
    }
  })
  
  pathway_input_pr<-eventReactive(input$pathway_analysis_pr,{
    withProgress(message = 'Pathway enrichment is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    enrichment_output_test(dep_pr(), input$pathway_database_pr)
    pathway_results<- test_gsea_mod(dep_pr(), databases=input$pathway_database_pr, contrasts = TRUE)
    plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_1_pr,
                                  databases=input$pathway_database_pr, nrow = 3, term_size = 8) + 
      aes(stringr::str_wrap(Term, 30)) +
      xlab(NULL)
    pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
    return(pathway_list)
    
  })
  
  
  #### Interactive UI
  output$significantBox_pr <- renderInfoBox({
    num_total <- dep_pr() %>%
      nrow()
    num_signif <- dep_pr() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
    frac <- num_signif / num_total
    
    info_box <- 		infoBox("Significant proteins",
                          paste0(num_signif,
                                 " out of ",
                                 num_total),
                          paste0(signif(frac * 100, digits = 3),
                                 "% of proteins differentially expressed across all conditions"),
                          icon = icon("stats", lib = "glyphicon"),
                          color = "olive",
                          # fill = TRUE,
                          width = 4)
    
    return(info_box)
  })
  
  
  ##### Get results dataframe from Summarizedexperiment object
  data_result_pr<-reactive({
    if(length(unique(exp_design_input_1()$condition)) <= 2){
      get_results_proteins(dep_pr(),FALSE)
    } else {
      get_results_proteins(dep_pr(),TRUE)
    }
  })
  
  
  #### Data table
  output$contents_pr <- DT::renderDataTable({
    withProgress(message = 'Result table calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    df<- data_result_pr()
    return(df)
  },
  options = list(scrollX = TRUE,
                 autoWidth=TRUE,
                 columnDefs= list(list(width = '400px', targets = c(-1))))
  )
  
  ## Deselect all rows button
  proxy <- dataTableProxy("contents_pr")
  
  observeEvent(input$clear_pr,{
    proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$original_pr,{
    output$contents_pr <- DT::renderDataTable({
      df<- data_result_pr()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })
  
  protein_name_brush_pr<- reactive({
    #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    protein_tmp<-brushedPoints(volcano_df_pr(), input$protein_brush_pr, 
                               xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  protein_name_click_pr<- reactive({
    protein_tmp<-nearPoints(volcano_df_pr(), input$protein_click_pr, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
    #xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  
  
  ## Select rows dynamically
  observeEvent(input$protein_brush_pr,{
    output$contents_pr <- DT::renderDataTable({
      df<- data_result_pr()[data_result_pr()[["Gene Name"]] %in% protein_name_brush_pr(), ]
      return(df)
    },
    options = list(scrollX= TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    proteins_selected<-data_result_pr()[data_result_pr()[["Gene Name"]] %in% protein_name_brush_pr(), ] #
    # get all rows selected
    ## convert contrast to x and padj to y
    diff_proteins <- grep(paste(input$volcano_cntrst_pr, "_log2", sep = ""),
                          colnames(proteins_selected))
    if(input$p_adj=="FALSE"){
      padj_proteins <- grep(paste(input$volcano_cntrst_pr, "_p.val", sep = ""),
                            colnames(proteins_selected))
    }
    else{
      padj_proteins <- grep(paste(input$volcano_cntrst_pr, "_p.adj", sep = ""),
                            colnames(proteins_selected))
    }
    df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                             y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                             name = proteins_selected$`Gene Name`)
    #print(df_protein)
    if(length(unique(exp_design_input_1()$condition)) <= 2) {
      p<-plot_volcano_new(dep_pr(),
                          input$volcano_cntrst_pr,
                          FALSE,
                          input$check_names_pr,
                          input$p_adj_pr)
    } else {
      p<-plot_volcano_new(dep_pr(),
                          input$volcano_cntrst_pr,
                          input$check_anova_pr,
                          input$check_names_pr,
                          input$p_adj_pr)
    }
    
    # 
    # p<-plot_volcano_new(dep_pr(),
    #                     input$volcano_cntrst_pr,
    #                     input$check_anova_pr,
    #                     input$check_names_pr,
    #                     input$p_adj_pr)
    
    p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
      ggrepel::geom_text_repel(data = df_protein,
                               aes(x, y, label = name),
                               size = 4,
                               box.padding = unit(0.1, 'lines'),
                               point.padding = unit(0.1, 'lines'),
                               segment.size = 0.5)
    
    output$volcano_pr <- renderPlot({
      withProgress(message = 'Volcano Plot calculations are in progress',
                   detail = 'Please wait for a while', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      p
    })
    return(p)
  })
  
  observeEvent(input$resetPlot_pr,{
    session$resetBrush("protein_brush_pr")
    brush <<- NULL
    
    output$contents_pr <- DT::renderDataTable({
      df<- data_result_pr()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    output$volcano_pr <- renderPlot({
      volcano_input_pr()
    })
  })
  
  observeEvent(input$protein_click_pr,{
    output$contents_pr <- DT::renderDataTable({
      df<- data_result_pr()[data_result_pr()[["Gene Name"]] %in% protein_name_click_pr(), ]
      return(df)
    },
    options = list(scrollX= TRUE)
    )
  })
  
  ## Render Result Plots
  output$pca_plot_pr<-renderPlot({
    pca_input_pr()
  })
  output$heatmap_pr<-renderPlot({
    withProgress(message = 'Heatmap rendering is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    heatmap_input_pr()
  })
  
  output$volcano_pr <- renderPlot({
    withProgress(message = 'Volcano Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_pr_rows_selected)){
      volcano_input_pr()
    }
    else if(!is.null(input$volcano_cntrst_pr)){
      volcano_input_selected_pr()
    } # else close
  })
  
  output$protein_plot_pr<-renderPlot({
    if(!is.null(input$contents_pr_rows_selected)){
      protein_input_pr()
    }
  })
  
  
  ### QC Outputs
  output$sample_corr_pr <-renderPlot({
    correlation_input_pr()
  })
  
  output$sample_cvs_pr <- renderPlot({
    cvs_input_pr()
  })
  
  output$norm_pr <- renderPlot({
    norm_input_pr()
  })
  
  output$missval_pr <- renderPlot({
    missval_input_pr()
  })
  
  output$detect_pr <- renderPlot({
    detect_input_pr()
  })
  
  output$imputation_pr <- renderPlot({
    imputation_input_pr()
  })
  
  # output$p_hist <- renderPlot({
  #   p_hist_input_pr()
  # })
  
  output$numbers_pr <- renderPlot({
    numbers_input_pr()
  })
  
  output$coverage_pr <- renderPlot({
    coverage_input_pr()
  })
  
  ## Enrichment Outputs
  output$go_enrichment_pr<-renderPlot({
    go_input_pr()$plot_go
  })
  
  output$pathway_enrichment_pr<-renderPlot({
    pathway_input_pr()$plot_pa
  })
  
  ##### Download Functions
  datasetInput_pr <- reactive({
    switch(input$dataset_pr,
           "Results" = get_results_proteins(dep_pr()),
           "Full dataset" = get_df_wide(dep_pr()))
  })
  
  output$downloadData_pr <- downloadHandler(
    filename = function() { paste(input$dataset_pr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(datasetInput_pr(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ### === Cluster Download ==== ####
  
  individual_cluster_pr <- reactive({
    cluster_number <- input$cluster_number_pr
    cluster_all <- heatmap_input_pr()
    data_result_pr()[cluster_all[[cluster_number]],]
  })
  
  
  
  output$downloadCluster_pr <- downloadHandler(
    filename = function() { paste("Cluster_info_",input$cluster_number_pr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(individual_cluster_pr(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  output$downloadVolcano_pr <- downloadHandler(
    filename = function() {
      paste0("Volcano_", input$volcano_cntrst_pr, ".pdf")
    },
    content = function(file) {
      pdf(file)
      print(volcano_input_selected_pr())
      dev.off()
    }
  )
  
  
  ## Protein plot download
  output$downloadProtein_pr <- downloadHandler(
    filename = function() {
      paste0(input$type_pr,".pdf")
    },
    content = function(file) {
      pdf(file)
      print(protein_input_pr())
      dev.off()
    }
  )
  
  ###### ==== DOWNLOAD GO TABLE ==== ####
  output$downloadGO_pr <- downloadHandler(
    filename = function() { paste("GO_enrichment_",input$go_database_pr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(go_input_pr()$go_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
  output$downloadPA_pr <- downloadHandler(
    filename = function() { paste("Pathway_enrichment_",input$pathway_database_pr, ".csv", sep = "") }, 
    ## use = instead of <-
    content = function(file) {
      write.table(pathway_input_pr()$pa_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  
  #####===== Download Report (proteinGroup)=====##### 
  output$downloadReport_pr <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "Phospho-Analyst(proteinGroup)report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "ProteinGroup_report.Rmd")
      file.copy("ProteinGroup_report.Rmd", tempReport, overwrite = TRUE)
      
      sig_proteins<-dep_pr() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()
      
      tested_contrasts<- gsub("_significant", "", 
                              colnames(SummarizedExperiment::rowData(dep_pr()))[grep("_significant", 
                                                                                  colnames(SummarizedExperiment::rowData(dep_pr())))])
      pg_width<- ncol(normalised_data_pr()) / 2.5
      # Set up parameters to pass to Rmd document
      params <- list(data = processed_data_pr,
                     alpha = input$p,
                     lfc = input$lfc,
                     num_signif= sig_proteins,
                     pg_width = pg_width,
                     tested_contrasts= tested_contrasts,
                     numbers_input= numbers_input_pr,
                     detect_input = detect_input_pr,
                     imputation_input = imputation_input_pr,
                     missval_input = missval_input_pr,
                     p_hist_input = p_hist_input_pr,
                     pca_input = pca_input_pr,
                     coverage_input= coverage_input_pr,
                     correlation_input =correlation_input_pr,
                     heatmap_input = heatmap_input_pr,
                     cvs_input= cvs_input_pr,
                     dep = dep_pr
      )
      
      # Knit the document, passing in the `params` list
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  #### Comparison page logic ========== #############
  ####======= Render Functions
  
  output$volcano_comp <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("volcano_comp",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$selected_gene = renderUI({
    if (!is.null(gene_names())){
      df <- gene_names()
      selectInput('selected_gene', 
                  label='Select one or more Gene names', 
                  choices = as.list(df$Gene.names),
                  selected = NULL,
                  multiple = TRUE) 
    }
  })
  
  # Reactive components
  phospho_df <- reactive({
    if (!is.null(input$volcano_comp)){
      # phospho_row <- rowData(dep()) %>% as.data.frame()
      phospho_row <- data_result() %>% mutate(rowname = Phosphosite) %>% as.data.frame()
      phospho_row <- column_to_rownames(phospho_row, 'rowname')
      phospho_intensity <- assay(dep())  %>% as.data.frame()
      phospho_df <- merge(phospho_row, phospho_intensity, by = 0) # Merge data according to row names
      if(length(unique(exp_design_input()$condition)) <= 2) {
        col_selected <- c(colnames(phospho_intensity),'Phosphosite','Gene.names',
                          paste(input$volcano_comp, "_log2 fold change", sep = ""),
                          paste(input$volcano_comp, "_p.val", sep = ""))
      } else {
        if (input$check_anova_comp =="FALSE") {
          col_selected <- c(colnames(phospho_intensity),'Phosphosite','Gene.names',
                            paste(input$volcano_comp, "_log2 fold change", sep = ""),
                            paste(input$volcano_comp, "_p.val", sep = ""))
        } else {
          col_selected <- c(colnames(phospho_intensity),'Phosphosite','Gene.names',
                            paste(input$volcano_comp, "_log2 fold change", sep = ""),
                            "ANOVA_p.val")
        }
      }
      phospho_df_1 <- phospho_df %>% 
        subset(select = col_selected) %>% 
        dplyr::rename(phospho_id = Phosphosite, phospho_diff = paste(input$volcano_comp, "_log2 fold change", sep = ""))
      return(phospho_df_1)
      print(head(phospho_df_1)) # test
    }
  })
  
  protein_df <- reactive({
    if (!is.null(input$volcano_comp)){
      # protein_row <- rowData(dep_pr()) %>% as.data.frame()
      data_df <- data_result_pr()
      colnames(data_df)[1]<-c("Gene.names")
      protein_row <- data_df %>% mutate(rowname = Gene.names) %>% as.data.frame()
      protein_row <- column_to_rownames(protein_row, 'rowname')
      protein_intensity <- assay(dep_pr()) %>% as.data.frame()
      protein_df <- merge(protein_row, protein_intensity, by = 0) # Merge data according to row names
      print(colnames(protein_df))
      if(length(unique(exp_design_input_1()$condition)) <= 2) {
        col_selected <- c(colnames(protein_intensity),"Protein ID",'Gene.names',
                          paste(input$volcano_comp, "_log2 fold change", sep = ""),
                          paste(input$volcano_comp, "_p.val", sep = ""))
      } else {
        if (input$check_anova_comp =="FALSE") {
          col_selected <- c(colnames(protein_intensity),"Protein ID",'Gene.names',
                            paste(input$volcano_comp, "_log2 fold change", sep = ""),
                            paste(input$volcano_comp, "_p.val", sep = ""))
        } else {
          col_selected <- c(colnames(protein_intensity),"Protein ID",'Gene.names',
                            paste(input$volcano_comp, "_log2 fold change", sep = ""),
                            "ANOVA_p.val")
        }
      }
      protein_df_1 <- protein_df %>% 
        subset(select = col_selected) %>% 
        dplyr::rename(protein_diff = paste(input$volcano_comp, "_log2 fold change", sep = "") )
      return(protein_df_1)
      cat(head(protein_df_1)) # test
    }
    
  })
  combined_df <- reactive({
    if (!is.null(phospho_df()) & !is.null(protein_df())){
      df <- phospho_df() %>%
        left_join(., protein_df(), by = "Gene.names")
      # df$protein_diff[is.na(df$protein_diff)] <- 0
      df$normalized_diff <- df$phospho_diff - df$protein_diff
      # get index of the p.val
      pval <- grep("_p.val.x",colnames(df))
      df$p_values <- as.numeric(df[, pval])
      
      df <- df %>%
        mutate(p_value_desc = case_when(phospho_diff > 1 & p_values < 0.05 ~ 'Up',
                                        phospho_diff < -1 & p_values < 0.05 ~ 'Down',
                                        TRUE ~ 'Not Sig'))
      df <- df %>% 
        mutate(n_p_value_desc = case_when(normalized_diff > 1 & p_values < 0.05 ~ 'Up',
                                          normalized_diff < -1 & p_values < 0.05 ~ 'Down',
                                          TRUE ~ 'Not Sig'))
      return(df)
      cat(head(df)) # test
    }
  })
  
  # gene names for selection input
  gene_names <- reactive({
    if (!is.null(phospho_df_long()) & !is.null(protein_df_long())){
      gene_names_list <- phospho_df_long() %>%
        inner_join(., protein_df_long(), by = "Gene.names") %>% 
        select('Gene.names') %>%
        unique()
      return(gene_names_list)
    }
  })
  
  # phosphosite and protein data 
  phospho_df_long <- reactive({
    if (!is.null(phospho_df()) & !is.null(protein_df())){
      combined_df <- phospho_df() %>%
        inner_join(., protein_df(), by = "Gene.names")
      exp_design <- exp_design_input()
      phospho_cols <- colnames(combined_df[grep('.x', colnames(combined_df))])
      phospho_cols_1 <- exp_design$label
      phospho_df_11 <- subset(combined_df, select = c(phospho_cols,'phospho_id','Gene.names','phospho_diff' ))
      names(phospho_df_11) <- c(phospho_cols_1,"p.val",'phospho_id','Gene.names','phospho_diff')
      phospho_df_22 <- phospho_df_11 %>% rownames_to_column() %>% 
        gather(label, intensity, -rowname,-"p.val",-"phospho_id", -"Gene.names", -"phospho_diff")
      phospho_df_33 <- phospho_df_22 %>%
        left_join(., exp_design, by = "label")
      return(phospho_df_33)
      cat(head(phospho_df_33)) # test
    }
  })
  
  protein_df_long <- reactive({
    if (!is.null(phospho_df()) & !is.null(protein_df())){
      combined_df <- phospho_df() %>%
        inner_join(., protein_df(), by = "Gene.names")
      exp_design <- exp_design_input_1()
      protein_cols <- colnames(combined_df[grep('.y', colnames(combined_df))])
      protein_cols_1 <- exp_design$label
      protein_df_11 <- subset(combined_df, select = c(protein_cols,"Protein ID",'Gene.names','protein_diff'))
      names(protein_df_11) <- c(protein_cols_1,"p.val","protein_id",'Gene.names','protein_diff')
      protein_df_22 <- protein_df_11 %>% rownames_to_column() %>% 
        gather(label, intensity, -rowname,-"p.val", -"protein_id", -"Gene.names",-"protein_diff")
      protein_df_33 <- protein_df_22 %>%
        left_join(., exp_design, by = "label")
      return(protein_df_33)
      cat(head(protein_df_33)) # test
    }
  })
  
  # interactive plots
  combined_inter <- reactive({
    if (!is.null(phospho_df_long()) & !is.null(protein_df_long())){
      if (is.null(input$selected_gene)){
        phospho_df <- phospho_df_long() %>% dplyr::filter(Gene.names == gene_names()$Gene.names[1])
        protein_df <- protein_df_long() %>% dplyr::filter(Gene.names == gene_names()$Gene.names[1])
      }
      else {
        phospho_df <- phospho_df_long() %>% dplyr::filter(Gene.names %in% input$selected_gene)
        protein_df <- protein_df_long() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      }
      
      # phospho_df <- phospho_df_long() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      # protein_df <- protein_df_long() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      
      p1 <- phospho_df %>% 
        unique() %>%
        ggplot(aes(x = condition, y = intensity)) +
        geom_point(aes(color = factor(replicate)),
                   size = 3) +
        geom_line(aes(group= factor(replicate), color= factor(replicate))) +
        scale_colour_discrete(name  ="Replicate") + 
        ylim(19,28) + labs(title = 'Phospho site', x = '') +
        facet_grid(. ~phospho_id + Gene.names) +
        theme(axis.text.x = element_text(angle = 45)) 
      
      p2 <- protein_df %>% 
        select(-'rowname') %>%
        unique() %>%
        ggplot(aes(x = condition, y = intensity)) +
        geom_point(aes(color = factor(replicate)),
                   size = 3) +
        geom_line(aes(group= factor(replicate), color= factor(replicate))) +
        scale_colour_discrete(name  ="Replicate") + 
        ylim(19,28) + labs(title = 'Protein', x = '') +
        facet_grid(. ~protein_id + Gene.names) +
        theme(axis.text.x = element_text(angle = 45)) 
      
      ggarrange(p2, 
                p1 + 
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank()), 
                widths=c(1,4), common.legend = TRUE, legend = 'right') 
    }
  })
  
  # # volcano plots input
  # volcano_phospho <- reactive({
  #   df <- combined_df() 
  #   print(colnames(df))  # test
  #   if(length(unique(exp_design_input()$condition)) <= 2) {
  #     pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
  #   } else {
  #     if (input$check_anova_comp =="FALSE") {
  #       pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
  #     } else {
  #       pval <- grep("ANOVA_p.val.x",colnames(df))
  #     }
  #   }
  #   print(pval) # test
  #   name1 <- gsub("_vs_.*", "", input$volcano_comp)
  #   name2 <- gsub(".*_vs_", "", input$volcano_comp)
  #   
  #   df <- df %>% dplyr::mutate(p_values = as.numeric(df[, pval]))
  #   
  #   df %>% ggplot(aes(x = phospho_diff, y = -log10(p_values))) +
  #     geom_point(aes(color = p_value_desc)) +
  #     geom_text(data = data.frame(), aes(x = c(Inf, -Inf),
  #                                        y = c(-Inf, -Inf),
  #                                        hjust = c(1, 0),
  #                                        vjust = c(-1, -1),
  #                                        label = c(name1, name2),
  #                                        size = 5,
  #                                        fontface = "bold")) +
  #     labs(x = 'Phosphosite log fold change',
  #          y = '-log10(p-value)') +
  #     geom_text_repel(
  #       data=df %>% filter(p_value_desc != 'Not Sig'), # Filter data first
  #       aes(label=phospho_id),
  #       nudge_x = 0.5, nudge_y = 0,
  #       size = 4) +
  #     theme(plot.title = element_text(hjust = 0.5)) +
  #     scale_color_manual(values = c("#003399", "#999999", "#b30000")) 
  # })
  # 
  # volcano_phospho_2 <- reactive({
  #   df <- combined_df() 
  #   if(length(unique(exp_design_input()$condition)) <= 2) {
  #     pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
  #   } else {
  #     if (input$check_anova_comp =="FALSE") {
  #       pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
  #     } else {
  #       pval <- grep("ANOVA_p.val.x",colnames(df))
  #     }
  #   }
  #   
  #   name1 <- gsub("_vs_.*", "", input$volcano_comp)
  #   name2 <- gsub(".*_vs_", "", input$volcano_comp)
  #   
  #   df$p_values <- as.numeric(df[, pval])
  #   
  #   df %>% 
  #     # filter(!is.na(protein_diff)) %>%
  #     ggplot(aes(x = normalized_diff, y = -log10(p_values))) +
  #     geom_point(aes(color = n_p_value_desc)) +
  #     geom_text(data = data.frame(), aes(x = c(Inf, -Inf),
  #                                        y = c(-Inf, -Inf),
  #                                        hjust = c(1, 0),
  #                                        vjust = c(-1, -1),
  #                                        label = c(name1, name2),
  #                                        size = 5,
  #                                        fontface = "bold")) +
  #     labs(x = 'Phosphosite log fold change',
  #          y = '-log10(p-value)') +
  #     geom_text_repel(
  #       data=df %>% filter(n_p_value_desc != 'Not Sig'), # Filter data first
  #       aes(label=phospho_id),
  #       nudge_x = 0.5, nudge_y = 0,
  #       size = 4)+
  #     scale_color_manual(values = c("#003399", "#999999", "#b30000"))
  # })
  
  # Phosphosite and protein log fold change scatter plot input
  scatter_plot <- reactive({
    if(!is.null(input$volcano_comp) & !is.null(combined_df())){
      df <- combined_df()
      df %>% filter(!is.na(protein_diff))  %>% 
        ggplot(aes(x=phospho_diff, y=protein_diff)) + 
        geom_point(size = 2, alpha = 0.8) +
        geom_text_repel( 
          data=df %>% filter(phospho_diff > 5 | phospho_diff < -5|
                               protein_diff>2| protein_diff < -2), # Filter data first
          aes(label=phospho_id),
          nudge_x = 0.5, nudge_y = 0,
          size = 4) + 
        labs(title = paste(input$volcano_comp,'Comparison between phospho and protein log fold change', sep = "\n"),
             x = 'Phosphosite log fold change', y = 'Protein log fold change') +
        theme(plot.title = element_text(hjust = 0.5))
    }
  })
  
  # Output plots
  # output$volcano_phospho <-renderPlot({
  #   withProgress(message = 'Volcano Plot calculations are in progress',
  #                detail = 'Please wait for a while', value = 0, {
  #                  for (i in 1:15) {
  #                    incProgress(1/15)
  #                    Sys.sleep(0.25)
  #                  }
  #                })
  #   volcano_phospho()
  # })
  # output$volcano_phospho_2 <-renderPlot({
  #   withProgress(message = 'Volcano Plot calculations are in progress',
  #                detail = 'Please wait for a while', value = 0, {
  #                  for (i in 1:15) {
  #                    incProgress(1/15)
  #                    Sys.sleep(0.25)
  #                  }
  #                })
  #   volcano_phospho_2()
  # })
  
  
  # Output combined QC plots
  output$pca_plot_c <- renderPlot({
    ggarrange(pca_input() + labs(title = 'Phospho'), 
              pca_input_pr() + labs(title = "Protein"), 
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })
  
  output$scatter_plot <- renderPlot({
    scatter_plot()
  })
  
  output$sample_corr_c1 <- renderPlot({
    correlation_input()
  })
  
  output$sample_corr_c2 <- renderPlot({
    correlation_input_pr()
  })
  
  output$sample_cvs_c <- renderPlot({
    ggarrange(cvs_input() + labs(title = 'Phospho'), 
              cvs_input_pr() + labs(title = "Protein"),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })    
  
  output$numbers_c <- renderPlot({
    ggarrange(numbers_input() + labs(title = 'Phospho'), 
              numbers_input_pr() + labs(title = "Protein"),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })
  
  output$coverage_c <- renderPlot({
    ggarrange(coverage_input() + labs(title = 'Phospho'), 
              coverage_input_pr() + labs(title = "Protein")) 
  })
  
  output$norm_c <- renderPlot({
    ggarrange(norm_input() + labs(title = 'Phospho'), 
              norm_input_pr() + labs(title = "Protein") ,
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })
  
  
  output$missval_c1 <- renderPlot({
    missval_input()
  })
  
  output$missval_c2 <- renderPlot({
    missval_input_pr()
  })
  
  output$imputation_c <- renderPlot({
    ggarrange(imputation_input() + labs(title = 'Phospho'), 
              imputation_input_pr() + labs(title = "Protein"),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })
  
  # Output interactive plots
  output$combined_inter <- renderPlot({
    combined_inter()
  })
  
  output$combined_point <- renderPlot({
    if (!is.null(phospho_df())){
      exp_design<- exp_design_input()
      conditions <- exp_design$condition %>% unique()
      # phospho_df <- phospho_df() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      
      if (is.null(input$selected_gene)){
        phospho_df <- phospho_df() %>% dplyr::filter(Gene.names == gene_names()$Gene.names[1])
      }
      else {
        phospho_df <- phospho_df() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      }
      phospho_df_1 <- phospho_df
      print(colnames(phospho_df_1)) # test
      for (i in 1:length(conditions)) {
        condition <- conditions[i]
        pattern <- paste(condition,"[[:digit:]]",sep = '_')
        phospho_df_1[paste0('mean',sep = "_",condition)] <- rowMeans(phospho_df_1 %>% select(grep(pattern, colnames(phospho_df_1))), na.rm = TRUE)
      }
      
      phospho_df_2 <- phospho_df_1 %>% 
        # dplyr::filter(Gene.names %in% input$selected_gene) %>% 
        select(Gene.names,phospho_id, grep('mean', colnames(phospho_df_1))) %>%
        pivot_longer(names_to = "group", values_to = "mean_intensity", cols = starts_with('mean'))
      
      phospho_df_2 %>% ggplot(aes(x = phospho_id, y = group)) +
        geom_point(aes(size = mean_intensity, color = group)) +
        facet_grid(rows = vars(Gene.names), scales = "free_y", space = "free_y") +
        coord_flip() +
        theme(axis.text.x = element_text(angle = 45)) 
    }
  })

  
  #### Normalized page logic ========== #############
  ####======= Render Functions
  
  output$volcano_cntrst_nr <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep_nr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("volcano_cntrst_nr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  ##comparisons
  output$contrast_nr <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep_nr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_nr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$contrast_1_nr <- renderUI({
    if (!is.null(comparisons())) {
      df <- SummarizedExperiment::rowData(dep_nr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_1_nr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$downloadTable_nr <- renderUI({
    if(!is.null(dep_nr())){
      selectizeInput("dataset_nr",
                     "Choose a dataset to save" ,
                     c("Results","Original_matrix",
                       "Imputed_matrix",
                       "Full_dataset"))
    }
  })
  
  output$downloadButton_nr <- renderUI({
    if(!is.null(dep_nr())){
      downloadButton('downloadData_nr', 'Save')
    }
  })
  
  output$downloadZip_nr <- renderUI({
    if(!is.null(dep_nr())){
      downloadButton('downloadZip1_nr', 'Download result plots')
    }
  })
  output$downloadreport_nr <- renderUI({
    if(!is.null(dep_nr())){
      downloadButton('downloadReport_nr', 'Download Report')
    }
  })
  
  output$downloadPlots <- renderUI({
    if(!is.null(dep_nr())){
      downloadButton('downloadPlots1_nr', 'Download Plots')
    }
  })
  
  # subtraction on raw data and pre-processing
  normalized_phospho_data <- reactive({
    if(!is.null (exp_design_input() )){
      exp_design<-reactive({exp_design_input()})
    }
    phospho_pre <- cleaned_data()
    protein_pre <- cleaned_data_pr()
    protein_pre <-  protein_pre %>% 
      select("Majority.protein.IDs", grep("LFQ.", colnames(protein_pre))) # select the intensity columns
    
    # save phospho intensity colnames for use
    intensity_names <- colnames(phospho_pre)[grep("Intensity.", colnames(phospho_pre))]
    # change protein intensity colnames 
    search_protein <- paste("LFQ.intensity.", exp_design()$label, sep = "") %>% unique()
    replace_protein <- paste(exp_design()$condition, exp_design()$replicate,sep = "_") %>% unique()
    colnames(protein_pre)[colnames(protein_pre) %in% search_protein] <- replace_protein[match(colnames(protein_pre), search_protein, nomatch = 0)]
    
    # change phospho intensity colnames
    search_phospho <- paste("Intensity.", exp_design()$label, sep = "") %>% unique()
    replace_phospho <- paste("Intensity",exp_design()$condition, exp_design()$replicate,sep = "_") %>% unique()
    colnames(phospho_pre)[colnames(phospho_pre) %in% search_phospho] <- replace_phospho[match(colnames(phospho_pre), search_phospho, nomatch = 0)]
    
    # get condition names
    conditions <- exp_design()$condition %>% unique()
    conditions
    
    # find median value in protein intensity
    for (i in 1:length(conditions)) {
      condition <- conditions[i]
      pattern <- paste(condition,"[[:digit:]]",sep = '_')
      protein_pre[paste0('median',sep = "_",condition)] <- rowMedians(as.matrix(protein_pre %>% select(grep(pattern, colnames(protein_pre)))), na.rm = TRUE)
    }
    
    protein_median <- protein_pre %>% 
      select("Majority.protein.IDs",grep("median", colnames(protein_pre)))
    
    # join two raw data 
    phospho_protein <- phospho_pre %>% left_join(., protein_median, by = c("Protein" = "Majority.protein.IDs"))
    # change NA median protein value to zero
    phospho_protein <- mutate_at(phospho_protein, grep("median_", colnames(phospho_protein)), ~replace(., is.na(.), 0))
    
    # use each phosphosite intensity value to subtract median protein intensity of a same group
    phospho_protein_1 <- phospho_protein
    for (i in 1:length(conditions)){
      condition <- conditions[i]
      intensity_col <- grep(paste0("Intensity", sep = "_",condition), colnames(phospho_protein_1)) %>% as.vector()
      median_col <- grep(paste0('median',sep = "_",condition), colnames(phospho_protein_1)) 
      for (j in intensity_col){
        phospho_protein_1[j] <- phospho_protein_1[j] - phospho_protein_1[median_col]
      }
    }
    
    # remove median columns
    phospho_protein_2 <- phospho_protein_1 %>% 
      select(-grep("median_", colnames(phospho_protein_1)))
    
    # rename phospho intensity names
    colnames(phospho_protein_2)[grep("Intensity_", colnames(phospho_protein_2))] <- intensity_names
    
    # Convert dataframe to SE object
    data_unique_names <- DEP::make_unique(phospho_protein_2, 'name','ID', delim = ";")
    intensity_ints <- grep("^Intensity.", colnames(data_unique_names))
    data_se <- make_se(data_unique_names, intensity_ints, exp_design())
    
    # Check number of replicates
    if(max(exp_design()$replicate)<3){
      threshold<-0
    } else  if(max(exp_design()$replicate)==3){
      threshold<-1
    } else if(max(exp_design()$replicate)<6 ){
      threshold<-2
    } else if (max(exp_design()$replicate)>=6){
      threshold<-trunc(max(exp_design()$replicate)/2)
    }
    
    filter_missval(data_se,thr = threshold)
    
  })
  
  unimputed_table_nr<-reactive({
    temp<-assay(normalized_phospho_data())
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"original_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) 
    #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  imputed_data_nr<-reactive({
    DEP::impute(normalized_phospho_data(),input$imputation)
  })
  
  normalised_data_nr<-reactive({
    normalize_vsn(imputed_data_nr())
  })
  
  imputed_table_nr<-reactive({
    temp<-assay(imputed_data_nr())
    #tibble::rownames_to_column(temp,var = "ProteinID")
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"imputed_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  diff_all_nr<-reactive({
    test_diff(normalised_data_nr(),type = 'all')
  })
  
  dep_nr<-reactive({
    if(input$fdr_correction=="BH"){
      diff_all<-test_limma(normalised_data_nr(),type='all', paired = input$paired)
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    else{
      diff_all<-test_diff(normalised_data_nr(),type='all')
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    
    if(length(unique(exp_design_input()$condition)) <= 2){
      return(diff_all_rej)
      
    }
    else if(length(unique(exp_design_input()$condition)) >= 3){
      anova_dep <- diff_all
      # get assay data
      intensity <- assay(anova_dep)
      exp_design <- exp_design_input()
      exp_design_rename<-exp_design
      exp_design_rename$label<-paste(exp_design_rename$condition, exp_design_rename$replicate, sep = "_")
      
      # reshape intensity columns
      data_reshape<-reshape2::melt(intensity,value.name = "intensity", variable.name = "label")
      colnames(data_reshape)<-c("uid", "label", "intensity")
      
      # Join the table
      data_experiment<-left_join(data_reshape, exp_design_rename, by="label")
      
      # apply anova function
      anova<-data_experiment %>%
        group_by(`uid`) %>%
        do(anova_function(.)) %>% dplyr::select(p.value) %>%
        ungroup()
      
      colnames(anova)<-c("name", "anova_p.val")
      
      # add anova p.value to row data
      rowData(anova_dep) <- merge(rowData(anova_dep), anova, by = 'name', sort = FALSE)
      anova_dep_rej <- DEP::add_rejections(anova_dep,alpha = input$p, lfc= input$lfc)
      
      # calculate adjusted anova p.value to data
      anova_adj <-anova
      anova_adj$anova_p.adj <- p.adjust(anova_adj$anova_p.val,method = input$fdr_correction)
      anova_adj <- anova_adj %>% select(-anova_p.val)
      
      # add adjusted anova p.value to row data
      rowData(anova_dep_rej) <- merge(rowData(anova_dep_rej), anova_adj, by = 'name', sort = FALSE)
      return(anova_dep_rej)
    }
    
  })
  
  comparisons_nr<-reactive ({
    temp<-capture.output(DEP::test_diff(normalised_data_nr(),type='all'),type = "message")
    gsub(".*: ","",temp)
    ## Split conditions into character vector
    unlist(strsplit(temp,","))
    ## Remove leading and trailing spaces
    trimws(temp)
  })
  
  ## Results plot inputs in Normalized page
  
  ## PCA Plot
  pca_input_nr<-eventReactive(input$analyze ,{ 
    if(input$analyze==0 ){
      return()
    }
    if (num_total_nr()<=500){
      if(length(levels(as.factor(colData(dep_nr())$replicate))) <= 6){
        pca_plot<-DEP::plot_pca(dep_nr(), n=num_total_nr(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
      else{
        pca_plot<-DEP::plot_pca(dep_nr(), n=num_total_nr(), point_size = 4, indicate = "condition") 
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
    }
    else{
      if(length(levels(as.factor(colData(dep_nr())$replicate))) <= 6){
        pca_plot<-DEP::plot_pca(dep_nr(),point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
      else{
        #pca_label<-SummarizedExperiment::colData(dep_nr())$replicate
        pca_plot<-DEP::plot_pca(dep_nr(), point_size = 4, indicate = "condition")
        #pca_plot<-pca_plot + geom_point()
        pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                                      size = 4,
                                                      box.padding = unit(0.1, 'lines'),
                                                      point.padding = unit(0.1, 'lines'),
                                                      segment.size = 0.5)
        pca_plot<-pca_plot + labs(title = "PCA Plot")
        return(pca_plot)
      }
    }
    
  })
  
  ### Heatmap Differentially expressed proteins
  heatmap_input_nr<-eventReactive(input$analyze ,{ 
    if(input$analyze==0 ){
      return()
    }
    get_cluster_heatmap(dep_nr(),
                        type="centered",kmeans = TRUE,
                        k=input$k_number, col_limit = 6,
                        indicate = "condition"
    )
  })
  
  ### Volcano Plot
  volcano_input_nr <- reactive({
    if(!is.null(input$volcano_cntrst_nr)) {
      if(length(unique(exp_design_input()$condition)) <= 2) {
        plot_volcano_new(dep_nr(),
                         input$volcano_cntrst_nr,
                         FALSE,
                         input$check_names_nr,
                         input$p_adj_nr)
      } else {
        plot_volcano_new(dep_nr(),
                         input$volcano_cntrst_nr,
                         input$check_anova_nr,
                         input$check_names_nr,
                         input$p_adj_nr)
      }
      
    }
  })
  
  volcano_df_nr<- reactive({
    if(!is.null(input$volcano_cntrst_nr)) {
      get_volcano_df(dep_nr(),
                     input$volcano_cntrst_nr) 
      
    }
  })
  
  
  volcano_input_selected_nr<-reactive({
    if(!is.null(input$volcano_cntrst_nr)){
      
      if (!is.null(input$contents_nr_rows_selected)){
        proteins_selected<-data_result_nr()[c(input$contents_nr_rows_selected),]## get all rows selected
      }
      else if(!is.null(input$protein_brush_nr)){
        proteins_selected<-data_result_nr()[data_result_nr()[["Phosphosite"]] %in% protein_name_brush_nr(), ] 
      }
      print(proteins_selected)
      ## convert contrast to x and padj to y
      diff_proteins <- grep(paste(input$volcano_cntrst_nr, "_log2", sep = ""),
                            colnames(proteins_selected))
      if(input$p_adj=="FALSE"){
        padj_proteins <- grep(paste(input$volcano_cntrst_nr, "_p.val", sep = ""),
                              colnames(proteins_selected))
      }
      else{
        padj_proteins <- grep(paste(input$volcano_cntrst, "_p.adj", sep = ""),
                              colnames(proteins_selected))
      }
      
      df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                               y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                               name = proteins_selected$Phosphosite)
      # print(df_protein)
      if(length(unique(exp_design_input()$condition)) <= 2){
        p<-plot_volcano_new(dep_nr(),
                            input$volcano_cntrst_nr,
                            FALSE,
                            input$check_names_nr,
                            input$p_adj_nr)
        
      } else {
        p<-plot_volcano_new(dep_nr(),
                            input$volcano_cntrst_nr,
                            input$check_anova_nr,
                            input$check_names_nr,
                            input$p_adj_nr)
      }
      p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
        ggrepel::geom_text_repel(data = df_protein,
                                 aes(x, y, label = name),
                                 size = 4,
                                 box.padding = unit(0.1, 'lines'),
                                 point.padding = unit(0.1, 'lines'),
                                 segment.size = 0.5)## use the dataframe to plot points
      
    }
  })
  
  protein_input_nr<-reactive({ 
    
    protein_selected  <- data_result_nr()[input$contents_nr_rows_selected,1]
    protein_selected <-as.character(protein_selected)
    if(length(levels(as.factor(colData(dep_nr())$replicate))) <= 8){
      plot_protein(dep_nr(), protein_selected, as.character(input$type))
    }
    else{
      protein_plot<-plot_protein(dep_nr(), protein_selected, as.character(input$type))
      protein_plot + scale_color_brewer(palette = "Paired")
    }
    
  })
  
  ## QC Inputs
  norm_input_nr <- reactive({
    plot_normalization(normalized_phospho_data(),
                       normalised_data_nr())
  })
  
  missval_input_nr <- reactive({
    plot_missval(normalized_phospho_data())
  })
  
  detect_input_nr <- reactive({
    plot_detect(normalized_phospho_data())
  })
  
  imputation_input_nr <- reactive({
    plot_imputation(normalized_phospho_data(),
                    diff_all_nr())
  })
  
  p_hist_input_nr <- reactive({
    plot_p_hist(dep_nr())
  })
  
  numbers_input_nr <- reactive({
    plot_numbers(normalised_data_nr()) +
      labs(title= "Phosphosites per sample", y = "Number of phosphosites")
  })
  
  coverage_input_nr <- reactive({
    plot_coverage(normalised_data_nr())+
      labs(title= "Phosphosites per sample", y = "Number of phosphosites")
  })
  
  correlation_input_nr<-reactive({
    plot_cor(dep_nr(),significant = FALSE)
  })
  
  cvs_input_nr<-reactive({
    plot_cvs(dep_nr())
  })
  
  num_total_nr<-reactive({
    dep_nr() %>%
      nrow()
  }) 
  
  ## Enrichment inputs
  
  go_input_nr<-eventReactive(input$go_analysis_nr,{
    withProgress(message = 'Gene ontology enrichment is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    
    if(!is.null(input$contrast_nr)){
      enrichment_output_test(dep_nr(), input$go_database_nr)
      go_results<- test_gsea_mod_phospho(dep_nr(), databases = input$go_database_nr, contrasts = TRUE)
      null_enrichment_test(go_results)
      plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast_nr,
                                databases = input$go_database_nr, nrow = 2, term_size = 8) + aes(stringr::str_wrap(Term, 60)) +
        xlab(NULL)
      go_list<-list("go_result"=go_results, "plot_go"=plot_go)
      return(go_list)
    }
  })
  
  pathway_input_nr<-eventReactive(input$pathway_analysis_nr,{
    progress_indicator("Pathway Analysis is running....")
    enrichment_output_test(dep_nr(), input$pathway_database_nr)
    pathway_results<- test_gsea_mod_phospho(dep_nr(), databases=input$pathway_database_nr, contrasts = TRUE)
    null_enrichment_test(pathway_results)
    plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_1_nr,
                                  databases=input$pathway_database_nr, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
      xlab(NULL)
    pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
    return(pathway_list)
  })
  
  #### Interactive UI (Normalized page)
  output$significantBox_nr <- renderInfoBox({
    num_total <- dep_nr() %>%
      nrow()
    num_signif <- dep_nr() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
    frac <- num_signif / num_total
    
    info_box <- 		infoBox("Significant phosphosites",
                          paste0(num_signif,
                                 " out of ",
                                 num_total),
                          paste0(signif(frac * 100, digits = 3),
                                 "% of phosphosites differentially expressed across all conditions"),
                          icon = icon("stats", lib = "glyphicon"),
                          color = "olive",
                          # fill = TRUE,
                          width = 4)
    
    return(info_box)
  })
  
  
  ##### Get results dataframe from Summarizedexperiment object
  data_result_nr<-reactive({
    if(length(unique(exp_design_input()$condition)) <= 2){
      get_results_phospho(dep_nr(),FALSE) %>% dplyr::select (-Residue.Both,-peptide.sequence,-Protein)
    } else {
      get_results_phospho(dep_nr(),TRUE) %>% dplyr::select (-Residue.Both,-peptide.sequence,-Protein)
    }
  })
  
  
  #### Data table
  output$contents_nr <- DT::renderDataTable({
    withProgress(message = 'Result table calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    df<- data_result_nr()
    return(df)
  },
  options = list(scrollX = TRUE,
                 autoWidth=TRUE,
                 columnDefs= list(list(width = '400px', targets = c(-1))))
  )
  
  ## Deselect all rows button
  proxy <- dataTableProxy("contents_nr")
  
  observeEvent(input$clear_nr,{
    proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$original_nr,{
    output$contents_nr <- DT::renderDataTable({
      df<- data_result_nr()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })
  
  protein_name_brush_nr<- reactive({
    #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    protein_tmp<-brushedPoints(volcano_df_nr(), input$protein_brush_nr, 
                               xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  protein_name_click_nr<- reactive({
    protein_tmp<-nearPoints(volcano_df_nr(), input$protein_click_nr, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
    #xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  
  
  ## Select rows dynamically
  observeEvent(input$protein_brush_nr,{
    output$contents_nr <- DT::renderDataTable({
      df<- data_result_nr()[data_result_nr()[["Phosphosite"]] %in% protein_name_brush_nr(), ]
      return(df)
    },
    options = list(scrollX= TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    proteins_selected<-data_result_nr()[data_result_nr()[["Phosphosite"]] %in% protein_name_brush_nr(), ] #
    # get all rows selected
    ## convert contrast to x and padj to y
    diff_proteins <- grep(paste(input$volcano_cntrst_nr, "_log2", sep = ""),
                          colnames(proteins_selected))
    if(input$p_adj=="FALSE"){
      padj_proteins <- grep(paste(input$volcano_cntrst_nr, "_p.val", sep = ""),
                            colnames(proteins_selected))
    }
    else{
      padj_proteins <- grep(paste(input$volcano_cntrst_nr, "_p.adj", sep = ""),
                            colnames(proteins_selected))
    }
    df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                             y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                             name = proteins_selected$Phosphosite)
    #print(df_protein)
    if(length(unique(exp_design_input()$condition)) <= 2) {
      p<-plot_volcano_new(dep_nr(),
                          input$volcano_cntrst_nr,
                          FALSE,
                          input$check_names_nr,
                          input$p_adj_nr)
    } else {
      p<-plot_volcano_new(dep_nr(),
                          input$volcano_cntrst_nr,
                          input$check_anova_nr,
                          input$check_names_nr,
                          input$p_adj_nr)
    }
    
    p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
      ggrepel::geom_text_repel(data = df_protein,
                               aes(x, y, label = name),
                               size = 4,
                               box.padding = unit(0.1, 'lines'),
                               point.padding = unit(0.1, 'lines'),
                               segment.size = 0.5)
    
    output$volcano_nr <- renderPlot({
      withProgress(message = 'Volcano Plot calculations are in progress',
                   detail = 'Please wait for a while', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      p
    })
    return(p)
  })
  
  observeEvent(input$resetPlot_nr,{
    session$resetBrush("protein_brush_nr")
    brush <<- NULL
    
    output$contents_nr <- DT::renderDataTable({
      df<- data_result_nr()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    output$volcano_nr <- renderPlot({
      volcano_input_nr()
    })
  })
  
  observeEvent(input$protein_click_nr,{
    output$contents_nr <- DT::renderDataTable({
      df<- data_result_nr()[data_result_nr()[["Phosphosite"]] %in% protein_name_click_nr(), ]
      return(df)
    },
    options = list(scrollX= TRUE)
    )
  })
  
  ## Render Result Plots
  output$pca_plot_nr<-renderPlot({
    pca_input_nr()
  })
  output$heatmap_nr<-renderPlot({
    withProgress(message = 'Heatmap rendering is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    heatmap_input_nr()
  })
  
  output$volcano_nr <- renderPlot({
    withProgress(message = 'Volcano Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_nr_rows_selected)){
      volcano_input_nr()
    }
    else if(!is.null(input$volcano_cntrst_nr)){
      volcano_input_selected_nr()
    } # else close
  })
  
  output$protein_plot_nr<-renderPlot({
    if(!is.null(input$contents_nr_rows_selected)){
      protein_input_nr()
    }
  })
  
  
  ### QC Outputs
  output$sample_corr_nr <-renderPlot({
    correlation_input_nr()
  })
  
  output$sample_cvs_nr <- renderPlot({
    cvs_input_nr()
  })
  
  output$norm_nr <- renderPlot({
    norm_input_nr()
  })
  
  output$missval_nr <- renderPlot({
    missval_input_nr()
  })
  
  output$detect_nr <- renderPlot({
    detect_input_nr()
  })
  
  output$imputation_nr <- renderPlot({
    imputation_input_nr()
  })
  
  # output$p_hist <- renderPlot({
  #   p_hist_input_nr()
  # })
  
  output$numbers_nr <- renderPlot({
    numbers_input_nr()
  })
  
  output$coverage_nr <- renderPlot({
    coverage_input_nr()
  })
  
  ## Enrichment Outputs
  output$go_enrichment_nr<-renderPlot({
    go_input_nr()$plot_go
  })
  
  output$pathway_enrichment_nr<-renderPlot({
    pathway_input_nr()$plot_pa
  })
  
  ##### Download Functions
  datasetInput_nr <- reactive({
    switch(input$dataset_nr,
           "Results" = get_results_proteins(dep_nr()),
           "Full dataset" = get_df_wide(dep_nr()))
  })
  
  output$downloadData_nr <- downloadHandler(
    filename = function() { paste(input$dataset_nr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(datasetInput_nr(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ### === Cluster Download ==== ####
  
  individual_cluster_nr <- reactive({
    cluster_number <- input$cluster_number_nr
    cluster_all <- heatmap_input_nr()
    data_result_nr()[cluster_all[[cluster_number]],]
  })
  
  
  
  output$downloadCluster_nr <- downloadHandler(
    filename = function() { paste("Cluster_info_",input$cluster_number_nr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(individual_cluster_nr(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  output$downloadVolcano_nr <- downloadHandler(
    filename = function() {
      paste0("Volcano_", input$volcano_cntrst_nr, ".pdf")
    },
    content = function(file) {
      pdf(file)
      print(volcano_input_selected_nr())
      dev.off()
    }
  )
  
  
  ## Protein plot download
  output$downloadProtein_nr <- downloadHandler(
    filename = function() {
      paste0(input$type_nr,".pdf")
    },
    content = function(file) {
      pdf(file)
      print(protein_input_nr())
      dev.off()
    }
  )
  
  ###### ==== DOWNLOAD GO TABLE ==== ####
  output$downloadGO_nr <- downloadHandler(
    filename = function() { paste("GO_enrichment_",input$go_database_nr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(go_input_nr()$go_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
  output$downloadPA_nr <- downloadHandler(
    filename = function() { paste("Pathway_enrichment_",input$pathway_database_nr, ".csv", sep = "") }, 
    ## use = instead of <-
    content = function(file) {
      write.table(pathway_input_nr()$pa_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  
  
  #####===== Download Report (normalized phosphosite)=====#####
  output$downloadReport_nr <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "Phospho-Analyst(normalization)report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "Normalised_phosphosite_report.Rmd")
      file.copy("Normalised_phosphosite_report.Rmd", tempReport, overwrite = TRUE)
      
      sig_proteins<-dep_nr() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()
      
      tested_contrasts<- gsub("_significant", "", 
                              colnames(SummarizedExperiment::rowData(dep_nr()))[grep("_significant", 
                                                                                  colnames(SummarizedExperiment::rowData(dep_nr())))])
      pg_width<- ncol(normalised_data_nr()) / 2.5
      # Set up parameters to pass to Rmd document
      params <- list(data = normalized_phospho_data,
                     alpha = input$p,
                     lfc = input$lfc,
                     num_signif= sig_proteins,
                     pg_width = pg_width,
                     tested_contrasts= tested_contrasts,
                     numbers_input= numbers_input_nr,
                     detect_input = detect_input_nr,
                     imputation_input = imputation_input_nr,
                     missval_input = missval_input_nr,
                     p_hist_input = p_hist_input_nr,
                     pca_input = pca_input_nr,
                     coverage_input= coverage_input_nr,
                     correlation_input =correlation_input_nr,
                     heatmap_input = heatmap_input_nr,
                     cvs_input= cvs_input_nr,
                     dep = dep_nr
      )
      
      # Knit the document, passing in the `params` list
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  #### Demo logic (Phosphosite)========== #############
  
  ####======= Render Functions
  output$volcano_cntrst_dm <- renderUI({
    if (!is.null(comparisons_dm())) {
      df <- SummarizedExperiment::rowData(dep_dm())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("volcano_cntrst_dm",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  ##comparisons
  output$contrast_dm <- renderUI({
    if (!is.null(comparisons_dm())) {
      df <- SummarizedExperiment::rowData(dep_dm())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_dm",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$contrast_1_dm <- renderUI({
    if (!is.null(comparisons_dm())) {
      df <- SummarizedExperiment::rowData(dep_dm())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_1_dm",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$downloadTable_dm <- renderUI({
    if(!is.null(dep_dm())){
      selectizeInput("dataset_dm",
                     "Choose a dataset to save" ,
                     c("Results","Original_matrix",
                       "Imputed_matrix",
                       "Full_dataset"))
    }
  })
  
  output$downloadButton_dm <- renderUI({
    if(!is.null(dep_dm())){
      downloadButton('downloadData_dm', 'Save')
    }
  })
  
  output$downloadZip_dm <- renderUI({
    if(!is.null(dep_dm())){
      downloadButton('downloadZip1_dm', 'Download result plots')
    }
  })
  output$downloadreport_dm <- renderUI({
    if(!is.null(dep_dm())){
      downloadButton('downloadReport_dm', 'Download Report')
    }
  })
  
  output$downloadPlots <- renderUI({
    if(!is.null(dep_dm())){
      downloadButton('downloadPlots1_dm', 'Download Plots')
    }
  })
  
  # load demo data
  env_dm<-reactive({
    LoadToEnvironment("data/phosphosite_demo_data.RData", env = globalenv())
  })
  
  processed_data_dm<-reactive({
    env_dm()[["data_missval"]]
  })
  
  
  unimputed_table_dm<-reactive({
    temp<-assay(processed_data_dm())
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"original_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) 
    #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  imputed_data_dm<-reactive({
    DEP::impute(processed_data_dm(),input$imputation)
  })
  
  normalised_data_dm<-reactive({
    normalize_vsn(imputed_data_dm())
  })
  
  imputed_table_dm<-reactive({
    temp<-assay(imputed_data_dm())
    #tibble::rownames_to_column(temp,var = "ProteinID")
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"imputed_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  diff_all_dm<-reactive({
    test_diff(normalised_data_dm(),type = 'all')
  })
  
  dep_dm<-reactive({
    env_dm()[["data_dep"]]
  })
  
  comparisons_dm<-reactive({
    comparisons<-gsub("_p.adj", "", 
                      colnames(SummarizedExperiment::rowData(dep_dm()))
                      [grep("p.adj", colnames(SummarizedExperiment::rowData(dep_dm())))])
  })
  
  ## Results plot inputs in Normalized page
  
  ## PCA Plot
  pca_label_dm<-reactive({
    pca_lable<-levels(as.factor(colData(dep_dm())$replicate))
    print(pca_label)
  })
  
  pca_input_dm<-reactive({
    if (num_total_dm()<=500){
      if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 6){
        pca_plot<- DEP::plot_pca(dep_dm(), n=num_total_dm(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
      else{
        pca_plot<-DEP::plot_pca(dep_dm(), n=num_total_dm(), point_size = 4, indicate = "condition")
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
    }
    else{
      if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 6){
        pca_plot<-DEP::plot_pca(dep_dm(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }else{
        pca_label<-SummarizedExperiment::colData(dep_dm())$replicate
        pca_plot<-DEP::plot_pca(dep_dm(), point_size = 4, indicate = "condition")
        #pca_plot<-pca_plot + geom_point()
        pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                                      size = 4,
                                                      box.padding = unit(0.1, 'lines'),
                                                      point.padding = unit(0.1, 'lines'),
                                                      segment.size = 0.5)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
    }
    
  })
  
  ### Heatmap Differentially expressed proteins
  heatmap_input_dm<-reactive({ 
    get_cluster_heatmap(dep_dm(),
                        type="centered",kmeans = TRUE,
                        k=6, col_limit = 6,
                        indicate = "condition"
    )
  })
  
  ### Volcano Plot
  volcano_input_dm <- reactive({
    if(!is.null(input$volcano_cntrst_dm)) {
      plot_volcano_new(dep_dm(),
                       input$volcano_cntrst_dm,
                       input$check_anova_dm,
                       input$check_names_dm,
                       input$p_adj_dm)
      
      
    }
  })
  
  volcano_df_dm<- reactive({
    if(!is.null(input$volcano_cntrst_dm)) {
      get_volcano_df(dep_dm(),
                     input$volcano_cntrst_dm) 
      
    }
  })
  
  
  volcano_input_selected_dm<-reactive({
    if(!is.null(input$volcano_cntrst_dm)){
      
      if (!is.null(input$contents_dm_rows_selected)){
        proteins_selected<-data_result_dm()[c(input$contents_dm_rows_selected),]## get all rows selected
      }
      else if(!is.null(input$protein_brush_dm)){
        proteins_selected<-data_result_dm()[data_result_dm()[["Phosphosite"]] %in% protein_name_brush_dm(), ] 
      }
      print(proteins_selected)
      ## convert contrast to x and padj to y
      diff_proteins <- grep(paste(input$volcano_cntrst_dm, "_log2", sep = ""),
                            colnames(proteins_selected))
      if(input$p_adj=="FALSE"){
        padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.val", sep = ""),
                              colnames(proteins_selected))
      }
      else{
        padj_proteins <- grep(paste(input$volcano_cntrst, "_p.adj", sep = ""),
                              colnames(proteins_selected))
      }
      
      df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                               y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                               name = proteins_selected$Phosphosite)
      # print(df_protein)
      p<-plot_volcano_new(dep_dm(),
                          input$volcano_cntrst_dm,
                          input$check_anova_dm,
                          input$check_names_dm,
                          input$p_adj_dm)
      
      p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
        ggrepel::geom_text_repel(data = df_protein,
                                 aes(x, y, label = name),
                                 size = 4,
                                 box.padding = unit(0.1, 'lines'),
                                 point.padding = unit(0.1, 'lines'),
                                 segment.size = 0.5)## use the dataframe to plot points
      
    }
  })
  
  protein_input_dm<-reactive({ 
    
    protein_selected  <- data_result_dm()[input$contents_dm_rows_selected,1]
    protein_selected <-as.character(protein_selected)
    if(length(levels(as.factor(colData(dep_dm())$replicate))) <= 8){
      plot_protein(dep_dm(), protein_selected, as.character(input$type))
    }
    else{
      protein_plot<-plot_protein(dep_dm(), protein_selected, as.character(input$type))
      protein_plot + scale_color_brewer(palette = "Paired")
    }
    
  })
  
  ## QC Inputs
  norm_input_dm <- reactive({
    plot_normalization(processed_data_dm(),
                       normalised_data_dm())
  })
  
  missval_input_dm <- reactive({
    plot_missval(processed_data_dm())
  })
  
  detect_input_dm <- reactive({
    plot_detect(processed_data_dm())
  })
  
  imputation_input_dm <- reactive({
    plot_imputation(processed_data_dm(),
                    diff_all_dm())
  })
  
  p_hist_input_dm <- reactive({
    plot_p_hist(dep_dm())
  })
  
  numbers_input_dm <- reactive({
    plot_numbers(processed_data_dm()) +
      labs(title= "Phosphosites per sample", y = "Number of phosphosites")
  })
  
  coverage_input_dm <- reactive({
    plot_coverage(processed_data_dm())+
      labs(title= "Phosphosites per sample", y = "Number of phosphosites")
  })
  
  correlation_input_dm<-reactive({
    plot_cor(dep_dm(),significant = FALSE)
  })
  
  cvs_input_dm<-reactive({
    plot_cvs(dep_dm())
  })
  
  num_total_dm<-reactive({
    dep_dm() %>%
      nrow()
  }) 
  
  ## Enrichment inputs
  
  go_input_dm<-eventReactive(input$go_analysis_dm,{
    withProgress(message = 'Gene ontology enrichment is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    
    if(!is.null(input$contrast_dm)){
      enrichment_output_test(dep_dm(), input$go_database_dm)
      go_results<- test_gsea_mod_phospho(dep_dm(), databases = input$go_database_dm, contrasts = TRUE)
      null_enrichment_test(go_results)
      if (input$go_database_dm == "KEGG" | input$go_database_dm == "Reactome"){
        plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_1_dm,
                                      databases=input$pathway_database_dm, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
          xlab(NULL)
      }
      else{
        plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast_dm,
                                  databases = input$go_database_dm, nrow = 2, term_size = 8) + aes(stringr::str_wrap(Term, 60)) +
          xlab(NULL)
      }
      go_list<-list("go_result"=go_results, "plot_go"=plot_go)
      return(go_list)
    }
  })
  
  KSEA_input_dm<-eventReactive(input$KSEA_analysis_dm,{
    progress_indicator("Kinase-Substrate Analysis is running....")
    
    result_df <- get_results_phospho(dep_dm(),FALSE)
    print(input$contrast_1_dm)  #test
    col_selected <- c('Protein','Gene.names','peptide.sequence', 'Residue.Both',
                      paste(input$contrast_1_dm, "_p.val", sep = ""),
                      paste(input$contrast_1_dm, "_log2 fold change", sep = ""))
    print(col_selected) #test
    
    # select required columns and rename them
    column_names <- c('Protein','Gene','Peptide','Residue.Both','p','FC')
    PX <- result_df %>% dplyr::select (col_selected)
    names(PX) <- column_names
    KSData <- KSEAapp::KSData 
    
    # Generate a summary bar plot using the KSEA.Barplot() function
    plot_KSEA_dm <- KSEAapp::KSEA.Barplot(KSData, PX, NetworKIN=TRUE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.01, export=FALSE)
    return(plot_KSEA_dm)
  })

  # pathway_input_dm<-eventReactive(input$pathway_analysis_dm,{
  #   progress_indicator("Pathway Analysis is running....")
  #   enrichment_output_test(dep_dm(), input$pathway_database_dm)
  #   pathway_results<- test_gsea_mod_phospho(dep_dm(), databases=input$pathway_database_dm, contrasts = TRUE)
  #   null_enrichment_test(pathway_results)
  #   plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_1_dm,
  #                                 databases=input$pathway_database_dm, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
  #     xlab(NULL)
  #   pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
  #   return(pathway_list)
  # })
  
  #### Interactive UI (Normalized page)
  output$significantBox_dm <- renderInfoBox({
    num_total <- dep_dm() %>%
      nrow()
    num_signif <- dep_dm() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
    frac <- num_signif / num_total
    
    info_box <- 		infoBox("Significant phosphosites",
                          paste0(num_signif,
                                 " out of ",
                                 num_total),
                          paste0(signif(frac * 100, digits = 3),
                                 "% of phosphosites differentially expressed across all conditions"),
                          icon = icon("stats", lib = "glyphicon"),
                          color = "olive",
                          # fill = TRUE,
                          width = 4)
    
    return(info_box)
  })
  
  
  ##### Get results dataframe from Summarizedexperiment object
  data_result_dm<-reactive({
      get_results_phospho(dep_dm(),TRUE) %>% dplyr::select (-Residue.Both,-peptide.sequence,-Protein)
  })
  
  
  #### Data table
  output$contents_dm <- DT::renderDataTable({
    withProgress(message = 'Result table calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    df<- data_result_dm()
    return(df)
  },
  options = list(scrollX = TRUE,
                 autoWidth=TRUE,
                 columnDefs= list(list(width = '400px', targets = c(-1))))
  )
  
  ## Deselect all rows button
  proxy <- dataTableProxy("contents_dm")
  
  observeEvent(input$clear_dm,{
    proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$original_dm,{
    output$contents_dm <- DT::renderDataTable({
      df<- data_result_dm()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })
  
  protein_name_brush_dm<- reactive({
    #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    protein_tmp<-brushedPoints(volcano_df_dm(), input$protein_brush_dm, 
                               xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  protein_name_click_dm<- reactive({
    protein_tmp<-nearPoints(volcano_df_dm(), input$protein_click_dm, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
    #xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  
  
  ## Select rows dynamically
  observeEvent(input$protein_brush_dm,{
    output$contents_dm <- DT::renderDataTable({
      df<- data_result_dm()[data_result_dm()[["Phosphosite"]] %in% protein_name_brush_dm(), ]
      return(df)
    },
    options = list(scrollX= TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    proteins_selected<-data_result_dm()[data_result_dm()[["Phosphosite"]] %in% protein_name_brush_dm(), ] #
    # get all rows selected
    ## convert contrast to x and padj to y
    diff_proteins <- grep(paste(input$volcano_cntrst_dm, "_log2", sep = ""),
                          colnames(proteins_selected))
    if(input$p_adj=="FALSE"){
      padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.val", sep = ""),
                            colnames(proteins_selected))
    }
    else{
      padj_proteins <- grep(paste(input$volcano_cntrst_dm, "_p.adj", sep = ""),
                            colnames(proteins_selected))
    }
    df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                             y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                             name = proteins_selected$Phosphosite)
    
      p<-plot_volcano_new(dep_dm(),
                          input$volcano_cntrst_dm,
                          input$check_anova_dm,
                          input$check_names_dm,
                          input$p_adj_dm)
    
    
    p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
      ggrepel::geom_text_repel(data = df_protein,
                               aes(x, y, label = name),
                               size = 4,
                               box.padding = unit(0.1, 'lines'),
                               point.padding = unit(0.1, 'lines'),
                               segment.size = 0.5)
    
    output$volcano_dm <- renderPlot({
      withProgress(message = 'Volcano Plot calculations are in progress',
                   detail = 'Please wait for a while', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      p
    })
    return(p)
  })
  
  observeEvent(input$resetPlot_dm,{
    session$resetBrush("protein_brush_dm")
    brush <<- NULL
    
    output$contents_dm <- DT::renderDataTable({
      df<- data_result_dm()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    output$volcano_dm <- renderPlot({
      volcano_input_dm()
    })
  })
  
  observeEvent(input$protein_click_dm,{
    output$contents_dm <- DT::renderDataTable({
      df<- data_result_dm()[data_result_dm()[["Phosphosite"]] %in% protein_name_click_dm(), ]
      return(df)
    },
    options = list(scrollX= TRUE)
    )
  })
  
  ## Render Result Plots
  output$pca_plot_dm<-renderPlot({
    pca_input_dm()
  })
  output$heatmap_dm<-renderPlot({
    withProgress(message = 'Heatmap rendering is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    heatmap_input_dm()
  })
  
  output$volcano_dm <- renderPlot({
    withProgress(message = 'Volcano Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_dm_rows_selected)){
      volcano_input_dm()
    }
    else if(!is.null(input$volcano_cntrst_dm)){
      volcano_input_selected_dm()
    } # else close
  })
  
  output$protein_plot_dm<-renderPlot({
    if(!is.null(input$contents_dm_rows_selected)){
      protein_input_dm()
    }
  })
  
  
  ### QC Outputs
  output$sample_corr_dm <-renderPlot({
    correlation_input_dm()
  })
  
  output$sample_cvs_dm <- renderPlot({
    cvs_input_dm()
  })
  
  output$norm_dm <- renderPlot({
    norm_input_dm()
  })
  
  output$missval_dm <- renderPlot({
    missval_input_dm()
  })
  
  output$detect_dm <- renderPlot({
    detect_input_dm()
  })
  
  output$imputation_dm <- renderPlot({
    imputation_input_dm()
  })
  
  # output$p_hist_dm <- renderPlot({
  #   p_hist_input_dm()
  # })
  
  output$numbers_dm <- renderPlot({
    numbers_input_dm()
  })
  
  output$coverage_dm <- renderPlot({
    coverage_input_dm()
  })
  
  ## Enrichment Outputs
  output$go_enrichment_dm<-renderPlot({
    go_input_dm()$plot_go
  })
  
  output$KSEA_enrichment_dm<-renderPlot({
    KSEA_input_dm()
  })
  
  ##### Download Functions
  datasetInput_dm <- reactive({
    switch(input$dataset_dm,
           "Results" = get_results_proteins(dep_dm()),
           "Full dataset" = get_df_wide(dep_dm()))
  })
  
  output$downloadData_dm <- downloadHandler(
    filename = function() { paste(input$dataset_dm, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(datasetInput_dm(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ### === Cluster Download ==== ####
  
  individual_cluster_dm <- reactive({
    cluster_number <- input$cluster_number_dm
    cluster_all <- heatmap_input_dm()
    data_result_dm()[cluster_all[[cluster_number]],]
  })
  
  
  
  output$downloadCluster_dm <- downloadHandler(
    filename = function() { paste("Cluster_info_",input$cluster_number_dm, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(individual_cluster_dm(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  output$downloadVolcano_dm <- downloadHandler(
    filename = function() {
      paste0("Volcano_", input$volcano_cntrst_dm, ".pdf")
    },
    content = function(file) {
      pdf(file)
      print(volcano_input_selected_dm())
      dev.off()
    }
  )
  
  
  ## Protein plot download
  output$downloadProtein_dm <- downloadHandler(
    filename = function() {
      paste0(input$type_dm,".pdf")
    },
    content = function(file) {
      pdf(file)
      print(protein_input_dm())
      dev.off()
    }
  )
  
  ###### ==== DOWNLOAD GO TABLE ==== ####
  output$downloadGO_dm <- downloadHandler(
    filename = function() { paste("GO_enrichment_",input$go_database_dm, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(go_input_dm()$go_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
  output$downloadPA_dm <- downloadHandler(
    filename = function() { paste("Pathway_enrichment_",input$pathway_database_dm, ".csv", sep = "") }, 
    ## use = instead of <-
    content = function(file) {
      write.table(pathway_input_dm()$pa_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  
  
  #####===== Download Report (demo phosphosite)=====#####
  output$downloadReport_dm <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "Phospho-Analyst(normalization)report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "Normalised_phosphosite_report.Rmd")
      file.copy("Normalised_phosphosite_report.Rmd", tempReport, overwrite = TRUE)
      
      sig_proteins<-dep_dm() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()
      
      tested_contrasts<- gsub("_significant", "", 
                              colnames(SummarizedExperiment::rowData(dep_dm()))[grep("_significant", 
                                                                                     colnames(SummarizedExperiment::rowData(dep_dm())))])
      pg_width<- ncol(normalised_data_dm()) / 2.5
      # Set up parameters to pass to Rmd document
      params <- list(data = normalized_phospho_data,
                     alpha = input$p,
                     lfc = input$lfc,
                     num_signif= sig_proteins,
                     pg_width = pg_width,
                     tested_contrasts= tested_contrasts,
                     numbers_input= numbers_input_dm,
                     detect_input = detect_input_dm,
                     imputation_input = imputation_input_dm,
                     missval_input = missval_input_dm,
                     p_hist_input = p_hist_input_dm,
                     pca_input = pca_input_dm,
                     coverage_input= coverage_input_dm,
                     correlation_input =correlation_input_dm,
                     heatmap_input = heatmap_input_dm,
                     cvs_input= cvs_input_dm,
                     dep = dep_dm
      )
      
      # Knit the document, passing in the `params` list
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  #### Demo logic (Protein Group)========== #############
  
  ####======= Render Functions
  output$volcano_cntrst_dm_pr <- renderUI({
    if (!is.null(comparisons_dm_pr())) {
      df <- SummarizedExperiment::rowData(dep_dm_pr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("volcano_cntrst_dm_pr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  ##comparisons
  output$contrast_dm_pr <- renderUI({
    if (!is.null(comparisons_dm_pr())) {
      df <- SummarizedExperiment::rowData(dep_dm_pr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_dm_pr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$contrast_1_dm_pr <- renderUI({
    if (!is.null(comparisons_dm_pr())) {
      df <- SummarizedExperiment::rowData(dep_dm_pr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_1_dm_pr",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })
  
  output$downloadTable_dm_pr <- renderUI({
    if(!is.null(dep_dm_pr())){
      selectizeInput("dataset_dm_pr",
                     "Choose a dataset to save" ,
                     c("Results","Original_matrix",
                       "Imputed_matrix",
                       "Full_dataset"))
    }
  })
  
  output$downloadButton_dm_pr <- renderUI({
    if(!is.null(dep_dm_pr())){
      downloadButton('downloadData_dm_pr', 'Save')
    }
  })
  
  output$downloadZip_dm_pr <- renderUI({
    if(!is.null(dep_dm_pr())){
      downloadButton('downloadZip1_dm_pr', 'Download result plots')
    }
  })
  output$downloadreport_dm_pr <- renderUI({
    if(!is.null(dep_dm_pr())){
      downloadButton('downloadReport_dm_pr', 'Download Report')
    }
  })
  
  output$downloadPlots <- renderUI({
    if(!is.null(dep_dm_pr())){
      downloadButton('downloadPlots1_dm_pr', 'Download Plots')
    }
  })
  
  # load demo data
  env_dm_pr<-reactive({
    LoadToEnvironment("data/proteinGroup_demo_data.RData", env = globalenv())
  })
  
  processed_data_dm_pr<-reactive({
    env_dm_pr()[["data_missval"]]
  })
  
  
  unimputed_table_dm_pr<-reactive({
    temp<-assay(processed_data_dm_pr())
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"original_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) 
    #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  imputed_data_dm_pr<-reactive({
    DEP::impute(processed_data_dm_pr(),input$imputation)
  })
  
  normalised_data_dm_pr<-reactive({
    normalize_vsn(imputed_data_dm_pr())
  })
  
  imputed_table_dm_pr<-reactive({
    temp<-assay(imputed_data_dm_pr())
    #tibble::rownames_to_column(temp,var = "ProteinID")
    temp1<-2^(temp)
    colnames(temp1)<-paste(colnames(temp1),"imputed_intensity",sep="_")
    temp1<-cbind(ProteinID=rownames(temp1),temp1) #temp1$ProteinID<-rownames(temp1)
    return(as.data.frame(temp1))
  })
  
  diff_all_dm_pr<-reactive({
    test_diff(normalised_data_dm_pr(),type = 'all')
  })
  
  dep_dm_pr<-reactive({
    env_dm_pr()[["data_dep"]]
  })
  
  comparisons_dm_pr<-reactive({
    comparisons<-gsub("_p.adj", "", 
                      colnames(SummarizedExperiment::rowData(dep_dm_pr()))
                      [grep("p.adj", colnames(SummarizedExperiment::rowData(dep_dm_pr())))])
  })
  
  ## Results plot inputs in Normalized page
  
  ## PCA Plot
  pca_label_dm_pr<-reactive({
    pca_lable<-levels(as.factor(colData(dep_dm_pr())$replicate))
    print(pca_label)
  })
  
  pca_input_dm_pr<-reactive({
    if (num_total_dm_pr()<=500){
      if(length(levels(as.factor(colData(dep_dm_pr())$replicate))) <= 6){
        pca_plot<- DEP::plot_pca(dep_dm_pr(), n=num_total_dm_pr(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
      else{
        pca_plot<-DEP::plot_pca(dep_dm_pr(), n=num_total_dm_pr(), point_size = 4, indicate = "condition")
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
    }
    else{
      if(length(levels(as.factor(colData(dep_dm_pr())$replicate))) <= 6){
        pca_plot<-DEP::plot_pca(dep_dm_pr(), point_size = 4)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }else{
        pca_label<-SummarizedExperiment::colData(dep_dm_pr())$replicate
        pca_plot<-DEP::plot_pca(dep_dm_pr(), point_size = 4, indicate = "condition")
        #pca_plot<-pca_plot + geom_point()
        pca_plot<-pca_plot + ggrepel::geom_text_repel(aes(label=factor(rowname)),
                                                      size = 4,
                                                      box.padding = unit(0.1, 'lines'),
                                                      point.padding = unit(0.1, 'lines'),
                                                      segment.size = 0.5)
        pca_plot<-pca_plot + labs(title = "PCA plot")
        return(pca_plot)
      }
    }
    
  })
  
  ### Heatmap Differentially expressed proteins
  heatmap_input_dm_pr<-reactive({ 
    get_cluster_heatmap(dep_dm_pr(),
                        type="centered",kmeans = TRUE,
                        k=6, col_limit = 6,
                        indicate = "condition"
    )
  })
  
  ### Volcano Plot
  volcano_input_dm_pr <- reactive({
    if(!is.null(input$volcano_cntrst_dm_pr)) {
      plot_volcano_new(dep_dm_pr(),
                         input$volcano_cntrst_dm_pr,
                         input$check_anova_dm_pr,
                         input$check_names_dm_pr,
                         input$p_adj_dm_pr)
    }
  })
  
  volcano_df_dm_pr<- reactive({
    if(!is.null(input$volcano_cntrst_dm_pr)) {
      get_volcano_df(dep_dm_pr(),
                     input$volcano_cntrst_dm_pr) 
      
    }
  })
  
  
  volcano_input_selected_dm_pr<-reactive({
    if(!is.null(input$volcano_cntrst_dm_pr)){
      
      if (!is.null(input$contents_dm_pr_rows_selected)){
        proteins_selected<-data_result_dm_pr()[c(input$contents_dm_pr_rows_selected),]## get all rows selected
      }
      else if(!is.null(input$protein_brush_dm_pr)){
        proteins_selected<-data_result_dm_pr()[data_result_dm_pr()[["Gene Name"]] %in% protein_name_brush_dm_pr(), ] 
      }
      print(proteins_selected)
      ## convert contrast to x and padj to y
      diff_proteins <- grep(paste(input$volcano_cntrst_dm_pr, "_log2", sep = ""),
                            colnames(proteins_selected))
      if(input$p_adj=="FALSE"){
        padj_proteins <- grep(paste(input$volcano_cntrst_dm_pr, "_p.val", sep = ""),
                              colnames(proteins_selected))
      }
      else{
        padj_proteins <- grep(paste(input$volcano_cntrst, "_p.adj", sep = ""),
                              colnames(proteins_selected))
      }
      
      df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                               y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                               name = proteins_selected$`Gene Name`)
      # print(df_protein)
      p<-plot_volcano_new(dep_dm_pr(),
                            input$volcano_cntrst_dm_pr,
                            input$check_anova_dm_pr,
                            input$check_names_dm_pr,
                            input$p_adj_dm_pr)
      
      p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
        ggrepel::geom_text_repel(data = df_protein,
                                 aes(x, y, label = name),
                                 size = 4,
                                 box.padding = unit(0.1, 'lines'),
                                 point.padding = unit(0.1, 'lines'),
                                 segment.size = 0.5)## use the dataframe to plot points
      
    }
  })
  
  protein_input_dm_pr<-reactive({ 
    
    protein_selected  <- data_result_dm_pr()[input$contents_dm_pr_rows_selected,1]
    protein_selected <-as.character(protein_selected)
    if(length(levels(as.factor(colData(dep_dm_pr())$replicate))) <= 8){
      plot_protein(dep_dm_pr(), protein_selected, as.character(input$type))
    }
    else{
      protein_plot<-plot_protein(dep_dm_pr(), protein_selected, as.character(input$type))
      protein_plot + scale_color_brewer(palette = "Paired")
    }
    
  })
  
  ## QC Inputs
  norm_input_dm_pr <- reactive({
    plot_normalization(processed_data_dm_pr(),
                       normalised_data_dm_pr())
  })
  
  missval_input_dm_pr <- reactive({
    plot_missval(processed_data_dm_pr())
  })
  
  detect_input_dm_pr <- reactive({
    plot_detect(processed_data_dm_pr())
  })
  
  imputation_input_dm_pr <- reactive({
    plot_imputation(processed_data_dm_pr(),
                    diff_all_dm_pr())
  })
  
  p_hist_input_dm_pr <- reactive({
    plot_p_hist(dep_dm_pr())
  })
  
  numbers_input_dm_pr <- reactive({
    plot_numbers(processed_data_dm_pr()) 
  })
  
  coverage_input_dm_pr <- reactive({
    plot_coverage(processed_data_dm_pr())
  })
  
  correlation_input_dm_pr<-reactive({
    plot_cor(dep_dm_pr(),significant = FALSE)
  })
  
  cvs_input_dm_pr<-reactive({
    plot_cvs(dep_dm_pr())
  })
  
  num_total_dm_pr<-reactive({
    dep_dm_pr() %>%
      nrow()
  }) 
  
  ## Enrichment inputs
  
  go_input_dm_pr<-eventReactive(input$go_analysis_dm_pr,{
    withProgress(message = 'Gene ontology enrichment is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    
    if(!is.null(input$contrast_dm_pr)){
      enrichment_output_test(dep_dm_pr(), input$go_database_dm_pr)
      go_results<- test_gsea_mod(dep_dm_pr(), databases = input$go_database_dm_pr, contrasts = TRUE)
      null_enrichment_test(go_results)
      plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast_dm_pr,
                                databases = input$go_database_dm_pr, nrow = 2, term_size = 8) + aes(stringr::str_wrap(Term, 60)) +
        xlab(NULL)
      go_list<-list("go_result"=go_results, "plot_go"=plot_go)
      return(go_list)
    }
  })
  
  pathway_input_dm_pr<-eventReactive(input$pathway_analysis_dm_pr,{
    progress_indicator("Pathway Analysis is running....")
    enrichment_output_test(dep_dm_pr(), input$pathway_database_dm_pr)
    pathway_results<- test_gsea_mod(dep_dm_pr(), databases=input$pathway_database_dm_pr, contrasts = TRUE)
    null_enrichment_test(pathway_results)
    plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_1_dm_pr,
                                  databases=input$pathway_database_dm_pr, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
      xlab(NULL)
    pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
    return(pathway_list)
  })
  
  #### Interactive UI (Normalized page)
  output$significantBox_dm_pr <- renderInfoBox({
    num_total <- dep_dm_pr() %>%
      nrow()
    num_signif <- dep_dm_pr() %>%
      .[SummarizedExperiment::rowData(.)$significant, ] %>%
      nrow()
    frac <- num_signif / num_total
    
    info_box <- 		infoBox("Significant Proteins",
                          paste0(num_signif,
                                 " out of ",
                                 num_total),
                          paste0(signif(frac * 100, digits = 3),
                                 "% of proteins differentially expressed across all conditions"),
                          icon = icon("stats", lib = "glyphicon"),
                          color = "olive",
                          # fill = TRUE,
                          width = 4)
    
    return(info_box)
  })
  
  
  ##### Get results dataframe from Summarizedexperiment object
  data_result_dm_pr<-reactive({
    get_results_proteins(dep_dm_pr(),TRUE)
  })
  
  
  #### Data table
  output$contents_dm_pr <- DT::renderDataTable({
    withProgress(message = 'Result table calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    df<- data_result_dm_pr()
    return(df)
  },
  options = list(scrollX = TRUE,
                 autoWidth=TRUE,
                 columnDefs= list(list(width = '400px', targets = c(-1))))
  )
  
  ## Deselect all rows button
  proxy <- dataTableProxy("contents_dm_pr")
  
  observeEvent(input$clear_dm_pr,{
    proxy %>% selectRows(NULL)
  })
  
  observeEvent(input$original_dm_pr,{
    output$contents_dm_pr <- DT::renderDataTable({
      df<- data_result_dm_pr()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
  })
  
  protein_name_brush_dm_pr<- reactive({
    #protein_tmp<-nearPoints(volcano_df(), input$protein_click, maxpoints = 1)
    protein_tmp<-brushedPoints(volcano_df_dm_pr(), input$protein_brush_dm_pr, 
                               xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  protein_name_click_dm_pr<- reactive({
    protein_tmp<-nearPoints(volcano_df_dm_pr(), input$protein_click_dm_pr, maxpoints = 1)
    # protein_tmp<-brushedPoints(volcano_df(), input$protein_brush, 
    #xvar = "diff", yvar = "p_values")
    protein_selected<-protein_tmp$name
  }) 
  
  
  ## Select rows dynamically
  observeEvent(input$protein_brush_dm_pr,{
    output$contents_dm_pr <- DT::renderDataTable({
      df<- data_result_dm_pr()[data_result_dm_pr()[["Gene Name"]] %in% protein_name_brush_dm_pr(), ]
      return(df)
    },
    options = list(scrollX= TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    proteins_selected<-data_result_dm_pr()[data_result_dm_pr()[["Gene Name"]] %in% protein_name_brush_dm_pr(), ] #
    # get all rows selected
    ## convert contrast to x and padj to y
    diff_proteins <- grep(paste(input$volcano_cntrst_dm_pr, "_log2", sep = ""),
                          colnames(proteins_selected))
    if(input$p_adj=="FALSE"){
      padj_proteins <- grep(paste(input$volcano_cntrst_dm_pr, "_p.val", sep = ""),
                            colnames(proteins_selected))
    }
    else{
      padj_proteins <- grep(paste(input$volcano_cntrst_dm_pr, "_p.adj", sep = ""),
                            colnames(proteins_selected))
    }
    df_protein <- data.frame(x = proteins_selected[, diff_proteins],
                             y = -log10(as.numeric(proteins_selected[, padj_proteins])),#)#,
                             name = proteins_selected$`Gene Name`)
    
    p<-plot_volcano_new(dep_dm_pr(),
                        input$volcano_cntrst_dm_pr,
                        input$check_anova_dm_pr,
                        input$check_names_dm_pr,
                        input$p_adj_dm_pr)
    
    
    p<- p + geom_point(data = df_protein, aes(x, y), color = "maroon", size= 3) +
      ggrepel::geom_text_repel(data = df_protein,
                               aes(x, y, label = name),
                               size = 4,
                               box.padding = unit(0.1, 'lines'),
                               point.padding = unit(0.1, 'lines'),
                               segment.size = 0.5)
    
    output$volcano_dm_pr <- renderPlot({
      withProgress(message = 'Volcano Plot calculations are in progress',
                   detail = 'Please wait for a while', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      p
    })
    return(p)
  })
  
  observeEvent(input$resetPlot_dm_pr,{
    session$resetBrush("protein_brush_dm_pr")
    brush <<- NULL
    
    output$contents_dm_pr <- DT::renderDataTable({
      df<- data_result_dm_pr()
      return(df)
    },
    options = list(scrollX = TRUE,
                   autoWidth=TRUE,
                   columnDefs= list(list(width = '400px', targets = c(-1))))
    )
    
    output$volcano_dm_pr <- renderPlot({
      volcano_input_dm_pr()
    })
  })
  
  observeEvent(input$protein_click_dm_pr,{
    output$contents_dm_pr <- DT::renderDataTable({
      df<- data_result_dm_pr()[data_result_dm_pr()[["Gene Name"]] %in% protein_name_click_dm_pr(), ]
      return(df)
    },
    options = list(scrollX= TRUE)
    )
  })
  
  ## Render Result Plots
  output$pca_plot_dm_pr<-renderPlot({
    pca_input_dm_pr()
  })
  output$heatmap_dm_pr<-renderPlot({
    withProgress(message = 'Heatmap rendering is in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    heatmap_input_dm_pr()
  })
  
  output$volcano_dm_pr <- renderPlot({
    withProgress(message = 'Volcano Plot calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    if(is.null(input$contents_dm_pr_rows_selected)){
      volcano_input_dm_pr()
    }
    else if(!is.null(input$volcano_cntrst_dm_pr)){
      volcano_input_selected_dm_pr()
    } # else close
  })
  
  output$protein_plot_dm_pr<-renderPlot({
    if(!is.null(input$contents_dm_pr_rows_selected)){
      protein_input_dm_pr()
    }
  })
  
  
  ### QC Outputs
  output$sample_corr_dm_pr <-renderPlot({
    correlation_input_dm_pr()
  })
  
  output$sample_cvs_dm_pr <- renderPlot({
    cvs_input_dm_pr()
  })
  
  output$norm_dm_pr <- renderPlot({
    norm_input_dm_pr()
  })
  
  output$missval_dm_pr <- renderPlot({
    missval_input_dm_pr()
  })
  
  output$detect_dm_pr <- renderPlot({
    detect_input_dm_pr()
  })
  
  output$imputation_dm_pr <- renderPlot({
    imputation_input_dm_pr()
  })
  
  # output$p_hist <- renderPlot({
  #   p_hist_input_dm_pr()
  # })
  
  output$numbers_dm_pr <- renderPlot({
    numbers_input_dm_pr()
  })
  
  output$coverage_dm_pr <- renderPlot({
    coverage_input_dm_pr()
  })
  
  ## Enrichment Outputs
  output$go_enrichment_dm_pr<-renderPlot({
    go_input_dm_pr()$plot_go
  })
  
  output$pathway_enrichment_dm_pr<-renderPlot({
    pathway_input_dm_pr()$plot_pa
  })
  
  ##### Download Functions
  datasetInput_dm_pr <- reactive({
    switch(input$dataset_dm_pr,
           "Results" = get_results_proteins(dep_dm_pr()),
           "Full dataset" = get_df_wide(dep_dm_pr()))
  })
  
  output$downloadData_dm_pr <- downloadHandler(
    filename = function() { paste(input$dataset_dm_pr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(datasetInput_dm_pr(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ### === Cluster Download ==== ####
  
  individual_cluster_dm_pr <- reactive({
    cluster_number <- input$cluster_number_dm_pr
    cluster_all <- heatmap_input_dm_pr()
    data_result_dm_pr()[cluster_all[[cluster_number]],]
  })
  
  
  
  output$downloadCluster_dm_pr <- downloadHandler(
    filename = function() { paste("Cluster_info_",input$cluster_number_dm_pr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(individual_cluster_dm_pr(),
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  output$downloadVolcano_dm_pr <- downloadHandler(
    filename = function() {
      paste0("Volcano_", input$volcano_cntrst_dm_pr, ".pdf")
    },
    content = function(file) {
      pdf(file)
      print(volcano_input_selected_dm_pr())
      dev.off()
    }
  )
  
  
  ## Protein plot download
  output$downloadProtein_dm_pr <- downloadHandler(
    filename = function() {
      paste0(input$type_dm_pr,".pdf")
    },
    content = function(file) {
      pdf(file)
      print(protein_input_dm_pr())
      dev.off()
    }
  )
  
  ###### ==== DOWNLOAD GO TABLE ==== ####
  output$downloadGO_dm_pr <- downloadHandler(
    filename = function() { paste("GO_enrichment_",input$go_database_dm_pr, ".csv", sep = "") }, ## use = instead of <-
    content = function(file) {
      write.table(go_input_dm_pr()$go_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  ###### ==== DOWNLOAD PATHWAY TABLE ==== ####
  output$downloadPA_dm_pr <- downloadHandler(
    filename = function() { paste("Pathway_enrichment_",input$pathway_database_dm_pr, ".csv", sep = "") }, 
    ## use = instead of <-
    content = function(file) {
      write.table(pathway_input_dm_pr()$pa_result,
                  file,
                  col.names = TRUE,
                  row.names = FALSE,
                  sep =",") }
  )
  
  
  
  #####===== Download Report (demo protein group)=====#####
  output$downloadReport_dm_pr <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "Phospho-Analyst(demo)report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "ProteinGroup_demo_report.Rmd")
      file.copy("ProteinGroup_demo_report.Rmd", tempReport, overwrite = TRUE)
      
      sig_proteins<-dep_dm_pr() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()
      
      tested_contrasts<- gsub("_significant", "", 
                              colnames(SummarizedExperiment::rowData(dep_dm_pr()))[grep("_significant", 
                                                                                     colnames(SummarizedExperiment::rowData(dep_dm_pr())))])
      pg_width<- ncol(normalised_data_dm_pr()) / 2.5
      # Set up parameters to pass to Rmd document
      params <- list(data = normalized_phospho_data,
                     alpha = input$p,
                     lfc = input$lfc,
                     num_signif= sig_proteins,
                     pg_width = pg_width,
                     tested_contrasts= tested_contrasts,
                     numbers_input= numbers_input_dm_pr,
                     detect_input = detect_input_dm_pr,
                     imputation_input = imputation_input_dm_pr,
                     missval_input = missval_input_dm_pr,
                     p_hist_input = p_hist_input_dm_pr,
                     pca_input = pca_input_dm_pr,
                     coverage_input= coverage_input_dm_pr,
                     correlation_input =correlation_input_dm_pr,
                     heatmap_input = heatmap_input_dm_pr,
                     cvs_input= cvs_input_dm_pr,
                     dep = dep_dm_pr
      )
      
      # Knit the document, passing in the `params` list
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  #### Comparison demo page logic ========== #############
  ####======= Render Functions
  
  # load demo data
  env_dm_c<-reactive({
    LoadToEnvironment("data/exp_demo_data.RData", env = globalenv())
  })
  
  exp_design_demo<-reactive({
    env_dm_c()[["exp_demo"]]
  })

  output$volcano_comp_dm <- renderUI({
    if (!is.null(comparisons_dm())) {
      df <- SummarizedExperiment::rowData(dep_dm())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("volcano_comp_dm",
                     "Comparison",
                     choices = gsub("_significant", "", colnames(df)[cols]))
    }
  })

  output$selected_gene_dm = renderUI({
    if (!is.null(gene_names_dm())){
      selectInput('selected_gene_dm',
                  label='Select one or more Gene names',
                  choices = as.list(gene_names_dm()$Gene.names),
                  selected = NULL,
                  multiple = TRUE)
    }
  })

  # Reactive components 
  phospho_df_dm <- reactive({
    if (!is.null(input$volcano_comp_dm)){
      # phospho_row <- rowData(dep()) %>% as.data.frame()
      phospho_row <- data_result_dm() %>% mutate(rowname = Phosphosite) %>% as.data.frame()
      phospho_row <- column_to_rownames(phospho_row, 'rowname')
      phospho_intensity <- assay(dep_dm())  %>% as.data.frame()
      phospho_df_dm <- merge(phospho_row, phospho_intensity, by = 0) # Merge data according to row names
      
      col_selected <- c(colnames(phospho_intensity),'Phosphosite','Gene.names',
                        paste(input$volcano_comp_dm, "_log2 fold change", sep = ""),
                        paste(input$volcano_comp_dm, "_p.val", sep = ""))
      
      phospho_df_1_dm <- phospho_df_dm %>%
        subset(select = col_selected) %>%
        dplyr::rename(phospho_id = Phosphosite, phospho_diff = paste(input$volcano_comp_dm, "_log2 fold change", sep = ""))
      return(phospho_df_1_dm)
    }
  })

  protein_df_dm <- reactive({
    if (!is.null(input$volcano_comp_dm)){
      # protein_row <- rowData(dep_pr()) %>% as.data.frame()
      data_df <- data_result_dm_pr()
      colnames(data_df)[1]<-c("Gene.names")
      protein_row <- data_df %>% mutate(rowname = Gene.names) %>% as.data.frame()
      protein_row <- column_to_rownames(protein_row, 'rowname')
      protein_intensity <- assay(dep_dm_pr()) %>% as.data.frame()
      protein_df_dm <- merge(protein_row, protein_intensity, by = 0) # Merge data according to row names
      print(colnames(protein_df_dm))
      
      col_selected <- c(colnames(protein_intensity),"Protein ID",'Gene.names',
                        paste(input$volcano_comp_dm, "_log2 fold change", sep = ""),
                        paste(input$volcano_comp_dm, "_p.val", sep = ""))
      
      protein_df_1_dm <- protein_df_dm %>%
        subset(select = col_selected) %>%
        dplyr::rename(protein_diff = paste(input$volcano_comp_dm, "_log2 fold change", sep = "") )
      return(protein_df_1_dm)
    }
  })

  combined_df_dm <- reactive({
    if (!is.null(phospho_df_dm()) & !is.null(protein_df_dm())){
      df <- phospho_df_dm() %>%
        left_join(., protein_df_dm(), by = "Gene.names")
      # df$protein_diff[is.na(df$protein_diff)] <- 0
      df$normalized_diff <- df$phospho_diff - df$protein_diff
      # get index of the p.val
      pval <- grep("_p.val.x",colnames(df))
      df$p_values <- as.numeric(df[, pval])
      
      df <- df %>%
        mutate(p_value_desc = case_when(phospho_diff > 1 & p_values < 0.05 ~ 'Up',
                                        phospho_diff < -1 & p_values < 0.05 ~ 'Down',
                                        TRUE ~ 'Not Sig'))
      df <- df %>%
        mutate(n_p_value_desc = case_when(normalized_diff > 1 & p_values < 0.05 ~ 'Up',
                                          normalized_diff < -1 & p_values < 0.05 ~ 'Down',
                                          TRUE ~ 'Not Sig'))
      return(df)
    }
  })

  # gene names for selection input
  gene_names_dm <- reactive({
    if (!is.null(phospho_df_long_dm()) & !is.null(protein_df_long_dm())){
      gene_names <- phospho_df_long_dm() %>%
        inner_join(., protein_df_long_dm(), by = "Gene.names") %>%
        select('Gene.names') %>%
        unique()
    }
  })

  # phosphosite and protein data
  phospho_df_long_dm <- reactive({
    if (!is.null(phospho_df_dm()) & !is.null(protein_df_dm())){
      combined_df <- phospho_df_dm() %>%
        inner_join(., protein_df_dm(), by = "Gene.names")
      exp_design <- exp_design_demo()
      phospho_cols <- colnames(combined_df[grep('.x', colnames(combined_df))])
      phospho_cols_1 <- exp_design$label
      phospho_df_11 <- subset(combined_df, select = c(phospho_cols,'phospho_id','Gene.names','phospho_diff' ))
      names(phospho_df_11) <- c(phospho_cols_1,"p.val",'phospho_id','Gene.names','phospho_diff')
      phospho_df_22 <- phospho_df_11 %>% rownames_to_column() %>%
        gather(label, intensity, -rowname,-"p.val",-"phospho_id", -"Gene.names", -"phospho_diff")
      phospho_df_33 <- phospho_df_22 %>%
        left_join(., exp_design, by = "label")
      return(phospho_df_33)
    }
  })

  protein_df_long_dm <- reactive({
    if (!is.null(phospho_df_dm()) & !is.null(protein_df_dm())){
      combined_df <- phospho_df_dm() %>%
        inner_join(., protein_df_dm(), by = "Gene.names")
      exp_design <- exp_design_demo()
      protein_cols <- colnames(combined_df[grep('.y', colnames(combined_df))])
      protein_cols_1 <- exp_design$label
      protein_df_11 <- subset(combined_df, select = c(protein_cols,"Protein ID",'Gene.names','protein_diff'))
      names(protein_df_11) <- c(protein_cols_1,"p.val","protein_id",'Gene.names','protein_diff')
      protein_df_22 <- protein_df_11 %>% rownames_to_column() %>%
        gather(label, intensity, -rowname,-"p.val", -"protein_id", -"Gene.names",-"protein_diff")
      protein_df_33 <- protein_df_22 %>%
        left_join(., exp_design, by = "label")
      return(protein_df_33)
      cat(head(protein_df_33)) # test
    }
  })

  # interactive plots
  combined_inter_dm <- reactive({
    if (!is.null(phospho_df_long_dm()) & !is.null(protein_df_long_dm())){
      if (is.null(input$selected_gene_dm)){
        phospho_df <- phospho_df_long_dm() %>% dplyr::filter(Gene.names == gene_names_dm()$Gene.names[1])
        protein_df <- protein_df_long_dm() %>% dplyr::filter(Gene.names == gene_names_dm()$Gene.names[1])
      }
      else {
        phospho_df <- phospho_df_long_dm() %>% dplyr::filter(Gene.names %in% input$selected_gene_dm)
        protein_df <- protein_df_long_dm() %>% dplyr::filter(Gene.names %in% input$selected_gene_dm)
      }
      
      # phospho_df <- phospho_df_long() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      # protein_df <- protein_df_long() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      
      p1 <- phospho_df %>%
        unique() %>%
        ggplot(aes(x = condition, y = intensity)) +
        geom_point(aes(color = factor(replicate)),
                   size = 3) +
        geom_line(aes(group= factor(replicate), color= factor(replicate))) +
        scale_colour_discrete(name  ="Replicate") +
        ylim(19,28) + labs(title = 'Phospho site', x = '') +
        facet_grid(. ~phospho_id + Gene.names) +
        theme(axis.text.x = element_text(angle = 45))
      
      p2 <- protein_df %>%
        select(-'rowname') %>%
        unique() %>%
        ggplot(aes(x = condition, y = intensity)) +
        geom_point(aes(color = factor(replicate)),
                   size = 3) +
        geom_line(aes(group= factor(replicate), color= factor(replicate))) +
        scale_colour_discrete(name  ="Replicate") +
        ylim(19,28) + labs(title = 'Protein', x = '') +
        facet_grid(. ~protein_id + Gene.names) +
        theme(axis.text.x = element_text(angle = 45))
      
      ggarrange(p2,
                p1 +
                  theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank()),
                widths=c(1,4), common.legend = TRUE, legend = 'right')
    }
  })

  # Phosphosite and protein log fold change scatter plot input
  scatter_plot_dm <- reactive({
    if(!is.null(input$volcano_comp_dm) & !is.null(combined_df_dm())){
      df <- combined_df_dm()
      df %>% filter(!is.na(protein_diff))  %>%
        ggplot(aes(x=phospho_diff, y=protein_diff)) +
        geom_point(size = 2, alpha = 0.8) +
        geom_text_repel(
          data=df %>% filter(phospho_diff > 5 | phospho_diff < -5|
                               protein_diff>2| protein_diff < -2), # Filter data first
          aes(label=phospho_id),
          nudge_x = 0.5, nudge_y = 0,
          size = 4) +
        labs(title = paste(input$volcano_comp_dm,'Comparison between phospho and protein log fold change', sep = "\n"),
             x = 'Phosphosite log fold change', y = 'Protein log fold change') +
        theme(plot.title = element_text(hjust = 0.5))
    }
  })

  # Output combined QC plots
  output$pca_plot_c_dm <- renderPlot({
    ggarrange(pca_input_dm() + labs(title = 'Phospho'),
              pca_input_dm_pr() + labs(title = "Protein"),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })

  output$scatter_plot_dm <- renderPlot({
    scatter_plot_dm()
  })

  output$sample_corr_c1_dm <- renderPlot({
    correlation_input_dm()
  })

  output$sample_corr_c2_dm <- renderPlot({
    correlation_input_dm_pr()
  })

  output$sample_cvs_c_dm <- renderPlot({
    ggarrange(cvs_input_dm() + labs(title = 'Phospho'),
              cvs_input_dm_pr() + labs(title = "Protein"),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })

  output$numbers_c_dm <- renderPlot({
    ggarrange(numbers_input_dm() + labs(title = 'Phospho'),
              numbers_input_dm_pr() + labs(title = "Protein"),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })

  output$coverage_c_dm <- renderPlot({
    ggarrange(coverage_input_dm() + labs(title = 'Phospho'),
              coverage_input_dm_pr() + labs(title = "Protein"))
  })

  output$norm_c_dm <- renderPlot({
    ggarrange(norm_input_dm() + labs(title = 'Phospho'),
              norm_input_dm_pr() + labs(title = "Protein") ,
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })


  output$missval_c1_dm <- renderPlot({
    missval_input_dm()
  })

  output$missval_c2_dm <- renderPlot({
    missval_input_dm_pr()
  })

  output$imputation_c_dm <- renderPlot({
    ggarrange(imputation_input_dm() + labs(title = 'Phospho'),
              imputation_input_dm_pr() + labs(title = "Protein"),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })

  # Output interactive plots
  output$combined_inter_dm <- renderPlot({
    combined_inter_dm()
  })

  output$combined_point_dm <- renderPlot({
    if(!is.null(phospho_df_dm())){
      exp_design<- exp_design_demo()
      conditions <- exp_design$condition %>% unique()
      # phospho_df <- phospho_df() %>% dplyr::filter(Gene.names %in% input$selected_gene)
      
      if (is.null(input$selected_gene_dm)){
        phospho_df <- phospho_df_dm() %>% dplyr::filter(Gene.names == gene_names_dm()$Gene.names[1])
      }
      else {
        phospho_df <- phospho_df_dm() %>% dplyr::filter(Gene.names %in% input$selected_gene_dm)
      }
      phospho_df_1 <- phospho_df
      print(colnames(phospho_df_1)) # test
      for (i in 1:length(conditions)) {
        condition <- conditions[i]
        pattern <- paste(condition,"[[:digit:]]",sep = '_')
        phospho_df_1[paste0('mean',sep = "_",condition)] <- rowMeans(phospho_df_1 %>% select(grep(pattern, colnames(phospho_df_1))), na.rm = TRUE)
      }
      
      phospho_df_2 <- phospho_df_1 %>%
        # dplyr::filter(Gene.names %in% input$selected_gene) %>%
        select(Gene.names,phospho_id, grep('mean', colnames(phospho_df_1))) %>%
        pivot_longer(names_to = "group", values_to = "mean_intensity", cols = starts_with('mean'))
      
      phospho_df_2 %>% ggplot(aes(x = phospho_id, y = group)) +
        geom_point(aes(size = mean_intensity, color = group)) +
        facet_grid(rows = vars(Gene.names), scales = "free_y", space = "free_y") +
        coord_flip() +
        theme(axis.text.x = element_text(angle = 45))
    }
  })
  
  # # used for save demo data
  # observeEvent(input$analyze ,{
  #   if(input$analyze==0 ){
  #     return()
  #   }
  # 
  #   data_missval <- processed_data()
  #   data_dep <- dep()
  #   save(data_missval, data_dep, file = "phosphosite_demo_data.RData")
  # 
  # })
  # 
  # observeEvent(input$analyze ,{
  #   if(input$analyze==0 ){
  #     return()
  #   }
  # 
  #   data_missval <- processed_data_pr()
  #   data_dep <- dep_pr()
  #   save(data_missval, data_dep, file = "proteinGroup_demo_data.RData")
  # 
  # })
  # 
  # observeEvent(input$analyze ,{
  #   if(input$analyze==0 ){
  #     return()
  #   }
  #   
  #   exp_demo <- exp_design_input()
  #   save(exp_demo, file = "exp_demo_data.RData")
  #   
  # })

  
}
