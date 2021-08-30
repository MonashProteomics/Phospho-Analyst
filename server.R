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
  
  observeEvent(input$analyze ,{ 
    if (is.null(input$file1)){
      hideTab(inputId = "panel_list", target = "PhosphoPage")
      hideTab(inputId = "panel_list", target = "ComparisonPage")
    }
    
    if (is.null(input$file2)){
      hideTab(inputId = "panel_list", target = "ProteinPage")
      hideTab(inputId = "panel_list", target = "ComparisonPage")
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
               timer = 10000) # timer in miliseconds (10 sec)
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
      selectizeInput("contrast",
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
  
  
  ## Read input files on shiny server
  ## NOTE: have to use reactive framework, otherwise throws out error
  # observeEvent(input$analyze,{
  #   maxquant_data<-reactive({
  #     inFile<-input$file1
  #     if(is.null(inFile))
  #       return(NULL)
  #     read.table(inFile$datapath,
  #                header = TRUE,
  #                fill= TRUE, # to fill any missing data
  #                sep = "\t"
  #     )
  #   })
  # })
  
  ## make reactive elements
  phospho_data_input<-reactive({NULL})
  protein_data_input<-reactive({NULL})
  
  exp_design_input<-reactive({NULL})
  # exp_design_example<-reactive({NULL})
  # maxquant_data_example<-reactive({NULL})
  
  phospho_data_input<-eventReactive(input$analyze,{
    inFile<-input$file1
    if(is.null(inFile))
      return(NULL)
    temp_data<-read.table(inFile$datapath,
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
    temp_data1<-read.table(inFile$datapath,
                           header = TRUE,
                           fill= TRUE, # to fill any missing data
                           sep = "\t"
    )
    #validate(maxquant_input_test(temp_data))
    return(temp_data1)
  })
  
  exp_design_input<-eventReactive(input$analyze,{
    inFile<-input$file3
    if (is.null(inFile))
      return(NULL)
    temp_df<-read.table(inFile$datapath,
                        header = TRUE,
                        sep="\t",
                        stringsAsFactors = FALSE)
    exp_design_test(temp_df)
    temp_df$label<-as.character(temp_df$label)
    temp_df$condition<-trimws(temp_df$condition, which = "left")
    return(temp_df)
  })
  
  
  
  
  
  ### Reactive components
  processed_data<- reactive({
    ## check which dataset
    if(!is.null (phospho_data_input() )){
      phospho_data <- reactive({phospho_data_input()})
    }
    
    if(!is.null (exp_design_input() )){
      exp_design<-reactive({exp_design_input()})
    }
    
    
    message(exp_design())

    if(grepl('+',phospho_data()$Reverse)){
      filtered_data<-dplyr::filter(phospho_data(),Reverse!="+")
    }
    else{filtered_data<-phospho_data()}
    if(grepl('+',filtered_data$Potential.contaminant)){
      filtered_data<-dplyr::filter(filtered_data,Potential.contaminant!="+")
    }

    filtered_data<-ids_test(filtered_data)

    
    # 
    # filtered_data <- phospho_data() %>% dplyr::filter(Reverse != '+') %>%
    #   dplyr::filter(Potential.contaminant != '+') %>%
    #   dplyr::select(- `Reverse`,-`Potential.contaminant`)
    # filtered_data<-ids_test(filtered_data)
    
    # 5. Expand Site table
    # get all intensity columns
    intensity <- grep("^Intensity.+|Intensity", colnames(filtered_data)) 
    # get the required intensity columns
    intensity_cols <- grep("^Intensity.+___\\d", colnames(filtered_data))
    intensity_names <- colnames( filtered_data[,intensity_cols])
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
    
    # 9. Convert the data into SummarisedExperiment object.
    data_pre <- data_pre %>% 
      mutate(name = paste(Gene.names,Positions.within.proteins, Multiplicity,sep = '_'))
    data_pre <- data_pre %>% 
      mutate(ID = paste(id,Multiplicity, sep = '_'))
    
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
  
  normalised_data<-reactive({
    print(dim(assay(processed_data())))
    # cat('Number of non-numeric Score.diff:',count(is.numeric(rowData(processed_data())$Score.diff) == FALSE),'\n') # test
    # cat('Non-numeric index:',which(is.na(as.numeric(rowData(processed_data())$Score.diff))),'\n')  # test
    
    normalize_vsn(processed_data())
  })
  
  imputed_data<-reactive({
    
    DEP::impute(normalised_data(),input$imputation)
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
    test_diff(imputed_data(),type = 'all')
  })
  
  dep<-reactive({
    # cat('Number of non-numeric Score.diff:',count(is.numeric(rowData(imputed_data())$Score.diff) == FALSE),'\n') # test
    # cat('Non-numeric index:',which(is.na(as.numeric(rowData(imputed_data())$Score.diff))),'\n')  # test
    
    if(input$fdr_correction=="BH"){
      diff_all<-test_limma(imputed_data(),type='all', paired = input$paired)
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    else{
      diff_all<-test_diff(imputed_data(),type='all')
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    
    if(length(unique(exp_design_input()$condition)) <= 2){
      return(diff_all_rej)
      
    }
    else if(length(unique(exp_design_input()$condition)) >= 3){
      anova_dep <- diff_all_rej
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
        do(anova_function(.))  %>% dplyr::select(p.value) %>%
        ungroup()
      
      colnames(anova)<-c("name", "anova_p.val")
      
      # add anova value to row data
      rowData(anova_dep) <- merge(rowData(anova_dep), anova, by = 'name', sort = FALSE)
      return(anova_dep)
    }
    
  })

  comparisons<-reactive ({
    temp<-capture.output(DEP::test_diff(imputed_data(),type='all'),type = "message")
    gsub(".*: ","",temp)
    ## Split conditions into character vector
    unlist(strsplit(temp,","))
    ## Remove leading and trailing spaces
    trimws(temp)
  })
  
  ## Select point on volcano plot
  # protein_graph_selected<- reactive({
  #   protein_row<-nearPoints(data_result(), input$protein_click,
  #                           maxpoints = 1)
  #  # as.character(protein_row$name)
  # })
  # 
  
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
    plot_imputation(normalised_data(),
                    diff_all())
  })
  
  p_hist_input <- reactive({
    plot_p_hist(dep())
  })
  
  numbers_input <- reactive({
    plot_numbers(normalised_data()) +
      labs(title= "Phosphosites per sample", y = "Number of phosphosites")
  })
  
  coverage_input <- reactive({
    plot_coverage(normalised_data())+
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
      go_results<- test_gsea_mod(dep(), databases = input$go_database, contrasts = TRUE)
      null_enrichment_test(go_results)
      plot_go<- plot_enrichment(go_results, number = 5, alpha = 0.05, contrasts =input$contrast,
                                databases = input$go_database, nrow = 2, term_size = 8) + aes(stringr::str_wrap(Term, 60)) +
        xlab(NULL)
      go_list<-list("go_result"=go_results, "plot_go"=plot_go)
      return(go_list)
    }
  })
  
  pathway_input<-eventReactive(input$pathway_analysis,{
    progress_indicator("Pathway Analysis is running....")
    enrichment_output_test(dep(), input$pathway_database)
    pathway_results<- test_gsea_mod(dep(), databases=input$pathway_database, contrasts = TRUE)
    null_enrichment_test(pathway_results)
    plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_1,
                                  databases=input$pathway_database, nrow = 3, term_size = 8) + aes(stringr::str_wrap(Term, 30)) +
      xlab(NULL)
    pathway_list<-list("pa_result"=pathway_results, "plot_pa"=plot_pathway)
    return(pathway_list)
  })
  
  
  #### Interactive UI
  output$significantBox <- renderInfoBox({
    num_total <- dep() %>%
      nrow()
    num_signif <- dep() %>%
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
  data_result<-reactive({
    if(length(unique(exp_design_input()$condition)) <= 2){
      get_results_phospho(dep(),FALSE)
    } else {
      get_results_phospho(dep(),TRUE)
    }
  })
  
  
  #### Data table
  output$contents <- DT::renderDataTable({
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
  #   observeEvent(input$protein_brush,{
  # output$protein_info<-renderPrint({
  # #  protein_selected()
  #   #nearPoints(rowData(dep()), input$protein_click, maxpoints = 1)
  #   brushedPoints(volcano_df(), input$protein_brush, 
  #               xvar = "diff", yvar = "p_values")
  #  # head(volcano_df())
  #   #input$protein_click
  #  # str(input$protein_hover)
  # })
  #  })
  
  ## Select rows dynamically
  
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
    # p<-plot_volcano_new(dep(),
    #                     input$volcano_cntrst,
    #                     input$check_anova,
    #                     input$check_names,
    #                     input$p_adj)
    
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
  
  output$p_hist <- renderPlot({
    p_hist_input()
  })
  
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
  
  output$pathway_enrichment<-renderPlot({
    pathway_input()$plot_pa
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
  
  #####===== Download Report =====#####
  output$downloadReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "LFQ-Analyst_report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "LFQ_report.Rmd")
      file.copy("LFQ_report.Rmd", tempReport, overwrite = TRUE)
      
      sig_proteins<-dep() %>%
        .[SummarizedExperiment::rowData(.)$significant, ] %>%
        nrow()
      
      tested_contrasts<- gsub("_p.adj", "", 
                              colnames(SummarizedExperiment::rowData(dep()))[grep("p.adj", 
                                                                                  colnames(SummarizedExperiment::rowData(dep())))])
      pg_width<- ncol(imputed_data()) / 2.5
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
  
  
  
  #### protein page logic ========== #############
  
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
  
  output$contrast_pr_1 <- renderUI({
    if (!is.null(comparisons_pr())) {
      df <- SummarizedExperiment::rowData(dep_pr())
      cols <- grep("_significant$",colnames(df))
      selectizeInput("contrast_pr_1",
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
  
  
  dep_pr<-reactive({
    if(input$fdr_correction=="BH"){
      diff_all<-test_limma(imputed_data_pr(),type='all', paired = input$paired)
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    else{
      diff_all<-test_diff(imputed_data_pr(),type='all')
      diff_all_rej <- add_rejections(diff_all,alpha = input$p, lfc= input$lfc)
    }
    
    if(length(unique(exp_design_input()$condition)) <= 2){
      return(diff_all_rej)
      
    }
    else if(length(unique(exp_design_input()$condition)) >= 3){
      anova_dep <- diff_all_rej
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
        do(anova_function(.))  %>% dplyr::select(p.value) %>%
        ungroup()
      
      colnames(anova)<-c("name", "anova_p.val")
      
      # add anova value to row data
      rowData(anova_dep) <- merge(rowData(anova_dep), anova, by = 'name', sort = FALSE)
      return(anova_dep)
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
      if(length(unique(exp_design_input()$condition)) <= 2){
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
      if (!is.null(input$contents_pr_rows_selected)){
        proteins_selected<-data_result_pr()[c(input$contents_pr_rows_selected),]## get all rows selected
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
      if(length(unique(exp_design_input()$condition)) <= 2){
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
  
  processed_data_pr<-reactive({
    ## check which dataset
    if(!is.null (protein_data_input() )){
      protein_data <- reactive({protein_data_input()})
    }
    
    if(!is.null (exp_design_input() )){
      exp_design<-reactive({exp_design_input()})
    }
    
    
    message(exp_design())
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
    
    # # test if there has no Gene names column
    # if(is.null(filtered_data$Gene.names)){
    #   filtered_data <- filtered_data %>% 
    #     dplyr::mutate(Gene.names = filtered_data$Fasta.headers %>% str_extract_all('GN=[:alnum:]+',simplify = FALSE) %>% gsub("GN=","",.))
    # }
    # if(is.null(filtered_data$Protein.names)){
    #   filtered_data <- filtered_data %>% 
    #     dplyr::mutate(Protein.names = filtered_data$Fasta.headers %>% 
    #                     str_extract_all('^.*HUMAN.*;',simplify = FALSE) %>% sub('^.*HUMAN ',"",.) %>% gsub(" OS=.*","",.))
    # }
    # 
    data_unique<- DEP::make_unique(filtered_data,"Gene.names","Protein.IDs",delim=";")
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
    
    if(length(unique(exp_design_input()$condition)) <= 2){
      return(diff_all_rej)
      
    }
    else if(length(unique(exp_design_input()$condition)) >= 3){
      anova_dep <- diff_all_rej
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
        do(anova_function(.))  %>% dplyr::select(p.value) %>%
        ungroup()
      
      colnames(anova)<-c("name", "anova_p.val")
      
      # add anova value to row data
      rowData(anova_dep) <- merge(rowData(anova_dep), anova, by = 'name', sort = FALSE)
      return(anova_dep)
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
    plot_imputation(normalised_data_pr(),
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
    plot_pathway<-plot_enrichment(pathway_results, number = 5, alpha = 0.05, contrasts =input$contrast_pr_1,
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
    if(length(unique(exp_design_input()$condition)) <= 2){
      get_results_proteins(dep_pr(),FALSE)
    } else {
      get_results_proteins(dep_pr(),TRUE)
    }
  })
  
  
  #### Data table
  output$contents_pr <- DT::renderDataTable({
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
    if(length(unique(exp_design_input()$condition)) <= 2) {
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
  
  output$p_hist <- renderPlot({
    p_hist_input_pr()
  })
  
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
  
  
  
  #####===== Download Report =====#####
  output$downloadReport_pr <- downloadHandler(
    
    filename = "LFQ-Analyst_report.pdf",
    content = function(file) {
      file.copy("www/LFQ-Analyst_report.pdf",file)   }
  )
  
  # Comparison page plots
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
  
  phospho_df <- reactive({
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
  })
  
  protein_df <- reactive({
    # protein_row <- rowData(dep_pr()) %>% as.data.frame()
    data_df <- data_result_pr()
    colnames(data_df)[1]<-c("Gene.names")
    protein_row <- data_df %>% mutate(rowname = Gene.names) %>% as.data.frame()
    protein_row <- column_to_rownames(protein_row, 'rowname')
    protein_intensity <- assay(dep_pr()) %>% as.data.frame()
    protein_df <- merge(protein_row, protein_intensity, by = 0) # Merge data according to row names
    print(colnames(protein_df))
    if(length(unique(exp_design_input()$condition)) <= 2) {
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
  })
  combined_df <- reactive({
    df <- phospho_df() %>%
      left_join(., protein_df(), by = "Gene.names")
    df$protein_diff[is.na(df$protein_diff)] <- 0
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
  })
  
  # gene names for selection input
  gene_names <- reactive({
    gene_names <- phospho_df() %>%
      inner_join(., protein_df(), by = "Gene.names") %>% 
      select('Gene.names') %>%
      unique()
  })
  
  # phosphosite and protein data 
  phospho_df_long <- reactive({
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
  })
  
  protein_df_long <- reactive({
    combined_df <- phospho_df() %>%
      inner_join(., protein_df(), by = "Gene.names")
    exp_design <- exp_design_input()
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
  })
  
  # interactive plots
  combined_inter <- reactive({
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
    
  })
  
  # volcano plots
  volcano_phospho <- reactive({
    df <- combined_df() 
    print(colnames(df))  # test
    if(length(unique(exp_design_input()$condition)) <= 2) {
      pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
    } else {
      if (input$check_anova_comp =="FALSE") {
        pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
      } else {
        pval <- grep("ANOVA_p.val.x",colnames(df))
      }
    }
    print(pval) # test
    df <- df %>% dplyr::mutate(p_values = as.numeric(df[, pval]))
    
    df %>% ggplot(aes(x = phospho_diff, y = -log10(p_values))) +
      geom_point(aes(color = p_value_desc)) +
      labs(x = 'Phosphosite log fold change',
           y = '-log10(p-value)') +
      geom_text_repel(
        data=df %>% filter(p_value_desc != 'Not Sig'), # Filter data first
        aes(label=phospho_id),
        nudge_x = 0.5, nudge_y = 0,
        size = 4) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#003399", "#999999", "#b30000")) 
  })
  
  volcano_phospho_2 <- reactive({
    df <- combined_df() 
    if(length(unique(exp_design_input()$condition)) <= 2) {
      pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
    } else {
      if (input$check_anova_comp =="FALSE") {
        pval <- grep(paste(input$volcano_comp, "_p.val.x", sep = ""),colnames(df))
      } else {
        pval <- grep("ANOVA_p.val.x",colnames(df))
      }
    }
    
    df$p_values <- as.numeric(df[, pval])
    
    df %>% ggplot(aes(x = phospho_diff, y = -log10(p_values))) +
      geom_point(aes(color = n_p_value_desc)) +
      labs(x = 'Phosphosite log fold change',
           y = '-log10(p-value)') +
      geom_text_repel(
        data=df %>% filter(n_p_value_desc != 'Not Sig'), # Filter data first
        aes(label=phospho_id),
        nudge_x = 0.5, nudge_y = 0,
        size = 4)+
      scale_color_manual(values = c("#003399", "#999999", "#b30000"))
  })
  
  scatter_plot <- reactive({
    df <- combined_df()
    df %>% 
      ggplot(aes(x=phospho_diff, y=protein_diff)) + 
      geom_point(size = 2, alpha = 0.8) +
      geom_text_repel( 
        data=df %>% filter(phospho_diff > 5 | phospho_diff < -5|
                             protein_diff>2| protein_diff < -2), # Filter data first
        aes(label=phospho_id),
        nudge_x = 0.5, nudge_y = 0,
        size = 4) + 
      labs(title = 'Comparison between phospho and protein log fold change',
           x = 'Phosphosite log fold change', y = 'Protein log fold change') +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$volcano_phospho <-renderPlot({
    volcano_phospho()
  })
  output$volcano_phospho_2 <-renderPlot({
    volcano_phospho_2()
  })
  # combined qc plots
  output$pca_plot_c <- renderPlot({
    ggarrange(pca_input() + labs(title = 'Phospho'), 
              pca_input_pr() + labs(title = "Protein"), 
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })
  
  output$scatter_plot <- renderPlot({
    scatter_plot()
  })
  
  # output$sample_corr_c <- renderPlot({
  #     plot_grid(correlation_input(), correlation_input_pr())
  # })
  
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
  
  # output$missval_c <- renderPlot({
  #     plot_grid(missval_input() + labs(title = 'Phospho'), 
  #               missval_input_pr() + labs(title = "Protein"))
  # })
  
  output$missval_c1 <- renderPlot({
    missval_input()
  })
  
  output$missval_c2 <- renderPlot({
    missval_input_pr()
  })
  
  output$imputation_c <- renderPlot({
    ggarrange(imputation_input() + labs(title = 'Phospho'), 
              imputation_input_pr(),
              widths=c(1,1), common.legend = TRUE, legend = 'right')
  })
  
  output$combined_inter <- renderPlot({
    combined_inter()
  })
  
  output$combined_point <- renderPlot({
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
  })
  
  output$selected_gene = renderUI({
    selectInput('selected_gene', 
                label='Select one or more Gene names', 
                choices = as.list(gene_names()$Gene.names),
                selected = NULL,
                multiple = TRUE) 
  })
  
  
  
}
