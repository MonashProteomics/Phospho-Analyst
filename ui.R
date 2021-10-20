# Create the theme
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#434C5E"
  ),
  adminlte_sidebar(
    width = "250px",
    dark_bg = "#bec8da",
    dark_hover_bg = "#81A1C1",
    dark_color = "#2E3440"
  ),
  adminlte_global(
    content_bg = "#FFF",
    box_bg = "#ffffff", 
    info_box_bg = "#ffffff"
  )
)

# Define UI for data upload app ----
ui <- function(request){
  shinyUI(
    dashboardPage(
      dashboardHeader(title = "Phospho-Analyst"),
      dashboardSidebar(
        useShinyalert(),
        sidebarMenu(
          id="tabs_selected",
          convertMenuItem(
            menuItem('Home', icon=icon("home"), selected = TRUE, tabName = "home"), tabName = 'home'),
          convertMenuItem(menuItem("Analysis",  tabName="analysis", icon=icon("flask"),
                                   fileInput('file1',
                                             p('Upload MaxQuant Phospho STY.txt', style = 'color:#2E3440'),
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv')),
                                   fileInput('file2',
                                             p('Upload MaxQuant ProteinGroup.txt', style = 'color:#2E3440'),
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv')),
                                   fileInput('file3',
                                             p('Upload Experimental Design Matrix', style = 'color:#2E3440'),
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv')),
                                   tags$hr(),
                                   menuItem("Advanced Options",tabName="advanced", icon = icon("cogs"), 
                                            numericInput("p", 
                                                         p("Adjusted p-value cutoff", style = 'color:#2E3440'),
                                                         min = 0.0001, max = 0.1, value = 0.05),
                                            numericInput("lfc",
                                                         p("Log2 fold change cutoff", style = 'color:#2E3440'),
                                                         min = 0, max = 10, value = 1),
                                            checkboxInput("paired",
                                                          p("Paired test", style = 'color:#2E3440'), FALSE),
                                            
                                            prettyRadioButtons("imputation",
                                                               p("Imputation type", style = 'color:#2E3440'),
                                                               choices = c("Perseus-type"="man", MsCoreUtils::imputeMethods())[1:9],
                                                               selected = "man"),
                                            
                                            prettyRadioButtons("fdr_correction",
                                                               p("Type of FDR correction", style = 'color:#2E3440'),
                                                               choices =  c("Benjamini Hochberg"="BH",
                                                                            "t-statistics-based"="fdrtool"
                                                               ), selected= "BH"),
                                            checkboxInput("single_peptide",
                                                          p("Include single peptide identifications", style = 'color:#2E3440'), FALSE),
                                            numericInput("k_number",
                                                         p("Number of clusters in heatmap", style = 'color:#2E3440'),
                                                         min = 1, max = 20, value = 6)
                                   ),
                                   tags$hr(),
                                   actionButton("analyze", "Start Analysis"),
                                   tags$hr(),
                                   p(a("Example Phosphoproteomics data", target= "_blank",
                                       href="data/Phospho (STY)Sites_example.txt", 
                                       download="Phospho (STY)Sites_example.txt")),
                                   p(a("Example ProteinGroup data", target= "_blank",
                                       href="data/proteinGroups_test_phospho_example.txt", 
                                       download="proteinGroups_test_phospho_example.txt")),
                                   p(a("Example Experimental Design file", target= "_blank",
                                       href="data/experimental_design_example.txt", 
                                       download="experimental_design_example.txt"))
          ), tabName = 'analysis'),
          convertMenuItem(
            menuItem('User Guide', icon=icon("question"), 
                     #href = "https://monashbioinformaticsplatform.github.io/LFQ-Analyst/", 
                     tabName = "info"), tabName = 'info' )
        )
      ), # sidebar close
      dashboardBody(
        use_theme(mytheme),
        useShinyjs(), #imp to use shinyjs functions
        tags$head(includeScript("google_analytics.js")),
        
        tabItems(
          tabItem(tabName = "home",
                  fluidRow( 
                    box(
                      title = "Overview",
                      h3("Phospho-Analyst: An easy-to-use interactive web-platform to analyze and visualize phosphoproteomics data 
                     preprocessed with MaxQuant."),
                      p("Phospho-Analyst is an easy-to-use, interactive web application developed to perform 
                  differential expression analysis with “one click” and to visualize label-free quantitative phosphoproteomics 
                  datasets preprocessed with MaxQuant.  Phospho-Analyst provides a wealth of user-analytic features 
                  and offers numerous publication-quality result output graphics and tables to facilitate statistical 
                  and exploratory analysis of label-free quantitative datasets. "), 
                      br(),
                      HTML('<center><img src="./Phospho_analyst.svg" width="600px"></center>'),
                      br(),
                      h4("Sidebar tabs"),
                      tags$ul(
                        tags$li(tags$b("Analysis: "),"perform your own analysis"), 
                        # tags$li(tags$b("Demo: "),"familiarise yourself with Phospho-Analyst by browsing through pre-analysed results"), 
                        tags$li(tags$b("User Guide: "), "download an in-depth manual") 
                      ),
                      width = 12,
                      solidHeader = TRUE,
                      #status = "info"
                      status = "primary"
                    )#box 1 closed
                    
                  ) #fluidrow close
          ), # home tab close
          tabItem(tabName = "analysis",
                  div(id="quickstart_info",
                      fluidPage(
                        box(
                          title = "Getting Started",
                          h3(tags$b(span("Quick Start", style="text-decoration:underline"))),
                          tags$ul(
                            tags$li("Upload your ", tags$b("Phospho SYT Sites.txt "), "generated by MaxQuant."),
                            tags$li("Upload your ", tags$b(" Protein Group.txt "),"generated by MaxQuant. "), 
                            tags$li("Upload your ", tags$b(" experimental design "),"table. "), 
                            "  (",tags$b("Phospho SYT Sites.txt "),'and',tags$b(" Protein Group.txt "),'could be optional to upload)',
                            tags$li(tags$b("Optional: "),"Adjust the p-value cut-off, the log2 fold change cut-off, 
                                                                         the imputation type, FDR correction method and/or number of clusters in heatmap in the",
                                    tags$b("Advanced Options")),
                            tags$li("Press ", tags$b("'Start Analysis' ")), 
                            tags$li(tags$b("Hint: "), " Use the ", tags$b("User Guide ")," tab for a detailed explanation of inputs,
                                                                         advanced options and outputs"), 
                            tags$li(tags$b("Note: "), " The experimental design file is not the" , 
                                    tags$b("'mqpar.xml' "),"file 
                                                                 from MaxQuant. Use the example file template provided.")
                          ),
                          br(),
                          HTML('<center><img src="./Phospho_analyst.svg" width="500px"></center>'),
                          width = 12,
                          solidHeader = TRUE,
                          status = "primary"
                        )
                      )
                  ), # QUICKSTART INFO CLOSE
                  shinyjs::hidden(
                    div(id = 'panels',
                        box(width = 12,
                            tabsetPanel(id = 'panel_list',
                                        tabPanel("Phosphosite",
                                                 fluidRow(
                                                   column(4,
                                                          column(7,uiOutput("downloadTable"),offset = 1),
                                                          column(4,uiOutput("downloadButton"))
                                                   ),
                                                   column(5,
                                                          infoBoxOutput("significantBox",width = 12)
                                                   ),
                                                   column(3,
                                                          uiOutput("downloadreport")
                                                   )
                                                   
                                                 ), # fluidRow closed
                                                 
                                                 # box(
                                                 #   column(7,uiOutput("downloadTable"),offset = 1),
                                                 #   column(4,uiOutput("downloadButton")), # make the button on same line
                                                 #   width = 4),
                                                 # infoBoxOutput("significantBox",width = 4),
                                                 # box(
                                                 #   column(4,uiOutput("downloadreport")), # offset for dist between buttons
                                                 #   #tags$br(),
                                                 #   #column(5,uiOutput('downloadPlots')),
                                                 #   width = 4)),
                                                 # align save button
                                                 tags$style(type='text/css', "#downloadTable { width:100%; margin-top: 15px;}"),
                                                 tags$style(type='text/css', "#downloadButton { width:100%; margin-top: 40px;}"), 
                                                 tags$style(type='text/css', "#downloadreport { width:100%; vertical-align- middle; margin-top: 40px;
                                                                    margin-bottom: 25px;}"),
                                                 tags$br(), # Blank lines
                                                 tags$br(),
                                                 ## Data table and result plots box
                                                 fluidRow(
                                                   box(
                                                     title = "Results Table",
                                                     DT::dataTableOutput("contents"),
                                                     #  actionButton("clear", "Deselect Rows"),
                                                     actionButton("original", "Refresh Table"),
                                                     width = 6,
                                                     status = "success",
                                                     #color=""
                                                     solidHeader = TRUE
                                                   ),
                                                   # column(
                                                   box(
                                                     width= 6,
                                                     collapsible = TRUE,
                                                     #status="primary",
                                                     #solidHeader=TRUE,
                                                     tabBox(
                                                       title = "Result Plots",
                                                       width = 12,
                                                       tabPanel(title = "Volcano plot",
                                                                fluidRow(
                                                                  box(width = 12,
                                                                    column(8,uiOutput("volcano_cntrst")),
                                                                    
                                                                    column(4,
                                                                           prettyCheckbox("check_anova",
                                                                                          "Apply ANOVA",
                                                                                          value = FALSE),
                                                                           prettyCheckbox("check_names",
                                                                                          "Display names",
                                                                                          value = FALSE),
                                                                           prettyCheckbox("p_adj",
                                                                                          "Adjusted p values",
                                                                                          value = FALSE)
                                                                           )
                                                                  ),
                                                                  tags$p("Select phosphosite from Results Table to highlight on the plot OR 
                                                  drag the mouse on plot to show expression of phosphosite in Table, ANOVA function only worked 
                                                                         for more than two groups.")
                                                                  #Add text line
                                                                  # tags$p("OR"),
                                                                  #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
                                                                ),
                                                                
                                                                fluidRow(
                                                                  plotOutput("volcano", height = 600,
                                                                             # hover = "protein_hover"),
                                                                             #),
                                                                             # click = "protein_click"),
                                                                             brush = "protein_brush",
                                                                             click = "protein_click"),
                                                                  downloadButton('downloadVolcano', 'Save Highlighted Plot'),
                                                                  actionButton("resetPlot", "Clear Selection")
                                                                  #)),
                                                                )),
                                                       tabPanel(title= "Heatmap",
                                                                fluidRow(
                                                                  plotOutput("heatmap", height = 600)
                                                                ),
                                                                fluidRow(
                                                                  box(numericInput("cluster_number",
                                                                                   "Cluster to download",
                                                                                   min=1, max=6, value = 1), width = 6),
                                                                  box(downloadButton('downloadCluster',"Save Cluster"),
                                                                      downloadButton('download_hm_svg', "Save svg"),
                                                                      width = 5)
                                                                )
                                                       ),
                                                       tabPanel(title = "Individual Plot",
                                                                fluidRow(
                                                                  box(prettyRadioButtons("type",
                                                                                         "Plot type",
                                                                                         choices = c("Box Plot"= "boxplot",
                                                                                                     "Violin Plot"="violin", 
                                                                                                     "Interaction Plot"= "interaction",
                                                                                                     "Intensity Plot"="dot"
                                                                                         ),
                                                                                         selected = "boxplot", 
                                                                                         inline = TRUE),
                                                                      width = 12
                                                                  ),
                                                                  tags$p("Select one or more rows from Results Table to plot individual 
                                                  phosphosite intesities across conditions and replicates")
                                                                ),
                                                                fluidRow(
                                                                  plotOutput("protein_plot"),
                                                                  downloadButton('downloadProtein', 'Download Plot')
                                                                )
                                                       )
                                                       # verbatimTextOutput("protein_info"))
                                                     )
                                                   ) # box or column end
                                                   
                                                   
                                                   
                                                 ), # result plot colsed
                                                 
                                                 ## QC Box
                                                 fluidRow(
                                                   div(id="qc_tab",
                                                       column(
                                                         width=6,
                                                         tabBox(title = "QC Plots", width = 12,
                                                                tabPanel(title = "PCA Plot",
                                                                         plotOutput("pca_plot", height=600),
                                                                         downloadButton('download_pca_svg', "Save svg")
                                                                ),
                                                                tabPanel(title="Sample Correlation",
                                                                         plotOutput("sample_corr", height = 600),
                                                                         downloadButton('download_corr_svg', "Save svg")
                                                                ),
                                                                tabPanel(title= "Sample CVs",
                                                                         plotOutput("sample_cvs", height = 600),
                                                                         downloadButton('download_cvs_svg', "Save svg")
                                                                ),
                                                                tabPanel(title = "Phosphosite Numbers",
                                                                         plotOutput("numbers", height = 600),
                                                                         downloadButton('download_num_svg', "Save svg")
                                                                ),
                                                                
                                                                tabPanel(title = "Sample coverage",
                                                                         plotOutput("coverage", height = 600),
                                                                         downloadButton('download_cov_svg', "Save svg")
                                                                ),
                                                                tabPanel(title = "Normalization",
                                                                         plotOutput("norm", height = 600),
                                                                         downloadButton('download_norm_svg', "Save svg")
                                                                ),
                                                                # tabPanel(title = "Missing values - Quant",
                                                                #          plotOutput("detect", height = 600)
                                                                # ),
                                                                tabPanel(title = "Missing values - Heatmap",
                                                                         plotOutput("missval", height = 600),
                                                                         downloadButton('download_missval_svg', "Save svg")
                                                                ),
                                                                tabPanel(title = "Imputation",
                                                                         plotOutput("imputation", height = 600),
                                                                         downloadButton('download_imp_svg', "Save svg")
                                                                )#,
                                                                # tabPanel(title = "p-value Histogram",
                                                                #          plotOutput("p_hist", height = 600)
                                                                # )
                                                         ) # Tab box close
                                                       ),
                                                       column(
                                                         width=6,
                                                         tabBox(title = "Enrichment", width = 12,
                                                                tabPanel(title="Gene Ontology",
                                                                         box(uiOutput("contrast"), width = 6),
                                                                         box(
                                                                           selectInput("go_database", "GO database:",
                                                                                       c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                                         "Cellular Component"="GO_Cellular_Component_2017b",
                                                                                         "Biological Process"="GO_Biological_Process_2017b")),
                                                                           width= 5),
                                                                         actionButton("go_analysis", "Run Enrichment"),
                                                                         plotOutput("go_enrichment"),
                                                                         downloadButton('downloadGO', 'Download Table')
                                                                         
                                                                ),
                                                                tabPanel(title= "Pathway enrichment",
                                                                         box(uiOutput("contrast_1"), width = 6),
                                                                         box(
                                                                           selectInput("pathway_database", "Pathway database:",
                                                                                       c("KEGG"="KEGG_2016",
                                                                                         "Reactome"="Reactome_2016")),
                                                                           width= 5),
                                                                         actionButton("pathway_analysis", "Run Enrichment"),
                                                                         plotOutput("pathway_enrichment"),
                                                                         downloadButton('downloadPA', 'Download Table')
                                                                )
                                                                
                                                         ) # Tab box close
                                                       )
                                                   )) # fluidrow qc close
                                                 
                                                 
                                        ),
                                        tabPanel("ProteinGroup",
                                                 fluidRow(
                                                   column(4,
                                                          column(7,uiOutput("downloadTable_pr"),offset = 1),
                                                          column(4,uiOutput("downloadButton_pr"))
                                                   ),
                                                   column(5,
                                                          infoBoxOutput("significantBox_pr",width = 12)
                                                   ),
                                                   column(3,
                                                          uiOutput("downloadreport_pr")
                                                   )
                                                   # box(
                                                   #   column(6,uiOutput("downloadTable_pr"),offset = 1),
                                                   #   column(4,uiOutput("downloadButton_pr")), # make the button on same line
                                                   #   width = 4),
                                                   # infoBoxOutput("significantBox_pr",width = 4),
                                                   # box(
                                                   #   column(5,uiOutput("downloadreport_pr")), # offset for dist between buttons
                                                   #   #tags$br(),
                                                   #   #column(5,uiOutput('downloadPlots')),
                                                   #   width = 4)
                                                 ),
                                                 # align save button
                                                 tags$style(type='text/css', "#downloadTable_pr { width:100%; margin-top: 15px;}"),
                                                 tags$style(type='text/css', "#downloadButton_pr { width:100%; margin-top: 40px;}"), 
                                                 tags$style(type='text/css', "#downloadreport_pr { width:100%; vertical-align- middle; margin-top: 40px;
                                                                    margin-bottom: 25px;}"),
                                                 tags$br(), # Blank lines
                                                 tags$br(),
                                                 ## Data table and result plots box
                                                 fluidRow(
                                                   box(
                                                     title = "Results Table",
                                                     DT::dataTableOutput("contents_pr"),
                                                     #  actionButton("clear", "Deselect Rows"),
                                                     actionButton("original_pr", "Refresh Table"),
                                                     width = 6,height = 800,
                                                     status = "success",
                                                     #color=""
                                                     solidHeader = TRUE
                                                   ),
                                                   # column(
                                                   box(
                                                     width= 6,
                                                     collapsible = TRUE,
                                                     #status="primary",
                                                     #solidHeader=TRUE,
                                                     tabBox(
                                                       title = "Result Plots",
                                                       width = 12,
                                                       tabPanel(title = "Volcano plot",
                                                                fluidRow(
                                                                  box(width = 12,
                                                                      column(8,uiOutput("volcano_cntrst_pr")),
                                                                      
                                                                      column(4,
                                                                             prettyCheckbox("check_anova_pr",
                                                                                            "Apply ANOVA",
                                                                                            value = FALSE),
                                                                             prettyCheckbox("check_names_pr",
                                                                                            "Display names",
                                                                                            value = FALSE),
                                                                             prettyCheckbox("p_adj_pr",
                                                                                            "Adjusted p values",
                                                                                            value = FALSE)
                                                                      )
                                                                  ),
                                                                  tags$p("Select phosphosite from Results Table to highlight on the plot OR 
                                                  drag the mouse on plot to show expression of phosphosite in Table, ANOVA function only worked 
                                                                         for more than two groups.")
                                                                  #Add text line
                                                                  # tags$p("OR"),
                                                                  #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
                                                                ),
                                                                
                                                                fluidRow(
                                                                  plotOutput("volcano_pr", height = 600,
                                                                             # hover = "protein_hover"),
                                                                             #),
                                                                             # click = "protein_click"),
                                                                             brush = "protein_brush_pr",
                                                                             click = "protein_click_pr"),
                                                                  downloadButton('downloadVolcano_pr', 'Save Highlighted Plot'),
                                                                  actionButton("resetPlot_pr", "Clear Selection")
                                                                  #)),
                                                                )),
                                                       tabPanel(title= "Heatmap",
                                                                fluidRow(
                                                                  plotOutput("heatmap_pr", height = 600)
                                                                ),
                                                                fluidRow(
                                                                  box(numericInput("cluster_number_pr",
                                                                                   "Cluster to download",
                                                                                   min=1, max=6, value = 1), width = 6),
                                                                  box(downloadButton('downloadCluster_pr',"Save Cluster"),
                                                                      downloadButton('download_hm_svg_pr', "Save svg"),
                                                                      width = 5)
                                                                )
                                                       ),
                                                       tabPanel(title = "Individual Plot",
                                                                fluidRow(
                                                                  box(prettyRadioButtons("type_pr",
                                                                                         "Plot type",
                                                                                         choices = c("Box Plot"= "boxplot",
                                                                                                     "Violin Plot"="violin", 
                                                                                                     "Interaction Plot"= "interaction",
                                                                                                     "Intensity Plot"="dot"
                                                                                         ),
                                                                                         selected = "boxplot", 
                                                                                         inline = TRUE),
                                                                      width = 12
                                                                  ),
                                                                  tags$p("Select one or more rows from Results Table to plot individual 
                                                  protein intesities across conditions and replicates")
                                                                ),
                                                                fluidRow(
                                                                  plotOutput("protein_plot_pr"),
                                                                  downloadButton('downloadProtein_pr', 'Download Plot')
                                                                )
                                                       )
                                                       # verbatimTextOutput("protein_info"))
                                                     )
                                                   ) # box or column end
                                                   
                                                   
                                                   
                                                 ), # result plot colsed
                                                 ## QC Box
                                                 fluidRow(
                                                   div(id="qc_tab",
                                                       column(
                                                         width=6,
                                                         tabBox(title = "QC Plots", width = 12,
                                                                tabPanel(title = "PCA Plot",
                                                                         plotOutput("pca_plot_pr", height=600),
                                                                         downloadButton('download_pca_svg_pr', "Save svg")
                                                                ),
                                                                tabPanel(title="Sample Correlation",
                                                                         plotOutput("sample_corr_pr", height = 600),
                                                                         downloadButton('download_corr_svg_pr', "Save svg")
                                                                ),
                                                                tabPanel(title= "Sample CVs",
                                                                         plotOutput("sample_cvs_pr", height = 600),
                                                                         downloadButton('download_cvs_svg_pr', "Save svg")
                                                                ),
                                                                tabPanel(title = "Protein Numbers",
                                                                         plotOutput("numbers_pr", height = 600),
                                                                         downloadButton('download_num_svg_pr', "Save svg")
                                                                ),
                                                                
                                                                tabPanel(title = "Sample coverage",
                                                                         plotOutput("coverage_pr", height = 600),
                                                                         downloadButton('download_cov_svg_pr', "Save svg")
                                                                ),
                                                                tabPanel(title = "Normalization",
                                                                         plotOutput("norm_pr", height = 600),
                                                                         downloadButton('download_norm_svg_pr', "Save svg")
                                                                ),
                                                                # tabPanel(title = "Missing values - Quant",
                                                                #          plotOutput("detect", height = 600)
                                                                # ),
                                                                tabPanel(title = "Missing values - Heatmap",
                                                                         plotOutput("missval_pr", height = 600),
                                                                         downloadButton('download_missval_svg_pr', "Save svg")
                                                                ),
                                                                tabPanel(title = "Imputation",
                                                                         plotOutput("imputation_pr", height = 600),
                                                                         downloadButton('download_imp_svg_pr', "Save svg")
                                                                )#,
                                                                # tabPanel(title = "p-value Histogram",
                                                                #          plotOutput("p_hist", height = 600)
                                                                # )
                                                         ) # Tab box close
                                                       ),
                                                       column(
                                                         width=6,
                                                         tabBox(title = "Enrichment", width = 12,
                                                                tabPanel(title="Gene Ontology",
                                                                         box(uiOutput("contrast_pr"), width = 6),
                                                                         box(
                                                                           selectInput("go_database_pr", "GO database:",
                                                                                       c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                                         "Cellular Component"="GO_Cellular_Component_2017b",
                                                                                         "Biological Process"="GO_Biological_Process_2017b")),
                                                                           width= 5),
                                                                         actionButton("go_analysis_pr", "Run Enrichment"),
                                                                         plotOutput("go_enrichment_pr"),
                                                                         downloadButton('downloadGO_pr', 'Download Table')
                                                                         
                                                                ),
                                                                tabPanel(title= "Pathway enrichment",
                                                                         box(uiOutput("contrast_1_pr"), width = 5),
                                                                         box(
                                                                           selectInput("pathway_database_pr", "Pathway database:",
                                                                                       c("KEGG"="KEGG_2016",
                                                                                         "Reactome"="Reactome_2016")),
                                                                           width= 5),
                                                                         actionButton("pathway_analysis_pr", "Run Enrichment"),
                                                                         plotOutput("pathway_enrichment_pr"),
                                                                         downloadButton('downloadPA_pr', 'Download Table')
                                                                )
                                                                
                                                         ) # Tab box close
                                                       )
                                                   )) # fluidrow qc close
                                        ),
                                        
                                        tabPanel("Comparison",
                                                 fluidRow(width = 12,
                                                          box(width = 12,
                                                              title = 'Comparison of two levels log fold change',
                                                              status = "primary",
                                                              solidHeader = TRUE,
                                                              column(12,
                                                                column(6,uiOutput("volcano_comp")),
                                                                column(6,
                                                                       br(),
                                                                       br(),
                                                                       prettyCheckbox("check_anova_comp",
                                                                                        "Apply ANOVA",
                                                                                        value = FALSE)),
                                                              ),
                                                              column(12, plotOutput("scatter_plot", height=600))
                                                              # box(width = 6,'Unnormalized Phospho volcano',
                                                              #     plotOutput('volcano_phospho')),
                                                              # box(width = 6,'Normalized Phospho volcano by protein',
                                                              #     plotOutput('volcano_phospho_2'))
                                                              
                                                          )
                                                 ), # fluid colsed
                                                 fluidRow(width = 12,
                                                          box(width = 12,
                                                              title = 'Whole experiment-level comparison',
                                                              status = "primary",
                                                              solidHeader = TRUE,
                                                              tabBox(title = "QC Plots", width = 12,
                                                                     tabPanel(title = "PCA Plot",
                                                                              plotOutput("pca_plot_c", height=600)
                                                                     ),
                                                                     # tabPanel(title = "Scatter plot",
                                                                     #          plotOutput("scatter_plot", height=600)
                                                                     # ),
                                                                     tabPanel(title="Sample Correlation",
                                                                              # plotOutput("sample_corr_c", height = 600)
                                                                              fluidRow(height=600,
                                                                                       column(6,'Phosphosite',plotOutput("sample_corr_c1", height = 600)),
                                                                                       column(6,'Protein',plotOutput("sample_corr_c2", height = 600))
                                                                              )
                                                                              
                                                                     ),
                                                                     tabPanel(title= "Sample CVs",
                                                                              plotOutput("sample_cvs_c", height = 600)
                                                                     ),
                                                                     tabPanel(title = "Numbers",
                                                                              plotOutput("numbers_c", height = 600)
                                                                     ),
                                                                     
                                                                     tabPanel(title = "Sample coverage",
                                                                              plotOutput("coverage_c", height = 600)
                                                                     ),
                                                                     tabPanel(title = "Normalization",
                                                                              plotOutput("norm_c", height = 600)
                                                                     ),
                                                                     tabPanel(title = "Missing values - Heatmap",
                                                                              # plotOutput("missval_c", height = 600)
                                                                              fluidRow(height=600,
                                                                                       column(6,'Phosphosite',plotOutput("missval_c1", height = 600)),
                                                                                       column(6,'Protein',plotOutput("missval_c2", height = 600))
                                                                              ) 
                                                                     ),
                                                                     tabPanel(title = "Imputation",
                                                                              plotOutput("imputation_c", height = 600)
                                                                     )
                                                              ) # Tab box close
                                                          )
                                                 ),
                                                 fluidRow(width = 12,
                                                          box(width = 12,
                                                              status = "primary",
                                                              solidHeader = TRUE,
                                                              title = 'Individual Protein plots',
                                                              uiOutput('selected_gene'),
                                                              fluidRow(height = 600,
                                                                       column(width = 8,'Interaction plot',
                                                                              plotOutput("combined_inter", height = 600)),
                                                                       column(width = 4,'Points plot',
                                                                              plotOutput("combined_point", height = 600)
                                                                       )
                                                              )
                                                              
                                                              
                                                          )
                                                 )
                                        ),
                                        tabPanel("Normalization",
                                                 fluidRow(
                                                   column(4,
                                                          column(7,uiOutput("downloadTable_nr"),offset = 1),
                                                          column(4,uiOutput("downloadButton_nr"))
                                                   ),
                                                   column(5,
                                                          infoBoxOutput("significantBox_nr",width = 12)
                                                   ),
                                                   column(3,
                                                          uiOutput("downloadreport_nr")
                                                   )
                                                   # box(
                                                   #   column(6,uiOutput("downloadTable_nr"),offset = 1),
                                                   #   column(4,uiOutput("downloadButton_nr")), # make the button on same line
                                                   #   width = 4),
                                                   # infoBoxOutput("significantBox_nr",width = 4),
                                                   # box(
                                                   #   column(5,uiOutput("downloadreport_nr")), # offset for dist between buttons
                                                   #   #tags$br(),
                                                   #   #column(5,uiOutput('downloadPlots')),
                                                   #   width = 4)
                                                 ),
                                                 # align save button
                                                 tags$style(type='text/css', "#downloadTable_nr { width:100%; margin-top: 15px;}"),
                                                 tags$style(type='text/css', "#downloadButton_nr { width:100%; margin-top: 40px;}"), 
                                                 tags$style(type='text/css', "#downloadreport_nr { width:100%; vertical-align- middle; margin-top: 40px;
                                                                    margin-bottom: 25px;}"),
                                                 tags$br(), # Blank lines
                                                 tags$br(),
                                                 ## Data table and result plots box
                                                 
                                                 
                                                 fluidRow(
                                                   box(
                                                     title = "Results Table",
                                                     DT::dataTableOutput("contents_nr"),
                                                     #  actionButton("clear", "Deselect Rows"),
                                                     actionButton("original_nr", "Refresh Table"),
                                                     width = 6,height = 800,
                                                     status = "success",
                                                     #color=""
                                                     solidHeader = TRUE
                                                   ),
                                                   # column(
                                                   box(
                                                     width= 6,
                                                     collapsible = TRUE,
                                                     #status="primary",
                                                     #solidHeader=TRUE,
                                                     tabBox(
                                                       title = "Result Plots",
                                                       width = 12,
                                                       tabPanel(title = "Volcano plot",
                                                                fluidRow(
                                                                  box(width = 12,
                                                                      column(8,uiOutput("volcano_cntrst_nr")),
                                                                      
                                                                      column(4,
                                                                             prettyCheckbox("check_anova_nr",
                                                                                            "Apply ANOVA",
                                                                                            value = FALSE),
                                                                             prettyCheckbox("check_names_nr",
                                                                                            "Display names",
                                                                                            value = FALSE),
                                                                             prettyCheckbox("p_adj_nr",
                                                                                            "Adjusted p values",
                                                                                            value = FALSE)
                                                                      )
                                                                  ),
                                                                  tags$p("Select phosphosite from Results Table to highlight on the plot OR 
                                                  drag the mouse on plot to show expression of phosphosite in Table, ANOVA function only worked 
                                                                         for more than two groups.")
                                                                  #Add text line
                                                                  # tags$p("OR"),
                                                                  #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
                                                                ),
                                                                
                                                                fluidRow(
                                                                  plotOutput("volcano_nr", height = 600,
                                                                             # hover = "protein_hover"),
                                                                             #),
                                                                             # click = "protein_click"),
                                                                             brush = "protein_brush_nr",
                                                                             click = "protein_click_nr"),
                                                                  downloadButton('downloadVolcano_nr', 'Save Highlighted Plot'),
                                                                  actionButton("resetPlot_nr", "Clear Selection")
                                                                  #)),
                                                                )),
                                                       tabPanel(title= "Heatmap",
                                                                fluidRow(
                                                                  plotOutput("heatmap_nr", height = 600)
                                                                ),
                                                                fluidRow(
                                                                  box(numericInput("cluster_number_nr",
                                                                                   "Cluster to download",
                                                                                   min=1, max=6, value = 1), width = 6),
                                                                  box(downloadButton('downloadCluster_nr',"Save Cluster"),
                                                                      downloadButton('download_hm_svg_nr', "Save svg"),
                                                                      width = 5)
                                                                )
                                                       ),
                                                       tabPanel(title = "Individual Plot",
                                                                fluidRow(
                                                                  box(prettyRadioButtons("type_nr",
                                                                                         "Plot type",
                                                                                         choices = c("Box Plot"= "boxplot",
                                                                                                     "Violin Plot"="violin", 
                                                                                                     "Interaction Plot"= "interaction",
                                                                                                     "Intensity Plot"="dot"
                                                                                         ),
                                                                                         selected = "boxplot", 
                                                                                         inline = TRUE),
                                                                      width = 12
                                                                  ),
                                                                  tags$p("Select one or more rows from Results Table to plot individual 
                                                  protein intesities across conditions and replicates")
                                                                ),
                                                                fluidRow(
                                                                  plotOutput("protein_plot_nr"),
                                                                  downloadButton('downloadProtein_nr', 'Download Plot')
                                                                )
                                                       )
                                                       # verbatimTextOutput("protein_info"))
                                                     )
                                                   ) # box or column end
                                                   
                                                   
                                                   
                                                 ), # result plot colsed
                                                 
                                                 
                                                 ## QC Box
                                                 fluidRow(
                                                   div(id="qc_tab",
                                                       column(
                                                         width=6,
                                                         tabBox(title = "QC Plots", width = 12,
                                                                tabPanel(title = "PCA Plot",
                                                                         plotOutput("pca_plot_nr", height=600),
                                                                         downloadButton('download_pca_svg_nr', "Save svg")
                                                                ),
                                                                tabPanel(title="Sample Correlation",
                                                                         plotOutput("sample_corr_nr", height = 600),
                                                                         downloadButton('download_corr_svg_nr', "Save svg")
                                                                ),
                                                                tabPanel(title= "Sample CVs",
                                                                         plotOutput("sample_cvs_nr", height = 600),
                                                                         downloadButton('download_cvs_svg_nr', "Save svg")
                                                                ),
                                                                tabPanel(title = "Protein Numbers",
                                                                         plotOutput("numbers_nr", height = 600),
                                                                         downloadButton('download_num_svg_nr', "Save svg")
                                                                ),
                                                                
                                                                tabPanel(title = "Sample coverage",
                                                                         plotOutput("coverage_nr", height = 600),
                                                                         downloadButton('download_cov_svg_nr', "Save svg")
                                                                ),
                                                                tabPanel(title = "Normalization",
                                                                         plotOutput("norm_nr", height = 600),
                                                                         downloadButton('download_norm_svg_nr', "Save svg")
                                                                ),
                                                                # tabPanel(title = "Missing values - Quant",
                                                                #          plotOutput("detect", height = 600)
                                                                # ),
                                                                tabPanel(title = "Missing values - Heatmap",
                                                                         plotOutput("missval_nr", height = 600),
                                                                         downloadButton('download_missval_svg_nr', "Save svg")
                                                                ),
                                                                tabPanel(title = "Imputation",
                                                                         plotOutput("imputation_nr", height = 600),
                                                                         downloadButton('download_imp_svg_nr', "Save svg")
                                                                )#,
                                                                # tabPanel(title = "p-value Histogram",
                                                                #          plotOutput("p_hist", height = 600)
                                                                # )
                                                         ) # Tab box close
                                                       ),
                                                       column(
                                                         width=6,
                                                         tabBox(title = "Enrichment", width = 12,
                                                                tabPanel(title="Gene Ontology",
                                                                         box(uiOutput("contrast_nr"), width = 6),
                                                                         box(
                                                                           selectInput("go_database_nr", "GO database:",
                                                                                       c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                                         "Cellular Component"="GO_Cellular_Component_2017b",
                                                                                         "Biological Process"="GO_Biological_Process_2017b")),
                                                                           width= 5),
                                                                         actionButton("go_analysis_nr", "Run Enrichment"),
                                                                         plotOutput("go_enrichment_nr"),
                                                                         downloadButton('downloadGO_nr', 'Download Table')
                                                                         
                                                                ),
                                                                tabPanel(title= "Pathway enrichment",
                                                                         box(uiOutput("contrast_1_nr"), width = 6),
                                                                         box(
                                                                           selectInput("pathway_database_nr", "Pathway database:",
                                                                                       c("KEGG"="KEGG_2016",
                                                                                         "Reactome"="Reactome_2016")),
                                                                           width= 5),
                                                                         actionButton("pathway_analysis_nr", "Run Enrichment"),
                                                                         plotOutput("pathway_enrichment_nr"),
                                                                         downloadButton('downloadPA_nr', 'Download Table')
                                                                )
                                                                
                                                         ) # Tab box close
                                                       )
                                                   )) # fluidrow qc close
                                                 
                                                 
                                                 
                                                 )) # tabPanel list closed
                        ) # panelBox closed
                    ) 
                    
                  )),
          tabItem(tabName = "info",
                  fluidRow( 
                    box(
                      title = "User Guide",
                      h3("Phospho-Analyst: Manual"),
                      #                    div(p(HTML(paste0('A detail online user manual can be accessed ',
                      # 			a(href = 'https://monashbioinformaticsplatform.github.io/Phospho-Analyst/', 
                      #                                       target='_blank', 'here'))))),
                      div(p(HTML(paste0("A detailed user manual can be accessed in the near future."
                                        # a(href = './Phospho-Analyst_manual.pdf', 
                                          # target='_blank', tags$b("here."))
                                        )))),  
                      h4("Contact Us"),
                      p("For any feedback or question regarding Phospho-Analyst, please contact the 
			  Monash Proteomics and Metabolomics Facility:"),
                      tags$ul(
                        tags$li("Anup Shah: anup.shah(at)monash.edu"),
                        tags$li("Ralf Schittenhelm: ralf.schittenhelm(at)monash.edu")
                      ),
                      
#                       h4("How to Cite Phospho-Analyst?"),
#                       
#                       div(p(HTML(paste0("Please Cite: Shah AD, Goode RJA, Huang C, Powell DR, Schittenhelm RB. 
# 		Phospho-Analyst: An easy-to-use interactive web-platform to analyze and 
# 		visualize proteomics data preprocessed with MaxQuant. DOI:",
#                                         a(href = 'https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00496', 
#                                           target='_blank', tags$b("0.1021/acs.jproteome.9b00496")))))),
                      
                      
                      h4("News and Updates"),
                      
                      tags$ul(
                        
                      ),   
                      width = 12,
                      solidHeader = TRUE,
                      status = "primary"
                    ) #includeMarkdown("www/Info.md")
                  )
          )# info tab close
        ),
        tags$footer(
          tags$p("Supported by: Monash Proteomics and Metabolomics Facility & Monash Bioinformatics Platform, 
         Monash University"),
          align = "right"))  # Dasbboardbody close 
    )  #Dashboard page close
  )#Shiny UI Close
}