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
                                             p('Upload MaxQuant Phospho (STY)Sites.txt', style = 'color:#2E3440'),
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv')),
                                   tags$hr(),
                                   fileInput('file3',
                                             p('Phosphosite Experimental Design', 'Either upload a text file:',style = 'color:#2E3440'),
                                             # p('Either upload a text file', style = 'color:#2E3440'),
                                             accept=c('text/csv',
                                                      'text/comma-separated-values,text/plain',
                                                      '.csv')),
                                   # editable table for phosphosite data
                                   tags$strong('Or modify the following template:', style = 'margin-left: 15px;color:#2E3440'),
                                   actionButton("showTable", "Template", icon = icon("table")),
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
                                            
                                            prettyRadioButtons("normalisation",
                                                               p("Normalisation type", style = 'color:#2E3440'), 
                                                               choices = c("vsn", "median", "median subraction" = "median_sub"),
                                                               selected = "vsn"),
                                            
                                            prettyRadioButtons("imputation",
                                                               p("Imputation type", style = 'color:#2E3440'),
                                                               choices = c("Perseus-type"="man", MsCoreUtils::imputeMethods())[1:9],
                                                               selected = "man"),
                                            
                                            prettyRadioButtons("fdr_correction",
                                                               p("Type of FDR correction", style = 'color:#2E3440'),
                                                               choices =  c("Benjamini Hochberg"="BH",
                                                                            "t-statistics-based"="fdrtool"
                                                               ), selected= "BH"),
                                            numericInput("k_number",
                                                         p("Number of clusters in heatmap", style = 'color:#2E3440'),
                                                         min = 1, max = 20, value = 6)
                                   ),
                                   tags$hr(style = 'border-top: 4px double #2E3440'),
                                   menuItem(strong("Total Proteome (Optional)",style = 'color:#2E3440'),tabName="proteinGroup",icon = icon("file-upload"),
                                            fileInput('file2',
                                                      p('Upload MaxQuant ProteinGroup.txt', style = 'color:#2E3440'),
                                                      accept=c('text/csv',
                                                               'text/comma-separated-values,text/plain',
                                                               '.csv')),
                                            tags$hr(),
                                            fileInput('file4',
                                                      p('Protein Experimental Design',  tags$br(),'Either upload a text file:', style = 'color:#2E3440'),
                                                      # p('Either upload a text file', style = 'color:#2E3440'),
                                                      accept=c('text/csv',
                                                               'text/comma-separated-values,text/plain',
                                                               '.csv')),
                                            # editable table for protein data
                                            tags$strong('Or modify the following template:', style = 'margin-left: 15px;color:#2E3440'),
                                            actionButton("showTable_pr", "Template", icon = icon("table")),
                                            tags$hr(),
                                            menuItem("Advanced Options",tabName="advanced_pr", icon = icon("cogs"), 
                                                     numericInput("p_pr", 
                                                                  p("Adjusted p-value cutoff", style = 'color:#2E3440'),
                                                                  min = 0.0001, max = 0.1, value = 0.05),
                                                     numericInput("lfc_pr",
                                                                  p("Log2 fold change cutoff", style = 'color:#2E3440'),
                                                                  min = 0, max = 10, value = 1),
                                                     checkboxInput("paired_pr",
                                                                   p("Paired test", style = 'color:#2E3440'), FALSE),
                                                     
                                                     prettyRadioButtons("imputation_pr",
                                                                        p("Imputation type", style = 'color:#2E3440'),
                                                                        choices = c("Perseus-type"="man", MsCoreUtils::imputeMethods())[1:9],
                                                                        selected = "man"),
                                                     
                                                     prettyRadioButtons("fdr_correction_pr",
                                                                        p("Type of FDR correction", style = 'color:#2E3440'),
                                                                        choices =  c("Benjamini Hochberg"="BH",
                                                                                     "t-statistics-based"="fdrtool"
                                                                        ), selected= "BH"),
                                                     checkboxInput("single_peptide_pr",
                                                                   p("Include single peptide identifications", style = 'color:#2E3440'), FALSE),
                                                     numericInput("k_number_pr",
                                                                  p("Number of clusters in heatmap", style = 'color:#2E3440'),
                                                                  min = 1, max = 20, value = 6)
                                            )),
                                   
                                   tags$hr(style = 'border-top: 4px double #2E3440'),
                                   actionButton("analyze", "Start Analysis"),
                                   br()
                                   # tags$hr(),
                                   # p(a("Example Phosphoproteomics data", target= "_blank",
                                   #     href="data/Phospho (STY)Sites_example.txt", 
                                   #     download="Phospho (STY)Sites_example.txt")),
                                   # p(a("Example ProteinGroup data", target= "_blank",
                                   #     href="data/proteinGroups_example.txt", 
                                   #     download="proteinGroups_example.txt")),
                                   # p(a("Example Experimental Design file", target= "_blank",
                                   #     href="data/experimental_design_example.txt", 
                                   #     download="experimental_design_example.txt"))
          ), tabName = 'analysis'),
          convertMenuItem(
            menuItem('Demo', icon=icon("eye"), tabName = "demo"), tabName = "demo"),
          convertMenuItem(
            menuItem('Example Datasets', icon=icon("file-alt"), tabName = "datasets"), tabName = "datasets"),
          convertMenuItem(
            menuItem('User Guide', icon=icon("question"), 
                     # href = "https://monashbioinformaticsplatform.github.io/LFQ-Analyst/",
                     tabName = "info"), tabName = 'info' )
        )
      ), # sidebar close
      
      ####=== Dashboard Body ===####
      dashboardBody(
        use_theme(mytheme),
        useShinyjs(), #imp to use shinyjs functions
        tags$head(includeScript("google_analytics.js"),
                  tags$style(type='text/css')),
        
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
                      HTML('<center><img src="./Phospho_Analyst.svg" width="600px"></center>'),
                      br(),
                      h4("Sidebar tabs"),
                      tags$ul(
                        tags$li(tags$b("Analysis: "),"perform your own analysis"), 
                        tags$li(tags$b("Demo: "),"familiarise yourself with Phospho-Analyst by browsing through pre-analysed results"),
                        tags$li(tags$b("Datasets: "),"download example dataset files"),
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
                            tags$li("Upload your ", tags$b("Phospho (STY)Sites.txt "), "generated by MaxQuant."),
                            tags$li(tags$b("Phosphosite experimental design: "), "either upload a text file or modify the template.",
                                    bsModal("modal", "Interactive phosphosite experiment design table", "showTable", size = "large",
                                            rHandsontableOutput("exp_phospho"),
                                            downloadButton('download_exp', 'Download'),
                                            actionButton("save_exp","Save"),
                                            br(),
                                            tags$b(textOutput("save_message", container = span)),
                                            br(),
                                            br(),
                                            p("(Hint: Drag the mouse to copy the same input)", style = 'color: red;')
                                    )
                            ), 
                            
                            tags$li("Upload your ", tags$b(" Protein Group.txt (Optional)"),"generated by MaxQuant. "), 
                            tags$li(tags$b("Protein experimental design (Optional): "), "either upload a text file or modify the template.",
                                    bsModal("modal_pr", "Interactive protein experiment design table", "showTable_pr", size = "large",
                                            rHandsontableOutput("exp_protein"),
                                            downloadButton('download_exp_pr', 'Download'),
                                            actionButton("save_exp_pr","Save"),
                                            br(),
                                            tags$b(textOutput("save_message_pr", container = span)),
                                            br(),
                                            br(),
                                            p("(Hint: Drag the mouse to copy the same input)", style = 'color: red;')
                                    )
                            ), 
                            
                            # "  (",tags$b("Phospho (STY)Sites.txt "),'and',tags$b(" Protein Group.txt "),'could be optional to upload)',
                            tags$li(tags$b("Optional: "),"Adjust the p-value cut-off, the log2 fold change cut-off, 
                                                                         the imputation type, FDR correction method and/or number of clusters in heatmap in the",
                                    tags$b("Advanced Options.")),
                            tags$li("Press ", tags$b('"Start Analysis".')), 
                            tags$li(tags$b("Hint: "), " Use the ", tags$b("User Guide ")," tab for a detailed explanation of inputs,
                                                                         advanced options and outputs."), 
                            tags$li(tags$b("Note: "), " The experimental design file is not the" , 
                                    tags$b("'mqpar.xml' "),"file 
                                                                 from MaxQuant. ", " Use the ", tags$b("Datasets ")," tab to download example files provided.")
                          ),
                          br(),
                          HTML('<center><img src="./Phospho_Analyst.svg" width="500px"></center>'),
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
                                                  drag the mouse on plot to show expression of phosphosites in Table, ANOVA function only worked 
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
                                                                tabPanel(title="Gene Ontology/ Pathway",
                                                                         box(uiOutput("contrast"), width = 6),
                                                                         box(
                                                                           selectInput("go_database", "Database:",
                                                                                       c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                                         "Cellular Component"="GO_Cellular_Component_2017b",
                                                                                         "Biological Process"="GO_Biological_Process_2017b",
                                                                                         "KEGG"="KEGG_2016",
                                                                                         "Reactome"="Reactome_2016")),
                                                                           width= 5),
                                                                         actionButton("go_analysis", "Run Enrichment"),
                                                                         plotOutput("go_enrichment"),
                                                                         downloadButton('downloadGO', 'Download Table')
                                                                         
                                                                ),
                                                                tabPanel(title= "Kinase-Substrate enrichment ",
                                                                         box(width = 12,
                                                                             column(12, 
                                                                                    uiOutput("contrast_1")),
                                                                             column(6,
                                                                                    numericInput("m.cutoff",
                                                                                                 p("substrate count cutoff", style = 'color:#2E3440'),
                                                                                                 min = 0, value = 5)),
                                                                             column(6,
                                                                                    numericInput("p.cutoff",
                                                                                                 p("p-value cutoff", style = 'color:#2E3440'),
                                                                                                 min = 0, max = 1, value = 0.05))
                                                                         ),
                                                                         actionButton("KSEA_analysis", "Run Enrichment"),
                                                                         plotOutput("KSEA_enrichment"),
                                                                         downloadButton('downloadKSEA', 'Download Table')
                                                                )
                                                                
                                                                # tabPanel(title= "Pathway enrichment",
                                                                #          box(uiOutput("contrast_1"), width = 6),
                                                                #          box(
                                                                #            selectInput("pathway_database", "Pathway database:",
                                                                #                        c("KEGG"="KEGG_2016",
                                                                #                          "Reactome"="Reactome_2016")),
                                                                #            width= 5),
                                                                #          actionButton("pathway_analysis", "Run Enrichment"),
                                                                #          plotOutput("pathway_enrichment"),
                                                                #          downloadButton('downloadPA', 'Download Table')
                                                                # )
                                                                # 
                                                         ) # Tab box close
                                                       )
                                                   )) # fluidrow qc close
                                                 
                                                 
                                        ),
                                        tabPanel("Protein",
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
                                                     width = 6,
                                                     # height = 800,
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
                                                                  tags$p("Select protein from Results Table to highlight on the plot OR 
                                                  drag the mouse on plot to show expression of proteins in Table, ANOVA function only worked 
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
                                                                            uiOutput("downloadreport_comp")
                                                                     ),
                                                              ),
                                                              column(12, plotOutput("scatter_plot", height=600)),
                                                              downloadButton('download_scatter_svg', "Save svg")  
                                                          )
                                                 ), # fluid colsed
                                                 tags$style(type='text/css', "#downloadreport_comp { width:100%; vertical-align- middle; margin-top: 5px;}"),
                                                 fluidRow(width = 12,
                                                          box(width = 12,
                                                              title = 'Whole experiment-level comparison',
                                                              status = "primary",
                                                              solidHeader = TRUE,
                                                              tabBox(title = "QC Plots", width = 12,
                                                                     tabPanel(title = "PCA Plot",
                                                                              plotOutput("pca_plot_c", height=600),
                                                                              downloadButton('download_pca_svg_c', "Save svg")
                                                                     ),
                                                                     # tabPanel(title = "Scatter plot",
                                                                     #          plotOutput("scatter_plot", height=600)
                                                                     # ),
                                                                     tabPanel(title="Sample Correlation",
                                                                              # plotOutput("sample_corr_c", height = 600)
                                                                              fluidRow(height=600,
                                                                                       column(6,'Phosphosite',
                                                                                              plotOutput("sample_corr_c1", height = 600),
                                                                                              downloadButton('download_corr_svg_c', "Save svg")),
                                                                                       column(6,'Protein',
                                                                                              plotOutput("sample_corr_c2", height = 600),
                                                                                              downloadButton('download_corr_svg_c_1', "Save svg"))
                                                                              )
                                                                              
                                                                     ),
                                                                     tabPanel(title= "Sample CVs",
                                                                              plotOutput("sample_cvs_c", height = 600),
                                                                              downloadButton('download_cvs_svg_c', "Save svg")
                                                                     ),
                                                                     tabPanel(title = "Numbers",
                                                                              plotOutput("numbers_c", height = 600),
                                                                              downloadButton('download_num_svg_c', "Save svg")
                                                                     ),
                                                                     
                                                                     tabPanel(title = "Sample coverage",
                                                                              plotOutput("coverage_c", height = 600),  
                                                                              downloadButton('download_cov_svg_c', "Save svg")
                                                                     ),
                                                                     tabPanel(title = "Normalization",
                                                                              plotOutput("norm_c", height = 600),
                                                                              downloadButton('download_norm_svg_c', "Save svg")
                                                                     ),
                                                                     tabPanel(title = "Missing values - Heatmap",
                                                                              # plotOutput("missval_c", height = 600)
                                                                              fluidRow(height=600,
                                                                                       column(6,'Phosphosite',
                                                                                              plotOutput("missval_c1", height = 600),
                                                                                              downloadButton('download_missval_svg_c', "Save svg")),
                                                                                       column(6,'Protein',
                                                                                              plotOutput("missval_c2", height = 600),
                                                                                              downloadButton('download_missval_svg_c_1', "Save svg"))
                                                                              ) 
                                                                     ),
                                                                     tabPanel(title = "Imputation",
                                                                              plotOutput("imputation_c", height = 600),
                                                                              downloadButton('download_imp_svg_c', "Save svg") 
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
                                                              tabBox(width = 12,
                                                                     tabPanel(title= "Interaction plot",
                                                                              plotOutput("combined_inter", height = 600),
                                                                              downloadButton('download_inter_svg_c', "Save svg")),
                                                                     tabPanel(title= "Bubble plot",
                                                                              plotOutput("combined_point", height = 600),
                                                                              downloadButton('download_point_svg_c', "Save svg"))
                                                              )
                                                          )
                                                 )
                                        ),
                                        tabPanel("Phosphosite(corrected)",
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
                                                     width = 6,
                                                     # height = 800,
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
                                                  drag the mouse on plot to show expression of phosphosites in Table, ANOVA function only worked 
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
                                                                tabPanel(title = "Normalization (normal vs corrected)",
                                                                         plotOutput("norm_nr", height = 600),
                                                                         downloadButton('download_norm_svg_nr', "Save svg")
                                                                ),
                                                                tabPanel(title = "Imputation (normal vs corrected)",
                                                                         plotOutput("imputation_nr", height = 600),
                                                                         downloadButton('download_imp_svg_nr', "Save svg")
                                                                )
                                                                # tabPanel(title = "Phosphosite Numbers",
                                                                #          plotOutput("numbers_nr", height = 600),
                                                                #          downloadButton('download_num_svg_nr', "Save svg")
                                                                # ),
                                                                
                                                                # tabPanel(title = "Sample coverage",
                                                                #          plotOutput("coverage_nr", height = 600),
                                                                #          downloadButton('download_cov_svg_nr', "Save svg")
                                                                # ),
                                                                
                                                                # tabPanel(title = "Missing values - Quant",
                                                                #          plotOutput("detect", height = 600)
                                                                # ),
                                                                # tabPanel(title = "Missing values - Heatmap",
                                                                #          plotOutput("missval_nr", height = 600),
                                                                #          downloadButton('download_missval_svg_nr', "Save svg")
                                                                # ),
                                                                
                                                                #,
                                                                # tabPanel(title = "p-value Histogram",
                                                                #          plotOutput("p_hist", height = 600)
                                                                # )
                                                         ) # Tab box close
                                                       ),
                                                       column(
                                                         width=6,
                                                         tabBox(title = "Enrichment", width = 12,
                                                                tabPanel(title="Gene Ontology/ Pathway",
                                                                         box(uiOutput("contrast_nr"), width = 6),
                                                                         box(
                                                                           selectInput("go_database_nr", "GO database:",
                                                                                       c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                                         "Cellular Component"="GO_Cellular_Component_2017b",
                                                                                         "Biological Process"="GO_Biological_Process_2017b",
                                                                                         "KEGG"="KEGG_2016",
                                                                                         "Reactome"="Reactome_2016")),
                                                                           width= 5),
                                                                         actionButton("go_analysis_nr", "Run Enrichment"),
                                                                         plotOutput("go_enrichment_nr"),
                                                                         downloadButton('downloadGO_nr', 'Download Table')
                                                                         
                                                                ),
                                                                tabPanel(title= "Kinase-Substrate enrichment ",
                                                                         box(width = 12,
                                                                             column(12, 
                                                                                    uiOutput("contrast_1_nr")),
                                                                             column(6,
                                                                                    numericInput("m.cutoff_nr",
                                                                                                 p("substrate count cutoff", style = 'color:#2E3440'),
                                                                                                 min = 0, value = 5)),
                                                                             column(6,
                                                                                    numericInput("p.cutoff_nr",
                                                                                                 p("p-value cutoff", style = 'color:#2E3440'),
                                                                                                 min = 0, max = 1, value = 0.05))
                                                                         ),
                                                                         actionButton("KSEA_analysis_nr", "Run Enrichment"),
                                                                         plotOutput("KSEA_enrichment_nr"),
                                                                         downloadButton('downloadKSEA_nr', 'Download Table')
                                                                )
                                                                # tabPanel(title= "Pathway enrichment",
                                                                #          box(uiOutput("contrast_1_nr"), width = 6),
                                                                #          box(
                                                                #            selectInput("pathway_database_nr", "Pathway database:",
                                                                #                        c("KEGG"="KEGG_2016",
                                                                #                          "Reactome"="Reactome_2016")),
                                                                #            width= 5),
                                                                #          actionButton("pathway_analysis_nr", "Run Enrichment"),
                                                                #          plotOutput("pathway_enrichment_nr"),
                                                                #          downloadButton('downloadPA_nr', 'Download Table')
                                                                # )
                                                                
                                                         ) # Tab box close
                                                       )
                                                   )) # fluidrow qc close
                                                 
                                                 
                                                 
                                        )) # tabPanel list closed
                        ) # panelBox closed
                    ) 
                    
                  )),
          
          #### demo ui ####
          tabItem(tabName = "demo",
                  tabsetPanel(id = 'panel_list_dm',
                              tabPanel("Phosphosite",
                                       fluidRow(
                                         column(4,
                                                column(7,uiOutput("downloadTable_dm"),offset = 1),
                                                column(4,uiOutput("downloadButton_dm"))
                                         ),
                                         column(5,
                                                infoBoxOutput("significantBox_dm",width = 12)
                                         ),
                                         column(3,
                                                uiOutput("downloadreport_dm")
                                         )
                                       ),
                                       # align save button
                                       tags$style(type='text/css', "#downloadTable_dm { width:100%; margin-top: 15px;}"),
                                       tags$style(type='text/css', "#downloadButton_dm { width:100%; margin-top: 40px;}"), 
                                       tags$style(type='text/css', "#downloadreport_dm { width:100%; vertical-align- middle; margin-top: 40px;
                                                                    margin-bottom: 25px;}"),
                                       tags$br(), # Blank lines
                                       tags$br(),
                                       ## Data table and result plots box
                                       
                                       
                                       fluidRow(
                                         box(
                                           title = "Results Table",
                                           DT::dataTableOutput("contents_dm"),
                                           #  actionButton("clear", "Deselect Rows"),
                                           actionButton("original_dm", "Refresh Table"),
                                           width = 6,
                                           # height = 800,
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
                                                            column(8,uiOutput("volcano_cntrst_dm")),
                                                            
                                                            column(4,
                                                                   prettyCheckbox("check_anova_dm",
                                                                                  "Apply ANOVA",
                                                                                  value = FALSE),
                                                                   prettyCheckbox("check_names_dm",
                                                                                  "Display names",
                                                                                  value = FALSE),
                                                                   prettyCheckbox("p_adj_dm",
                                                                                  "Adjusted p values",
                                                                                  value = FALSE)
                                                            )
                                                        ),
                                                        tags$p("Select phosphosite from Results Table to highlight on the plot OR 
                                                  drag the mouse on plot to show expression of phosphosites in Table, ANOVA function only worked 
                                                                         for more than two groups.")
                                                        #Add text line
                                                        # tags$p("OR"),
                                                        #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
                                                      ),
                                                      
                                                      fluidRow(
                                                        plotOutput("volcano_dm", height = 600,
                                                                   # hover = "protein_hover"),
                                                                   #),
                                                                   # click = "protein_click"),
                                                                   brush = "protein_brush_dm",
                                                                   click = "protein_click_dm"),
                                                        downloadButton('downloadVolcano_dm', 'Save Highlighted Plot'),
                                                        actionButton("resetPlot_dm", "Clear Selection")
                                                        #)),
                                                      )),
                                             tabPanel(title= "Heatmap",
                                                      fluidRow(
                                                        plotOutput("heatmap_dm", height = 600)
                                                      ),
                                                      fluidRow(
                                                        box(numericInput("cluster_number_dm",
                                                                         "Cluster to download",
                                                                         min=1, max=6, value = 1), width = 6),
                                                        box(downloadButton('downloadCluster_dm',"Save Cluster"),
                                                            downloadButton('download_hm_svg_dm', "Save svg"),
                                                            width = 5)
                                                      )
                                             ),
                                             tabPanel(title = "Individual Plot",
                                                      fluidRow(
                                                        box(prettyRadioButtons("type_dm",
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
                                                        plotOutput("protein_plot_dm"),
                                                        downloadButton('downloadProtein_dm', 'Download Plot')
                                                      )
                                             )
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
                                                               plotOutput("pca_plot_dm", height=600),
                                                               downloadButton('download_pca_svg_dm', "Save svg")
                                                      ),
                                                      tabPanel(title="Sample Correlation",
                                                               plotOutput("sample_corr_dm", height = 600),
                                                               downloadButton('download_corr_svg_dm', "Save svg")
                                                      ),
                                                      tabPanel(title= "Sample CVs",
                                                               plotOutput("sample_cvs_dm", height = 600),
                                                               downloadButton('download_cvs_svg_dm', "Save svg")
                                                      ),
                                                      tabPanel(title = "Phosphosite Numbers",
                                                               plotOutput("numbers_dm", height = 600),
                                                               downloadButton('download_num_svg_dm', "Save svg")
                                                      ),
                                                      
                                                      tabPanel(title = "Sample coverage",
                                                               plotOutput("coverage_dm", height = 600),
                                                               downloadButton('download_cov_svg_dm', "Save svg")
                                                      ),
                                                      tabPanel(title = "Normalization",
                                                               plotOutput("norm_dm", height = 600),
                                                               downloadButton('download_norm_svg_dm', "Save svg")
                                                      ),
                                                      # tabPanel(title = "Missing values - Quant",
                                                      #          plotOutput("detect", height = 600)
                                                      # ),
                                                      tabPanel(title = "Missing values - Heatmap",
                                                               plotOutput("missval_dm", height = 600),
                                                               downloadButton('download_missval_svg_dm', "Save svg")
                                                      ),
                                                      tabPanel(title = "Imputation",
                                                               plotOutput("imputation_dm", height = 600),
                                                               downloadButton('download_imp_svg_dm', "Save svg")
                                                      )#,
                                                      # tabPanel(title = "p-value Histogram",
                                                      #          plotOutput("p_hist", height = 600)
                                                      # )
                                               ) # Tab box close
                                             ),
                                             column(
                                               width=6,
                                               tabBox(title = "Enrichment", width = 12,
                                                      tabPanel(title="Gene Ontology/ Pathway",
                                                               box(uiOutput("contrast_dm"), width = 6),
                                                               box(
                                                                 selectInput("go_database_dm", "Database:",
                                                                             c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                               "Cellular Component"="GO_Cellular_Component_2017b",
                                                                               "Biological Process"="GO_Biological_Process_2017b",
                                                                               "KEGG"="KEGG_2016",
                                                                               "Reactome"="Reactome_2016")),
                                                                 width= 5),
                                                               actionButton("go_analysis_dm", "Run Enrichment"),
                                                               plotOutput("go_enrichment_dm"),
                                                               downloadButton('downloadGO_dm', 'Download Table')
                                                               
                                                      ),
                                                      tabPanel(title= "Kinase-Substrate enrichment ",
                                                               box(width = 12,
                                                                   column(12, 
                                                                          uiOutput("contrast_1_dm")),
                                                                   column(6,
                                                                          numericInput("m.cutoff_dm",
                                                                                       p("substrate count cutoff", style = 'color:#2E3440'),
                                                                                       min = 0, value = 5)),
                                                                   column(6,
                                                                          numericInput("p.cutoff_dm",
                                                                                       p("p-value cutoff", style = 'color:#2E3440'),
                                                                                       min = 0, max = 1, value = 0.05))
                                                               ),
                                                               actionButton("KSEA_analysis_dm", "Run Enrichment"),
                                                               plotOutput("KSEA_enrichment_dm"),
                                                               downloadButton('downloadKSEA_dm', 'Download Table')
                                                               
                                                               # tabPanel(title= "Pathway enrichment",
                                                               #          box(uiOutput("contrast_1_dm"), width = 6),
                                                               #          box(
                                                               #            selectInput("pathway_database_dm", "Pathway database:",
                                                               #                        c("KEGG"="KEGG_2016",
                                                               #                          "Reactome"="Reactome_2016")),
                                                               #            width= 5),
                                                               #          actionButton("pathway_analysis_dm", "Run Enrichment"),
                                                               #          plotOutput("pathway_enrichment_dm"),
                                                               #          downloadButton('downloadPA_dm', 'Download Table')
                                                      )
                                                      
                                               ) # Tab box close
                                             )
                                         )) # fluidrow qc close
                              ), # phosphosite demo panel close
                              
                              tabPanel("Protein",
                                       fluidRow(
                                         column(4,
                                                column(7,uiOutput("downloadTable_dm_pr"),offset = 1),
                                                column(4,uiOutput("downloadButton_dm_pr"))
                                         ),
                                         column(5,
                                                infoBoxOutput("significantBox_dm_pr",width = 12)
                                         ),
                                         column(3,
                                                uiOutput("downloadreport_dm_pr")
                                         )
                                       ),
                                       # align save button
                                       tags$style(type='text/css', "#downloadTable_dm_pr { width:100%; margin-top: 15px;}"),
                                       tags$style(type='text/css', "#downloadButton_dm_pr { width:100%; margin-top: 40px;}"), 
                                       tags$style(type='text/css', "#downloadreport_dm_pr { width:100%; vertical-align- middle; margin-top: 40px;
                                                                    margin-bottom: 25px;}"),
                                       tags$br(), # Blank lines
                                       tags$br(),
                                       ## Data table and result plots box
                                       
                                       
                                       fluidRow(
                                         box(
                                           title = "Results Table",
                                           DT::dataTableOutput("contents_dm_pr"),
                                           #  actionButton("clear", "Deselect Rows"),
                                           actionButton("original_dm_pr", "Refresh Table"),
                                           width = 6,
                                           # height = 800,
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
                                                            column(8,uiOutput("volcano_cntrst_dm_pr")),
                                                            
                                                            column(4,
                                                                   prettyCheckbox("check_anova_dm_pr",
                                                                                  "Apply ANOVA",
                                                                                  value = FALSE),
                                                                   prettyCheckbox("check_names_dm_pr",
                                                                                  "Display names",
                                                                                  value = FALSE),
                                                                   prettyCheckbox("p_adj_dm_pr",
                                                                                  "Adjusted p values",
                                                                                  value = FALSE)
                                                            )
                                                        ),
                                                        tags$p("Select protein from Results Table to highlight on the plot OR 
                                                  drag the mouse on plot to show expression of proteins in Table, ANOVA function only worked 
                                                                         for more than two groups.")
                                                        #Add text line
                                                        # tags$p("OR"),
                                                        #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
                                                      ),
                                                      
                                                      fluidRow(
                                                        plotOutput("volcano_dm_pr", height = 600,
                                                                   # hover = "protein_hover"),
                                                                   #),
                                                                   # click = "protein_click"),
                                                                   brush = "protein_brush_dm_pr",
                                                                   click = "protein_click_dm_pr"),
                                                        downloadButton('downloadVolcano_dm_pr', 'Save Highlighted Plot'),
                                                        actionButton("resetPlot_dm_pr", "Clear Selection")
                                                        #)),
                                                      )),
                                             tabPanel(title= "Heatmap",
                                                      fluidRow(
                                                        plotOutput("heatmap_dm_pr", height = 600)
                                                      ),
                                                      fluidRow(
                                                        box(numericInput("cluster_number_dm_pr",
                                                                         "Cluster to download",
                                                                         min=1, max=6, value = 1), width = 6),
                                                        box(downloadButton('downloadCluster_dm_pr',"Save Cluster"),
                                                            downloadButton('download_hm_svg_dm_pr', "Save svg"),
                                                            width = 5)
                                                      )
                                             ),
                                             tabPanel(title = "Individual Plot",
                                                      fluidRow(
                                                        box(prettyRadioButtons("type_dm_pr",
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
                                                        plotOutput("protein_plot_dm_pr"),
                                                        downloadButton('downloadProtein_dm_pr', 'Download Plot')
                                                      )
                                             )
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
                                                               plotOutput("pca_plot_dm_pr", height=600),
                                                               downloadButton('download_pca_svg_dm_pr', "Save svg")
                                                      ),
                                                      tabPanel(title="Sample Correlation",
                                                               plotOutput("sample_corr_dm_pr", height = 600),
                                                               downloadButton('download_corr_svg_dm_pr', "Save svg")
                                                      ),
                                                      tabPanel(title= "Sample CVs",
                                                               plotOutput("sample_cvs_dm_pr", height = 600),
                                                               downloadButton('download_cvs_svg_dm_pr', "Save svg")
                                                      ),
                                                      tabPanel(title = "Protein Numbers",
                                                               plotOutput("numbers_dm_pr", height = 600),
                                                               downloadButton('download_num_svg_dm_pr', "Save svg")
                                                      ),
                                                      
                                                      tabPanel(title = "Sample coverage",
                                                               plotOutput("coverage_dm_pr", height = 600),
                                                               downloadButton('download_cov_svg_dm_pr', "Save svg")
                                                      ),
                                                      tabPanel(title = "Normalization",
                                                               plotOutput("norm_dm_pr", height = 600),
                                                               downloadButton('download_norm_svg_dm_pr', "Save svg")
                                                      ),
                                                      # tabPanel(title = "Missing values - Quant",
                                                      #          plotOutput("detect", height = 600)
                                                      # ),
                                                      tabPanel(title = "Missing values - Heatmap",
                                                               plotOutput("missval_dm_pr", height = 600),
                                                               downloadButton('download_missval_svg_dm_pr', "Save svg")
                                                      ),
                                                      tabPanel(title = "Imputation",
                                                               plotOutput("imputation_dm_pr", height = 600),
                                                               downloadButton('download_imp_svg_dm_pr', "Save svg")
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
                                                               box(uiOutput("contrast_dm_pr"), width = 6),
                                                               box(
                                                                 selectInput("go_database_dm_pr", "GO database:",
                                                                             c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                               "Cellular Component"="GO_Cellular_Component_2017b",
                                                                               "Biological Process"="GO_Biological_Process_2017b")),
                                                                 width= 5),
                                                               actionButton("go_analysis_dm_pr", "Run Enrichment"),
                                                               plotOutput("go_enrichment_dm_pr"),
                                                               downloadButton('downloadGO_dm_pr', 'Download Table')
                                                               
                                                      ),
                                                      tabPanel(title= "Pathway enrichment",
                                                               box(uiOutput("contrast_1_dm_pr"), width = 6),
                                                               box(
                                                                 selectInput("pathway_database_dm_pr", "Pathway database:",
                                                                             c("KEGG"="KEGG_2016",
                                                                               "Reactome"="Reactome_2016")),
                                                                 width= 5),
                                                               actionButton("pathway_analysis_dm_pr", "Run Enrichment"),
                                                               plotOutput("pathway_enrichment_dm_pr"),
                                                               downloadButton('downloadPA_dm_pr', 'Download Table')
                                                      )
                                                      
                                               ) # Tab box close
                                             )
                                         )) # fluidrow qc close
                              ), # protein demo panel close
                              
                              tabPanel(
                                "Comparison",
                                fluidRow(width = 12,
                                         box(width = 12,
                                             title = 'Comparison of two levels log fold change',
                                             status = "primary",
                                             solidHeader = TRUE,
                                             column(12,
                                                    column(6,uiOutput("volcano_comp_dm")),
                                                    column(6,
                                                           br(),
                                                           uiOutput("downloadreport_comp_dm")
                                                           # prettyCheckbox("check_anova_comp_dm",
                                                           #                "Apply ANOVA",
                                                           #                value = FALSE)
                                                    ),
                                             ),
                                             column(12, plotOutput("scatter_plot_dm", height=600)),
                                             downloadButton('download_scatter_svg_dm', "Save svg")
                                         )
                                ), # fluid colsed
                                tags$style(type='text/css', "#downloadreport_comp_dm { width:100%; vertical-align- middle; margin-top: 5px;}"),
                                fluidRow(width = 12,
                                         box(width = 12,
                                             title = 'Whole experiment-level comparison',
                                             status = "primary",
                                             solidHeader = TRUE,
                                             tabBox(title = "QC Plots", width = 12,
                                                    tabPanel(title = "PCA Plot",
                                                             plotOutput("pca_plot_c_dm", height=600),
                                                             downloadButton('download_pca_svg_c_dm', "Save svg")
                                                    ),
                                                    tabPanel(title="Sample Correlation",
                                                             fluidRow(height=600,
                                                                      column(6,'Phosphosite',
                                                                             plotOutput("sample_corr_c1_dm", height = 600),
                                                                             downloadButton('download_corr_svg_c_dm', "Save svg")),
                                                                      column(6,'Protein',
                                                                             plotOutput("sample_corr_c2_dm", height = 600),
                                                                             downloadButton('download_corr_svg_c_dm_1', "Save svg"))
                                                             )
                                                             # downloadButton('download_corr_pdf_c_dm', "Save pdf")
                                                             
                                                    ),
                                                    tabPanel(title= "Sample CVs",
                                                             plotOutput("sample_cvs_c_dm", height = 600),
                                                             downloadButton('download_cvs_svg_c_dm', "Save svg")
                                                    ),
                                                    tabPanel(title = "Numbers",
                                                             plotOutput("numbers_c_dm", height = 600),
                                                             downloadButton('download_num_svg_c_dm', "Save svg")
                                                    ),
                                                    
                                                    tabPanel(title = "Sample coverage",
                                                             plotOutput("coverage_c_dm", height = 600),
                                                             downloadButton('download_cov_svg_c_dm', "Save svg")
                                                    ),
                                                    tabPanel(title = "Normalization",
                                                             plotOutput("norm_c_dm", height = 600),
                                                             downloadButton('download_norm_svg_c_dm', "Save svg")
                                                    ),
                                                    tabPanel(title = "Missing values - Heatmap",
                                                             fluidRow(height=600,
                                                                      column(6,'Phosphosite',
                                                                             plotOutput("missval_c1_dm", height = 600),
                                                                             downloadButton('download_missval_svg_c_dm', "Save svg")),
                                                                      column(6,'Protein',
                                                                             plotOutput("missval_c2_dm", height = 600),
                                                                             downloadButton('download_missval_svg_c_dm_1', "Save svg"))
                                                             )
                                                    ),
                                                    tabPanel(title = "Imputation",
                                                             plotOutput("imputation_c_dm", height = 600),
                                                             downloadButton('download_imp_svg_c_dm', "Save svg")   
                                                    )
                                             ) # Tab box close
                                         )
                                ),
                                fluidRow(width = 12,
                                         box(width = 12,
                                             status = "primary",
                                             solidHeader = TRUE,
                                             title = 'Individual Protein plots',
                                             uiOutput('selected_gene_dm'),
                                             tabBox(width = 12,
                                                    tabPanel(title= "Interaction plot",
                                                             plotOutput("combined_inter_dm", height = 600),
                                                             downloadButton('download_inter_svg_c_dm', "Save svg")),
                                                    tabPanel(title= "Bubble plot",
                                                             plotOutput("combined_point_dm", height = 600),
                                                             downloadButton('download_point_svg_c_dm', "Save svg"))
                                             )
                                         )
                                )
                                
                              ), # comparison demo page closed
                              
                              tabPanel(
                                "Phosphosite(corrected)",
                                fluidRow(
                                  column(4,
                                         column(7,uiOutput("downloadTable_dm_nr"),offset = 1),
                                         column(4,uiOutput("downloadButton_dm_nr"))
                                  ),
                                  column(5,
                                         infoBoxOutput("significantBox_dm_nr",width = 12)
                                  ),
                                  column(3,
                                         uiOutput("downloadreport_dm_nr")
                                  )
                                  # box(
                                  #   column(6,uiOutput("downloadTable_dm_nr"),offset = 1),
                                  #   column(4,uiOutput("downloadButton_dm_nr")), # make the button on same line
                                  #   width = 4),
                                  # infoBoxOutput("significantBox_dm_nr",width = 4),
                                  # box(
                                  #   column(5,uiOutput("downloadreport_dm_nr")), # offset for dist between buttons
                                  #   #tags$br(),
                                  #   #column(5,uiOutput('downloadPlots')),
                                  #   width = 4)
                                ),
                                # align save button
                                tags$style(type='text/css', "#downloadTable_dm_nr { width:100%; margin-top: 15px;}"),
                                tags$style(type='text/css', "#downloadButton_dm_nr { width:100%; margin-top: 40px;}"), 
                                tags$style(type='text/css', "#downloadreport_dm_nr { width:100%; vertical-align- middle; margin-top: 40px;
                                                                    margin-bottom: 25px;}"),
                                tags$br(), # Blank lines
                                tags$br(),
                                ## Data table and result plots box
                                
                                
                                fluidRow(
                                  box(
                                    title = "Results Table",
                                    DT::dataTableOutput("contents_dm_nr"),
                                    #  actionButton("clear", "Deselect Rows"),
                                    actionButton("original_dm_nr", "Refresh Table"),
                                    width = 6,
                                    # height = 800,
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
                                                     column(8,uiOutput("volcano_cntrst_dm_nr")),
                                                     
                                                     column(4,
                                                            prettyCheckbox("check_anova_dm_nr",
                                                                           "Apply ANOVA",
                                                                           value = FALSE),
                                                            prettyCheckbox("check_names_dm_nr",
                                                                           "Display names",
                                                                           value = FALSE),
                                                            prettyCheckbox("p_adj_dm_nr",
                                                                           "Adjusted p values",
                                                                           value = FALSE)
                                                     )
                                                 ),
                                                 tags$p("Select phosphosite from Results Table to highlight on the plot OR 
                                                  drag the mouse on plot to show expression of phosphosites in Table, ANOVA function only worked 
                                                                         for more than two groups.")
                                                 #Add text line
                                                 # tags$p("OR"),
                                                 #  tags$p("Drag the mouse on plot to show expression of proteins in Table") 
                                               ),
                                               
                                               fluidRow(
                                                 plotOutput("volcano_dm_nr", height = 600,
                                                            # hover = "protein_hover"),
                                                            #),
                                                            # click = "protein_click"),
                                                            brush = "protein_brush_dm_nr",
                                                            click = "protein_click_dm_nr"),
                                                 downloadButton('downloadVolcano_dm_nr', 'Save Highlighted Plot'),
                                                 actionButton("resetPlot_dm_nr", "Clear Selection")
                                                 #)),
                                               )),
                                      tabPanel(title= "Heatmap",
                                               fluidRow(
                                                 plotOutput("heatmap_dm_nr", height = 600)
                                               ),
                                               fluidRow(
                                                 box(numericInput("cluster_number_dm_nr",
                                                                  "Cluster to download",
                                                                  min=1, max=6, value = 1), width = 6),
                                                 box(downloadButton('downloadCluster_dm_nr',"Save Cluster"),
                                                     downloadButton('download_hm_svg_dm_nr', "Save svg"),
                                                     width = 5)
                                               )
                                      ),
                                      tabPanel(title = "Individual Plot",
                                               fluidRow(
                                                 box(prettyRadioButtons("type_dm_nr",
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
                                                 plotOutput("protein_plot_dm_nr"),
                                                 downloadButton('downloadProtein_dm_nr', 'Download Plot')
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
                                                        plotOutput("pca_plot_dm_nr", height=600),
                                                        downloadButton('download_pca_svg_dm_nr', "Save svg")
                                               ),
                                               tabPanel(title="Sample Correlation",
                                                        plotOutput("sample_corr_dm_nr", height = 600),
                                                        downloadButton('download_corr_svg_dm_nr', "Save svg")
                                               ),
                                               tabPanel(title= "Sample CVs",
                                                        plotOutput("sample_cvs_dm_nr", height = 600),
                                                        downloadButton('download_cvs_svg_dm_nr', "Save svg")
                                               ),
                                               
                                               tabPanel(title = "Normalization (normal vs corrected)",
                                                        plotOutput("norm_dm_nr", height = 600),
                                                        downloadButton('download_norm_svg_dm_nr', "Save svg")
                                               ),
                                               
                                               tabPanel(title = "Imputation (normal vs corrected)",
                                                        plotOutput("imputation_dm_nr", height = 600),
                                                        downloadButton('download_imp_svg_dm_nr', "Save svg")
                                               )
                                               # tabPanel(title = "Phosphosite Numbers",
                                               #          plotOutput("numbers_dm_nr", height = 600),
                                               #          downloadButton('download_num_svg_dm_nr', "Save svg")
                                               # ),
                                               # 
                                               # tabPanel(title = "Sample coverage",
                                               #          plotOutput("coverage_dm_nr", height = 600),
                                               #          downloadButton('download_cov_svg_dm_nr', "Save svg")
                                               # ),
                                               
                                               # tabPanel(title = "Missing values - Quant",
                                               #          plotOutput("detect", height = 600)
                                               # ),
                                               # tabPanel(title = "Missing values - Heatmap",
                                               #          plotOutput("missval_dm_nr", height = 600),
                                               #          downloadButton('download_missval_svg_dm_nr', "Save svg")
                                               # ),
                                               
                                               #,
                                               # tabPanel(title = "p-value Histogram",
                                               #          plotOutput("p_hist", height = 600)
                                               # )
                                        ) # Tab box close
                                      ),
                                      column(
                                        width=6,
                                        tabBox(title = "Enrichment", width = 12,
                                               tabPanel(title="Gene Ontology/ Pathway",
                                                        box(uiOutput("contrast_dm_nr"), width = 6),
                                                        box(
                                                          selectInput("go_database_dm_nr", "GO database:",
                                                                      c("Molecular Function"="GO_Molecular_Function_2017b",
                                                                        "Cellular Component"="GO_Cellular_Component_2017b",
                                                                        "Biological Process"="GO_Biological_Process_2017b",
                                                                        "KEGG"="KEGG_2016",
                                                                        "Reactome"="Reactome_2016")),
                                                          width= 5),
                                                        actionButton("go_analysis_dm_nr", "Run Enrichment"),
                                                        plotOutput("go_enrichment_dm_nr"),
                                                        downloadButton('downloadGO_dm_nr', 'Download Table')
                                                        
                                               ),
                                               tabPanel(title= "Kinase-Substrate enrichment ",
                                                        box(width = 12,
                                                            column(12, 
                                                                   uiOutput("contrast_1_dm_nr")),
                                                            column(6,
                                                                   numericInput("m.cutoff_dm_nr",
                                                                                p("substrate count cutoff", style = 'color:#2E3440'),
                                                                                min = 0, value = 5)),
                                                            column(6,
                                                                   numericInput("p.cutoff_dm_nr",
                                                                                p("p-value cutoff", style = 'color:#2E3440'),
                                                                                min = 0, max = 1, value = 0.05))
                                                        ),
                                                        actionButton("KSEA_analysis_dm_nr", "Run Enrichment"),
                                                        plotOutput("KSEA_enrichment_dm_nr"),
                                                        downloadButton('downloadKSEA_dm_nr', 'Download Table')
                                               )
                                               
                                        ) # Tab box close
                                      )
                                  )) # fluidrow qc close
                              ) # phosphosite(corrected) demo page closed
                              
                  ) # panel list close
                  
                  
                  
          ),
          
          #### example datasets ui ####
          tabItem(tabName = "datasets",
                  fluidRow( 
                    box(
                      title = "Example Datasets",
                      h3("Dowanloadable Example datasets:"),
                      #                    div(p(HTML(paste0('A detail online user manual can be accessed ',
                      # 			a(href = 'https://monashbioinformaticsplatform.github.io/Phospho-Analyst/', 
                      #                                       target='_blank', 'here'))))),
                      div(p(HTML(paste0("The datasets are the raw data files used in Demo Tab"
                                        # a(href = './Phospho-Analyst_manual.pdf', 
                                        # target='_blank', tags$b("here."))
                      )))),
                      tags$ul(
                        tags$li(a("Example Phosphoproteomics data", target= "_blank",
                                  href="data/Phospho (STY)Sites_example.txt", 
                                  download="Phospho (STY)Sites_example.txt")),
                        tags$li(a("Example ProteinGroup data", target= "_blank",
                                  href="data/proteinGroups_example.txt", 
                                  download="proteinGroups_example.txt")),
                        tags$li(a("Example Experimental Design file (Phosphosite)", target= "_blank",
                                  href="data/experimental_design_example.txt", 
                                  download="experimental_design_example.txt")),
                        tags$li(a("Example Experimental Design file (Proteome)", target= "_blank",
                                  href="data/experimental_design_protein_example.txt", 
                                  download="experimental_design_protein_example.txt"))
                      ),
                      width = 12,
                      solidHeader = TRUE,
                      status = "primary"
                    ) #includeMarkdown("www/Info.md")
                  )
          ),# example dataset tab close
          
          #### user guide ui ####
          tabItem(tabName = "info",
                  fluidRow( 
                    box(
                      title = "User Guide",
                      h3("Phospho-Analyst: Manual"),
                      #                    div(p(HTML(paste0('A detail online user manual can be accessed ',
                      # 			a(href = 'https://monashbioinformaticsplatform.github.io/Phospho-Analyst/', 
                      #                                       target='_blank', 'here'))))),
                      div(p(HTML(paste0('A detail online user manual can be accessed ',
                                        a(href = './Phospho-Analyst-manual.pdf',
                                        target='_blank', tags$b("here."))
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
                        tags$li(p(HTML(paste0("Related App: ", tags$b("LFQ-Analyst")," for analysing label-free quantitative proteomics dataset",
                                              a(href = 'https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst/',
                                                target='_blank', tags$b("https://bioinformatics.erc.monash.edu/apps/LFQ-Analyst"))))))
                        
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
