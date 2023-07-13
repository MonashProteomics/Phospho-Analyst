# module for saving plots
save_plot_left_ui <- function(id){
  ns <- NS(id)
  dropdownButton(
    inputId = ns("plot_dropdown"),
    label = "Save plot",
    circle = FALSE,
    width = "200px",
    up = TRUE,
    right = FALSE,
    numericInput(ns("plot_h"), "Height", value = 7),
    numericInput(ns("plot_w"), "Width", value = 7),
    radioButtons(ns("plot_format"), label = NULL, inline = T, 
                 choices = c("svg",  "png", "tiff")),
    downloadButton(ns('download_plot_svg'), "Download")
  )
}

# module for saving results plot
save_plot_right_ui <- function(id){
  ns <- NS(id)
  dropdownButton(
    inputId = ns("plot_dropdown"),
    label = "Save plot",
    circle = FALSE,
    width = "200px",
    up = TRUE,
    right = TRUE,
    numericInput(ns("plot_h"), "Height", value = 7),
    numericInput(ns("plot_w"), "Width", value = 7),
    radioButtons(ns("plot_format"), label = NULL, inline = T, 
                 choices = c("svg",  "png", "tiff")),
    downloadButton(ns('download_plot_svg'), "Download")
  )
}

save_plot_server <- function(id, plot_input,plot_option){
  moduleServer(
    id,
    function(input,output, session){
      ns <- session$ns
      
      plot_name <- reactive({
        name <- id %>% gsub("(_).*","",.)
        if (grepl("^volcano|^protein_plot|^abundance_comp",id)){
          plot_name <- paste0(name, "_plot_", plot_option(), ".", input$plot_format)
        } else {
          plot_name <- paste0(name, "_plot.", input$plot_format)
        }
        return(plot_name)
      })
      
      output$download_plot_svg<-downloadHandler(
        filename = function() {
          plot_name()
        },
        content = function(file) {
          if(input$plot_format=="svg"){
            svg(file, height = input$plot_h, width = input$plot_w)
          } else if (input$plot_format=="png"){
            png(file, height = input$plot_h, width = input$plot_w, units = "in", res = 300)
          } else {
            tiff(file, height = input$plot_h, width = input$plot_w, units = "in", res = 300)
          }
          print(plot_input())
          dev.off()
        }
      )
    }
  )
}
