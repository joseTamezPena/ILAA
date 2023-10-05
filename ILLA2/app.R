#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Install and load required packages if not already installed
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("heatmaply", quietly = TRUE)) {
  install.packages("heatmaply")
}

if (!requireNamespace("FRESA.CAD", quietly = TRUE)) {
  devtools::install_github("joseTamezPena/FRESA.CAD")
}


# Load required packages
library(shiny)
library(ggplot2)
library(dplyr)
library(heatmaply)
library(FRESA.CAD)

# Define UI
ui <- fluidPage(
  titlePanel("ILAA Calculator @0.80"),
  
  sidebarLayout(
    sidebarPanel(
      verbatimTextOutput("Message"),      
      fileInput("file", "Choose a CSV file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      downloadButton("downloadCSV", "Download ILAA Results"),
      downloadButton("downloadUPLTM", "Download the UPLTM")
    ),
    
    mainPanel(
      plotlyOutput("heatmapin"),
      plotlyOutput("heatmapout")
    )
  )
)

# Define server
server <- function(input, output) {
  
  output$Message <- renderText({
    "ILAA will compute a linear transformation 
    of the data to handle multicollineariy 
    issues in data sets.
    
     Add your tabular .csv data.
    On return you will have a transformed dataset.
    and the corresponding linear transformation."
  })
  
  # Load dataset
  dataset <- reactive({
    req(input$file)
    inFile <- input$file
    
    # Check if file is a CSV
    if (grepl("\\.csv$", inFile$name)) {
      read.csv(inFile$datapath)
    } else {
      return(NULL)
    }
  })
  
  # Plot data
  output$heatmapin <- renderPlotly({
    if (!is.null(dataset())) {
      heatmaply(percentize(dataset()), plot_method = "plotly",main = "Input")
    }
  })
  
  # Perform ILAA
  ILAA_result <- reactive({
    if (!is.null(dataset())) {
      illa <- ILAA(dataset())
      return(illa)
    }
  })
  
  # Generate scatter plot
  output$heatmapout <- renderPlotly({
    if (!is.null(ILAA_result())) {
      heatmaply(percentize(ILAA_result()), plot_method = "plotly",main = "Output")
    }
  })
  
  # Download ILAA results as CSV
  output$downloadCSV <- downloadHandler(
    filename = function() {
      paste("ILAA_results", ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(ILAA_result())) {
        write.csv(as.data.frame(ILAA_result()), file, row.names = FALSE)
        write.csv(as.data.frame(attr(ILAA_result(),"UPLTM")), paste(file,"UPLTM.csv",sep="_"), row.names = TRUE)
      }
    }
  )
  # Download UPLTM results as CSV
  output$downloadUPLTM <- downloadHandler(
    filename = function() {
      paste("ILAA_UPLTM", ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(ILAA_result())) {
        write.csv(as.data.frame(attr(ILAA_result(),"UPLTM")), file, row.names = TRUE)
      }
    }
  )
}

# Run the application
shinyApp(ui, server)
