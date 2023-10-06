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

if (!requireNamespace("Rfast", quietly = TRUE)) {
  install.packages("Rfast")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

if (!requireNamespace("FRESA.CAD", quietly = TRUE)) {
  install.packages("FRESA.CAD")
}


# Load required packages
library(shiny)
library(igraph)
library(ggplot2)
library(dplyr)
library(heatmaply)
library(FRESA.CAD)

# Define UI
ui <- fluidPage(
  titlePanel("UPLTM Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      htmlOutput("Message"),
      br(),
      uiOutput("tab"),
      uiOutput("tab2"),
      uiOutput("tab3"),
      br(),
      sliderInput("selected_Thr", "Target Maximum Correlation:", value = 80, min = 0, max = 99,step = 5),
      br(),
      fileInput("file", "Choose a CSV file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      downloadButton("downloadCSV", "Download the Transformed Data set"),
      downloadButton("downloadUPLTM", "Download the Transformation"),
    ),
    
    mainPanel(
      column(6,
        plotlyOutput("heatmapin"),
        plotOutput("Network")
      ),
      column(6,
        plotlyOutput("heatmapout"),
        dataTableOutput("Formulas")
      )
    )
  )
)

# Define server
server <- function(input, output) {
  output$Message <- renderText({
"<p>This app performs a Linear Transformation (UPLTM) on the dataset, addressing issues related to multicollinearity.</p>

<p>Follow these steps:</p>
<ol>
  <li>Select the desired maximum correlation (Pearson).</li>
  <li>Upload your tabular .csv data. Columns should represent variables, and rows should represent independent observations.</li>
</ol>

<p>Upon completion, you will receive:</p>
<ol>
  <li>The transformed dataset.</li>
  <li>The corresponding linear transformation.</li>
</ol>

<p>In the main panel, you'll find visualizations of:</p> 
<ol>
  <li>the input Pearson correlation matrix
  <li>the transformed Pearson correlation matrix.
  <li>the Variable Association Network.
  <li>the latent variable formula with its corresponding explained variance(R2).
</ol>

<p>Please note:<br>
The transformation is applied to continuous or ordinal variables with more than 4 categories<br>
Character columns, factors and binary variables will not be affected by the transformation<br>.
.</p>
<p>Contact: jose.tamezpena@tec.mx</p>

"    
 
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
  
  
  # Perform ILAA
  ILAA_result <- reactive({
    if (!is.null(dataset())) {
      thrvalue <- input$selected_Thr/100
      illa <- ILAA(dataset(),thr=thrvalue)
      return(illa)
    }
  })
  
  # Plot of correlated data
  output$heatmapin <- renderPlotly({
    if (!is.null(ILAA_result())) {
      datavars <- rownames(attr(ILAA_result(),"UPLTM"))
      heatmaply(cor(dataset()[,datavars]), plot_method = "plotly",main = "Pearson Correlation of the Input Data")
    }
  })
  
  # Generate heat map of output correlation  
    output$heatmapout <- renderPlotly({
      if (!is.null(ILAA_result())) {
        datavars <- colnames(attr(ILAA_result(),"UPLTM"))
        heatmaply(cor(ILAA_result()[,datavars]), plot_method = "plotly",main = "Pearson Correlation of Transformed Data")
      }
    })

    # Report the formulas  
    output$Formulas <- renderDataTable({
      if (!is.null(ILAA_result())) {
        R2 <- 1.0 - attr(ILAA_result(),"VarRatio")
        dta <- attr(getLatentCoefficients(ILAA_result()),"LatentCharFormulas")
        dtb <- cbind(Variable=names(dta),R2=R2[names(dta)],formula=dta)
        data.frame(dtb)
      }
    })
    
    # Report the formulas Network 
    output$Network <- renderPlot({
      if (!is.null(ILAA_result())) {
        VertexSize <-attr(ILAA_result(),"VarRatio")
        VertexSize <- 10*sqrt((VertexSize-min(VertexSize))/(max(VertexSize)-min(VertexSize)))
        transform <- 1*(attr(ILAA_result(),"UPLTM") != 0)
        colnames(transform) <- str_remove_all(colnames(transform),"La_")  # For network analysis
        names(VertexSize) <- str_remove_all(names(VertexSize),"La_")  # For network analysis
        tnames <- colnames(transform)
        
        rsum <- apply(transform,1,sum) + 0.01*VertexSize[tnames]
        csum <- apply(transform,2,sum) + 0.01*VertexSize[tnames]

        ntop <- min(10,length(rsum))
        topfeatures <- unique(c(names(rsum[order(-rsum)])[1:ntop],names(csum[order(-csum)])[1:ntop]))
        rtrans <- transform[topfeatures,]
        csum <- (apply(rtrans !=0,2,sum) > 1)
        rtrans <- rtrans[,csum]
        topfeatures <- unique(c(topfeatures,colnames(rtrans)))

        transform <- transform[topfeatures,topfeatures]

        VertexSize <- VertexSize[colnames(transform)]
        
        gr <- graph_from_adjacency_matrix(transform,mode = "directed",diag = FALSE,weighted=TRUE)
        gr$layout <- layout_with_fr
        
        fc <- cluster_optimal(gr)
        plot(fc, gr,
             edge.width = 2*E(gr)$weight,
             vertex.size=VertexSize,
             edge.arrow.size=0.65,
             edge.arrow.width=0.75,
             vertex.label.cex=(0.5+0.1*VertexSize),
             vertex.label.dist=0.5 + 0.2*VertexSize,
             main="Top Feature Association")
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
  url <- a("ILAA Tutorial", href="https://rpubs.com/J_Tamez/ILAA_Tutorial")
  url2 <- a("Validation Repository", href="https://github.com/joseTamezPena/LatentBiomarkers")
  url3 <- a("FRESA.CAD Repository", href="https://github.com/joseTamezPena/FRESA.CAD")
  output$tab <- renderUI({
    tagList("Tutorial:", url)
  })
  output$tab2 <- renderUI({
    tagList("Validation:", url2)
  })
  output$tab3 <- renderUI({
    tagList("FRESA.CAD:", url3)
  })
  
  
}

# Run the application
shinyApp(ui, server)
