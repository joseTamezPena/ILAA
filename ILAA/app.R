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

if (!requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS")
}

if (!requireNamespace("nlme", quietly = TRUE)) {
  install.packages("nlme")
}

if (!requireNamespace("mda", quietly = TRUE)) {
  install.packages("mda")
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

  tags$head(
    tags$meta(name = "description", content = "Estimate Linear Transformation Matrices"),
    tags$meta(name = "keywords", content = "Multicollinearity, Decorrelation, PCA, EFA, Whitening"),
    tags$meta(property = "og:title", content="Linear Decorrelation"),
    tags$meta(property = "og:description", content="Estimation of Linear Matrix to Address Multicolliearity")
  ),
    
  sidebarLayout(
    sidebarPanel(
      htmlOutput("Message"),
      br(),
      sliderInput("selected_Thr", "Target Maximum Correlation:", value = 50, min = 0, max = 99,step = 5),
      fileInput("file", "Choose a CSV file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      br(),
#      sliderInput("Bootstrap", "Number of Bootstraps", value = 0, min = 0, max = 120,step = 30),
      selectInput("corMeasure", "Select Correlation Measure", 
                  choices=c("pearson","spearman"),
                  selected="pearson"),
      br(),
      downloadButton("downloadCSV", "Download the Transformed Data set"),
      downloadButton("downloadUPLTM", "Download the Transformation"),
      br(),
      uiOutput("tab"),
      uiOutput("tab2"),
      uiOutput("tab3"),
      uiOutput("Twitter"),
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
"<p>This app estimates and performs a Linear Transformation (UPLTM) on the dataset, addressing issues related to multicollinearity.</p>
<p>Follow at least these two steps:</p>
<ol>
  <li>Select the target maximum correlation.</li>
  <li>Upload your tabular .csv data. Columns should represent variables, and rows should represent independent observations.</li>
</ol>
<p> If desired select:</p>
<ol>
<!--  <li>The number of bootstraps.</li> -->
  <li>The correlation measure: Pearson or Spearman.</li>
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
<p style='font-size:10px; '> Please note:<br>
The transformation is applied to continuous or ordinal variables with more than 4 categories.
Character columns, factors and binary variables will not be affected by the transformation.
.</p> "
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
#      illa <- ILAA(dataset(),thr=thrvalue,bootstrap=input$Bootstrap,method=input$corMeasure)
      illa <- ILAA(dataset(),thr=thrvalue,method=input$corMeasure)
#      illa <- ILAA(dataset(),thr=thrvalue)
      return(illa)
    }
  })
  
  # Plot of correlated data
  output$heatmapin <- renderPlotly({
    if (!is.null(ILAA_result())) {
      datavars <- rownames(attr(ILAA_result(),"UPLTM"))
      heatmaply(cor(dataset()[,datavars],method=input$corMeasure), 
                plot_method = "plotly",
                main = paste(input$corMeasure,"Correlation of the Input Data"))
    }
  })
  
  # Generate heat map of output correlation  
    output$heatmapout <- renderPlotly({
      if (!is.null(ILAA_result())) {
        datavars <- colnames(attr(ILAA_result(),"UPLTM"))
        heatmaply(cor(ILAA_result()[,datavars],method=input$corMeasure), 
                  plot_method = "plotly",
                  main = paste(input$corMeasure,"Correlation of Transformed Data"))
      }
    })

    # Report the formulas  
    output$Formulas <- renderDataTable({
      if (!is.null(ILAA_result())) {
        R2 <- 1.0 - attr(ILAA_result(),"VarRatio")
        R2 <- round(R2, digits = 3)
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
             edge.arrow.width=1.00,
             vertex.label.cex=(0.5 + 0.1*VertexSize),
             vertex.label.dist=0.75 + 0.2*VertexSize,
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
  url4 <- a("@jtamezpena", href="https://twitter.com/jtamezpena")
  output$tab <- renderUI({
    tagList("Tutorial:", url)
  })
  output$tab2 <- renderUI({
    tagList("Validation:", url2)
  })
  output$tab3 <- renderUI({
    tagList("FRESA.CAD:", url3)
  })
  output$Twitter <- renderUI({
    tagList("Twitter:", url4)
  })
  
  
}

# Run the application
shinyApp(ui, server)
