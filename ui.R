library(shiny)
options(java.parameters = "- Xmx1024m")
library(pca3d)
library(shinyRGL)
library(rgl)
library(d3heatmap)
library(htmlwidgets)

shinyUI(fluidPage(
  h1("PCA/HeatMap"),
  sidebarPanel(
    fileInput("inputfile","Upload Data File",accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    
    conditionalPanel(condition="input.tabs=='Heat Map'",checkboxInput("cluster","Cluster",value = FALSE)),
    conditionalPanel(condition="input.tabs=='Home'",selectInput("imputation_type","Imputation",choices = c("Zero","Stratified Average"))),
    conditionalPanel(condition="input.tabs!='Haaaaaome'",sliderInput("cf_level","Confidence Level",min = 0,max=1,step = 0.01,value = 0.05)),
    conditionalPanel(condition="input.tabs=='PCA 3D'",numericInput("width","Width (IN)",value=4)),
    conditionalPanel(condition="input.tabs=='Heat Map'",checkboxInput("pro","Protein Name",value=T)),
    conditionalPanel(condition="input.tabs=='Heat Map'",checkboxInput("gene","Gene Name",value=T)),
    conditionalPanel(condition="input.tabs=='PCA 3D'",numericInput("resol","Resolution (DPI)",value=600)),
    conditionalPanel(condition="input.tabs=='Heat Map'",numericInput("heatmap_width", "Width (in)",value = 4)),
    conditionalPanel(condition="input.tabs=='Heat Map'",numericInput("heatmap_height", "Height (in)",value = 4)),
    conditionalPanel(condition="input.tabs=='Heat Map'",numericInput("heatmap_resolution", "Resolution (dpi)",value = 300)),
    conditionalPanel(condition="input.tabs=='PCA 3D'",selectInput("type_of_download",label = "Format",choices = c("PNG","TIFF","JPEG","PDF","SVG","PS"))),
    conditionalPanel(condition="input.tabs=='Heat Map'",selectInput("type_of_download_new",label = "Format",choices = c("PNG","TIFF","JPEG","PDF"))),
    conditionalPanel(condition="input.tabs!='Home'",textInput("change_wd","Download Directory",value = paste0(getwd(),"/"))),
    conditionalPanel(condition="input.tabs!='Home'",textInput("file_name",label = "File Name",value = Sys.Date())),
    conditionalPanel(condition="input.tabs=='PCA 3D'",actionButton("download_link_pca3d","Download")),
    conditionalPanel(condition="input.tabs=='Heat Map'",actionButton("download_heatmap", "Download"))
  ,width = 3),

  mainPanel(
    tabsetPanel(id="tabs",
    tabPanel("Home",tags$b(h4(textOutput("Preview"))),tableOutput("Datafile_show")),
    tabPanel("PCA 3D", plotOutput("NORMAL3d")),
    tabPanel("Heat Map",
             tags$div(id = "plotDiv", style = "min-height: 100px; min-width: 100px; border: 2px solid; padding: 20px; resize: both; overflow: auto;",d3heatmapOutput("heatmap", height = "400px", width = "auto")),
             tags$script('document.getElementById("plotDiv").onmouseup = function(e) {$(window).trigger("resize"); document.getElementById("heatmap").style.height = (parseFloat(this.style.height) - 30) + \'px\';};')
             
    ))
  )
  )
)
  