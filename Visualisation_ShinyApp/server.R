

require(shiny)
require(igraph)
require(visNetwork)
require(org.Hs.eg.db)

source("serverFunc__enrichLibNet.R")
source("serverFunc__targetEffSafNet.R")



# Define server logic to plot  ----
server <- function(input, output, session) {
  
  
  serverFunc__enrichLibNet(input, output, session)


  serverFunc__targetEffSafNet(input, output, session)
  
  

  
  # output$OUT_output_contents <- renderText({ 
  #   paste(class(output))
  # })
  
}