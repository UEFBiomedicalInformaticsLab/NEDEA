

require(shiny)
require(visNetwork)


source("serverFunc__enrichLib_net.R")




# Define server logic to plot  ----
server <- function(input, output, session) {
  
  
  serverFunc__enrichLib_net(input, output, session)

  
  
  # output$OUT_output_contents <- renderText({ 
  #   paste(class(output))
  # })
  
}