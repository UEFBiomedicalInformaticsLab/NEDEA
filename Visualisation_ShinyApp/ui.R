
require(shiny)
require(shinydashboard)
require(visNetwork)

source("uiTabItem__home.R")
source("uiTabItem__enrichLibNet.R")
source("uiTabItem__targetEffSafNet.R")


ui <- dashboardPage(
  skin = "black", 
  
  
  
  # Title of the app
  dashboardHeader(title = "DrugComb explorer"),
  
  
  # Navigation panel
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Enrichment library gene network", tabName = "enrichLibNet", icon = icon("circle-nodes")),
      menuItem("Target feature interaction", tabName = "targetEffSafNet", icon = icon("circle-nodes"))
    )
  ),
  
  
  dashboardBody(
    
    tabItems(
      
      uiTabItem__home,
      uiTabItem__enrichLibNet,
      uiTabItem__targetEffSafNet
      

      
      
    ) # end of tabItems
  ) # end of dashboardBody
) # end of dashboardPage