
require(shiny)
require(shinydashboard)
require(visNetwork)


# Define UI 
ui <- dashboardPage(
  skin = "black", 
  
  # Title of the app
  dashboardHeader(title = "DrugComb explorer"),
  
  
  # Navigation panel
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Enrichment library gene network", tabName = "enrichLib_net", icon = icon("circle-nodes")),
      menuItem("Target feature interaction", tabName = "target_feature_ixn", icon = icon("circle-nodes"))
    )
  ),
  
  # Content of the main panel
  dashboardBody(
    
    tabItems(
      tabItem(tabName = "home", h3("Home page")),
      
      
      tabItem(tabName = "enrichLib_net", 
              h3("Network of the genes in the enrichment library genes"),
              fluidPage(
                
                # Define the menu ----
                fluidRow(
                  style = "border: 1px solid black; padding: 5px;",
                  
                  # Row 1 menu ----
                  fluidRow(
                    column(selectInput(inputId = "disease", 
                                       label = "Disease", 
                                       choices = c("Breast Cancer" = "BreastCancer",
                                                   "Kidney Cancer" = "KidneyCancer",
                                                   "Lung Cancer" = "LungCancer",
                                                   "Ovary Cancer" = "OvaryCancer",
                                                   "Prostate Cancer" = "ProstateCancer",
                                                   "Skin Cancer" = "SkinCancer")), 
                           width = 4),
                    
                    column(selectInput(inputId = "input_network_name", 
                                       label = "Base network", 
                                       choices = c("STRING")), 
                           width = 4),
                    
                    column(selectInput(inputId = "feature_type", 
                                       label = "Feature type", 
                                       choices = c("Efficacy" = "efficacy",
                                                   "Safety" = "safety",
                                                   "KEGG pathways" = "kegg",
                                                   "SMPDB [Drug Metabolism]" = "smpdbDrugMet",
                                                   "SMPDB [Drug Action]" = "smpdbDrugAct",
                                                   "Miscellaneous" = "misc")), 
                           width = 4)
                  ), 
                  
                  # Row 2 menu ----
                  fluidRow(
                    column(selectizeInput("enrich_lib_select", 
                                          label = "Select enrichment library", 
                                          choices = NULL, 
                                          selected = NULL,
                                          width = "90%"), 
                           width = 12)
                    
                  )
                ),
                
                
                
                fluidRow(style = "border: 1px solid black; padding: 5px; margin-top: 10px; height: 900px",
                         
                         uiOutput("OUT_input_network_size"),
                         
                         visNetworkOutput("OUT_vis_subnet", width = "100%", height = "100%")
                         
                         
                ))
              ),

      
      tabItem(tabName = "target_feature_ixn", h3("Network of the genes selected for FGSEA and the enrichment library genes"))
      
    )
    
    
    
    
  )
)
