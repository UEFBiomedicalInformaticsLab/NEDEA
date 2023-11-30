uiTabItem__targetEffSafNet <- tabItem(
  
  tabName = "targetEffSafNet",
  
  h3("Network of the genes selected for FGSEA and the enrichment library genes"),
  
  fluidPage(
    
    
    # Start menu panel ----
    fluidRow(style = "border: 1px solid black; padding: 5px;",
             
             # Start row 1 menu ----
             fluidRow(
               
               column(selectInput(inputId = "targetEffSafNet_disease",
                                  label = "Disease",
                                  choices = c(
                                    "Breast Cancer" = "BreastCancer",
                                    "Kidney Cancer" = "KidneyCancer",
                                    "Lung Cancer" = "LungCancer",
                                    "Ovary Cancer" = "OvaryCancer",
                                    "Prostate Cancer" = "ProstateCancer",
                                    "Skin Cancer" = "SkinCancer"),
                                  selected = "Breast Cancer"
               ),
               width = 6),
               
               column(selectInput(
                 inputId = "targetEffSafNet_inputNetworkName",
                 label = "Base network",
                 choices = c("STRING"),
                 selected = "STRING"
               ),
               width = 6)
             ),
             # End row 1 menu ----
             
             # Start row 2 menu ----
             fluidRow(
               
               column(
                 shinyjs::useShinyjs(),
                 checkboxInput("targetEffSafNet_showEfficacy",
                               "Show efficacy genes",
                               value = FALSE,
                               width = NULL),
                 
                 selectizeInput("targetEffSafNet_efficacyLibSelect",
                                label = "Select efficacy gene set",
                                choices = NULL,
                                selected = NULL,
                                width = "90%"),
                 width = 5),
               
               
               column(
                 checkboxInput("targetEffSafNet_showSafety",
                               "Show safety genes",
                               value = FALSE,
                               width = NULL),
                 
                 selectizeInput("targetEffSafNet_safetyLibSelect",
                                label = "Select safety gene set",
                                choices = NULL,
                                selected = NULL,
                                width = "90%"),
                 width = 5)
               
             ),
             # End row 2 menu ----
             
             
             # Start row 3 menu ----
             fluidRow(
               
               column(selectizeInput("targetEffSafNet_drugComb",
                                     label = "Select the drug combination",
                                     choices = NULL,
                                     selected = NULL,
                                     width = "90%"),
                      width = 6),
               
               
               column(selectInput("targetEffSafNet_drugTargetType",
                                  label = "Select the drug target type",
                                  choices = c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR"), 
                                  selected = "known"),
                      width = 3),
               
               column(checkboxInput("targetEffSafNet_showRWRgenes",
                                    "Show genes selected by RWR",
                                    value = FALSE,
                                    width = NULL),
                      width = 3)
             )
             # End row 3 menu ----
             
    ),
    # End menu panel ----
    
    # Start display panel ----
    fluidRow(style = "border: 1px solid black; padding: 5px; margin-top: 10px; height: 1000px",
             
             uiOutput({ "OUT_targetEffSafNet_networkProp" }),
             visNetworkOutput("OUT_targetEffSafNet_visSubnet", width = "100%", height = "100%")
             
    )# End of display panel ----
    
  ) # End of fluidPage
) # End of tabItem










