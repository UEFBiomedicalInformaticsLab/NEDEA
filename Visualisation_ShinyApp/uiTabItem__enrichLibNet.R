

uiTabItem__enrichLibNet <- tabItem(

  tabName = "enrichLibNet",

  h3("Network of the genes in the enrichment library genes"),

  fluidPage(

    # Define the menu ----
    fluidRow(style = "border: 1px solid black; padding: 5px;",

             # Start of row 1 menu ----
             fluidRow(


               column(selectInput(inputId = "enrichLibNet_disease",
                                  label = "Disease",
                                  choices = c("Breast Cancer" = "BreastCancer",
                                              "Kidney Cancer" = "KidneyCancer",
                                              "Lung Cancer" = "LungCancer",
                                              "Ovary Cancer" = "OvaryCancer",
                                              "Prostate Cancer" = "ProstateCancer",
                                              "Skin Cancer" = "SkinCancer")
               ),
               width = 4),

               column(selectInput(inputId = "enrichLibNet_inputNetworkName",
                                  label = "Base network",
                                  choices = c("STRING")
               ),
               width = 4),

               column(selectInput(inputId = "enrichLibNet_featureType",
                                  label = "Feature type",
                                  choices = c("Efficacy" = "efficacy",
                                              "Safety" = "safety",
                                              "KEGG pathways" = "kegg",
                                              "SMPDB [Drug Metabolism]" = "smpdbDrugMet",
                                              "SMPDB [Drug Action]" = "smpdbDrugAct",
                                              "Miscellaneous" = "misc")
               ),
               width = 4)
             ),
             # End of Row 1 menu ----

             # Start of row 2 menu ----
             fluidRow(
               column(selectizeInput("enrichLibNet_enrichLibSelect",
                                     label = "Select enrichment library",
                                     choices = NULL,
                                     selected = NULL,
                                     width = "90%"
               ),
               width = 12)
             )
             # End of two 2 menu

    ),
    # End of menu panel ----

    # Start of display panel ----
    fluidRow(style = "border: 1px solid black; padding: 5px; margin-top: 10px; height: 900px",

             uiOutput("OUT_enrichLibNet_inputNetworkProp"),

             visNetworkOutput("OUT_enrichLibNet_visSubnet", width = "100%", height = "100%")
    )
    # End of display panel ----

  ) # End of fluidPage
) # End of tabItem