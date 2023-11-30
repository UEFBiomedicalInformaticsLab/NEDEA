# Code for tab to show network from the enrichment library genes


require(igraph)
require(tidyverse)

serverFunc__enrichLibNet <- function(input, output, session) {
  
  
  observe({
    # Read the enrichment library ----
    switch(input$enrichLibNet_featureType,
           "efficacy" = {enrichment_lib <- readRDS(paste0("../InputFiles/Enrichment_analysis_libraries/Disease2Gene_", input$enrichLibNet_disease, "_lib.rds"))},
           "safety" = {enrichment_lib <- readRDS("../InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")},
           "kegg" = {enrichment_lib <- readRDS("../InputFiles/Enrichment_analysis_libraries/CHG_keggPath2Gene_lib.rds")},
           "smpdbDrugMet" = {enrichment_lib <- readRDS("../InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"); enrichment_lib <- enrichment_lib$`Drug Metabolism`},
           "smpdbDrugAct" = {enrichment_lib <- readRDS("../InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"); enrichment_lib <- enrichment_lib$`Drug Action`},
           "misc" = {enrichment_lib <- readRDS("../InputFiles/Enrichment_analysis_libraries/miscellaneous_gene_lib.rds")}
    )
    
    # Create the drop down menu to select enrichment library
    updateSelectizeInput(session, "enrichLibNet_enrichLibSelect", choices = names(enrichment_lib), server = TRUE)
    
    
    # Read the network
    switch(input$enrichLibNet_inputNetworkName,
           "STRING" = {inputNetwork <- readRDS("../InputFiles/Networks/STRING_PPI_Net.rds")})
    
    
    
    observeEvent(input$enrichLibNet_enrichLibSelect, {
      
      # Extract the sub network consisting of the selected enrichment library genes
      subnet <- induced.subgraph(graph = inputNetwork, vids = V(inputNetwork)[V(inputNetwork)$name %in% enrichment_lib[[input$enrichLibNet_enrichLibSelect]]])
      
      
      # Print the size of the network
      output$OUT_enrichLibNet_inputNetworkProp <- renderUI({
        
        text <- paste0(
          "<p>",
          "<span style = 'color:blue'>Network properties::</span>  ",
          " <b>Vertices</b> = ", vcount(subnet),
          " <b>Edges</b> = ", ecount(subnet),
          " <b>Density</b> = ", round(edge_density(subnet), 3),
          " <b>Number of components</b> = ", count_components(subnet),
          " </p>"
        )
        
        HTML(text)
        
      })
      
      
      # Show the network
      output$OUT_enrichLibNet_visSubnet <- renderVisNetwork({
        visIgraph(subnet,
                  physics = FALSE,
                  smooth = FALSE,
                  randomSeed = 5081) %>%
          visIgraphLayout(layout = "layout_nicely") %>%
          visOptions(width = "100%",
                     height = "100%",
                     highlightNearest = list("hover" = TRUE),
                     clickToUse = TRUE)
      })
    })
  })
}
