# Code to show the network of drug targets with the efficacy and safety library



serverFunc__targetEffSafNet <- function(input, output, session){
  
  
  
  
  rv <- reactiveValues(
    enrichment_lib_efficacy = NULL,
    enrichment_lib_safety = NULL,
    target_set = NULL,
    drugCombs_targets = NULL
  )
  
  
  # Read the network
  # switch(input$enrichLibNet_inputNetworkName,
  #        "STRING" = {inputNetwork <- readRDS("../InputFiles/Networks/STRING_PPI_Net.rds")})
  
  inputNetwork <- reactive({
    switch(input$enrichLibNet_inputNetworkName,
           "STRING" = {readRDS("../InputFiles/Networks/STRING_PPI_Net.rds")})
  })
  
  
  # Select the efficacy library
  observeEvent(input$targetEffSafNet_showEfficacy, {
    if(input$targetEffSafNet_showEfficacy == FALSE) {
      shinyjs::disable("targetEffSafNet_efficacyLibSelect") 
    } else {
      shinyjs::enable("targetEffSafNet_efficacyLibSelect")
    }
  }, ignoreNULL = TRUE)
  
  observeEvent(input$targetEffSafNet_disease, {
    rv$enrichment_lib_efficacy <- readRDS(paste0("../InputFiles/Enrichment_analysis_libraries/Disease2Gene_", input$targetEffSafNet_disease, "_lib.rds"))
    updateSelectizeInput(session, "targetEffSafNet_efficacyLibSelect", choices = names(enrichment_lib_efficacy), server = TRUE)
  })
  
  
  
  
  
  # Select the safety library
  observeEvent(input$targetEffSafNet_showSafety  , {
    if(input$targetEffSafNet_showSafety == FALSE) {
      shinyjs::disable("targetEffSafNet_safetyLibSelect") 
    } else {
      shinyjs::enable("targetEffSafNet_safetyLibSelect")
    }
  }, ignoreNULL = TRUE)
  
  observeEvent(input$targetEffSafNet_disease, {
    rv$enrichment_lib_safety <- readRDS("../InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
    updateSelectizeInput(session, "targetEffSafNet_safetyLibSelect", choices = names(enrichment_lib_safety), server = TRUE)
  })
  
  
  
  # Read the drug combinations and take the input
  
  observeEvent(input$targetEffSafNet_disease, {
    rv$drugCombs_targets <- readRDS(paste0("../InputFiles/Drug_combination_targets/drugCombs_targets_extended_", input$targetEffSafNet_disease, ".rds"))
    rv$drugCombs_targets$comb_name <- paste(rv$drugCombs_targets$Drug1_DrugBank_id,
                                            rv$drugCombs_targets$Drug2_DrugBank_id,
                                         sep = "_")
    updateSelectizeInput(session, "targetEffSafNet_drugComb", choices = rv$drugCombs_targets$comb_name, server = TRUE)
    
  })
  
  
  
  
  # Get the gene list
  
  observe({
    
    if(input$targetEffSafNet_showEfficacy){
      rv$enrichment_lib_efficacy <- lapply(rv$enrichment_lib_efficacy, function(x){x[x %in% V(inputNetwork)$name]})
      rv$enrichment_lib_efficacy <- rv$enrichment_lib_efficacy[[input$targetEffSafNet_efficacyLibSelect]]
    }
    
    if(input$targetEffSafNet_showSafety){
      rv$enrichment_lib_safety <- lapply(rv$enrichment_lib_safety, function(x){x[x %in% V(inputNetwork)$name]})
      rv$enrichment_lib_safety <- rv$enrichment_lib_safety[[input$targetEffSafNet_safetyLibSelect]]
    }
    
    
    switch(input$targetEffSafNet_drugTargetType,
           "known" = { drug_target_col <- c("drugTarget_geneSymbol") },
           "PS" = { drug_target_col <- c("ext_PS_targets") },
           "SIGNOR" = { drug_target_col <- c("ext_SIGNOR_targets") },
           "NPA" = { drug_target_col <- c("ext_NPA_targets") },
           "RI" = { drug_target_col <- c("ext_RI_targets") },
           "KEGG" = { drug_target_col <- c("ext_KEGG_targets") })
    
    rv$target_set <- rv$drugCombs_targets[rv$drugCombs_targets$comb_name %in% input$targetEffSafNet_drugComb, drug_target_col, drop = TRUE]
    rv$target_set <- unlist(strsplit(rv$target_set, ","))
    suppressMessages(rv$target_set <- select(org.Hs.eg.db,
                                                     keys = rv$target_set,
                                                     columns = "ENSEMBL",
                                                     keytype = "SYMBOL"))
    rv$target_set <- rv$target_set$ENSEMBL
    
    
    
  })
  
  
  
  
  
  output$OUT_targetEffSafNet_text <- renderText({ paste0("Disease: ", input$targetEffSafNet_disease,
                                                         "Drug comb: ", input$targetEffSafNet_drugComb,
                                                         "Targets: ", class(rv$target_set) #paste(target_set, collapse = ", ")
  )
  })
  
}

