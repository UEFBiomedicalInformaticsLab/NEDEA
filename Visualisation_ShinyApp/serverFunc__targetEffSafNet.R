# Code to show the network of drug targets with the efficacy and safety library



serverFunc__targetEffSafNet <- function(input, output, session){
  
  
  rv <- reactiveValues(inputNetwork = NULL,
                       enrichment_lib_efficacy = NULL,
                       enrichment_lib_safety = NULL,
                       efficacy_genes = NULL,
                       safety_genes = NULL,
                       drugCombs_targets = NULL,
                       drug_target_col = NULL,
                       target_set = NULL, 
                       nodeLabels = NULL)
  
  
  
  # Read the network
  observe({
    rv$inputNetwork <- switch(input$enrichLibNet_inputNetworkName,
                              "STRING" = {readRDS("../InputFiles/Networks/STRING_PPI_Net.rds")})
  })
  
  
  # Select the efficacy library
  observeEvent(input$targetEffSafNet_showEfficacy, {
    if(input$targetEffSafNet_showEfficacy == FALSE) {
      shinyjs::disable("targetEffSafNet_efficacyLibSelect") 
      rv$efficacy_genes <- NULL
    } else {
      shinyjs::enable("targetEffSafNet_efficacyLibSelect")
    }
  }, ignoreNULL = TRUE)
  
  observeEvent(input$targetEffSafNet_disease, {
    rv$enrichment_lib_efficacy <- readRDS(paste0("../InputFiles/Enrichment_analysis_libraries/Disease2Gene_", input$targetEffSafNet_disease, "_lib.rds"))
    updateSelectizeInput(session, "targetEffSafNet_efficacyLibSelect", 
                         choices = names(rv$enrichment_lib_efficacy), 
                         selected = names(rv$enrichment_lib_efficacy)[1],
                         server = TRUE)
  })
  
  
  
  # Select the safety library
  observeEvent(input$targetEffSafNet_showSafety  , {
    if(input$targetEffSafNet_showSafety == FALSE) {
      shinyjs::disable("targetEffSafNet_safetyLibSelect") 
      rv$safety_genes <- NULL
    } else {
      shinyjs::enable("targetEffSafNet_safetyLibSelect")
    }
  }, ignoreNULL = TRUE)
  
  observeEvent(input$targetEffSafNet_disease, {
    rv$enrichment_lib_safety <- readRDS("../InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
    updateSelectizeInput(session, "targetEffSafNet_safetyLibSelect", 
                         choices = names(rv$enrichment_lib_safety), 
                         selected = names(rv$enrichment_lib_safety)[1],
                         server = TRUE)
  })
  
  
  # Read the drug combinations and take the input
  
  observeEvent(input$targetEffSafNet_disease, {
    rv$drugCombs_targets <- readRDS(paste0("../InputFiles/Drug_combination_targets/drugCombs_targets_extended_", input$targetEffSafNet_disease, ".rds"))
    rv$drugCombs_targets$comb_name <- paste(rv$drugCombs_targets$Drug1_DrugBank_id,
                                            rv$drugCombs_targets$Drug2_DrugBank_id,
                                            sep = "_")
    updateSelectizeInput(session, "targetEffSafNet_drugComb", 
                         choices = rv$drugCombs_targets$comb_name, 
                         selected = NULL,
                         server = TRUE)
  })
  
  
  # Get the gene list
  observeEvent(input$targetEffSafNet_efficacyLibSelect, {
    if(input$targetEffSafNet_showEfficacy){
      rv$efficacy_genes <- lapply(rv$enrichment_lib_efficacy, function(x){x[x %in% V(rv$inputNetwork)$name]})
      rv$efficacy_genes <- rv$efficacy_genes[[input$targetEffSafNet_efficacyLibSelect]]
    }
  }, ignoreNULL = FALSE)
  
  
  observeEvent(input$targetEffSafNet_safetyLibSelect, {
    if(input$targetEffSafNet_showSafety){
      rv$safety_genes <- lapply(rv$enrichment_lib_safety, function(x){x[x %in% V(rv$inputNetwork)$name]})
      rv$safety_genes <- rv$safety_genes[[input$targetEffSafNet_safetyLibSelect]]
    }
  }, ignoreNULL = FALSE)
  
  
  
  
  observe({
    
    rv$drug_target_col <- switch(input$targetEffSafNet_drugTargetType,
                                 "known" = {c("drugTarget_geneSymbol") },
                                 "PS" = {c("ext_PS_targets") },
                                 "SIGNOR" = {c("ext_SIGNOR_targets") },
                                 "NPA" = {c("ext_NPA_targets") },
                                 "RI" = {c("ext_RI_targets") },
                                 "KEGG" = {c("ext_KEGG_targets") })
  })
  
  
  observeEvent({input$targetEffSafNet_drugComb
    input$targetEffSafNet_drugTargetType
  }, 
  {
    if(input$targetEffSafNet_drugComb != ""){
      rv$target_set <- rv$drugCombs_targets[rv$drugCombs_targets$comb_name %in% input$targetEffSafNet_drugComb, c(rv$drug_target_col), drop = TRUE]
      rv$target_set <- unlist(strsplit(rv$target_set, ","))
      suppressMessages(rv$target_set <- AnnotationDbi::select(org.Hs.eg.db,
                                                              keys = rv$target_set,
                                                              columns = "ENSEMBL",
                                                              keytype = "SYMBOL"))
      rv$target_set <- rv$target_set$ENSEMBL
      rv$target_set <- rv$target_set[rv$target_set %in% V(rv$inputNetwork)$name]
    }
  })
  
  
  
  
  
  
  # observeEvent(input$targetEffSafNet_showRWRgenes,
  #              {
  #               if(input$targetEffSafNet_showRWRgenes){
  #                 rwr_data <- readRDS(file = rwr_result_file_path)
  #                 rwr_data_select <- rwr_data[, drug_comb_name]
  #                 rwr_threshold <- func_RWR_threshold(rwr_data_select[!names(rwr_data_select) %in% target_set])
  #                 rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold$ELB], decreasing = TRUE)
  #               }
  #                 
  #              })
  
  
  
  
  
  
  observeEvent(
    
    {input$targetEffSafNet_drugComb
      input$targetEffSafNet_drugTargetType
      req(rv$target_set)
      input$targetEffSafNet_efficacyLibSelect
      input$targetEffSafNet_safetyLibSelect}, 
    
    {  
      rv$nodeLabels <- as.data.frame(matrix(NA, 
                                            nrow = vcount(rv$inputNetwork), 
                                            ncol = 5, 
                                            dimnames = list(c(V(rv$inputNetwork)$name), 
                                                            c("Efficacy", "Safety", "Targets", "RWR_genes", "RWR_probs")
                                            )
      ))
      

      
      rv$nodeLabels[rv$target_set, ]$Targets <- "target"
      
      if( input$targetEffSafNet_showEfficacy & length(rv$efficacy_genes) > 0 ){ rv$nodeLabels[rv$efficacy_genes, ]$Efficacy <- "efficacy" }
      if( input$targetEffSafNet_showSafety & length(rv$safety_genes) > 0 ){ rv$nodeLabels[rv$safety_genes, ]$Safety <- "safety" }
      

      
      rv$nodeLabels <- rv$nodeLabels %>%
        unite(col = label, Efficacy, Safety, Targets, RWR_genes,
              sep = "_", remove = FALSE, na.rm = TRUE)
      
      
      if(nrow(rv$nodeLabels[grep("_", rv$nodeLabels$label), ]) > 0 ){
        rv$nodeLabels[grep("_", rv$nodeLabels$label), ]$label <- "overlap"
      }
      rv$nodeLabels[rv$nodeLabels$label == "", ]$label <- "others"
      
      
      
      rv$nodeLabels$color <- unlist(sapply(rv$nodeLabels$label, 
                                           FUN = function(x){switch(x, 
                                                                    "efficacy" = "green", 
                                                                    "safety" = "red",
                                                                    "overlap" = "orange",
                                                                    "rwrGenes" = "cyan", 
                                                                    "target" = "blue",
                                                                    "others" = "lightblue"
                                           )}, 
                                           USE.NAMES = FALSE))
      
      
      # Map the node attributes to the graph
      for(attr in colnames(rv$nodeLabels)){
        
        print(paste(attr, length(rv$nodeLabels[,attr]), sep = "___"))
        rv$inputNetwork <- set_vertex_attr(graph = rv$inputNetwork,
                                           name = attr,
                                           index = V(rv$inputNetwork)[V(rv$inputNetwork)$name %in% row.names(rv$nodeLabels)],
                                           value = rv$nodeLabels[,attr])
      }
      
      
      

      
      
      select_nodes <- rv$nodeLabels[rv$nodeLabels$label != "others", ]
      subnet <- induced_subgraph(graph = rv$inputNetwork, 
                                 vids = V(rv$inputNetwork)[V(rv$inputNetwork)$name %in% row.names(select_nodes)])
      
      

      
      
      # Show the network
      output$OUT_targetEffSafNet_visSubnet <- renderVisNetwork({
        visIgraph(subnet, 
                  physics = FALSE, 
                  smooth = FALSE, 
                  randomSeed = 5081) %>%
          visIgraphLayout(layout = "layout_nicely") %>% 
          visOptions(width = 1920, 
                     height = 800,
                     highlightNearest = list("hover" = TRUE),
                     clickToUse = TRUE)
      })
      
      
      
      
    })
  
  
  
  
  
  
  
  output$OUT_targetEffSafNet_networkProp <- renderUI({
    
    text <- paste0("<p>", 
                   "<br>  <b>Diease:</b>  &nbsp", input$enrichLibNet_disease, 
                   "      <b>Network:</b>  &nbsp", input$enrichLibNet_inputNetworkName, 
                   "<br>  <b>Efficacy:</b>  &nbsp", input$targetEffSafNet_efficacyLibSelect, " (genes = ", length(rv$efficacy_genes), ")",
                   "<br>  <b>Safety:</b>  &nbsp", input$targetEffSafNet_safetyLibSelect, " (genes = ", length(rv$safety_genes), ")",
                   "<br>  <b>Drug combination:</b>  &nbsp", input$targetEffSafNet_drugComb,
                   "      <b>Target:</b>  &nbsp", input$targetEffSafNet_drugTargetType, " (genes = ", length(rv$target_set), ")",
                   "</p>"
    )
    
    HTML(text)
    
  })
  
  
  
  
  
}

