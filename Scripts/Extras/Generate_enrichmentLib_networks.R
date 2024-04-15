set.seed(5081)



# Script to visualise the networks formed by the enrichment analysis libs



# Load libraries
library(igraph)
library(visNetwork)


# Read the network from which to extract the sub-networks
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))



enrichment_lib_list <- c()

# Compile efficacy libraries
diseases <- c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")
for(disease in diseases){
  enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  names(enrichment_lib) <- paste0("[Efficacy_", disease, "] ", names(enrichment_lib))
  enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)
}


# Compile safety library
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
names(enrichment_lib) <- paste0("[Safety] ", names(enrichment_lib))
enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)


# Compile KEGG library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/CHG_keggPath2Gene_lib.rds"))
names(enrichment_lib) <- paste0("[KEGG] ", names(enrichment_lib))
enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)


# Compile SMPDB (Drug Metabolism) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Metabolism`
names(enrichment_lib) <- paste0("[SMPDB_DrugMetabolism] ", names(enrichment_lib))
enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)


# Compile SMPDB (Drug Action) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Action`
names(enrichment_lib) <- paste0("[SMPDB_DrugAction] ", names(enrichment_lib))
enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)


# Compile miscellaneous library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/miscellaneous_gene_lib.rds"))
names(enrichment_lib) <- paste0("[miscellaneous] ", names(enrichment_lib))
enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)

enrichment_lib_list <- lapply(enrichment_lib_list, unique)

rm(enrichment_lib)


# Generate the graphs
if(!dir.exists("OutputFiles/Network_visualisation/Enrichment_lib_net/")){
  dir.create("OutputFiles/Network_visualisation/Enrichment_lib_net/", recursive = TRUE)
}
setwd("OutputFiles/Network_visualisation/Enrichment_lib_net/")


for(lib_name in names(enrichment_lib_list)){
  
  
  enrichment_lib <- enrichment_lib_list[[lib_name]]
  subnet <- induced.subgraph(graph = input_network, vids = V(input_network)[V(input_network)$name %in% enrichment_lib])
  
  
  if(ecount(subnet) == 0){
    cat(paste0("\n- Graph not generated due to missing edges for ", lib_name, 
               "\n\t -- Nodes: ", vcount(subnet), " Edges: ", ecount(subnet)))
  }else{
    visIgraph(subnet, 
              physics = FALSE, 
              smooth = FALSE, 
              randomSeed = 5081) %>%
      visIgraphLayout(layout = "layout_nicely") %>% 
      visOptions(width = 1920, 
                 height = 1200,
                 highlightNearest = list("hover" = TRUE),
                 clickToUse = TRUE) %>%
      visSave(file = paste0(gsub("/", "", lib_name), 
                            ".html"), 
              selfcontained = TRUE)
  }
  
}



print(warnings())
