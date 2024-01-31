set.seed(5081)



# Script to calculate the similarity between enrichment gene sets


# Load libraries
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_Barabasi_metrics.R")
source("Scripts/Functions/Functions_parallelprocesses.R")



enrichment_lib_list <- c()

# Compile efficacy libraries
diseases <- c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")

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



# # Compile SMPDB (Drug Metabolsim) library
# enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
# enrichment_lib <- enrichment_lib$`Drug Metabolism`
# names(enrichment_lib) <- paste0("[SMPDB_DrugMetabolism] ", names(enrichment_lib))
# enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)



# # Compile SMPDB (Drug Action) library
# enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
# enrichment_lib <- enrichment_lib$`Drug Action`
# names(enrichment_lib) <- paste0("[SMPDB_DrugAction] ", names(enrichment_lib))
# enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)



# Compile miscellaneous library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/miscellaneous_gene_lib.rds"))
names(enrichment_lib) <- paste0("[miscellaneous] ", names(enrichment_lib))
enrichment_lib_list <- c(enrichment_lib_list, enrichment_lib)

enrichment_lib_list <- lapply(enrichment_lib_list, unique)

# Read the network on which to calculate the separation
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")


cl <- makeCluster(50)
registerDoParallel(cl) 


lib_term_similarity <- foreach(lib1=names(enrichment_lib_list), 
                               .combine = "rbind", 
                               .packages = c("igraph")) %:%
  foreach(lib2=names(enrichment_lib_list), 
          .combine = "rbind", 
          .packages = c("igraph")) %dopar% {
            
            # Get the gene list
            lib1_genes <- enrichment_lib_list[[lib1]]
            lib2_genes <- enrichment_lib_list[[lib2]]
            
            len_lib1_genes <- length(lib1_genes)
            len_lib2_genes <- length(lib2_genes)
            
            # Calculate Jaccard similarity 
            intersection <- intersect(lib1_genes, lib2_genes)
            union <- union(lib1_genes, lib2_genes)
            jaccard <- length(intersection) / length(union)
            
            # Keep genes in network
            lib1_genes <- lib1_genes[lib1_genes %in% V(input_network)$name]
            lib2_genes <- lib2_genes[lib2_genes %in% V(input_network)$name]
            
            # Calculate Barabasi proximity
            proximity_separation <- Barabasi_proximity_separation(input_network, lib1_genes, lib2_genes)
            proximity_closest <- Barabasi_proximity_closest(input_network, lib1_genes, lib2_genes)
            proximity_shortest<- Barabasi_proximity_shortest(input_network, lib1_genes, lib2_genes)
            proximity_centre <- Barabasi_proximity_centre(input_network, lib1_genes, lib2_genes)
            proximity_kernel <- Barabasi_proximity_kernel(input_network, lib1_genes, lib2_genes)
            
            # Calculate the distances between the genes in the network
            distance_mat <- distances(graph = input_network, v = lib1_genes, to = lib1_genes)
            
            
            res <- data.frame(GeneSet_1 = lib1,
                              GeneSet_2 = lib2,
                              
                              GeneSet_1_size = len_lib1_genes,
                              GeneSet_2_size = len_lib2_genes,
                              
                              Jaccard = jaccard,
                              
                              GeneSet_1_inNet_size = length(lib1_genes),
                              GeneSet_2_inNet_size = length(lib2_genes),
                              
                              proximity_separation = proximity_separation,
                              proximity_closest = proximity_closest,
                              proximity_shortest = proximity_shortest, 
                              proximity_centre = proximity_centre, 
                              proximity_kernel = proximity_kernel, 
                              
                              min_distance = min(distance_mat),
                              max_distance = max(distance_mat)
            )
            
          }


stopCluster(cl)
unregister_dopar()





# Export to file
if(!dir.exists("OutputFiles/Tables/")){
  dir.create("OutputFiles/Tables/", recursive = TRUE)
}
write.csv(lib_term_similarity, "OutputFiles/Tables/Enrichment_library_similarity.csv", row.names = FALSE)
saveRDS(lib_term_similarity, "OutputFiles/Tables/Enrichment_library_similarity.rds")


print(warnings())