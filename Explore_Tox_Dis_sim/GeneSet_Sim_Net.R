


# Calculate similarity between ADR and disease genes #



# Load libraries
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")


enrichment_lib_1 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level3_lib.rds")
names(enrichment_lib_1) <- paste0("[ADR] ", names(enrichment_lib_1))

enrichment_lib_2 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level4_lib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))

enrichment_lib_3 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
names(enrichment_lib_3) <- paste0("[DISEASE] ", names(enrichment_lib_3))

enrichment_lib_4 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
names(enrichment_lib_4) <- paste0("[DISEASE] ", names(enrichment_lib_4))

enrichment_lib_5 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/msigdb_keggPath2Gene_lib.rds")
names(enrichment_lib_5) <- paste0("[PATHWAY] ", names(enrichment_lib_5))

enrichment_lib_6 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
names(enrichment_lib_6) <- paste0("[DISEASE] ", names(enrichment_lib_6))

enrichment_lib_7 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
names(enrichment_lib_7) <- paste0("[DISEASE] ", names(enrichment_lib_7))

enrichment_lib_8 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
names(enrichment_lib_8) <- paste0("[DISEASE] ", names(enrichment_lib_8))



enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2, enrichment_lib_3,
                    enrichment_lib_4, enrichment_lib_5, enrichment_lib_6,
                    enrichment_lib_7, enrichment_lib_8)


cl <- makeCluster(35)
registerDoParallel(cl) 


similarity_df <- foreach(lib_1=names(enrichment_lib), 
                         .combine = "rbind",
                         .verbose = FALSE,
                         .errorhandling = "pass") %:%
  foreach(lib_2=names(enrichment_lib),
          .combine = "rbind",
          .verbose = FALSE,
          .errorhandling = "pass") %dopar% {
            
            if(lib_1 != lib_2){
              intersection <- intersect(enrichment_lib[[lib_1]], enrichment_lib[[lib_2]])
              union <- union(enrichment_lib[[lib_1]], enrichment_lib[[lib_2]])
              jaccard <- length(intersection) / length(union)
              
              tmp1 <-  data.frame(lib_1 = lib_1,
                                  lib_2 = lib_2, 
                                  lib_1_size = length(enrichment_lib[[lib_1]]),
                                  lib_2_size = length(enrichment_lib[[lib_2]]),
                                  similarity = jaccard)
              
              tmp1
            }

          }

stopCluster(cl)
unregister_dopar()


write.csv(similarity_df, "Explore_Tox_Dis_sim/GeneSet_similarity_Net.csv", row.names = FALSE)
