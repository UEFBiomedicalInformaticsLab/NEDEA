set.seed(5081)


# Script to calculate the similarity between gene set libraries



# Load libraries
library(foreach)
library(doParallel)
library(openxlsx)
source("Scripts/Functions/Functions_parallelprocesses.R")



cl <- makeCluster(6)
registerDoParallel(cl) 

disease_list = c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer", "LiverCancer")

similarity_df <- foreach(disease=disease_list, 
                         .verbose = TRUE,
                         .errorhandling = "stop",
                         .final = function(x){setNames(x, disease_list)}) %dopar% {
                           
                           enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
                           names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
                           
                           enrichment_lib_2 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
                           names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
                           
                           
                           enrichment_lib_3 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/CHG_keggPath2Gene_lib.rds"))
                           names(enrichment_lib_3) <- paste0("[KEGG] ", names(enrichment_lib_3))
                           
                           enrichment_lib_4 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
                           enrichment_lib_4 <- enrichment_lib_4$`Drug Metabolism`
                           names(enrichment_lib_4) <- paste0("[DrugMet] ", names(enrichment_lib_4))


                           enrichment_lib_5 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
                           enrichment_lib_5 <- enrichment_lib_5$`Drug Action`
                           names(enrichment_lib_5) <- paste0("[DrugAction] ", names(enrichment_lib_5))

                           enrichment_lib_6 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/miscellaneous_gene_lib.rds"))
                           names(enrichment_lib_6) <- paste0("[misc] ", names(enrichment_lib_6))
                           
                           
                           enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2, enrichment_lib_3,
                                               enrichment_lib_4, enrichment_lib_5, enrichment_lib_6)

                           
                           # Use Jaccard index to calculate similarity between terms
                           lib_term_similarity <- data.frame()
                           for (i in unique(names(enrichment_lib))) {
                             for (j in unique(names(enrichment_lib))) {
                               intersection <-
                                 intersect(enrichment_lib[[i]], enrichment_lib[[j]])
                               union <- union(enrichment_lib[[i]], enrichment_lib[[j]])
                               jaccard <- length(intersection) / length(union)
                               lib_term_similarity <- rbind(lib_term_similarity,
                                                            data.frame(
                                                              geneset_1 = i,
                                                              geneset_2 = j ,
                                                              Jaccard = jaccard))
                             }
                           }
                           # similarity_df[[disease]] <- lib_term_similarity
                           lib_term_similarity
                         }


stopCluster(cl)
unregister_dopar()


if(!dir.exists("OutputFiles/Tables/")){ dir.create("OutputFiles/Tables/", recursive = TRUE) }
write.xlsx(similarity_df, 
           "OutputFiles/Tables/Enrichment_lib_similarity.xlsx", overwrite = TRUE)


print(warnings())