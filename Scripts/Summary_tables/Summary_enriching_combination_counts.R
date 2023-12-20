set.seed(5081)



# Script to check the number of drug combinations that produce enrichment scores


# Load libraries
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")


cl <- makeCluster(10)
registerDoParallel(cl) 

summary_df <- foreach(drug_target_type=c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all"), 
                      .packages = c("tidyverse")) %:%
  foreach(disease=c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer"), 
          .packages = c("tidyverse")) %:%
  foreach(feature_type=c("efficacy", "safety", "combinedEfficacySafety", "kegg", "smpdbDrugMet", "smpdbDrugAct", "misc"), 
          .packages = c("tidyverse")) %dopar% {
            
            # print(paste(drug_target_type, disease, feature_type, sep = (" --- ")))
            
            # Read the FGSEA result
            if(feature_type %in% c("efficacy", "safety", "combinedEfficacySafety")){
              fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
              fgsea_result <- fgsea_result[[feature_type]]
            }
            
            if(feature_type %in% c("kegg", "smpdbDrugMet", "smpdbDrugAct")){
              fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_Pathway_", disease, "_", drug_target_type, ".rds"))
              fgsea_result <- fgsea_result[[feature_type]]
            }
            
            if(feature_type %in% c("misc")){
              fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_Miscellaneous_", disease, "_", drug_target_type, ".rds"))
              fgsea_result <- fgsea_result[[feature_type]]
            }
            
            
            fgsea_result <- as.data.frame(t(fgsea_result))
            fgsea_result <- rownames_to_column(fgsea_result, "comb_name")
            
            
            if(feature_type %in% c("combinedEfficacySafety")){
              tmp1 <- fgsea_result %>% 
                rowwise() %>% 
                mutate(all_total_NES = sum(c_across(starts_with("[ADR]") | starts_with("[DISEASE]"))),
                       all_total_DISEASE_NES = sum(c_across(starts_with("[DISEASE]"))),
                       all_total_ADR_NES = sum(c_across(starts_with("[ADR]")))) %>%
                select(starts_with("all_total_"))
              
              tmp2 <- colSums(tmp1 != 0)
            } else {
              tmp1 <- fgsea_result %>% 
                rowwise() %>% 
                mutate(all_total_NES = sum(c_across(!"comb_name"))
                ) %>%
                select(starts_with("all_total_"))
              
              tmp2 <- colSums(tmp1 != 0)
            }
            
            
            # Read the drug combination category
            drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
            drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
            drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("comb_name", "class_EffAdv")]
            
            
            # Keep NES for the combinations with EFF-ADV categorization
            fgsea_result <- fgsea_result[fgsea_result$comb_name %in% drugCombs_cat$comb_name,]
            
            
            if(feature_type %in% c("combinedEfficacySafety")){
              tmp3 <- fgsea_result %>% 
                rowwise() %>% 
                mutate(filtered_total_NES = sum(c_across(starts_with("[ADR]") | starts_with("[DISEASE]"))),
                       filtered_total_DISEASE_NES = sum(c_across(starts_with("[DISEASE]"))),
                       filtered_total_ADR_NES = sum(c_across(starts_with("[ADR]")))) %>%
                select(starts_with("filtered_total_"))
              
              tmp4 <- colSums(tmp3 != 0)
            } else {
              tmp3 <- fgsea_result %>% 
                rowwise() %>% 
                mutate(filtered_total_NES = sum(c_across(!"comb_name"))
                ) %>%
                select(starts_with("filtered_total_"))
              
              tmp4 <- colSums(tmp3 != 0)
            }
            
            
            # Check the number of drug combinations per class having NES
            tmp5 <- fgsea_result %>% select(-comb_name) %>% rowSums()  
            names(tmp5) <- fgsea_result$comb_name
            tmp5 <- as.data.frame(tmp5)
            tmp5 <- merge(drugCombs_cat, tmp5, by.x = "comb_name", by.y = 0)
            tmp5 <- tmp5[tmp5$tmp5 != 0, ]
            tmp5$class_EffAdv <- paste0(tmp5$class_EffAdv, "_with_NES")
            tmp5 <- table(tmp5$class_EffAdv)

            
            res_df <- data.frame(feature_type = feature_type,
                                 drug_target_type = drug_target_type,
                                 disease = disease,
                                 all_combs = nrow(tmp1),
                                 rbind(tmp2), 
                                 filtered_combs = nrow(tmp3),
                                 rbind(tmp4),
                                 rbind(tmp5))
            row.names(res_df) <- NULL
            res_df
            
            
          }

stopCluster(cl)
unregister_dopar()


summary_df <- unlist(unlist(summary_df, recursive = FALSE), recursive = FALSE)
summary_df <- plyr::rbind.fill(summary_df)


write.csv(summary_df, "OutputFiles/Tables/Summary_enriching_drugComb_count.csv", row.names = FALSE)



print(warnings())