set.seed(5081)



# Script to check the number of drug combinations that produce enrichment scores


# Load libraries
library(foreach)
library(doParallel)
library(tidyverse)
source("Scripts/Functions/Functions_parallelprocesses.R")


cl <- makeCluster(10)
registerDoParallel(cl) 

summary_df <- foreach(drug_target_type=c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all"), 
                      .packages = c("tidyverse")) %:%
  foreach(disease=c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer"), 
          .packages = c("tidyverse")) %dopar% {
            
            # print(paste(drug_target_type, disease, feature_type, sep = (" --- ")))
            
            fgsea_result <- readRDS(paste0("Extend_Eff_Saf/FGSEA_results/fgseaNES_extendedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
            fgsea_result <- fgsea_result$extendedEfficacySafety
            
            
            
            fgsea_result <- as.data.frame(t(fgsea_result))
            fgsea_result <- rownames_to_column(fgsea_result, "comb_name")
            
            tmp1 <- fgsea_result %>% 
              rowwise() %>% 
              mutate(all_total_NES = sum(c_across(!"comb_name"))
              ) %>%
              select(starts_with("all_total_"))
            
            tmp2 <- colSums(tmp1 != 0)
            
            
            # Read the drug combination category
            drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
            drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
            drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("comb_name", "class_EffAdv")]
            
            
            # Keep NES for the combinations with EFF-ADV categorization
            fgsea_result <- fgsea_result[fgsea_result$comb_name %in% drugCombs_cat$comb_name,]
            
            
            tmp3 <- fgsea_result %>% 
              rowwise() %>% 
              mutate(filtered_total_NES = sum(c_across(!"comb_name"))
              ) %>%
              select(starts_with("filtered_total_"))
            
            tmp4 <- colSums(tmp3 != 0)
            
            tmp5 <- data.frame(drug_target_type = drug_target_type,
                               disease = disease,
                               all_combs = nrow(tmp1),
                               rbind(tmp2), 
                               filtered_combs = nrow(tmp3),
                               rbind(tmp4))
            row.names(tmp5) <- NULL
            tmp5
            
            
          }

stopCluster(cl)
unregister_dopar()


summary_df <- unlist(summary_df, recursive = FALSE)
summary_df <- plyr::rbind.fill(summary_df)


# write.csv(summary_df, "OutputFiles/Tables/Summary_enriching_drugComb_count.csv", row.names = FALSE)



print(warnings())