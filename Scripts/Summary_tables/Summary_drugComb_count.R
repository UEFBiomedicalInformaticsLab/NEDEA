set.seed(5081)


# script to summarise the number of drug combinations by category


# Load libraries
library(tidyverse)


#####


summary_df <- data.frame()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
 
  drugCombs_cat_1 <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
  drugCombs_cat_1$comb_name <- paste(drugCombs_cat_1$Drug1_DrugBank_id, drugCombs_cat_1$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat_1 <- drugCombs_cat_1[, c("comb_name", "Syn_level", "class_synergyScore")]
  
  drugCombs_cat_2 <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
  drugCombs_cat_2$comb_name <- paste(drugCombs_cat_2$Drug1_DrugBank_id, drugCombs_cat_2$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat_2 <- drugCombs_cat_2[, c("comb_name", "ADR_status")]

  drugCombs_cat_3 <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_cat_3$comb_name <- paste(drugCombs_cat_3$Drug1_DrugBank_id, drugCombs_cat_3$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat_3 <- drugCombs_cat_3[, c("comb_name", "class_EffAdv")]
  

  drugCombs_cat <- reduce(list(drugCombs_cat_1, drugCombs_cat_2, drugCombs_cat_3), merge, by = "comb_name", all.x = TRUE)
  drugCombs_cat <- column_to_rownames(drugCombs_cat, "comb_name")
  

  # Replace the missing categories as unknown
  drugCombs_cat <- drugCombs_cat %>% 
    mutate_at(c("Syn_level", "class_synergyScore", "ADR_status", "class_EffAdv"), 
              ~replace_na(., "unknown"))
  
  
  drugCombs_cat <- apply(drugCombs_cat, 2, function(x){as.data.frame(table(x, useNA = "ifany"))})
  drugCombs_cat <- bind_rows(drugCombs_cat, .id = "category_type")
  drugCombs_cat$disease <- disease

  drugCombs_cat <- pivot_wider(drugCombs_cat, names_from = x, values_from = Freq)
  
  summary_df <- rbind(summary_df, drugCombs_cat)
}


summary_df <- summary_df[order(summary_df$category_type, decreasing = TRUE), ]
summary_df <- summary_df %>% select(c("category_type", "disease", "unknown", everything()))
summary_df$total_combs <- rowSums(summary_df[, -c(1:2)], na.rm = TRUE)
summary_df <- apply(summary_df, 2, as.character)
summary_df[is.na(summary_df)] <- ""


write.csv(summary_df, "OutputFiles/Tables/Summary_drugComb_count.csv", row.names = FALSE)


#####


print(warnings())