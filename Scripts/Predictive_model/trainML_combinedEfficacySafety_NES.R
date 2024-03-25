set.seed(5081)


# Load libraries
library(tidyverse)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_parallelprocesses.R")



# Define global options for this script 
# disease <- opt$disease
# data_balance_method <- opt$data_balance_method
# nproc <- opt$nproc

disease <- "KidneyCancer"
data_balance_method <- "none"
drug_target_type <- "known"
nproc <- 3



cat(paste0("\n\nTraining model for: ", disease, "\n"))
cat(paste0("\nData balance method: ", data_balance_method, "\n"))




# Read the train-test split
train_test_split <- readRDS(paste0("OutputFiles/ML_data_split/ML_dataSplit_", disease, ".rds"))


# Read the FGSEA result
fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
fgsea_result <- fgsea_result[["combinedEfficacySafety"]]








# Read the drug combination category
plot_col <- "class_EffAdv"
drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]

stat_data <- fgsea_result









# Add the drug combination categories to the stat data
stat_data <- as.data.frame(t(stat_data))
stat_data <- rownames_to_column(stat_data, "comb_name")

stat_data$category <- drugCombs_cat[,c(plot_col), drop = TRUE][match(stat_data$comb_name, drugCombs_cat$comb_name)]
stat_data <- stat_data %>% mutate_at("category", ~replace_na(., "unknown"))
stat_data <- stat_data[stat_data$category != "unknown", ]
stat_data$category <- as.factor(stat_data$category)

stat_data <- pivot_longer(data = stat_data, 
                          cols = colnames(stat_data)[!colnames(stat_data) %in% c("comb_name", "category")], 
                          cols_vary = "fastest", 
                          names_to = "feature", 
                          values_to = "value")

# Calculate the statistical difference
stat_res_final <- data.frame()
for(lib_name in unique(stat_data$feature)){
  
  stat_data_select <- stat_data[stat_data$feature == lib_name, ]
  
  
  if(grepl("^\\[ADR\\]", lib_name)){
    stat_res <- wilcox.test(x = stat_data_select$value[stat_data_select$category == "Eff"],
                            y = stat_data_select$value[stat_data_select$category == "Adv"], 
                            alternative = "less")
  }
  
  if(grepl("^\\[DISEASE\\]", lib_name)){
    stat_res <- wilcox.test(x = stat_data_select$value[stat_data_select$category == "Eff"],
                            y = stat_data_select$value[stat_data_select$category == "Adv"], 
                            alternative = "greater")
  }
  
  
  stat_res_final <- rbind(stat_res_final, data.frame("feature" = lib_name,
                                                     "W" = unname(stat_res$statistic),
                                                     "p_val" = stat_res$p.value))
}


stat_res_final <- stat_res_final[stat_res_final$p_val <= 0.0001, ]
stat_res_final <- stat_res_final[order(stat_res_final$p_val, decreasing = FALSE), ]


# # Select the feature
# safety_feature_select <- stat_res_final[grep("^\\[ADR\\]", stat_res_final$feature), ]
# if(nrow(safety_feature_select) > 0){
#   safety_feature_select <- safety_feature_select$feature[safety_feature_select$p_val == min(safety_feature_select$p_val)]
# }else{ safety_feature_select <- c() }
# 
# efficacy_feature_select <- stat_res_final[grep("^\\[DISEASE\\]", stat_res_final$feature), ]
# if(nrow(efficacy_feature_select) > 0){
#   efficacy_feature_select <- efficacy_feature_select$feature[efficacy_feature_select$p_val == min(efficacy_feature_select$p_val)]
# }else{ efficacy_feature_select <- c() }







# Sub-set the FGSEA results for only the selected results



fgsea_result_select <- fgsea_result[row.names(fgsea_result) %in% stat_res_final$feature, ]
fgsea_result_select <- as.data.frame(t(fgsea_result_select))
fgsea_result_select$category <-  drugCombs_cat[,c(plot_col), drop = TRUE][match(row.names(fgsea_result_select), drugCombs_cat$comb_name)]
fgsea_result_select <- fgsea_result_select[!is.na(fgsea_result_select$category), ]
























cl <- makeCluster(nproc)
registerDoParallel(cl) 

















print(warnings())