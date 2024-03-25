set.seed(5081)



# Use Wilcoxon's test to identify the features to use for model
# Only significant (at p <= 0.05) selected



# Load libraries
library(tidyverse)
library(unixtools)
library(optparse)



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known", 
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call. = FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call. = FALSE)
}



# Define global options for this script 
disease <- opt$disease
drug_target_type <- opt$drug_target_type



cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))



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

stat_res_final <- stat_res_final[stat_res_final$p_val <= 0.05, ]
stat_res_final <- stat_res_final[order(stat_res_final$p_val, decreasing = FALSE), ]



# Save to file 
if(!dir.exists("OutputFiles/Feature_selection/")){dir.create("OutputFiles/Feature_selection/", recursive = TRUE)}
write.csv(stat_res_final, paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_", drug_target_type, ".csv"), row.names = FALSE)



print(warnings())