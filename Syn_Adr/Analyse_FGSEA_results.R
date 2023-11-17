

disease <- "SkinCancer"
drug_target_type <- "KEGG"



# Read the FGSEA
fgsea_result <- readRDS(paste0("Syn_Adr/FGSEA_results/fgseaNES_allADR_", disease, "_", drug_target_type, ".rds"))
fgsea_result <- fgsea_result$allADR




# # For all combinations
# # Extract dataframes with only features which show enrichment
# 
# tmp1 <- fgsea_result[rowSums(fgsea_result != 0) > 0,, drop = FALSE]
# View(tmp1)
# 
# 
# # List the features and the number of the drugs combinations that show enrichment
# tmp2 <- data.frame(sort(apply(tmp1 != 0, 1, sum), decreasing = TRUE))
# colnames(tmp2) <- "count"
# View(tmp2)
# 
# 
# 
# 
# # Read the categories
# drugCombs_cat <- readRDS(file = paste0("Syn_Adr/Drug_combination_class/drugCombs_cat_4way_", disease, ".rds"))
# drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
# 
# 
# 
# 
# # For Ant_Adv
# select_drugCombs <- drugCombs_cat[drugCombs_cat$class_4way %in% "Ant_Adv", "comb_name"]
# fgsea_result_select <- fgsea_result[, select_drugCombs]
# tmp1 <- fgsea_result_select[rowSums(fgsea_result_select != 0) > 0,, drop = FALSE]
# View(tmp1)l
# tmp2 <- data.frame(sort(apply(tmp1 != 0, 1, sum), decreasing = TRUE))
# colnames(tmp2) <- "count"
# View(tmp2)
# 
# 
# 
# # For Ant_Eff
# select_drugCombs <- drugCombs_cat[drugCombs_cat$class_4way %in% "Ant_Eff", "comb_name"]
# fgsea_result_select <- fgsea_result[, select_drugCombs]
# tmp1 <- fgsea_result_select[rowSums(fgsea_result_select != 0) > 0,, drop = FALSE]
# View(tmp1)
# tmp2 <- data.frame(sort(apply(tmp1 != 0, 1, sum), decreasing = TRUE))
# colnames(tmp2) <- "count"
# View(tmp2)
# 
# 
# 
# # For Syn_Adv
# select_drugCombs <- drugCombs_cat[drugCombs_cat$class_4way %in% "Syn_Adv", "comb_name"]
# fgsea_result_select <- fgsea_result[, select_drugCombs]
# tmp1 <- fgsea_result_select[rowSums(fgsea_result_select != 0) > 0,, drop = FALSE]
# View(tmp1)
# tmp2 <- data.frame(sort(apply(tmp1 != 0, 1, sum), decreasing = TRUE))
# colnames(tmp2) <- "count"
# View(tmp2)
# 
# 
# 
# # For Syn_Eff
# select_drugCombs <- drugCombs_cat[drugCombs_cat$class_4way %in% "Syn_Eff", "comb_name"]
# fgsea_result_select <- fgsea_result[, select_drugCombs]
# tmp1 <- fgsea_result_select[rowSums(fgsea_result_select != 0) > 0,, drop = FALSE]
# View(tmp1)
# tmp2 <- data.frame(sort(apply(tmp1 != 0, 1, sum), decreasing = TRUE))
# colnames(tmp2) <- "count"
# View(tmp2)
















# Check the number of targets

drugs_with_NES <- colnames(fgsea_result[, colSums(fgsea_result != 0) > 0, drop = FALSE])
drugs_without_NES <- colnames(fgsea_result[, colSums(fgsea_result != 0) == 0, drop = FALSE])









plot_data_1 <- targets_drugs_with_NES %>% select("comb_name", contains("_count"), contains("_cnt"))
plot_data_1$NES <- "yes"
plot_data_2 <- targets_drugs_without_NES %>% select("comb_name", contains("_count"), contains("_cnt"))
plot_data_2$NES <- "no"


plot_data <- rbind(plot_data_1, plot_data_2)


plot_data <- plot_data %>% pivot_longer(cols = c(contains("count") | contains("_cnt")), 
                                        names_to = "target_type", 
                                        values_to = "count")



ggplot(plot_data, aes(x = count)) +
  geom_histogram() +
  facet_grid(rows = vars(plot_data$NES), cols = vars(plot_data$target_type))


ggplot(plot_data, aes(x = target_type, y = count, fill = NES)) +
  geom_violin() 
