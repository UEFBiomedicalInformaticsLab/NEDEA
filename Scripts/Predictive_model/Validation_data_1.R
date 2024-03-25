set.seed(5081)




# Generate validation data and check applicibility domain w.r.t known drug target FGSEA 




# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call. = FALSE)
}

if(!opt$disease %in% c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  print_help(opt_parser)
  stop("--disease must be one of the following: BreastCancer, KidneyCancer, LungCancer, OvaryCancer, ProstateCancer, SkinCancer", call. = FALSE)
}


# Define global options for this script
disease <- opt$disease


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))


# Read the synergy level of the  drug combination
drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
drugCombs_data <- drugCombs_data[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "Syn_level")]


# Read the DDI data
DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
DrugBank_ddi <- DrugBank_ddi[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", paste0("ADR_", disease))]


# Merge the data
valid_drugCombs_cat <- merge(DrugBank_ddi, drugCombs_data, 
                             by = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"),
                             all.y = TRUE)


# Assign the drug combination categories
valid_drugCombs_cat$class_EffAdv <- NA

valid_drugCombs_cat[valid_drugCombs_cat$Syn_level >= 2 & valid_drugCombs_cat$Syn_level < 3 & valid_drugCombs_cat[, paste0("ADR_", disease)] %in% "unknown",]$class_EffAdv <- "Eff"
valid_drugCombs_cat[valid_drugCombs_cat$Syn_level >= 2 & valid_drugCombs_cat$Syn_level < 3 & valid_drugCombs_cat[, paste0("ADR_", disease)] %in% "adr_positive",]$class_EffAdv <- "Adv"


# valid_drugCombs_cat <- valid_drugCombs_cat[!is.na(valid_drugCombs_cat$class_EffAdv), c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "class_EffAdv")]
valid_drugCombs_cat <- valid_drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "class_EffAdv")]
valid_drugCombs_cat <- valid_drugCombs_cat[!is.na(valid_drugCombs_cat$class_EffAdv), ]
valid_drugCombs_cat$comb_name <- paste(valid_drugCombs_cat$Drug1_DrugBank_id, valid_drugCombs_cat$Drug2_DrugBank_id, sep = "_")


# Read the training set drug combination 
train_drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
train_drugCombs_cat$comb_name <- paste(train_drugCombs_cat$Drug1_DrugBank_id, train_drugCombs_cat$Drug2_DrugBank_id, sep = "_")
train_drugCombs_cat <- train_drugCombs_cat[!is.na(train_drugCombs_cat$class_EffAdv), ]


# Check for overlapping drug combinations
remove_rows <- c()
for(i in 1:nrow(valid_drugCombs_cat)){
  drug1 <- valid_drugCombs_cat[i, "Drug1_DrugBank_id"]
  drug2 <- valid_drugCombs_cat[i, "Drug2_DrugBank_id"]
  
  tmp1 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug1 & train_drugCombs_cat$Drug2_DrugBank_id == drug2, ]
  tmp2 <- train_drugCombs_cat[train_drugCombs_cat$Drug1_DrugBank_id == drug2 & train_drugCombs_cat$Drug2_DrugBank_id == drug1, ]
  
  if(nrow(tmp1) > 0 | nrow(tmp2) > 0 ){
    remove_rows <- c(remove_rows, i)
  }
}

if(length(remove_rows) > 0){
  valid_drugCombs_cat <- valid_drugCombs_cat[-remove_rows,]
}


if(!dir.exists("InputFiles/Validation_data/")){
  dir.create("InputFiles/Validation_data/", recursive = TRUE)
}
saveRDS(valid_drugCombs_cat, file = paste0("InputFiles/Validation_data/Validation1_drugCombs_cat_effVadv_", disease, ".rds"))



###
# Check applicability domain


# Read the FGSEA results
fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_known.rds"))
fgsea_result <- fgsea_result[["combinedEfficacySafety"]]


# Read the important features and select the top to plot
selected_features <- read.csv(paste0("OutputFiles/Feature_selection/NES_EfficacySafety_selectedFeatures_", disease, "_known.csv"))

safety_feature_select <- selected_features[grep("^\\[ADR\\]", selected_features$feature), ]
if(nrow(safety_feature_select) > 0){
  safety_feature_select <- safety_feature_select$feature[safety_feature_select$p_val == min(safety_feature_select$p_val)][1]
}else{ safety_feature_select <- c() }


efficacy_feature_select <- selected_features[grep("^\\[DISEASE\\]", selected_features$feature), ]
if(nrow(efficacy_feature_select) > 0){
  efficacy_feature_select <- efficacy_feature_select$feature[efficacy_feature_select$p_val == min(efficacy_feature_select$p_val)][1]
}else{ efficacy_feature_select <- c() }


# Get the FGSEA results for the train and test data
train_fgsea_result <- fgsea_result[row.names(fgsea_result) %in% c(safety_feature_select, efficacy_feature_select), 
                                   colnames(fgsea_result) %in% train_drugCombs_cat$comb_name]

valid_fgsea_result <- fgsea_result[row.names(fgsea_result) %in% c(safety_feature_select, efficacy_feature_select), 
                                   colnames(fgsea_result) %in% valid_drugCombs_cat$comb_name]


train_fgsea_result <- as.data.frame(t(train_fgsea_result))
train_fgsea_result$category <- train_drugCombs_cat$class_EffAdv[match(row.names(train_fgsea_result), train_drugCombs_cat$comb_name)]
train_fgsea_result$data_group <- "Training"

valid_fgsea_result <- as.data.frame(t(valid_fgsea_result))
valid_fgsea_result$category <- valid_drugCombs_cat$class_EffAdv[match(row.names(valid_fgsea_result), valid_drugCombs_cat$comb_name)]
valid_fgsea_result$data_group <- "Validation"


# Plot

plot_data <- rbind(train_fgsea_result, valid_fgsea_result)

colnames(plot_data)[colnames(plot_data) %in% efficacy_feature_select] <- "F1"
colnames(plot_data)[colnames(plot_data) %in% safety_feature_select] <- "F2"

x_axis_label = str_wrap(efficacy_feature_select, 40)
y_axis_label = str_wrap(safety_feature_select, 40)



if(!dir.exists("OutputFiles/Plots/validation_data_applicibility_domain/")){
  dir.create("OutputFiles/Plots/validation_data_applicibility_domain/", recursive = TRUE)
}

tiff(paste0("OutputFiles/Plots/validation_data_applicibility_domain/plot_validation1_AD_", disease, ".tiff"),
     width = 7, height = 6,
     units = "cm", compression = "lzw", res = 1200)


ggplot() +
  geom_point(data = plot_data, 
             mapping = aes(x = F1, y = F2, color = category, shape = data_group),  
             size = 0.5, 
             stroke = 0.1) +
  scale_shape_manual(values = c("Training" = 1, "Validation" = 3)) + 
  scale_color_manual(values = c("Eff" = "#0000FF", "Adv" = "#FF0000")) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "right",
        legend.key = element_rect(fill = NA), 
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 2.5),
        legend.text = element_text(size = 2),
        legend.margin = margin(1,1,1,1),
        legend.spacing = unit(0, "cm")  ) +
  labs(title = disease,
       x = x_axis_label,
       y = y_axis_label,
       color = "Category", 
       shape = "Group")

dev.off()



print(warnings())