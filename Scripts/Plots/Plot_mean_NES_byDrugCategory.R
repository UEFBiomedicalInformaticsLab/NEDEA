set.seed(5081)



# Script to plot the mean NES from FGSEA for each drug category

# Load libraries
library(tidyverse)




# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known",
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character"),
  make_option(c("--feature_type"), type = "character", default = NULL,
              help = "The feature type to use for plotting. Possible values: efficacy, safety, combinedEfficacySafety. Default: NULL", metavar = "character"),
  make_option(c("--drug_category_type"), type = "character", default = NULL,
              help = "The drug category type that will be plotted. Possible values: pk, adr. Default: NULL", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
}

if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call.=FALSE)
}

if(!opt$feature_type %in% c("efficacy", "safety", "combinedEfficacySafety")){
  print_help(opt_parser)
  stop("--feature_type should be: efficacy, safety, combinedEfficacySafety", call.=FALSE)
}

if(is.null(opt$drug_category_type)){
  print_help(opt_parser)
  stop("--drug_category_type argument needed", call.=FALSE)
}

if(!opt$drug_category_type %in% c("pk", "adr")){
  print_help(opt_parser)
  stop("--drug_category_type should be: pk, adr", call.=FALSE)
}


# Define global options for this script
disease <- opt$disease
drug_target_type <- opt$drug_target_type
feature_type <- opt$feature_type
drug_category_type <- opt$drug_category_type


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))
cat(paste0("\nFeature type: ", feature_type, "\n\n"))


# Read the FGSEA result
if(feature_type %in% c("efficacy", "safety", "combinedEfficacySafety")){
  fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
  fgsea_result <- fgsea_result[[feature_type]]
}



# Read the extended drug target data
drugCombs_targets <- readRDS(paste0("InputFiles/Drug_targets/Drug_targets_extended_", disease, ".rds"))
drugCombs_targets$combination_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")


# Select the column to be plotted
switch(drug_category_type,
       "pk" = {plot_col <- "phamk"},
       "adr" = {plot_col <- paste0("ADR_", disease)})


# Plot
plot_data <- fgsea_result


# Find the features with top variance
feature_vars <- apply(plot_data, 1 ,var)
feature_vars <- sort(feature_vars, decreasing = TRUE)
feature_vars <- feature_vars[1:15]

plot_data <- plot_data[row.names(plot_data) %in% names(feature_vars), ]
plot_data <- as.data.frame(t(plot_data))
plot_data <- rownames_to_column(plot_data, "combination_name")


plot_data$category <- drugCombs_targets[,c(plot_col), drop = TRUE][match(plot_data$combination_name, drugCombs_targets$combination_name)]

plot_data <- pivot_longer(data = plot_data, 
                          cols = colnames(plot_data)[!colnames(plot_data) %in% c("combination_name", "category")], 
                          cols_vary = "fastest", 
                          names_to = "feature", 
                          values_to = "value")

plot_data <- plot_data %>% 
  group_by(feature, category) %>% 
  summarise(mean_value = mean(value), 
            se = sd(value) / sqrt(n()), 
            .groups = "drop")

plot_data <- plot_data[!is.na(plot_data$category), ]

plot_data$feature <- str_wrap(plot_data$feature, width = 25)


if(!dir.exists("OutputFiles/Plots/mean_NES_barPlot/")){
  dir.create("OutputFiles/Plots/mean_NES_barPlot/", recursive = TRUE)
}
tiff(paste0("OutputFiles/Plots/mean_NES_barPlot/meanNES_", disease, "_", drug_target_type, "_", feature_type, "_", drug_category_type, ".tiff"),
     width = 30, height = 20,
     units = "cm", compression = "lzw", res = 1200)


ggplot(plot_data, aes(x = category, y = mean_value, fill = category)) +
  geom_bar(width = 0.5, lwd = 0.1, stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_value - se, 
                    ymax = mean_value + se), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  facet_wrap(~ feature) +
  labs(y = "Mean Value +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 6, margin = margin(1,1,1,1)),
        text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 5),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))

dev.off()



print(warnings())