set.seed(5081)



# Script to plot

# Load libraries
library(tidyverse)



disease <- "BreastCancer"
drug_target_type <- "all"




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
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
}



# Define global options for this script 
disease <- opt$disease
drug_target_type <- opt$drug_target_type


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))


# Read the FGSEA result
fgsea_result <- readRDS(paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))


# Read the extended drug target data
drugCombs_targets <- readRDS(paste0("InputFiles/Drug_targets/Drug_targets_extended_", disease, ".rds"))
drugCombs_targets$combination_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")



# Plot the mean NES for efficacy gene sets w.r.t. pharmakokinetics info
plot_data <- fgsea_result$efficacy

# plot_data <- plot_data[-which(apply(plot_data, 1 ,var) < 0.01),] # Filter the features with low variance
plot_data <- plot_data[which(apply(plot_data, 1 ,var) >= 0.01),] # Filter the features with low variance

plot_data <- as.data.frame(t(plot_data))
plot_data <- rownames_to_column(plot_data, "combination_name")

plot_data$category <- drugCombs_targets$phamk[match(plot_data$combination_name, drugCombs_targets$combination_name)]
plot_data <- pivot_longer(data = plot_data, 
                          cols = colnames(plot_data)[grep("\\[DISEASE\\]", colnames(plot_data))], 
                          cols_vary = "fastest", 
                          names_to = "feature", 
                          values_to = "value")


plot_data <- plot_data %>% 
  group_by(feature, category) %>% 
  summarise(mean_value = mean(value), 
            se = sd(value) / sqrt(n()), 
            .groups = "drop")


plot_data <- plot_data[which(!is.na(plot_data$category)), ]


plot_data$feature <- str_wrap(plot_data$feature, width = 25)


if(!dir.exists("OutputFiles/Plots/meanNES_barPlot/")){
  dir.create("OutputFiles/Plots/meanNES_barPlot/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/meanNES_barPlot/meanNES_diseaseVSpharmk_", disease, "_", drug_target_type, ".tiff"),
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




# Plot the mean NES for safety gene sets w.r.t. ADR association info
plot_data <- fgsea_result$safety

plot_data <- plot_data[which(apply(plot_data, 1 ,var) >= 0.01),] # Filter the features with low variance

plot_data <- as.data.frame(t(plot_data))
plot_data <- rownames_to_column(plot_data, "combination_name")

plot_data$category <- drugCombs_targets[,c(paste0("ADR_", disease)), drop = TRUE][match(plot_data$combination_name, drugCombs_targets$combination_name)]
plot_data <- pivot_longer(data = plot_data, 
                          cols = colnames(plot_data)[grep("\\[ADR\\]", colnames(plot_data))], 
                          cols_vary = "fastest", 
                          names_to = "feature", 
                          values_to = "value")


plot_data <- plot_data %>% 
  group_by(feature, category) %>% 
  summarise(mean_value = mean(value), 
            se = sd(value) / sqrt(n()), 
            .groups = "drop")


plot_data <- plot_data[which(!is.na(plot_data$category)), ]


plot_data$feature <- str_wrap(plot_data$feature, width = 25)


if(!dir.exists("OutputFiles/Plots/meanNES_barPlot/")){
  dir.create("OutputFiles/Plots/meanNES_barPlot/", recursive = TRUE)
}  
tiff(paste0("OutputFiles/Plots/meanNES_barPlot/meanNES_safetyVSadr_", disease, "_", drug_target_type, ".tiff"),
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