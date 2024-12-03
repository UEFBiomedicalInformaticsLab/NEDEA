set.seed(5081)


# Script to plot the mean network proximity for each drug category


# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
library(ggpubr)

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known",
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character"),
  make_option(c("--feature_type"), type = "character", default = NULL,
              help = "The feature type to use for plotting. Possible values: efficacy, safety, combinedEfficacySafety, kegg, smpdbDrugMet, smpdbDrugAct, misc. Default: NULL", metavar = "character"),
  make_option(c("--drugComb_category_type"), type = "character", default = NULL,
              help = "The drug category type that will be plotted. Possible values: SL (synergy level), SS (class synergy score), ADR (disease specific ADR status), EA (effective vs adverse). Default: NULL", metavar = "character"),
  make_option(c("--proximity_type"), type = "character", default = NULL,
              help = "The network proximity type to be used. Possible values:  Centre, Closest, Kernel, Separation, Shortest. Default: Separation", metavar = "character"),
  make_option(c("--top9varying"), type = "logical", default = TRUE,
              help = "Plots the top nine features with most variance. Set FALSE to plot all features above variance 0.01. Default: TRUE", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$disease %in% c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  print_help(opt_parser)
  stop("--disease must be one of the following: BreastCancer, KidneyCancer, LungCancer, OvaryCancer, ProstateCancer, SkinCancer", call. = FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
}

if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call. = FALSE)
}

if(!opt$feature_type %in% c("efficacy", "safety", "combinedEfficacySafety", 
                            "kegg", "smpdbDrugMet", "smpdbDrugAct", "misc")){
  print_help(opt_parser)
  stop("--feature_type should be: efficacy, safety, combinedEfficacySafety, kegg, smpdbDrugMet, smpdbDrugAct, misc", call. = FALSE)
}

if(is.null(opt$drugComb_category_type)){
  print_help(opt_parser)
  stop("--drugComb_category_type argument needed", call.=FALSE)
}

if(!opt$drugComb_category_type %in% c("SL", "SS", "ADR", "EA")){
  print_help(opt_parser)
  stop("--drugComb_category_type should be: SL, SS, ADR, EA", call. = FALSE)
}

if(!opt$proximity_type %in% c("Centre", "Closest", "Kernel", "Separation", "Shortest")){
  print_help(opt_parser)
  stop("--proximity_type should be: Centre, Closest, Kernel, Separation, Shortest", call.=FALSE)
}



# Define global options for this script
disease <- opt$disease
drug_target_type <- opt$drug_target_type
feature_type <- opt$feature_type
drugComb_category_type <- opt$drugComb_category_type
proximity_type <- opt$proximity_type
top9varying <- opt$top9varying


#####


# Read the network proximity data result
if(feature_type %in% c("efficacy", "safety", "combinedEfficacySafety")){
  proximity_result <- readRDS(paste0("OutputFiles/Network_proximity_results/netProx_EfficacySafety_", disease, "_", drug_target_type, ".rds"))
  proximity_result <- proximity_result[[feature_type]][[proximity_type]]
  proximity_result <-  as.matrix(column_to_rownames(proximity_result, "lib_name"))
}

if(feature_type %in% c("kegg", "smpdbDrugMet", "smpdbDrugAct")){
  proximity_result <- readRDS(paste0("OutputFiles/Network_proximity_results/netProx_Pathway_", disease, "_", drug_target_type, ".rds"))
  proximity_result <- proximity_result[[feature_type]][[proximity_type]]
  proximity_result <-  as.matrix(column_to_rownames(proximity_result, "lib_name"))
}

if(feature_type %in% c("misc")){
  proximity_result <- readRDS(paste0("OutputFiles/Network_proximity_results/netProx_Miscellaneous_", disease, "_", drug_target_type, ".rds"))
  proximity_result <- proximity_result[[feature_type]][[proximity_type]]
  proximity_result <- as.matrix(column_to_rownames(proximity_result, "lib_name"))
}


#####


# Read the drug combination category
switch(drugComb_category_type,
       "SL" = {plot_col <- "Syn_level"},
       "SS" = {plot_col <- "class_synergyScore"},
       "ADR" = {plot_col <- "ADR_status"},
       "EA" = {plot_col <- "class_EffAdv"})


if(drugComb_category_type %in% c("SL", "SS")){
  drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]
}

if(drugComb_category_type %in% c("ADR")){
  drugCombs_cat <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]
}

if(drugComb_category_type %in% c("EA")){
  drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
  drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
  drugCombs_cat <- drugCombs_cat[, c("comb_name", plot_col)]
}
#####


# Get the plot data
plot_data <- proximity_result

# Find the features with top variance (executed only if top9varying is TRUE)
if(isTRUE(top9varying)){
  feature_vars <- apply(plot_data, 1, var)
  feature_vars <- sort(feature_vars, decreasing = TRUE)
  feature_vars <- feature_vars[1:9]
  plot_data <- plot_data[row.names(plot_data) %in% names(feature_vars), ,  drop = FALSE]
}

# Find the features with variance above 0.01 (executed only if top9varying is FALSE)
if(isFALSE(top9varying)){
  feature_vars <- apply(plot_data, 1, var)
  feature_vars <- feature_vars[feature_vars > 0.01]
  plot_data <- plot_data[row.names(plot_data) %in% names(feature_vars), , drop = FALSE]
  
  if(nrow(plot_data) == 0){
    write("No features with variance > 0.01. No output generated.", 
          paste0("OutputFiles/Plots/mean_NES_barPlot_all/meanNES_", disease, "_", drug_target_type, "_", feature_type, "_", drugComb_category_type, ".txt"))
    quit()
  }
}


# Add the drug combination categories to the plot data
plot_data <- as.data.frame(t(plot_data))
plot_data <- rownames_to_column(plot_data, "comb_name")

plot_data$category <- drugCombs_cat[,c(plot_col), drop = TRUE][match(plot_data$comb_name, drugCombs_cat$comb_name)]
plot_data <- plot_data %>% mutate_at("category", ~replace_na(., "unknown"))
plot_data$category <- as.factor(plot_data$category)

plot_data <- pivot_longer(data = plot_data, 
                          cols = colnames(plot_data)[!colnames(plot_data) %in% c("comb_name", "category")], 
                          cols_vary = "fastest", 
                          names_to = "feature", 
                          values_to = "value")


# Plot the proximity as box plot
if(isTRUE(top9varying)){
  if(!dir.exists("OutputFiles/Plots/Proximity_boxPlot_top9Var/")){
    dir.create("OutputFiles/Plots/Proximity_boxPlot_top9Var/", recursive = TRUE)
  }
  tiff(paste0("OutputFiles/Plots/Proximity_boxPlot_top9Var/Proximity", proximity_type, "_", disease, "_", drug_target_type, "_", feature_type, "_", drugComb_category_type, ".tiff"),
       width = 30,
       height = 15,
       units = "cm", compression = "lzw", res = 1200)
}

if(isFALSE(top9varying)){
  if(!dir.exists("OutputFiles/Plots/Proximity_boxPlot_all/")){
    dir.create("OutputFiles/Plots/Proximity_boxPlot_all/", recursive = TRUE)
  }
  tiff(paste0("OutputFiles/Plots/Proximity_boxPlot_all/Proximity", proximity_type, "_", disease, "_", drug_target_type, "_", feature_type, "_", drugComb_category_type, ".tiff"),
       width = 30, 
       height = length(unique(plot_data$feature))/6 + 10,
       units = "cm", compression = "lzw", res = 1200)
}


ggplot(plot_data, aes(x = category, y = value)) +
  geom_boxplot(width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3, 
               outlier.size = 1) +
  facet_wrap(~ feature, scales = "free_y", labeller = label_wrap_gen(width = 30, 
                                                                     multi_line = TRUE)) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "red") +
  geom_pwc(method = "wilcox_test",
           label = "p.signif", 
           hide.ns = TRUE, 
           vjust = 0.5, 
           color = "blue") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 6, margin = margin(1,1,1,1)),
        text = element_text(size = 8), 
        plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
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


#####


# Calculate the mean and SD for each category for each feature
plot_data <- plot_data %>% 
  group_by(feature, category) %>% 
  summarise(mean_value = mean(value), 
            se = sd(value) / sqrt(n()), 
            .groups = "drop")


# Plot 
if(isTRUE(top9varying)){
  if(!dir.exists("OutputFiles/Plots/mean_Proximity_barPlot_top9Var/")){
    dir.create("OutputFiles/Plots/mean_Proximity_barPlot_top9Var/", recursive = TRUE)
  }
  tiff(paste0("OutputFiles/Plots/mean_Proximity_barPlot_top9Var/meanProximity", proximity_type, "_", disease, "_", drug_target_type, "_", feature_type, "_", drugComb_category_type, ".tiff"),
       width = 30,
       height = 15,
       units = "cm", compression = "lzw", res = 1200)
}

if(isFALSE(top9varying)){
  if(!dir.exists("OutputFiles/Plots/mean_Proximity_barPlot_all/")){
    dir.create("OutputFiles/Plots/mean_Proximity_barPlot_all/", recursive = TRUE)
  }
  tiff(paste0("OutputFiles/Plots/mean_Proximity_barPlot_all/meanProximity", proximity_type, "_", disease, "_", drug_target_type, "_", feature_type, "_", drugComb_category_type, ".tiff"),
       width = 30, 
       height = length(unique(plot_data$feature))/6 + 10,
       units = "cm", compression = "lzw", res = 1200)
}


ggplot(plot_data, aes(x = category, y = mean_value)) +
  geom_bar(width = 0.5, lwd = 0.1, stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_value - se, 
                    ymax = mean_value + se), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  facet_wrap(~ feature, labeller = label_wrap_gen(width = 30, 
                                                  multi_line = TRUE)) +
  labs(title = paste0("Disease: ", disease, 
                      "; Drug target: ", drug_target_type,
                      "; Feature: ", feature_type,
                      "; Drug comb. cat.: ", drugComb_category_type,
                      "; top9varying: ", top9varying),
       y = "Mean Value +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 6, margin = margin(1,1,1,1)),
        text = element_text(size = 8), 
        plot.title = element_text(size = 8, hjust = 0.5, face = "plain"),
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


#####


print(warnings())