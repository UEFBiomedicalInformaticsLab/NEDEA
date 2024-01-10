set.seed(5081)



# Script to plot number of targets for the drug combinations that show enrichment for any features 



# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)
library(ggpubr)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Get arguments
option_list = list(
  make_option(c("--feature_type"), type = "character", default = NULL,
              help = "The feature type to use for plotting. Possible values: efficacy, safety, combinedEfficacySafety, kegg, smpdbDrugMet, smpdbDrugAct, misc. Default: NULL", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$feature_type)){
  print_help(opt_parser)
  stop("--feature_type argument needed", call. = FALSE)
}

if(!opt$feature_type %in% c("efficacy", "safety", "combinedEfficacySafety",
                            "kegg", "smpdbDrugMet", "smpdbDrugAct", "misc")){
  print_help(opt_parser)
  stop("--feature_type should be: efficacy, safety, combinedEfficacySafety, kegg, smpdbDrugMet, smpdbDrugAct, misc", call. = FALSE)
}



# Define global options for this script
feature_type <- opt$feature_type

cat("\n\nUsing the following parameters: ")
cat(paste0("\nFeature type: ", feature_type, " \n"))

plot_list <- list()

for(drug_target_type in c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR")){
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    
    # Select the column containing the drug target based on th user input
    switch(drug_target_type,
           "known" = { drug_target_col <- c("drugTarget_count") },
           "PS" = { drug_target_col <- c("ext_PS_tar_cnt") },
           "SIGNOR" = { drug_target_col <- c("ext_SIGNOR_tar_cnt") },
           "NPA" = { drug_target_col <- c("ext_NPA_tar_cnt") },
           "RI" = { drug_target_col <- c("ext_RI_tar_cnt") },
           "KEGG" = { drug_target_col <- c("ext_KEGG_tar_cnt") },
           "all" = { drug_target_col <- c("drugTarget_count", "ext_PS_tar_cnt", 
                                          "ext_SIGNOR_tar_cnt", "ext_NPA_tar_cnt", 
                                          "ext_RI_tar_cnt", "ext_KEGG_tar_cnt") })
    
    # Read the drug targets
    drugCombs_targets <- readRDS(paste0("InputFiles/Drug_combination_targets/drugCombs_targets_extended_", disease, ".rds"))
    drugCombs_targets$comb_name <- paste(drugCombs_targets$Drug1_DrugBank_id, drugCombs_targets$Drug2_DrugBank_id, sep = "_")
    drugCombs_targets <- drugCombs_targets[, c("comb_name", drug_target_col), drop = FALSE]
    
    
    # Read the drug combination category
    drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
    drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
    drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), c("comb_name", "class_EffAdv")]
    
    
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
    
    
    # Keep the FGSEA result for only the drug combinations with valid category
    fgsea_result <- fgsea_result[, colnames(fgsea_result) %in% drugCombs_cat$comb_name, drop = FALSE]
    
    
    
    # Extract the combinations that has significant NES for atleast one gene set
    drugs_with_NES <- colnames(fgsea_result[, colSums(fgsea_result != 0) > 0, drop = FALSE])
    drugs_without_NES <- colnames(fgsea_result[, colSums(fgsea_result != 0) == 0, drop = FALSE])
    
    
    # Extract the targets for the drug combinations
    targets_drugs_with_NES <- drugCombs_targets[drugCombs_targets$comb_name %in% drugs_with_NES, ]
    targets_drugs_without_NES <- drugCombs_targets[drugCombs_targets$comb_name %in% drugs_without_NES, ]
    
    
    
    # plot the counts
    plot_data_1 <- targets_drugs_with_NES %>% select("comb_name", contains("_count"), contains("_cnt"))
    plot_data_1$NES <- "yes"
    plot_data_2 <- targets_drugs_without_NES %>% select("comb_name", contains("_count"), contains("_cnt"))
    plot_data_2$NES <- "no"
    
    
    plot_data <- rbind(plot_data_1, plot_data_2)
    plot_data <- plot_data %>% pivot_longer(cols = c(contains("count") | contains("_cnt")), 
                                            names_to = "target_type", 
                                            values_to = "count")
    
    
    plot_data$target_type <- gsub("drugTarget_count", "known", plot_data$target_type)
    plot_data$target_type <- gsub("^ext_|_tar_cnt|", "", plot_data$target_type)
    
    plot_list[[drug_target_type]][[disease]] <- ggplot(plot_data, aes(x = target_type, y = count, fill = NES)) +
      geom_boxplot(width = 0.5, lwd = 0.1, 
                   outlier.shape = 3, 
                   outlier.size = 0.5) +
      geom_pwc(aes(group = NES),
               method = "wilcox_test",
               label = "p.signif", 
               hide.ns = TRUE, 
               vjust = 0.5, 
               color = "blue") +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            text = element_text(size = 4),
            plot.title = element_text(hjust = 0.5, size = 4),
            axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            legend.position = "bottom",
            legend.key = element_blank(),
            legend.key.size = unit(0.5, 'cm'),
            legend.title = element_text(size = 3),
            legend.text = element_text(size = 3),
            legend.margin = margin(1,1,1,1),
            legend.box.spacing = unit(0.1, 'cm'),
            legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
      labs(title = paste0("Disease: ", disease,
                          "\nTarget: ", drug_target_type),
           x = "Target type",
           y = "Number of drug targets",
           color = "NES")
    
  }}




if(!dir.exists("OutputFiles/Plots/target_count_enriching_combs/")){
  dir.create("OutputFiles/Plots/target_count_enriching_combs/", recursive = TRUE)
}


tiff(paste0("OutputFiles/Plots/target_count_enriching_combs/target_count_enriching_combs_", feature_type, ".tiff"),
     width = 30, height = 30,
     units = "cm", compression = "lzw", res = 1200)


ggarrange(plotlist = unlist(plot_list, recursive = FALSE), 
          ncol = 6, nrow = 6, 
          common.legend = TRUE, legend = "bottom")


dev.off()


print(warnings())