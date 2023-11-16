set.seed(5081)





# Script to plot the NES of top (by variance) features from FGSEA 



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
cat(paste0("\nFeature type: ", feature_type))



plot_list <- list()


for(drug_target_type in c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR","all")){
  for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
    
    
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
    
    
    plot_data <- fgsea_result
    
    
    # Find the features to plot
    if(feature_type %in% c("combinedEfficacySafety")){
      feature_vars <- apply(plot_data, 1, var)
      feature_vars <- sort(feature_vars, decreasing = TRUE)
      efficacy_feature_select <- names(feature_vars[grep("^\\[DISEASE\\]", names(feature_vars))][1])
      safety_feature_select <- names(feature_vars[grep("^\\[ADR\\]", names(feature_vars))][1])
      
      plot_data <- plot_data[row.names(plot_data) %in% c(efficacy_feature_select, safety_feature_select), ,  drop = FALSE]
      
      x_axis_label = str_wrap(efficacy_feature_select, 30)
      y_axis_label = str_wrap(safety_feature_select, 30)
    }

    if(feature_type %in% c("efficacy", "safety", "kegg", "smpdbDrugMet", "smpdbDrugAct", "misc")){
      feature_vars <- apply(plot_data, 1, var)
      feature_vars <- sort(feature_vars, decreasing = TRUE)
      feature_vars <- feature_vars[1:2]
      
      plot_data <- plot_data[row.names(plot_data) %in% names(feature_vars), ,  drop = FALSE]
      
      x_axis_label = str_wrap(names(feature_vars)[1], 30)
      y_axis_label = str_wrap(names(feature_vars)[2], 30)
    }
    

    plot_data <- as.data.frame(t(plot_data))
    colnames(plot_data) <- c("F1", "F2")
    
    plot_data <- merge(plot_data, drugCombs_cat, by.x = 0, by.y = "comb_name")
    
    plot_list[[drug_target_type]][[disease]] <- ggplot() +
      geom_point(data = plot_data, 
                 mapping = aes(x = F1, y = F2, color = class_EffAdv),  
                 size = 0.5,
                 shape = 3) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
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
           x = x_axis_label,
           y = y_axis_label,
           color = "Class")
  }
}


if(!dir.exists("OutputFiles/Plots/top_feature_NES_scatter/")){
  dir.create("OutputFiles/Plots/top_feature_NES_scatter/", recursive = TRUE)
}


tiff(paste0("OutputFiles/Plots/top_feature_NES_scatter/plot_NES_", feature_type, "_topFeature_byVar.tiff"),
     width = 30, height = 25,
     units = "cm", compression = "lzw", res = 1200)

ggarrange(plotlist = unlist(plot_list, recursive = FALSE), 
          ncol = 6, nrow = 7, 
          common.legend = TRUE, legend = "bottom")

dev.off()

print(warnings())