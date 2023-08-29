set.seed(5081)



# Script to plot the distribution of synergy scores from FIMM



# Load libraries
library(httr)
library(jsonlite)
httr::set_config(config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
source("Scripts/Functions/Functions_data_manipulation.R")

library(ggfortify)
library(ggpubr)
library(tidyverse)
library(arules)

source("Explore_FIMM/Process_DrugBank_DDI.R")




# Download drug combinations from FIMM DrugComb 
if(!dir.exists("Databases/FimmDrugComb/"))dir.create("Databases/FimmDrugComb/")
if(!file.exists("Databases/FimmDrugComb/summary_v_1_5.csv")){
  download.file(url = "https://drugcomb.fimm.fi/jing/summary_v_1_5.csv",
                destfile = "Databases/FimmDrugComb/summary_v_1_5.csv", method = "wget")
}


# Read the DrugComb data
FimmDrugCombs <- read.csv("Databases/FimmDrugComb/summary_v_1_5.csv", header = TRUE)
FimmDrugCombs <- FimmDrugCombs[FimmDrugCombs$drug_row != "NULL", ]
FimmDrugCombs <- FimmDrugCombs[FimmDrugCombs$drug_col != "NULL", ]


# Get cell line information
FimmDrugComb_cellLine <- GET("https://api.drugcomb.org/cell_lines")
FimmDrugComb_cellLine <- fromJSON(rawToChar(FimmDrugComb_cellLine$content))


# Download disease IDs from NCI Thesaurus (NCIt) 
if(!dir.exists("Databases/NCIt/"))dir.create("Databases/NCIt/")
if(!file.exists("Databases/NCIt/Thesaurus.txt")){
  download.file(url = "https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/22.12d_Release/Thesaurus_22.12d.FLAT.zip",
                destfile = "Databases/NCIt/Thesaurus_22.12d.FLAT.zip", method = "wget")
  unzip("Databases/NCIt/Thesaurus_22.12d.FLAT.zip", exdir = "Databases/NCIt/", file = "Thesaurus.txt")
}
NCIthesaurus <- read.table("Databases/NCIt/Thesaurus.txt", sep = "\t", header = FALSE, comment.char = "", fill = TRUE, quote = "")
colnames(NCIthesaurus) <- c("code", "concept IRI", "parents", "synonyms", "definition", "display name", "concept status", "semantic type")


# Filter drug combinations tested on cancer related cell lines
NCIthesaurus <- NCIthesaurus[NCIthesaurus$code %in% FimmDrugComb_cellLine$disease_id, ]
NCIthesaurus <- NCIthesaurus[grep("cancer|carcinoma|sarcoma|lymphoma|leukemia|melanoma", NCIthesaurus$synonyms, ignore.case = TRUE), ]

FimmDrugComb_cellLine <- FimmDrugComb_cellLine[FimmDrugComb_cellLine$disease_id %in% NCIthesaurus$code, ]

FimmDrugCombs <- FimmDrugCombs[FimmDrugCombs$cell_line_name %in% FimmDrugComb_cellLine$name, ]



# Map drugs to DrugBank drug ID
FimmDrugComb_drugs <- GET("https://api.drugcomb.org/drugs")
FimmDrugComb_drugs <- fromJSON(rawToChar(FimmDrugComb_drugs$content))


FimmDrugCombs$Drug1_DrugBank_drug_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugCombs$drug_row, FimmDrugComb_drugs$dname)]
FimmDrugCombs$Drug2_DrugBank_drug_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugCombs$drug_col, FimmDrugComb_drugs$dname)]


FimmDrugCombs <- FimmDrugCombs[(FimmDrugCombs$Drug1_DrugBank_drug_id != "NA"),]
FimmDrugCombs <- FimmDrugCombs[(FimmDrugCombs$Drug2_DrugBank_drug_id != "NA"),]

# FimmDrugCombs$log_synergy_hsa <- log(FimmDrugCombs$synergy_hsa + 1)



FimmDrugCombs <- FimmDrugCombs[FimmDrugCombs$tissue_name %in% c("breast", "kidney", "lung", "ovary", "prostate", "skin"), ]

# Get DDI for the drug combinations
DrugBank_drugInteractions <- read.csv("Databases/DrugBank/drug_drug_interactions.csv")

top_freq = find_top_k_frequent_sequences(DrugBank_drugInteractions$description, n = 4, 100)
risk_sev_res = find_top_k_sequences_with_keyword(DrugBank_drugInteractions$description, keyword = "The risk or severity", n = 6, k = 50)


FimmDrugCombs$DDI <- DrugBank_drugInteractions$description[match(paste(FimmDrugCombs$Drug1_DrugBank_drug_id, FimmDrugCombs$Drug2_DrugBank_drug_id), 
                                                                 paste(DrugBank_drugInteractions$drugbank.id, DrugBank_drugInteractions$parent_key))]
FimmDrugCombs$DDI <- DrugBank_drugInteractions$description[match(paste(FimmDrugCombs$Drug1_DrugBank_drug_id, FimmDrugCombs$Drug2_DrugBank_drug_id), 
                                                                 paste(DrugBank_drugInteractions$parent_key, DrugBank_drugInteractions$drugbank.id))]



FimmDrugCombs$test_eff_inc <- grepl("The therapeutic efficacy of", FimmDrugCombs$DDI) & grepl("can be increased when used in combination", FimmDrugCombs$DDI)
FimmDrugCombs$test_eff_dec = grepl("The therapeutic efficacy of", FimmDrugCombs$DDI) & grepl("can be decreased when used in combination", FimmDrugCombs$DDI)
FimmDrugCombs$test_anticinc = grepl("may increase the anticoagulant activities", FimmDrugCombs$DDI)
FimmDrugCombs$test_anticdec = grepl("may decrease the anticoagulant activities", FimmDrugCombs$DDI)
FimmDrugCombs$test_antihinc = grepl("may increase the antihypertensive activities", FimmDrugCombs$DDI)
FimmDrugCombs$test_antihdec = grepl("may decrease the antihypertensive activities", FimmDrugCombs$DDI)
FimmDrugCombs$test_nsd = grepl("(CNS depressant)", FimmDrugCombs$DDI)
FimmDrugCombs$test_scinc = grepl("The serum concentration of", FimmDrugCombs$DDI) & grepl("can be increased", FimmDrugCombs$DDI)
FimmDrugCombs$test_scdec = grepl("The serum concentration of", FimmDrugCombs$DDI) & grepl("can be decreased", FimmDrugCombs$DDI)
FimmDrugCombs$test_meinc = grepl("The metabolism of", FimmDrugCombs$DDI) & grepl("can be increased", FimmDrugCombs$DDI)
FimmDrugCombs$test_medec = grepl("The metabolism of", FimmDrugCombs$DDI) & grepl("can be decreased", FimmDrugCombs$DDI)
FimmDrugCombs$test_absde = grepl("a decrease in the absorption ", FimmDrugCombs$DDI)




list_tox_test = list()
for(i in 1:length(risk_sev_res$ngrams)){
  list_tox_test[[length(list_tox_test) + 1]] = grepl(risk_sev_res$ngrams[i], FimmDrugCombs$DDI)
}
names(list_tox_test) = risk_sev_res$ngrams

tmp <- data.frame(list_tox_test, check.names = FALSE)

FimmDrugCombs$test_risk_sev <- apply(tmp, 1, any)





plot_cols <- colnames(FimmDrugCombs)[grep("^test_", colnames(FimmDrugCombs))]


# Plot the overall distribution of synergy scores
for(plot_col_select in plot_cols){
  
  plot_list <- list()
  for(tissue in unique(FimmDrugCombs$tissue_name)){
    
    plot_data <- FimmDrugCombs[FimmDrugCombs$tissue_name == tissue, c("Drug1_DrugBank_drug_id", "drug_row", 
                                                                      "Drug2_DrugBank_drug_id", "drug_col", 
                                                                      "cell_line_name", "tissue_name", "synergy_hsa", plot_col_select)]
    colnames(plot_data) <- c("Drug1_DrugBank_drug_id", "drug_row", 
                             "Drug2_DrugBank_drug_id", "drug_col", 
                             "cell_line_name", "tissue_name", "synergy", "category")
    
    plot_list[[tissue]] <- ggplot() +
      geom_boxplot(data = plot_data, 
                   mapping = aes(x = cell_line_name, 
                                 y = synergy, 
                                 fill = category), 
                   outlier.shape = NA,
                   position = position_dodge(preserve = "single")) + 
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            strip.background = element_rect(color = "black", linewidth = 0.25,),
            strip.text = element_text(margin = margin(1,1,1,1)),
            text = element_text(size = 4), 
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            legend.position = "bottom",
            legend.key = element_blank(),
            legend.key.size = unit(0.5, 'cm'),
            legend.text = element_text(size = 4),
            legend.margin = margin(1,1,1,1),
            legend.box.spacing = unit(0.1, 'cm'),
            legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
      ggtitle(tissue) +
      labs(fill = plot_col_select) 
    
  }
  
  if(!dir.exists("Explore_FIMM/Overall_dist/")){
    dir.create("Explore_FIMM/Overall_dist/", recursive = TRUE)
  }  
  tiff(paste0("Explore_FIMM/Overall_dist/Overall_dist__", plot_col_select, ".tiff"),
       width = 40, height = 40,
       units = "cm", compression = "lzw", res = 1200)
  
  plot <- ggarrange(plotlist = plot_list,
                    ncol = 1,
                    common.legend = TRUE, legend = "bottom")
  print(plot)
  
  dev.off()
  
}








# Plot the cell line specific discrtized distribution
for(plot_col_select in plot_cols){
  
  for(tissue in unique(FimmDrugCombs$tissue_name)){
    
    
    FimmDrugCombs_select <- FimmDrugCombs[FimmDrugCombs$tissue_name == tissue, c("Drug1_DrugBank_drug_id", "drug_row", 
                                                                                 "Drug2_DrugBank_drug_id", "drug_col", 
                                                                                 "cell_line_name", "tissue_name", "synergy_hsa", plot_col_select)]
    colnames(FimmDrugCombs_select) <- c("Drug1_DrugBank_drug_id", "drug_row", 
                                        "Drug2_DrugBank_drug_id", "drug_col", 
                                        "cell_line_name", "tissue_name", "synergy", "category")
    
    
    plot_list <- list()
    for(cell_line_select in unique(FimmDrugCombs_select$cell_line_name)){
      
      
      plot_data <- FimmDrugCombs_select[FimmDrugCombs_select$cell_line_name == cell_line_select, ]
      
      if(length(unique((plot_data$synergy))) > 2){
        plot_data$synergy_level <- discretize(plot_data$synergy, breaks = 3, method = "interval", labels = c("Low", "Mid", "High"), infinity = TRUE)
      }else{
        plot_data$synergy_level <- NA
      }
      
      
      plot_data <- plot_data[plot_data$category == TRUE,]
      
      
      plot_list[[cell_line_select]] <- ggplot() +
        geom_boxplot(data = plot_data, 
                     mapping = aes(x = synergy_level, 
                                   y = synergy, 
                                   fill = synergy_level), 
                     outlier.shape = NA,
                     position = position_dodge(preserve = "single")) + 
        scale_colour_manual(values = c("Low" = "red", "Mid" = "blue", "High" = "green")) +
        theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
              panel.grid = element_blank(),
              panel.spacing = unit(0.1, "cm"),
              strip.background = element_rect(color = "black", linewidth = 0.25,),
              strip.text = element_text(margin = margin(1,1,1,1)),
              text = element_text(size = 4), 
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
              axis.ticks = element_line(colour = "black", linewidth = 0.2),
              legend.position = "bottom",
              legend.key = element_blank(),
              legend.key.size = unit(0.5, 'cm'),
              legend.text = element_text(size = 4),
              legend.margin = margin(1,1,1,1),
              legend.box.spacing = unit(0.1, 'cm'),
              legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
        ggtitle(cell_line_select) 
      
    }
    
    if(!dir.exists("Explore_FIMM/Cell_line_discretized/")){
      dir.create("Explore_FIMM/Cell_line_discretized/", recursive = TRUE)
    }  
    tiff(paste0("Explore_FIMM/Cell_line_discretized/Cell_line_discretized__", plot_col_select, "__", tissue, ".tiff"),
         width = 40, height = 40,
         units = "cm", compression = "lzw", res = 1200)
    
    plot <- ggarrange(plotlist = plot_list,
                      common.legend = TRUE, legend = "bottom")
    print(plot)
    
    dev.off()
    
    
  }
}




# Plot the tissue specific discrtized distribution

for(plot_col_select in plot_cols){
  
  for(tissue in unique(FimmDrugCombs$tissue_name)){
    
    FimmDrugCombs_select <- FimmDrugCombs[FimmDrugCombs$tissue_name == tissue, c("Drug1_DrugBank_drug_id", "drug_row", 
                                                                                 "Drug2_DrugBank_drug_id", "drug_col", 
                                                                                 "cell_line_name", "tissue_name", "synergy_hsa", plot_col_select)]
    colnames(FimmDrugCombs_select) <- c("Drug1_DrugBank_drug_id", "drug_row", 
                                        "Drug2_DrugBank_drug_id", "drug_col", 
                                        "cell_line_name", "tissue_name", "synergy", "category")
    
    
    plot_data <- FimmDrugCombs_select
    
    if(length(unique((plot_data$synergy))) > 2){
      plot_data$synergy_level <- discretize(plot_data$synergy, breaks = 3, method = "interval", labels = c("Low", "Mid", "High"), infinity = TRUE)
    }else{
      plot_data$synergy_level <- NA
    }
    
    
    plot <- ggplot() +
      geom_boxplot(data = plot_data, 
                   mapping = aes(x = synergy_level, 
                                 y = synergy, 
                                 fill = category), 
                   outlier.shape = NA,
                   position = position_dodge(preserve = "single")) + 
      # scale_colour_manual(values = c("Low" = "red", "Mid" = "blue", "High" = "green")) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
            panel.grid = element_blank(),
            panel.spacing = unit(0.1, "cm"),
            strip.background = element_rect(color = "black", linewidth = 0.25,),
            strip.text = element_text(margin = margin(1,1,1,1)),
            text = element_text(size = 4), 
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.ticks = element_line(colour = "black", linewidth = 0.2),
            legend.position = "bottom",
            legend.key = element_blank(),
            legend.key.size = unit(0.5, 'cm'),
            legend.text = element_text(size = 4),
            legend.margin = margin(1,1,1,1),
            legend.box.spacing = unit(0.1, 'cm'),
            legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
      labs(fill = plot_col_select) +
      facet_wrap(vars(cell_line_name))
    
    if(!dir.exists("Explore_FIMM/Tissue_discretized/")){
      dir.create("Explore_FIMM/Tissue_discretized/", recursive = TRUE)
    }  
    tiff(paste0("Explore_FIMM/Tissue_discretized/Tissue_discretized__", plot_col_select, "__", tissue, ".tiff"),
         width = 40, height = 40,
         units = "cm", compression = "lzw", res = 1200)
    
    print(plot)
    
    dev.off()
    
  }
  
}




res <- data.frame()
for(plot_col_select in plot_cols){
  
  for(tissue in unique(FimmDrugCombs$tissue_name)){
    
    
    FimmDrugCombs_select <- FimmDrugCombs[FimmDrugCombs$tissue_name == tissue, c("Drug1_DrugBank_drug_id", "drug_row", 
                                                                                 "Drug2_DrugBank_drug_id", "drug_col", 
                                                                                 "cell_line_name", "tissue_name", "synergy_hsa", plot_col_select)]
    colnames(FimmDrugCombs_select) <- c("Drug1_DrugBank_drug_id", "drug_row", 
                                        "Drug2_DrugBank_drug_id", "drug_col", 
                                        "cell_line_name", "tissue_name", "synergy", "category")
    
    
    plot_list <- list()
    for(cell_line_select in unique(FimmDrugCombs_select$cell_line_name)){
      
      
      plot_data <- FimmDrugCombs_select[FimmDrugCombs_select$cell_line_name == cell_line_select, ]
      
      res <- rbind(res, data.frame("Tissue" = tissue,
                                   "Cell line" = cell_line_select,
                                   "DDI category" = plot_col_select,
                                   "Total combinations" = nrow(plot_data),
                                   "True for DDI" = sum(plot_data$category),
                                   check.names = FALSE
      ))
      
      
    }
  }
}

write.csv(res, "Explore_FIMM/Summary.csv", row.names = FALSE)