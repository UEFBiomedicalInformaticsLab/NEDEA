set.seed(5081)



# Script to compile the library terms and sizes used for enrichment



# Load libraries
library(tidyverse)
library(openxlsx)

library_size <- list()


# Compile efficacy libraries
diseases <- c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")

for(disease in diseases){
  enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  size <- as.data.frame(lengths(enrichment_lib))
  colnames(size) <- "Size"
  size <- rownames_to_column(size, "Description")
  library_size[[paste0("Efficacy_", disease)]] <- size
}



# Compile safety library
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["Safety"]] <- size



# Compile KEGG library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/CHG_keggPath2Gene_lib.rds"))
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["KEGG"]] <- size



# Compile SMPDB (Drug Metabolsim) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Metabolism`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["SMPDB_DrugMetabolism"]] <- size



# Compile SMPDB (Drug Action) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Action`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["SMPDB_DrugAction"]] <- size



# Compile miscelleneous library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/miscellaneous_gene_lib.rds"))
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "Size"
size <- rownames_to_column(size, "Description")
library_size[["miscellaneous"]] <- size



# Plot the library sizes as violin plot

if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}
tiff("OutputFiles/Plots/Enrichment_library_size.tiff",
     width = 15, height = 5,
     units = "cm", compression = "lzw", res = 1200)

plot_data <- bind_rows(library_size, .id = "Source")

plot_data$Size[plot_data$Size >= 500] <- 500

plot_data$Source <- factor(plot_data$Source, 
                           levels = c("Efficacy_LungCancer", "Efficacy_BreastCancer", "Efficacy_ProstateCancer", 
                                      "Efficacy_OvaryCancer", "Efficacy_KidneyCancer", "Efficacy_SkinCancer", 
                                      "Safety", "KEGG", "SMPDB_DrugMetabolism", "SMPDB_DrugAction", "miscellaneous"))

plot_data_summary <- plot_data %>% group_by(Source) %>% summarise(n = n())

ggplot() +
  geom_violin(data = plot_data, 
              aes(x = Source, y = Size), 
              fill = "grey", lwd = 0.1) +
  geom_text(data = plot_data_summary,
            aes(x = Source, y = Inf, label = n),
            color = "red", size = 2, vjust = 2) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.2)) +
  labs(x = "Library",
       y = "Library size")

dev.off()


# Export to file
if(!dir.exists("OutputFiles/Tables/")){
  dir.create("OutputFiles/Tables/", recursive = TRUE)
}
write.xlsx(library_size, "OutputFiles/Tables/Enrichment_library_size.xlsx")



print(warnings())