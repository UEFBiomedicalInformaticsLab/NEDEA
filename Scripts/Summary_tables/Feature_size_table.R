set.seed(5081)



# Script to compile the library terms and sizes used for enrichment



# Load libraries
library(tidyverse)
library(openxlsx)
library(igraph)


library_size <- list()

# Read the network on which to execute RWR
input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net.rds")
cat(paste0("\n\nInput network size:: vertices = ", vcount(input_network), ", edges = ", ecount(input_network), "\n\n"))



# Compile efficacy libraries
diseases <- c("LungCancer", "BreastCancer", "ProstateCancer", "OvaryCancer", "KidneyCancer", "SkinCancer")

for(disease in diseases){
  enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  size <- as.data.frame(lengths(enrichment_lib))
  colnames(size) <- "all_size"
  size_1 <- rownames_to_column(size, "Description")
  enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
  size <- as.data.frame(lengths(enrichment_lib))
  colnames(size) <- "inNet_size"
  size_2 <- rownames_to_column(size, "Description")
  size <- merge(size_1, size_2)
  library_size[[paste0("Efficacy_", disease)]] <- size
}



# Compile safety library
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/curatedAdr2Gene_lib.rds")
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- merge(size_1, size_2)
library_size[["Safety"]] <- size



# Compile KEGG library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/CHG_keggPath2Gene_lib.rds"))
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- merge(size_1, size_2)
library_size[["KEGG"]] <- size



# Compile SMPDB (Drug Metabolsim) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Metabolism`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- merge(size_1, size_2)
library_size[["SMPDB_DrugMetabolism"]] <- size



# Compile SMPDB (Drug Action) library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Action`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- merge(size_1, size_2)
library_size[["SMPDB_DrugAction"]] <- size



# Compile miscelleneous library
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/miscellaneous_gene_lib.rds"))
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- merge(size_1, size_2)
library_size[["miscellaneous"]] <- size



# Plot the library sizes as violin plot

if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}
tiff("OutputFiles/Plots/Enrichment_library_size.tiff",
     width = 20, height = 10,
     units = "cm", compression = "lzw", res = 1200)

plot_data <- bind_rows(library_size, .id = "Source")

plot_data$all_size[plot_data$all_size >= 500] <- 500
plot_data$inNet_size[plot_data$inNet_size >= 500] <- 500

plot_data$Source <- factor(plot_data$Source, 
                           levels = c("Efficacy_LungCancer", "Efficacy_BreastCancer", "Efficacy_ProstateCancer", 
                                      "Efficacy_OvaryCancer", "Efficacy_KidneyCancer", "Efficacy_SkinCancer", 
                                      "Safety", "KEGG", "SMPDB_DrugMetabolism", "SMPDB_DrugAction", "miscellaneous"))



plot_data_summary <- plot_data %>% group_by(Source) %>% summarise(n = n())



plot_data <- pivot_longer(plot_data, 
                          cols = c("all_size", "inNet_size"), 
                          names_to = "type", 
                          values_to = "size")


ggplot() +
  geom_boxplot(data = plot_data, 
               aes(x = Source, y = size, fill = type), 
               width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3,
               outlier.size = 0.5) +
  geom_text(data = plot_data_summary,
            aes(x = Source, y = Inf, label = n),
            color = "red", size = 2, vjust = 2) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 3),
        legend.text = element_text(size = 3),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
  labs(x = "Library",
       y = "Library size",
       fill = "Type")

dev.off()


# Export to file
if(!dir.exists("OutputFiles/Tables/")){
  dir.create("OutputFiles/Tables/", recursive = TRUE)
}
write.xlsx(library_size, "OutputFiles/Tables/Enrichment_library_size.xlsx")



print(warnings())