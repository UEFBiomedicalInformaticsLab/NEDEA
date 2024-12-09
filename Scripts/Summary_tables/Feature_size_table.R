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
diseases <- c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")

for(disease in diseases){
  enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
  size <- as.data.frame(lengths(enrichment_lib))
  colnames(size) <- "all_size"
  size_1 <- rownames_to_column(size, "Description")
  enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
  size <- as.data.frame(lengths(enrichment_lib))
  colnames(size) <- "inNet_size"
  size_2 <- rownames_to_column(size, "Description")
  size <- reduce(list(size_1, size_2), merge, by = "Description", all = TRUE)
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
size <- reduce(list(size_1, size_2), merge, by = "Description", all = TRUE)
library_size[["Safety"]] <- size



# Compile KEGG library
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/CHG_keggPath2Gene_lib.rds")
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- reduce(list(size_1, size_2), merge, by = "Description", all = TRUE)
library_size[["KEGG"]] <- size



# Compile SMPDB (Drug Metabolsim) library
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds")
enrichment_lib <- enrichment_lib$`Drug Metabolism`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- reduce(list(size_1, size_2), merge, by = "Description", all = TRUE)
library_size[["SMPDB_DrugMetabolism"]] <- size



# Compile SMPDB (Drug Action) library
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds")
enrichment_lib <- enrichment_lib$`Drug Action`
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- reduce(list(size_1, size_2), merge, by = "Description", all = TRUE)
library_size[["SMPDB_DrugAction"]] <- size



# Compile miscelleneous library
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/miscellaneous_gene_lib.rds")
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "all_size"
size_1 <- rownames_to_column(size, "Description")
enrichment_lib <- lapply(enrichment_lib, function(x){x[x %in% V(input_network)$name]})
size <- as.data.frame(lengths(enrichment_lib))
colnames(size) <- "inNet_size"
size_2 <- rownames_to_column(size, "Description")
size <- reduce(list(size_1, size_2), merge, by = "Description", all = TRUE)
library_size[["miscellaneous"]] <- size


#####


# Plot the library sizes as box plot

if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}
tiff("OutputFiles/Plots/Enrichment_library_size.tiff",
     width = 20, height = 10,
     units = "cm", compression = "lzw", res = 1200)

plot_data <- bind_rows(library_size, .id = "Source")

plot_data$all_size[plot_data$all_size >= 2000] <- 2000
plot_data$inNet_size[plot_data$inNet_size >= 2000] <- 2000

plot_data$Source <- factor(plot_data$Source, 
                           levels = c("Efficacy_BreastCancer", "Efficacy_KidneyCancer", "Efficacy_LungCancer",  
                                      "Efficacy_OvaryCancer",  "Efficacy_ProstateCancer","Efficacy_SkinCancer", 
                                      "Safety", "KEGG", "SMPDB_DrugMetabolism", "SMPDB_DrugAction", "miscellaneous"))


plot_data_summary <- plot_data %>% group_by(Source) %>% summarise(n = n())



plot_data <- pivot_longer(plot_data, 
                          cols = c("all_size", "inNet_size"), 
                          names_to = "type", 
                          values_to = "size")

plot_data$type <- factor(plot_data$type, levels = c("all_size", "inNet_size"))

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


#####


# Shortened plot for publication

if(!dir.exists("OutputFiles/Plots_publication/")){
  dir.create("OutputFiles/Plots_publication/", recursive = TRUE)
}
tiff("OutputFiles/Plots_publication/Publication_Enrichment_library_size.tiff",
     width = 8, height = 4,
     units = "cm", compression = "lzw", res = 1200)



plot_data <- plot_data[plot_data$Source %in% c("Efficacy_BreastCancer", "Efficacy_KidneyCancer", "Efficacy_LungCancer", 
                                               "Efficacy_OvaryCancer", "Efficacy_ProstateCancer", "Efficacy_SkinCancer", 
                                               "Safety"), ]
plot_data$Source <- droplevels(plot_data$Source)

plot_data$Source <- factor(x = plot_data$Source, 
                           levels = c("Efficacy_BreastCancer", "Efficacy_KidneyCancer", "Efficacy_LungCancer", 
                                      "Efficacy_OvaryCancer", "Efficacy_ProstateCancer", "Efficacy_SkinCancer", 
                                      "Safety"), 
                           labels = c("BreastCancer", "KidneyCancer", "LungCancer", 
                                      "OvaryCancer", "ProstateCancer", "SkinCancer", 
                                      "Safety") )

plot_data_summary$Source <- gsub("^Efficacy_", "", plot_data_summary$Source)
plot_data_summary <- plot_data_summary[plot_data_summary$Source %in% plot_data$Source, ]

plot_data$size[plot_data$size >= 1000] <- 1000


ggplot() +
  geom_boxplot(data = plot_data, 
               aes(x = Source, y = size, fill = type), 
               width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3,
               outlier.size = 0.5,
               outlier.stroke = 0.1) +
  geom_text(data = plot_data_summary,
            aes(x = Source, y = Inf, label = n),
            color = "red", size = 2, vjust = 2) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        text = element_text(size = 4),
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 3),
        legend.text = element_text(size = 3),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.25)) +
  labs(x = "Library",
       y = "Library size",
       fill = "Type") 


dev.off()


#####


# Plot the library composition as bar plot

if(!dir.exists("OutputFiles/Plots_publication/")){
  dir.create("OutputFiles/Plots_publication/", recursive = TRUE)
}
tiff("OutputFiles/Plots_publication/Publication_Enrichment_library_list_size.tiff",
     width = 21, height = 29,
     units = "cm", compression = "lzw", res = 1200)

ggplot() +
  geom_bar(data = plot_data, 
           aes(x = Description, y = size, fill = type), 
           stat = "identity", 
           position = "dodge") +
  facet_wrap(.~Source, scales = "free", ncol = 1, labeller = labeller(.cols = function(x){ gsub("Cancer$", " Cancer", x) })   ) +
  scale_x_discrete(labels = function(x) scales::label_wrap(28)(x)) +
  scale_fill_manual(values = c("all_size" = "#FF9F00", "inNet_size" = "#007FFF"), labels = c("all_size" = "All", "inNet_size" = "In network")) +
  labs(x = "Gene sets", 
       y = "Library size",
       fill = "Number of genes:") +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = NULL),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", linewidth = 0.25,),
        strip.text = element_text(size = 4, margin = margin(1,1,1,1)),
        text = element_text(size = 4), 
        axis.title = element_text(size = 4), 
        axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust = 1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        legend.position = "bottom",
        legend.title = element_text(margin = margin(r = 2)),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.spacing.x = unit(0.1, "cm"),
        legend.text = element_text(size = 5, margin = margin(l = 1)),
        legend.margin = margin(1,1,1,1),
        legend.box.spacing = unit(0.1, 'cm'),
        legend.box.background = element_rect(colour = "black", linewidth = 0.1))

dev.off()



print(warnings())