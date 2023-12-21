set.seed(5081)



# Script to plot the overlap percentages as box plot


# Load libraries
library(tidyverse)


# Read the data
overlap_data <- list()
overlap_data[["target_enrichLib"]] <- readRDS("OutputFiles/Geneset_overlap_check/OverlapPerc_target_enrichLib.rds")
overlap_data[["target_extendedEnrichLib"]] <- readRDS("OutputFiles/Geneset_overlap_check/OverlapPerc_target_extendedEnrichLib.rds")
overlap_data[["rankedGene_enrichLib"]] <- readRDS("OutputFiles/Geneset_overlap_check/OverlapPerc_rankedGene_enrichLib.rds")
overlap_data[["rankedGene_extendedEnrichLib"]] <- readRDS("OutputFiles/Geneset_overlap_check/OverlapPerc_rankedGene_extendedEnrichLib.rds")


# Prepare the data for plotting
overlap_data <- unlist(unlist(overlap_data, recursive = FALSE), recursive = FALSE)
overlap_data <- bind_rows(overlap_data, .id = "name")
overlap_data <- separate(overlap_data, col = "name", 
                         into = c("Overlap_type", "Drug_target_type", "Disease"), 
                         sep = "\\.")


overlap_data <- pivot_longer(data = overlap_data, 
                            cols = c(starts_with("[DISEASE]"), starts_with("[ADR]")),
                            names_to = "Enrichment_lib", 
                            values_to = "Overlap_percentage", 
                            values_drop_na = TRUE)

overlap_data$Overlap_type <- factor(overlap_data$Overlap_type, 
                                       levels = c("target_enrichLib", "target_extendedEnrichLib", "rankedGene_enrichLib", "rankedGene_extendedEnrichLib"))
overlap_data$Drug_target_type <- factor(overlap_data$Drug_target_type, 
                                    levels = c("known", "KEGG", "NPA", "PS", "RI", "SIGNOR", "all"))

# Plot


if(!dir.exists("OutputFiles/Geneset_overlap_check/")){
  dir.create("OutputFiles/Geneset_overlap_check/", recursive = TRUE)
}


tiff("OutputFiles/Geneset_overlap_check/OverlapPerc_boxplot.tiff",
     width = 25, height = 30,
     units = "cm", compression = "lzw", res = 1200)




ggplot(overlap_data) +
  geom_boxplot(aes(x = Disease, y = Overlap_percentage, fill = Overlap_type), 
               width = 0.5, 
               lwd = 0.1,
               outlier.shape = 3,
               outlier.size = 0.5) + 
  facet_grid(rows = vars(Drug_target_type)) + 
  stat_summary(fun = mean, 
               show.legend = FALSE,
               mapping = aes(x = Disease, y = Overlap_percentage, fill = Overlap_type), 
               position = position_dodge2(width = 0.5,   
                                          preserve = "single"),
               geom = "point", size = 1, color = "red") +
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
        legend.box.background = element_rect(colour = "black", linewidth = 0.25)) 

dev.off()


print(warnings())