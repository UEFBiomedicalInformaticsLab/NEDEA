set.seed(5081)


# Drug combinations from FIMM DrugComb (https://drugcomb.fimm.fi/)


# Load libraries
library(httr)
library(jsonlite)
library(tidyverse)
httr::set_config(config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))


#####


# Download drug combinations from FIMM DrugComb 
if(!dir.exists("Databases/FimmDrugComb/"))dir.create("Databases/FimmDrugComb/", recursive = TRUE)
if(!file.exists("Databases/FimmDrugComb/summary_v_1_5.csv")){
  download.file(url = "https://drugcomb.fimm.fi/jing/summary_v_1_5.csv",
                destfile = "Databases/FimmDrugComb/summary_v_1_5.csv", method = "wget")
}


# Read the DrugComb data
FimmDrugComb_drugCombCat <- read.csv("Databases/FimmDrugComb/summary_v_1_5.csv", header = TRUE)
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$drug_row != "NULL", ]
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$drug_col != "NULL", ]
FimmDrugComb_drugCombCat$synergy_loewe[which(FimmDrugComb_drugCombCat$synergy_loewe == "\\N")] <- 0
FimmDrugComb_drugCombCat$synergy_loewe <- as.numeric(FimmDrugComb_drugCombCat$synergy_loewe)


# Count number of drug combinations
# FimmDrugComb_drugCombCat %>% select(drug_row, drug_col) %>% distinct() %>% count()

# Number of drugs
# length(unique(c(FimmDrugComb_drugCombCat$drug_row, FimmDrugComb_drugCombCat$drug_col)))


# Get cell line information
FimmDrugComb_cellLine <- GET("https://api.drugcomb.org/cell_lines")
FimmDrugComb_cellLine <- fromJSON(rawToChar(FimmDrugComb_cellLine$content))
write.csv(FimmDrugComb_cellLine, "Databases/FimmDrugComb/cell_lines.csv", quote = TRUE, row.names = FALSE)


#####


# Download disease IDs from NCI Thesaurus (NCIt) 
if(!dir.exists("Databases/NCIt/"))dir.create("Databases/NCIt/", recursive = TRUE)
if(!file.exists("Databases/NCIt/Thesaurus.txt")){
  download.file(url = "https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/archive/2023/23.06d_Release/Thesaurus.FLAT.zip",
                destfile = "Databases/NCIt/Thesaurus.FLAT.zip", method = "wget")
  unzip("Databases/NCIt/Thesaurus.FLAT.zip", exdir = "Databases/NCIt/", file = "Thesaurus.txt")
}
NCIthesaurus <- read.table("Databases/NCIt/Thesaurus.txt", sep = "\t", header = FALSE, comment.char = "", fill = TRUE, quote = "")
colnames(NCIthesaurus) <- c("code", "concept IRI", "parents", "synonyms", "definition", "display name", "concept status", "semantic type", "concept in subset")


# Filter drug combinations tested on cancer related cell lines
NCIthesaurus <- NCIthesaurus[NCIthesaurus$code %in% FimmDrugComb_cellLine$disease_id, ]
NCIthesaurus <- NCIthesaurus[grep("cancer|carcinoma|sarcoma|lymphoma|leukemia|melanoma", NCIthesaurus$synonyms, ignore.case = TRUE), ]

FimmDrugComb_cellLine <- FimmDrugComb_cellLine[FimmDrugComb_cellLine$disease_id %in% NCIthesaurus$code, ]

FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[FimmDrugComb_drugCombCat$cell_line_name %in% FimmDrugComb_cellLine$name, ]

used_cellLines <- unique(FimmDrugComb_drugCombCat[, c("cell_line_name", "tissue_name")])
write.csv(used_cellLines, "OutputFiles/Tables/FIMM_used_cell_lines.csv", row.names = FALSE)


#####


# * ZIP (Zero Interaction Potency):
# The ZIP model calculates the expected effect of two drugs if there were no interaction between them. 
# The observed combined effect is then compared to this expected effect.
#
# Synergy: Observed effect > Expected effect.
# Antagonism: Observed effect < Expected effect.
# Additivity: Observed effect = Expected effect.
#
# * Loewe Additivity:
# The Loewe model uses a combination index (CI) derived from the expected dose of each drug 
# to achieve a particular effect if used alone versus when used in combination.
#
# Synergy: Observed effect > Expected effect.
# Antagonism: Observed effect < Expected effect.
# Additivity: Observed effect = Expected effect.
#
# * Bliss Independence:
# Bliss independence calculates the expected effect assuming that two drugs act independently.
#
# Synergy: Observed effect > Expected effect.
# Antagonism: Observed effect < Expected effect.
# Additivity: Observed effect = Expected effect.
#
# * HSA (Highest Single Agent):
# HSA compares the effect of a drug combination to the effect of the most effective drug in the combination used alone.
# 
# Synergy: Combined effect > Effect of the most effective single agent.
# Antagonism: Combined effect < Effect of the most effective single agent.
# Additivity: Combined effect = Effect of the most effective single agent.
#
# * FIMM Scores: The FIMM group introduced some scores to quantify drug interactions:
#   - S-score (Synergy score): This score reflects the difference between the observed and expected responses. 
#     A score near zero indicates additivity. A positive score indicates synergy, 
#     and a negative score may indicate antagonism, although the interpretation of negative scores can sometimes be ambiguous.



# Aggregate scores across tissue for each drug combination
# i.e. the mean scores across cell lines are used

dcc_agg_all <- FimmDrugComb_drugCombCat %>%
  group_by(drug_row, drug_col, tissue_name)  %>%
  summarise(avgSynSM = mean(S_mean, na.rm = TRUE),
            avgSynZIP = mean(synergy_zip, na.rm = TRUE),
            avgSynLoewe = mean(synergy_loewe, na.rm = TRUE),
            avgSynHSA = mean(synergy_hsa, na.rm = TRUE),
            avgSynBliss = mean(synergy_bliss, na.rm=  TRUE), 
            .groups="drop")


# Plot the average synergy scores 
if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}


tiff("OutputFiles/Plots/Avg_SynergyScores_byCancer.tiff",
     width = 25, height = 10,
     units = "cm", compression = "lzw", res = 1200)


plot_data <- dcc_agg_all %>%
  pivot_longer(
    cols = starts_with("avgSyn"), 
    names_to = "SynergyType", 
    values_to = "Value"
  )

ggplot(plot_data, aes(x = tissue_name, y = Value)) +
  geom_boxplot(width = 0.5, 
               lwd = 0.2, 
               outlier.shape = NA, 
               show.legend = FALSE,
               fill = "grey") +
  facet_wrap(~ SynergyType, scales = "free_y") +
  labs(title = "Average Synergy Scores Across Cancer Types",
       x = "Cancer Type",
       y = "Average Synergy Value") +
  coord_cartesian(ylim = c(-100, 100)) +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        linewidth = 0.25, 
                                        linetype = NULL),
        panel.grid = element_line(colour = "grey", 
                                  linewidth = 0.05),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", 
                                        linewidth = 0.25,),
        strip.text = element_text(size = 5, 
                                  margin = margin(1,1,1,1)),
        text = element_text(size = 5), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 4, 
                                   angle = 45, 
                                   hjust = 1), 
        axis.ticks = element_line(colour = "black", 
                                  linewidth = 0.2))


dev.off()


#####


# Identify synergistic, antagonist and neutral drugs

# When you aggregate synergy scores across multiple cell lines of the same tumor tissue or cancer type, 
# you're essentially estimating the "typical" effect of a drug pair on that tumor type. 
# This can be useful for drawing broader conclusions, especially if there is substantial variability between cell lines of the same tissue type. 
# The mean of these average scores then becomes a representative metric for that drug pair's effectiveness against a particular tumor type.

# If you compute the expected mean synergistic score for a given cancer type (by averaging all the average scores for that cancer type), 
# you're essentially setting a "baseline" or "benchmark" of synergy for that cancer type.

# Given this setup, you can apply the same reasoning to classify drug interactions at this aggregated level:

# Synergy: Observed average effect > Expected mean effect for the cancer type.
# Antagonism: Observed average effect < Expected mean effect for the cancer type.
# Additivity: Observed average effect = Expected mean effect for the cancer type.

# Note
# If your mean synergy score is positive and an individual score is just below this mean, 
# this might indicate that the interaction is less synergistic than typical for that tumor type, 
# but not necessarily antagonistic. It might still be synergistic relative to a global baseline (like 0 or 1, depending on the scale of your scores), 
# but less so than the average for that tumor type.

cancer_stats <- dcc_agg_all %>%
  group_by(tissue_name) %>%
  summarise(meanSynS = mean(avgSynSM, na.rm = TRUE),
            sdSynS = sd(avgSynSM, na.rm = TRUE),
            meanSynZIP = mean(avgSynZIP, na.rm = TRUE),
            sdSynZIP = sd(avgSynZIP, na.rm = TRUE),
            meanSynLoewe = mean(avgSynLoewe, na.rm = TRUE),
            sdSynLoewe = sd(avgSynLoewe, na.rm = TRUE),
            meanSynHSA = mean(avgSynHSA, na.rm = TRUE),
            sdSynHSA = sd(avgSynHSA, na.rm = TRUE),
            meanSynBliss = mean(avgSynBliss, na.rm = TRUE),
            sdSynBliss = sd(avgSynBliss, na.rm = TRUE))

df_with_stats <- left_join(dcc_agg_all, cancer_stats, by = "tissue_name")

df_with_stats$cS <- rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cZIP <- rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cLoewe <- rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cHSA <- rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$cBliss <- rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'
df_with_stats$Syn_level <- rep(0, nrow(df_with_stats)) # 0 means 'NEUTRAL'


k <- 1
for(i in 1:nrow(df_with_stats)) {
  
  # S
  if(df_with_stats$avgSynSM[i] > df_with_stats$meanSynS[i] + k * df_with_stats$sdSynS[i] & df_with_stats$avgSynSM[i] > 0) df_with_stats$cS[i] = 1
  if(df_with_stats$avgSynSM[i] < df_with_stats$meanSynS[i] - k * df_with_stats$sdSynS[i] & df_with_stats$avgSynSM[i] < 0) df_with_stats$cS[i] = -1
  
  # ZIP
  if(df_with_stats$avgSynZIP[i] > df_with_stats$meanSynZIP[i] + k * df_with_stats$sdSynZIP[i] & df_with_stats$avgSynZIP[i] > 0) df_with_stats$cZIP[i] = 1
  if(df_with_stats$avgSynZIP[i] < df_with_stats$meanSynZIP[i] - k * df_with_stats$sdSynZIP[i] & df_with_stats$avgSynZIP[i] < 0) df_with_stats$cZIP[i] = -1
  
  # LOEWE
  if(df_with_stats$avgSynLoewe[i] > df_with_stats$meanSynLoewe[i] + k * df_with_stats$sdSynLoewe[i] & df_with_stats$avgSynLoewe[i] > 0) df_with_stats$cLoewe[i] = 1
  if(df_with_stats$avgSynLoewe[i] < df_with_stats$meanSynLoewe[i] - k * df_with_stats$sdSynLoewe[i] & df_with_stats$avgSynLoewe[i] < 0) df_with_stats$cLoewe[i] = -1
  
  # HSA
  if(df_with_stats$avgSynHSA[i] > df_with_stats$meanSynHSA[i] + k * df_with_stats$sdSynHSA[i] & df_with_stats$avgSynHSA[i] > 0) df_with_stats$cHSA[i] = 1
  if(df_with_stats$avgSynHSA[i] < df_with_stats$meanSynHSA[i] - k * df_with_stats$sdSynHSA[i] & df_with_stats$avgSynHSA[i] < 0) df_with_stats$cHSA[i] = -1
  
  # BLISS
  if(df_with_stats$avgSynBliss[i] > df_with_stats$meanSynBliss[i] + k * df_with_stats$sdSynBliss[i] & df_with_stats$avgSynBliss[i] > 0) df_with_stats$cBliss[i] = 1
  if(df_with_stats$avgSynBliss[i] < df_with_stats$meanSynBliss[i] - k * df_with_stats$sdSynBliss[i] & df_with_stats$avgSynBliss[i] < 0) df_with_stats$cBliss[i] = -1
  
  # consensus
  df_with_stats$Syn_level[i] = df_with_stats$cS[i] + df_with_stats$cZIP[i] + df_with_stats$cLoewe[i] + df_with_stats$cHSA[i] + df_with_stats$cBliss[i]
}



# Plot the average synergy scores as points along with class
if(!dir.exists("OutputFiles/Plots/")){
  dir.create("OutputFiles/Plots/", recursive = TRUE)
}

tiff("OutputFiles/Plots/Avg_SynergyScores_scatter_byCancer.tiff",
     width = 25, height = 10,
     units = "cm", compression = "lzw", res = 1200)

plot_data_1 <- df_with_stats %>%
  select(c("drug_row", "drug_col", "tissue_name"), starts_with("avgSyn")) %>%
  pivot_longer(
    cols = starts_with("avgSyn"), 
    names_to = "SynergyType", 
    values_to = "synergy_value"
  )
plot_data_1$SynergyType <- gsub("avgSyn", "", plot_data_1$SynergyType)
plot_data_1$SynergyType <- gsub("SM", "S", plot_data_1$SynergyType)


plot_data_2 <- df_with_stats %>%
  select(c("drug_row", "drug_col", "tissue_name"), starts_with("c")) %>%
  pivot_longer(
    cols = starts_with("c"), 
    names_to = "SynergyType", 
    values_to = "synergy_class"
  )
plot_data_2$SynergyType <- gsub("c", "", plot_data_2$SynergyType)



plot_data_3 <- df_with_stats %>%
  select(c("drug_row", "drug_col", "tissue_name", "Syn_level")) 

plot_data <- merge(plot_data_1, plot_data_2, by = c("drug_row", "drug_col", "tissue_name", "SynergyType"))
plot_data <- plot_data %>% left_join(plot_data_3, by = c("drug_row", "drug_col", "tissue_name"))

plot_data <- plot_data[plot_data$tissue_name %in% c("breast", "kidney", "lung", "ovary", "prostate", "skin"), ]

plot_data$synergy_class <- as.factor(plot_data$synergy_class)
plot_data$Syn_level <- as.factor(plot_data$Syn_level)

row.names(plot_data) <- NULL

ggplot(data = plot_data, mapping = aes(x = tissue_name, y = synergy_value, color = synergy_class)) +
  geom_point(size = 0.5, shape = 3) +
  facet_wrap(~ SynergyType, scales = "free_y") +
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 1,
    shape = 24,
    fill = "red"
  ) +
  scale_x_discrete(labels = c("breast" = "Breast", "kidney" = "Kidney", "lung" = "Lung", 
                              "ovary" = "Ovary", "prostate" = "Prostate", "skin" = "Skin")) +
  labs(x = "Cancer type",
       y = "Average synergy value",
       color = "Synergy class") +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        linewidth = 0.25, 
                                        linetype = NULL),
        panel.grid = element_line(colour = "grey", 
                                  linewidth = 0.05),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_rect(color = "black", 
                                        linewidth = 0.25,),
        strip.text = element_text(size = 8, 
                                  margin = margin(1,1,1,1)),
        text = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, 
                                   angle = 45, 
                                   hjust = 1), 
        axis.ticks = element_line(colour = "black", 
                                  linewidth = 0.2))
dev.off()


#####


# Map drugs to DrugBank drug ID
FimmDrugComb_drugs <- GET("https://api.drugcomb.org/drugs")
FimmDrugComb_drugs <- fromJSON(rawToChar(FimmDrugComb_drugs$content))

FimmDrugComb_drugCombCat <- df_with_stats

FimmDrugComb_drugCombCat$Drug1_DrugBank_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat$drug_row, FimmDrugComb_drugs$dname)]
FimmDrugComb_drugCombCat$Drug2_DrugBank_id <- FimmDrugComb_drugs$drugbank_id[match(FimmDrugComb_drugCombCat$drug_col, FimmDrugComb_drugs$dname)]

FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[(FimmDrugComb_drugCombCat$Drug1_DrugBank_id != "NA" & 
                                                        FimmDrugComb_drugCombCat$Drug2_DrugBank_id != "NA"),]
FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[(FimmDrugComb_drugCombCat$Drug1_DrugBank_id != "" & 
                                                        FimmDrugComb_drugCombCat$Drug2_DrugBank_id != ""),]

# Assign the categories
FimmDrugComb_drugCombCat$class_synergyScore <- rep("neutral", nrow(FimmDrugComb_drugCombCat))
FimmDrugComb_drugCombCat$class_synergyScore[which(FimmDrugComb_drugCombCat$Syn_level >= 3)] <- "synergism"
FimmDrugComb_drugCombCat$class_synergyScore[which(FimmDrugComb_drugCombCat$Syn_level <= -3)] <- "antagonism"


#####


# Read the drug type information
DrugBank_drug_type <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_type <- DrugBank_drug_type$drugs$general_information
DrugBank_drug_type <- DrugBank_drug_type[DrugBank_drug_type$type == "small molecule", ] # retain only small molecular drugs


FimmDrugComb_drugCombCat <- FimmDrugComb_drugCombCat[(FimmDrugComb_drugCombCat$Drug1_DrugBank_id %in% DrugBank_drug_type$primary_key) &
                                                       (FimmDrugComb_drugCombCat$Drug2_DrugBank_id %in% DrugBank_drug_type$primary_key), ]


if(!dir.exists("InputFiles/Reference_list/"))dir.create("InputFiles/Reference_list/", recursive = TRUE)
saveRDS(FimmDrugComb_drugCombCat, "InputFiles/Reference_list/FimmDrugComb_drugCombinations.rds")


#####


print(warnings())