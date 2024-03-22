set.seed(5081)



# Analysis of the significant features selected (NES)


# Load libraries
library(tidyverse)
library(openxlsx)


# Set the options
feature_type <- "combinedEfficacySafety"



# Read the enrichment library sizes
sheet_names <- getSheetNames("OutputFiles/Tables/Enrichment_library_size.xlsx")
sheet_names <- sheet_names[grepl("^Efficacy_|^Safety", sheet_names)]
library_size <- data.frame()
for(sheet in sheet_names){
  
  if(grepl("^Efficacy_", sheet)){
    tmp1 <- read.xlsx("OutputFiles/Tables/Enrichment_library_size.xlsx", sheet = sheet)
    tmp1$Description <- paste0("[DISEASE] ", tmp1$Description)
  }
  
  if(grepl("^Safety", sheet)){
    tmp1 <- read.xlsx("OutputFiles/Tables/Enrichment_library_size.xlsx", sheet = sheet)
    tmp1$Description <- paste0("[ADR] ", tmp1$Description)
  }
  
  library_size <- rbind(library_size, tmp1)
}

library_size$feature_source <-  gsub("\\[(.*)\\].*\\[(.*)\\]", paste0("\\1", "__", "\\2"), library_size$Description)


# Read the Wilcoxon test results
signif_feature <- read.csv(paste0("OutputFiles/Plots_publication/NES_", feature_type, "_boxplot_topSignifFeatures_EA/NES_", feature_type, "_WilcoxTest_rest.csv"))
signif_feature$feature_source <-  gsub("\\[(.*)\\].*\\[(.*)\\]", paste0("\\1", "__", "\\2"), signif_feature$feature)


# Add the library sizes
# signif_feature$lib_size <- library_size$inNet_size[match(gsub("^\\[ADR\\] |\\[DISEASE\\] ", "", signif_feature$feature), library_size$Description)]
signif_feature$lib_size <- library_size$inNet_size[match(signif_feature$feature, library_size$Description)]


# Filter the data 
signif_feature <- signif_feature[signif_feature$drug_target_type == "known", ]
signif_feature <- signif_feature[signif_feature$p_val <= 0.05, ]


# Count the number of significant features from each source
signif_feature_count <- signif_feature %>%
  group_by(drug_target_type, disease, feature_source) %>%  # Group by combination
  summarize(count = n()) 
signif_feature_count <- pivot_wider(data = signif_feature_count, 
                                    names_from = "feature_source", 
                                    values_from = "count")
View(signif_feature_count)
sort(apply(signif_feature_count[, -c(1:2)], 2, FUN = function(x){ sum(!is.na(x)) }))




tmp1 <- signif_feature[, c("disease", "feature", "W", "p_val")]
write.csv(tmp1, "Trash/Significant_features_NES.csv", row.names = FALSE)



print(warnings())