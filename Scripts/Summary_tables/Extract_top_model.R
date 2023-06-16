# Script to find the top models for each feature type across all cancers



# Load libraries
library(tidyverse)
library(openxlsx)



mean_accuracy_scores <- read.csv("OutputFiles/Tables/panCancer_meanModelAccuracy.csv", header = TRUE)

# imbalance <- "none"
# featureType <- "DrgDisAdr_BbsiProx_centre"
# disease <- "LungCancer"

top_model_list <- list()
for(imbalance in unique(mean_accuracy_scores$imbalance)){
  top_model_df <- data.frame()
  for(featureType in unique(mean_accuracy_scores$featureType)){
    for(disease in unique(mean_accuracy_scores$disease)){
      tmp1 <- mean_accuracy_scores[mean_accuracy_scores$imbalance == imbalance &
                                     mean_accuracy_scores$featureType == featureType &
                                     mean_accuracy_scores$disease == disease, ]
      tmp1 <- tmp1$model[tmp1$median_PRAUC_test == max(tmp1$median_PRAUC_test, na.rm = TRUE)]
      top_model_df <- rbind(top_model_df, data.frame("featureType" = featureType,
                                                     "disease" = disease,
                                                     "top_model" = tmp1))
    }
  }
  top_model_df <- reshape(top_model_df,
                          timevar = "disease",
                          idvar = "featureType",
                          v.names = "top_model",
                          direction = "wide")
  colnames(top_model_df) <- gsub("^top_model.", "", colnames(top_model_df))
  top_model_list[[imbalance]] <-    top_model_df
}

write.xlsx(top_model_list, "OutputFiles/Tables/panCancer_topModels_PRAUC.xlsx", overwrite = TRUE)



print(warnings())








