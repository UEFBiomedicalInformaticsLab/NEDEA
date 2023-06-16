set.seed(5081)



# Script to compare similarity between selected genes from RWR using dynamic quartile probability thresholding (version 5)



# Load libraries
library(optparse)
library(foreach)
library(doParallel)
library(tidyverse)
library(openxlsx)
library(ComplexHeatmap)
library(svglite)





# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
              help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer (if used).", call.=FALSE)
  }
}

# Define global options for this script 
disease <- opt$disease
nproc <- opt$nproc

cat(paste0("\n\nCalculating similarity for: ", disease, "\n\n"))



# Load RWR results
load(file = paste0("OutputFiles/Model_train/", disease, "/dNetRWR050_", disease, ".rda"))


effectiveComb_rwr <- drugCombs_rwr_res_final$effectiveCombinations
adverseComb_rwr <- drugCombs_rwr_res_final$adverseCombinations

# Check if all resuls have valied number of columns
effectiveComb_rwr <- effectiveComb_rwr[lapply(effectiveComb_rwr, ncol) == 7]
adverseComb_rwr <- adverseComb_rwr[lapply(adverseComb_rwr, ncol) == 7]
rm(drugCombs_rwr_res_final)


cat("\n\nNumbe rof drug combinations: ")
cat(paste0("\n - Effective combinations = ", length(effectiveComb_rwr)))
cat(paste0("\n - Adverse combinations = ", length(adverseComb_rwr), "\n\n"))


drugCombs <- c(effectiveComb_rwr, adverseComb_rwr)





# Calculate the Jaccard similarity between the top ranked genes
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl)


similarity_df <- foreach(i=names(drugCombs), .combine = "rbind") %:%
  foreach(j=names(drugCombs), .combine = "rbind") %dopar% {
    
    ranked_gene_i <- drugCombs[[i]]$union_seed
    names(ranked_gene_i) <- drugCombs[[i]]$node_name
    ranked_gene_i <- ranked_gene_i[order((ranked_gene_i), decreasing = TRUE)]
    ranked_gene_i <- ranked_gene_i[ranked_gene_i > quantile(unique(ranked_gene_i), prob = c(.9))]
    
    ranked_gene_j <- drugCombs[[j]]$union_seed
    names(ranked_gene_j) <- drugCombs[[j]]$node_name
    ranked_gene_j <- ranked_gene_j[order(ranked_gene_j, decreasing = TRUE)]
    ranked_gene_j <- ranked_gene_j[ranked_gene_j > quantile(unique(ranked_gene_j), prob = c(.9))]
    
    intersection <- intersect(names(ranked_gene_i), names(ranked_gene_j))
    union <- union(names(ranked_gene_i), names(ranked_gene_j))
    jaccard <- length(intersection) / length(union)
    
    tmp1 <-  data.frame(drugComb_1 = i,
                        drugComb_2 = j, 
                        similarity = jaccard)
    tmp1
  }
stopCluster(cl)



similarity_matrix <- reshape(similarity_df, 
                             idvar = "drugComb_1",
                             timevar = "drugComb_2",
                             v.names = "similarity",
                             direction = "wide")


colnames(similarity_matrix) <- gsub("similarity.", "", colnames(similarity_matrix))
colnames(similarity_matrix) <- gsub("drugComb_1", "drugComb", colnames(similarity_matrix))


tmp2 <- list("similarity_matrix" = similarity_matrix,
             "similarity_df" = similarity_df)

write.xlsx(tmp2, file = paste0("OutputFiles/Model_train/", disease, "/rwrGeneRank_probCut_", disease, ".xlsx"), overwrite = TRUE)





# Create heatmap of the similarity
row.names(similarity_matrix) <- NULL
similarity_matrix <- column_to_rownames(similarity_matrix, var = "drugComb")

svglite(paste0("OutputFiles/Model_train/", disease, "/rwrSimMat_probCut_", disease, ".svg"), width = 10, height = 10)
col_annot <- HeatmapAnnotation(Class = substring(colnames(similarity_matrix), 1, 3),
                               col = list(Class = c("Adv" = "Orange", "Eff" = "Green")))
Heatmap(as.matrix(similarity_matrix),
        name = paste("Jaccard", "similarity", sep = "\n"),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE, 
        top_annotation = col_annot)
dev.off()



print(warnings())