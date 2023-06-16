# Script to download CPTAC-3 lung cancer data from GDC
# Downloads only tumor data

# set temporary directory 
library(unixtools) 
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)} 
set.tempdir("tmp_dir/")

library(tidyverse)
library(GenomicDataCommons)
GenomicDataCommons::status()
if(!dir.exists("Databases/CPTAC/")){dir.create("Databases/CPTAC/", recursive = TRUE)} 
gdc_set_cache(paste0(getwd(), "/Databases/CPTAC/cache"))

# Check all possible project IDs
# available_values('files','cases.project.project_id')

# Create manifest for download
CPTAC_manifest <- files() %>%
  GenomicDataCommons::filter( cases.project.project_id == "CPTAC-3") %>% 
  GenomicDataCommons::filter( data_category == "transcriptome profiling" ) %>%
  GenomicDataCommons::filter( data_type == "Gene Expression Quantification" ) %>% 
  GenomicDataCommons::filter( experimental_strategy == "RNA-Seq" ) %>%                                                            
  GenomicDataCommons::filter( cases.primary_site == "bronchus and lung")  %>%
  GenomicDataCommons::filter( cases.samples.tissue_type == "tumor") %>%
  GenomicDataCommons::filter( data_format == "tsv")  %>%
  manifest() 

# Download data
file_names <- as.list(gdcdata(CPTAC_manifest$id, progress = TRUE))

# Retrieve clinical data for the samples
CPTAC_id_mapping <- files() %>%
  GenomicDataCommons::filter( cases.project.project_id == "CPTAC-3") %>% 
  GenomicDataCommons::filter( data_category == "transcriptome profiling" ) %>%
  GenomicDataCommons::filter( data_type == "Gene Expression Quantification" ) %>% 
  GenomicDataCommons::filter( experimental_strategy == "RNA-Seq" ) %>%                                                            
  GenomicDataCommons::filter( cases.primary_site == "bronchus and lung")  %>%
  GenomicDataCommons::filter( cases.samples.tissue_type == "tumor") %>%
  GenomicDataCommons::filter( data_format == "tsv")  %>%
  GenomicDataCommons::select(fields = c("cases.case_id", "cases.submitter_id", "cases.samples.tissue_type")) %>%
  results_all()
CPTAC_id_mapping <- do.call(rbind, CPTAC_id_mapping$cases)
CPTAC_id_mapping$tissue_type <- unlist(lapply(CPTAC_id_mapping$samples, unique), use.names = FALSE)
CPTAC_id_mapping$file_id <- row.names(CPTAC_id_mapping)
row.names(CPTAC_id_mapping) <- NULL
CPTAC_id_mapping$submitter_id <- gsub("-", "_", CPTAC_id_mapping$submitter_id)

duplicate_ids <- CPTAC_id_mapping$submitter_id[duplicated(CPTAC_id_mapping$submitter_id)]

for(i in duplicate_ids){
  tmp <- CPTAC_id_mapping[CPTAC_id_mapping$submitter_id == i, "submitter_id"]
  CPTAC_id_mapping[CPTAC_id_mapping$submitter_id == i, "submitter_id"] <- paste(tmp, 1:length(tmp), sep = "_")
}

CPTAC_clinical_data <- gdc_clinical(CPTAC_id_mapping$case_id) 


df_list <- list(CPTAC_clinical_data$main[, c("case_id", "disease_type", "primary_site")],
                CPTAC_id_mapping[, c("case_id", "submitter_id", "tissue_type", "file_id")],
                CPTAC_clinical_data$demographic[, c("case_id", "gender")],
                CPTAC_clinical_data$diagnoses[, c("case_id", "tumor_grade", "primary_diagnosis", "ajcc_pathologic_stage")])

CPTAC_clinical_data <- df_list %>% reduce(full_join, by='case_id')



write.csv(CPTAC_clinical_data, "Databases/CPTAC/CPTAC_LungCancer_clinical_data.csv", row.names = FALSE, quote = TRUE)


# Compile counts from all files to single matrix
counts_unstranded <- data.frame()
tpm_unstranded <- data.frame()
fpkm_unstranded <- data.frame()
fpkm_uq_unstranded <- data.frame()

for(names in names(file_names)){
  
  # Read table
  file_path <- file_names[[names]]
  submitter_id <- CPTAC_id_mapping[CPTAC_id_mapping$file_id == names, "submitter_id"]
  count_table <- read.table(file_path, sep = "\t", fill = TRUE, header = FALSE, skip = 6)
  colnames(count_table) <- c("gene_id",	"gene_name", "gene_type", 
                             paste(submitter_id, c("unstranded", "stranded_first", "stranded_second", 
                                                   "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"), sep = "__"))
  
  # Extract raw counts
  if(nrow(counts_unstranded) == 0){
    counts_unstranded <- count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "unstranded", sep = "__"))]
  }else{
    counts_unstranded <- merge(counts_unstranded, 
                               count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "unstranded", sep = "__"))], 
                               all = TRUE, by = c("gene_id", "gene_name", "gene_type"))
  }
  colnames(counts_unstranded) <- gsub(pattern = "__unstranded$", replacement = "", x = colnames(counts_unstranded))
  
  
  # Extract TPM normalised counts
  if(nrow(tpm_unstranded) == 0){
    tpm_unstranded <- count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "tpm_unstranded", sep = "__"))]
  }else{
    tpm_unstranded <- merge(tpm_unstranded, 
                            count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "tpm_unstranded", sep = "__"))], 
                            all = TRUE, by = c("gene_id", "gene_name", "gene_type"))
  }
  colnames(tpm_unstranded) <- gsub(pattern = "__tpm_unstranded$", replacement = "", x = colnames(tpm_unstranded))
  
  # Extract FPKM 
  if(nrow(fpkm_unstranded) == 0){
    fpkm_unstranded <- count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "fpkm_unstranded", sep = "__"))]
  }else{
    fpkm_unstranded <- merge(fpkm_unstranded, 
                             count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "fpkm_unstranded", sep = "__"))], 
                             all = TRUE, by = c("gene_id", "gene_name", "gene_type"))
  }
  colnames(fpkm_unstranded) <- gsub(pattern = "__fpkm_unstranded$", replacement = "", x = colnames(fpkm_unstranded))
  
  # Extract FPKM-UQ
  if(nrow(fpkm_uq_unstranded) == 0){
    fpkm_uq_unstranded <- count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "fpkm_uq_unstranded", sep = "__"))]
  }else{
    fpkm_uq_unstranded <- merge(fpkm_uq_unstranded, 
                                count_table[, c("gene_id", "gene_name", "gene_type", paste(submitter_id, "fpkm_uq_unstranded", sep = "__"))], 
                                all = TRUE, by = c("gene_id", "gene_name", "gene_type"))
  }
  colnames(fpkm_uq_unstranded) <- gsub(pattern = "__fpkm_uq_unstranded$", replacement = "", x = colnames(fpkm_uq_unstranded))
  }

write.csv(counts_unstranded, "Databases/CPTAC/CPTAC_LungCancer_RNAseq_counts.csv", row.names = FALSE, quote = TRUE)
write.csv(tpm_unstranded, "Databases/CPTAC/CPTAC_LungCancer_RNAseq_TPM.csv", row.names = FALSE, quote = TRUE)
write.csv(fpkm_unstranded, "Databases/CPTAC/CPTAC_LungCancer_RNAseq_FPKM.csv", row.names = FALSE, quote = TRUE)
write.csv(fpkm_uq_unstranded, "Databases/CPTAC/CPTAC_LungCancer_RNAseq_FPKMUQ.csv", row.names = FALSE, quote = TRUE)
