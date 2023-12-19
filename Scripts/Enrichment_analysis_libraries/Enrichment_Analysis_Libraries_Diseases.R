set.seed(5081)



# Libraries for enrichment analysis of diseases



# Load libraries 
library(unixtools)
library(org.Hs.eg.db)
library(tidyverse)
library(sparklyr)
library(sparklyr.nested)
source("Scripts/Functions/Functions_data_manipulation.R")
source("Scripts/Functions/Functions_ID_Conversion.R")



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")



## Enrichment of DisGeNET diseases -------------------------------------------------------------


# Retrieve Curated gene-disease associations from https://www.disgenet.org/
# The file contains gene-disease associations from UNIPROT, CGI, ClinGen, Genomics England, CTD (human subset), PsyGeNET, and Orphanet.
if(!dir.exists("Databases/DisGeNET/")){dir.create("Databases/DisGeNET/", recursive = TRUE)}
if(!file.exists("Databases/DisGeNET/curated_gene_disease_associations.tsv.gz")){
  download.file(url = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz",
                destfile = "Databases/DisGeNET/curated_gene_disease_associations.tsv.gz", method = "wget")
}


DisGeNET_data <- read.table(gzfile("Databases/DisGeNET/curated_gene_disease_associations.tsv.gz"), header = T, sep = "\t", quote = "")


# Add ensembl gene IDs
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# entrezId_2_ensemblId <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "entrezgene_accession", "description", "gene_biotype"), 
#                               filters = "entrezgene_id",
#                               values = DisGeNET_data$geneId,
#                               mart = ensembl)
# DisGeNET_data$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_gene_id[match(DisGeNET_data$geneId, entrezId_2_ensemblId$entrezgene_id)]
DisGeNET_data$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(DisGeNET_data$geneId, entrezId_2_ensemblId$gene_id)]
DisGeNET_data <- DisGeNET_data[!is.na(DisGeNET_data$ensembl_gene_id ), ]


### Enrichment of DisGeNET diseases 
# diseaseId_2_diseaseName <- unique(DisGeNET_data[DisGeNET_data$diseaseType == "disease",c("diseaseId", "diseaseName")])

DisGeNET_Gene_Disease <- na.exclude(DisGeNET_data[DisGeNET_data$diseaseType == "disease",c("ensembl_gene_id", "diseaseId", "diseaseName")])
DisGeNET_Gene_Disease$diseaseName <- paste0(DisGeNET_Gene_Disease$diseaseName, " (", DisGeNET_Gene_Disease$diseaseId, ")")

DisGeNET_Disease2Gene_lib <- split(x = DisGeNET_Gene_Disease$ensembl_gene_id, f = DisGeNET_Gene_Disease$diseaseName)
DisGeNET_Disease2Gene_lib <- DisGeNET_Disease2Gene_lib[lengths(DisGeNET_Disease2Gene_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(DisGeNET_Disease2Gene_lib, "InputFiles/Enrichment_analysis_libraries/DisGeNET_Disease2Gene_lib.rds")



### Enrichment of DisGeNET phenotype
DisGeNET_Gene_Phenotype <- na.exclude(DisGeNET_data[DisGeNET_data$diseaseType == "phenotype",c("ensembl_gene_id", "diseaseId", "diseaseName")])
DisGeNET_Gene_Phenotype$diseaseName <- paste0(DisGeNET_Gene_Phenotype$diseaseName, " (", DisGeNET_Gene_Phenotype$diseaseId, ")")

DisGeNET_Phenotype2Gene_lib <- split(x = DisGeNET_Gene_Phenotype$ensembl_gene_id, f = DisGeNET_Gene_Phenotype$diseaseName)
DisGeNET_Phenotype2Gene_lib <- DisGeNET_Phenotype2Gene_lib[lengths(DisGeNET_Phenotype2Gene_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(DisGeNET_Phenotype2Gene_lib, "InputFiles/Enrichment_analysis_libraries/DisGeNET_Phenotype2Gene_lib.rds")



### Enrichment of DisGeNET groups
DisGeNET_Gene_DiseaseGroup <- na.exclude(DisGeNET_data[DisGeNET_data$diseaseType == "group",c("ensembl_gene_id", "diseaseId", "diseaseName")])
DisGeNET_Gene_DiseaseGroup$diseaseName <- paste0(DisGeNET_Gene_DiseaseGroup$diseaseName, " (", DisGeNET_Gene_DiseaseGroup$diseaseId, ")")

DisGeNET_DiseaseGroup2Gene_lib <- split(x = DisGeNET_Gene_DiseaseGroup$ensembl_gene_id, f = DisGeNET_Gene_DiseaseGroup$diseaseName)
DisGeNET_DiseaseGroup2Gene_lib <- DisGeNET_DiseaseGroup2Gene_lib[lengths(DisGeNET_DiseaseGroup2Gene_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(DisGeNET_DiseaseGroup2Gene_lib, "InputFiles/Enrichment_analysis_libraries/DisGeNET_DiseaseGroup2Gene_lib.rds")





## Enrichment of OpenTargets diseases  -------------------------------------------------------------  
# Different types explained: https://platform-docs.opentargets.org/evidence


# Download data using wget in terminal
if(!dir.exists("Databases/OpenTargets/")){dir.create("Databases/OpenTargets/", recursive = TRUE)}

if(!dir.exists("Databases/OpenTargets/associationByDatatypeDirect")){
  system("wget --recursive --no-parent --no-host-directories -P Databases/OpenTargets/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/23.09/output/etl/parquet/associationByDatatypeDirect
")
}

if(!dir.exists("Databases/OpenTargets/diseases")){
  system("wget --recursive --no-parent --no-host-directories -P Databases/OpenTargets/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/23.09/output/etl/parquet/diseases")
}


sc <- spark_connect(master = "local")


# Read the disease info
OpenTargets_Diseases <- spark_read_parquet(sc, "OpenTargets_Diseases", "Databases/OpenTargets/diseases/")
OpenTargets_Diseases <- OpenTargets_Diseases %>% select(id, name, description)
OpenTargets_Diseases <- OpenTargets_Diseases %>% collect()

# Read the drug-gene associations
OpenTargets_Target_Disease <- spark_read_parquet(sc, "OpenTargets_Target_Disease", "Databases/OpenTargets/associationByDatatypeDirect/")


## Disease gene enrichment library based genetic association (GA) 
OpenTargets_Target_Disease_GA <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "genetic_association")

OpenTargets_Target_Disease_GA <- OpenTargets_Target_Disease_GA %>% collect()  
OpenTargets_Target_Disease_GA <- OpenTargets_Target_Disease_GA[OpenTargets_Target_Disease_GA$score > quantile(unique(OpenTargets_Target_Disease_GA$score), 0.9),]

OpenTargets_Target_Disease_GA$diseaseLabel  <- OpenTargets_Diseases$name[match(OpenTargets_Target_Disease_GA$diseaseId, OpenTargets_Diseases$id)]
OpenTargets_Target_Disease_GA$diseaseLabel <- paste0(OpenTargets_Target_Disease_GA$diseaseLabel, " (", OpenTargets_Target_Disease_GA$diseaseId, ")")
OpenTargets_Disease2Gene_GA_lib <- split(OpenTargets_Target_Disease_GA$targetId, OpenTargets_Target_Disease_GA$diseaseLabel)

OpenTargets_Disease2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[lengths(OpenTargets_Disease2Gene_GA_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_GA_lib, "InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_GA_lib.rds")


## Disease gene enrichment library based RNA expression score 
OpenTargets_Target_Disease_RNA <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "rna_expression") 

OpenTargets_Target_Disease_RNA <- OpenTargets_Target_Disease_RNA %>% collect()  
OpenTargets_Target_Disease_RNA <- OpenTargets_Target_Disease_RNA[OpenTargets_Target_Disease_RNA$score > quantile(unique(OpenTargets_Target_Disease_RNA$score), 0.9),]

OpenTargets_Target_Disease_RNA$diseaseLabel  <- OpenTargets_Diseases$name[match(OpenTargets_Target_Disease_RNA$diseaseId, OpenTargets_Diseases$id)]
OpenTargets_Target_Disease_RNA$diseaseLabel <- paste0(OpenTargets_Target_Disease_RNA$diseaseLabel, " (", OpenTargets_Target_Disease_RNA$diseaseId, ")")
OpenTargets_Disease2Gene_RNA_lib <- split(OpenTargets_Target_Disease_RNA$targetId, OpenTargets_Target_Disease_RNA$diseaseLabel)

OpenTargets_Disease2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[lengths(OpenTargets_Disease2Gene_RNA_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_RNA_lib, "InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_RNA_lib.rds")


## Disease gene enrichment library based literature 
OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "literature") 

OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease_lit %>% collect()  
OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease_lit[OpenTargets_Target_Disease_lit$score > quantile(unique(OpenTargets_Target_Disease_lit$score), 0.9),]

OpenTargets_Target_Disease_lit$diseaseLabel  <- OpenTargets_Diseases$name[match(OpenTargets_Target_Disease_lit$diseaseId, OpenTargets_Diseases$id)]
OpenTargets_Target_Disease_lit$diseaseLabel <- paste0(OpenTargets_Target_Disease_lit$diseaseLabel, " (", OpenTargets_Target_Disease_lit$diseaseId, ")")
OpenTargets_Disease2Gene_lit_lib <- split(OpenTargets_Target_Disease_lit$targetId, OpenTargets_Target_Disease_lit$diseaseLabel)

OpenTargets_Disease2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[lengths(OpenTargets_Disease2Gene_lit_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_lit_lib, "InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_lit_lib.rds")


## Disease gene enrichment library based on somatic mutation 
OpenTargets_Target_Disease_SM <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "somatic_mutation") 

OpenTargets_Target_Disease_SM <- OpenTargets_Target_Disease_SM %>% collect()  
OpenTargets_Target_Disease_SM <- OpenTargets_Target_Disease_SM[OpenTargets_Target_Disease_SM$score > quantile(unique(OpenTargets_Target_Disease_SM$score), 0.9),]

OpenTargets_Target_Disease_SM$diseaseLabel  <- OpenTargets_Diseases$name[match(OpenTargets_Target_Disease_SM$diseaseId, OpenTargets_Diseases$id)]
OpenTargets_Target_Disease_SM$diseaseLabel <- paste0(OpenTargets_Target_Disease_SM$diseaseLabel, " (", OpenTargets_Target_Disease_SM$diseaseId, ")")
OpenTargets_Disease2Gene_SM_lib <- split(OpenTargets_Target_Disease_SM$targetId, OpenTargets_Target_Disease_SM$diseaseLabel)

OpenTargets_Disease2Gene_SM_lib <- OpenTargets_Disease2Gene_SM_lib[lengths(OpenTargets_Disease2Gene_SM_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(OpenTargets_Disease2Gene_SM_lib, "InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_SM_lib.rds")


spark_disconnect(sc)





## Enrichment of cancer driver genes from intogen.org/ -------------------------------------------------------------



# Retrieve cancer driver gene list from https://intogen.org/download
# Download and unzip files from: 
# (a) https://intogen.org/download?file=IntOGen-Drivers-20200201.zip
# (b) https://intogen.org/download?file=IntOGen-Cohorts-20200201.zip

if(!dir.exists("Databases/Intogen")){dir.create("Databases/Intogen", recursive = TRUE)}

if(!file.exists("Databases/Intogen/Compendium_Cancer_Genes.tsv")){
  download.file(url = "https://www.intogen.org/download?file=IntOGen-Drivers-20230531.zip",
                destfile = "Databases/Intogen/IntOGen-Drivers-20230531.zip", method = "wget")
  unzip("Databases/Intogen/IntOGen-Drivers-20230531.zip", files = "2023-05-31_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", 
        exdir = "Databases/Intogen/", junkpaths = TRUE)
}

if(!file.exists("Databases/Intogen/cohorts.tsv")){
  download.file(url = "https://www.intogen.org/download?file=IntOGen-Cohorts-20230531.zip",
                destfile = "Databases/Intogen/IntOGen-Cohorts-20230531.zip", method = "wget")
  unzip("Databases/Intogen/IntOGen-Cohorts-20230531.zip", files = "2023-05-31_IntOGen-Cohorts/cohorts.tsv", 
        exdir = "Databases/Intogen/", junkpaths = TRUE)
}



# Read the cancer types in Intogen
Intogen_cancers <- read.table("Databases/Intogen/cohorts.tsv", sep = "\t", header = TRUE, check.names = FALSE, quote = "")
Intogen_cancers <- unique(Intogen_cancers[, c("CANCER", "CANCER_NAME")])


# Read disease gene association data
Intogen_data <- read.table("Databases/Intogen/Compendium_Cancer_Genes.tsv", sep = "\t", header = TRUE, check.names = FALSE)


# Add ensembl gene IDs
Intogen_data$ensembl_gene_id <- geneSymbol_2_ensemblId$ensembl_id[match(Intogen_data$SYMBOL, geneSymbol_2_ensemblId$symbol)]
Intogen_data <- Intogen_data[!is.na(Intogen_data$ensembl_gene_id),]

Intogen_data$CANCER_NAME <- Intogen_cancers$CANCER_NAME[match(Intogen_data$CANCER_TYPE, Intogen_cancers$CANCER)]
Intogen_data$CANCER_TYPE <- paste0(Intogen_data$CANCER_NAME, " (", Intogen_data$CANCER_TYPE, ")")

Intogen_Disease2Gene_lib <- split(Intogen_data$ensembl_gene_id, f = Intogen_data$CANCER_TYPE)
Intogen_Disease2Gene_lib <- lapply(Intogen_Disease2Gene_lib, unique)

Intogen_Disease2Gene_lib <- Intogen_Disease2Gene_lib[lengths(Intogen_Disease2Gene_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(Intogen_Disease2Gene_lib, "InputFiles/Enrichment_analysis_libraries/Intogen_Disease2Gene_lib.rds")






# Enrichment of Enrichr diseases ---------------------------------------------------------



## Disease gene enrichment library based on Disease_Signatures_from_GEO_2014 (Enrichr)

# Download libraries from Enrichr
if(!dir.exists("Databases/Enrichr")){dir.create("Databases/Enrichr", recursive = TRUE)}
if(!file.exists("Databases/Enrichr/Disease_Signatures_from_GEO_down_2014.txt")){
  download.file(url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Disease_Signatures_from_GEO_down_2014",
                destfile = "Databases/Enrichr/Disease_Signatures_from_GEO_down_2014.txt", method = "wget")
  download.file(url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Disease_Signatures_from_GEO_up_2014",
                destfile = "Databases/Enrichr/Disease_Signatures_from_GEO_up_2014.txt", method = "wget")
}


Enrichr_GeoDiseaseSignatures_Up <- read.table("Databases/Enrichr/Disease_Signatures_from_GEO_up_2014.txt", sep = "\t", quote = "", fill = TRUE)
result <- list()
for(i in 1:nrow(Enrichr_GeoDiseaseSignatures_Up)){
  tmp1 <- Enrichr_GeoDiseaseSignatures_Up[i, -c(1:2)]
  tmp2 <- as.character(apply(tmp1, 1, function(x)gsub(pattern = ",[0-9]+.[0-9]+", replacement = "", x = x)))
  tmp2 <- tmp2[!is.na(tmp2)]
  result[[i]] <- geneSymbol_2_ensemblId[geneSymbol_2_ensemblId$symbol %in% tmp2,]$ensembl_id 
  names(result)[i] <- Enrichr_GeoDiseaseSignatures_Up[i,1]
}
Enrichr_GeoDiseaseSignatures_Up <- result
Enrichr_GeoDiseaseSignatures_Up <- Enrichr_GeoDiseaseSignatures_Up[lengths(Enrichr_GeoDiseaseSignatures_Up) >= 5]

Enrichr_GeoDiseaseSignatures_Down <- read.table("Databases/Enrichr/Disease_Signatures_from_GEO_down_2014.txt", sep = "\t", quote = "", fill = TRUE)
result <- list()
for(i in 1:nrow(Enrichr_GeoDiseaseSignatures_Down)){
  tmp1 <- Enrichr_GeoDiseaseSignatures_Down[i, -c(1:2)]
  tmp2 <- as.character(apply(tmp1, 1, function(x)gsub(pattern = ",[0-9]+.[0-9]+", replacement = "", x = x)))
  tmp2 <- tmp2[!is.na(tmp2)]
  result[[i]] <- geneSymbol_2_ensemblId[geneSymbol_2_ensemblId$symbol %in% tmp2,]$ensembl_id 
  names(result)[i] <- Enrichr_GeoDiseaseSignatures_Down[i,1]
}
Enrichr_GeoDiseaseSignatures_Down <- result
Enrichr_GeoDiseaseSignatures_Down <- Enrichr_GeoDiseaseSignatures_Down[lengths(Enrichr_GeoDiseaseSignatures_Down) >= 5]

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- list()
Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up <- Enrichr_GeoDiseaseSignatures_Up
Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down <- Enrichr_GeoDiseaseSignatures_Down

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(Enrichr_Disease2Gene_GeoDiseaseSig_lib, "InputFiles/Enrichment_analysis_libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")





# Enrichment of PharmGKB diseases ---------------------------------------------------------
if(!dir.exists("Databases/PharmGKB/")){dir.create("Databases/PharmGKB/", recursive = TRUE)}
if(!file.exists("Databases/PharmGKB/relationships.tsv")){
  download.file(url = "https://api.pharmgkb.org/v1/download/file/data/relationships.zip",
                destfile = "Databases/PharmGKB/relationships.zip", method = "wget")
  unzip("Databases/PharmGKB/relationships.zip", files = "relationships.tsv", exdir = "Databases/PharmGKB/")
}

PharmGKB_data <- read.table("Databases/PharmGKB/relationships.tsv", sep = "\t", fill = TRUE, header = TRUE, quote = "")
PharmGKB_Gene_Disease <- PharmGKB_data[(PharmGKB_data$Entity1_type == "Gene" & PharmGKB_data$Entity2_type == "Disease" & PharmGKB_data$Association == "associated"),]

# Map PharmGKB gene ID to Ensembl gene ID
tmp <- func_pharmgkb_mapping(PharmGKB_Gene_Disease$Entity1_id, map_to = "ensembl_gene_id")
PharmGKB_Gene_Disease <- merge(PharmGKB_Gene_Disease,  tmp, by.x = c("Entity1_id"), by.y = "Entity_id", all.x = TRUE)
names(PharmGKB_Gene_Disease)[names(PharmGKB_Gene_Disease) == "ensembl_gene_id"] <- "Entity1_ensembl_gene_id"

PharmGKB_Gene_Disease <- PharmGKB_Gene_Disease[,c("Entity1_ensembl_gene_id", "Entity2_name")]
PharmGKB_Gene_Disease[PharmGKB_Gene_Disease == ""] <- NA
PharmGKB_Gene_Disease <- na.exclude(PharmGKB_Gene_Disease)

PharmGKB_Disease2Gene_lib <- split(x = PharmGKB_Gene_Disease$Entity1_ensembl_gene_id, f = PharmGKB_Gene_Disease$Entity2_name)

PharmGKB_Disease2Gene_lib <- PharmGKB_Disease2Gene_lib[lengths(PharmGKB_Disease2Gene_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(PharmGKB_Disease2Gene_lib, "InputFiles/Enrichment_analysis_libraries/PharmGKB_Disease2Gene_lib.rds")





## Enrichment of Comparative Toxicogenimics Database (CTD) diseases -------------------------------------------------------------

# Download and read the gene disease association file from CTD
if(!dir.exists("Databases/ComparativeToxicogenomicsDatabase/")){dir.create("Databases/ComparativeToxicogenomicsDatabase/", recursive = TRUE)}
if(!file.exists("Databases/ComparativeToxicogenomicsDatabase/CTD_genes_diseases.csv.gz")){
  download.file(url = "http://ctdbase.org/reports/CTD_genes_diseases.csv.gz",
                destfile = "Databases/ComparativeToxicogenomicsDatabase/CTD_genes_diseases.csv.gz", method = "wget")
}

CTD_Gene_Disease <- read.csv(gzfile("Databases/ComparativeToxicogenomicsDatabase/CTD_genes_diseases.csv.gz"), header = FALSE, fill = TRUE, skip = 29)
colnames(CTD_Gene_Disease) <- c("GeneSymbol", "GeneID", "DiseaseName", "DiseaseID", "DirectEvidence", "InferenceChemicalName", "InferenceScore", "OmimIDs", "PubMedIDs")

# Filter to keep only direct associations
CTD_Gene_Disease <- CTD_Gene_Disease[CTD_Gene_Disease$DirectEvidence %in% c("marker/mechanism", "therapeutic", "marker/mechanism|therapeutic"), ]


# Map entrez gene IDs to ensembl IDs
CTD_Gene_Disease$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(CTD_Gene_Disease$GeneID, entrezId_2_ensemblId$gene_id)]
CTD_Gene_Disease <- CTD_Gene_Disease[!is.na(CTD_Gene_Disease$ensembl_gene_id),]


CTD_Gene_Disease$DiseaseName <- paste0(CTD_Gene_Disease$DiseaseName, " (", CTD_Gene_Disease$DiseaseID, ")")
CTD_Disease2Gene_lib <- split(x = CTD_Gene_Disease$ensembl_gene_id, f = CTD_Gene_Disease$DiseaseName)

CTD_Disease2Gene_lib <- CTD_Disease2Gene_lib[lengths(CTD_Disease2Gene_lib) >= 5]

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(CTD_Disease2Gene_lib, "InputFiles/Enrichment_analysis_libraries/CTD_Disease2Gene_lib.rds")



print(warnings())