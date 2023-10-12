set.seed(5081)



# Libraries for enrichment analysis of pathways



# Load libraries
library(unixtools)
library(org.Hs.eg.db)
# library(msigdbr)
library(httr)
library(readxl)
library(tidyverse)




# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





# ## Enrichment of KEGG pathway -------------------------------------------------------------
# # Refer for categories and sub-categories: http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
# msigdb_Gene_keggPath <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
# msigdb_keggPath2Gene_lib <- split(x = msigdb_Gene_keggPath$ensembl_gene, f = msigdb_Gene_keggPath$gs_exact_source)
# 
# pathwayName_2_pathwayId <- unique(as.data.frame(msigdb_Gene_keggPath[, c("gs_name", "gs_exact_source")]))
# for(i in 1:length(msigdb_keggPath2Gene_lib)){
#   pathwayName <- paste0(pathwayName_2_pathwayId[pathwayName_2_pathwayId$gs_exact_source == names(msigdb_keggPath2Gene_lib)[i],]$gs_name, " (", 
#                         pathwayName_2_pathwayId[pathwayName_2_pathwayId$gs_exact_source == names(msigdb_keggPath2Gene_lib)[i],]$gs_exact_source, ")")
#   names(msigdb_keggPath2Gene_lib)[i] <- pathwayName
# }
# if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
# saveRDS(msigdb_keggPath2Gene_lib, "InputFiles/Enrichment_analysis_libraries/msigdb_keggPath2Gene_lib.rds")






# ## Enrichment of KEGG modules -------------------------------------------------------------
# res <- GET("http://rest.kegg.jp/link/module/hsa")
# kegg_gene_keggModule <- read.table(text = httr::content(res), col.names = c("entrez_gene_id", "kegg_module"))
# kegg_gene_keggModule$entrez_gene_id <- gsub(pattern = "^hsa:", x = kegg_gene_keggModule$entrez_gene_id, replacement = "")
# kegg_gene_keggModule$kegg_module <- gsub(pattern = "^md:", x = kegg_gene_keggModule$kegg_module, replacement = "")
# 
# kegg_gene_keggModule$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(kegg_gene_keggModule$entrez_gene_id, entrezId_2_ensemblId$gene_id)]
# kegg_module2gene_lib <- split(x = kegg_gene_keggModule$ensembl_gene_id, f = kegg_gene_keggModule$kegg_module)
# 
# # res <- GET("http://rest.kegg.jp/list/module") # Unable to retrieve the human specific module ids
# # moduleId_2_moduleName <- read.table(text = content(res), col.names = c("kegg_module_id", "kegg_module_name"), sep = "\t")
# # for(i in 1:length(kegg_module2gene_lib)){
# #   moduleName <- paste0(moduleId_2_moduleName[moduleId_2_moduleName$kegg_module_id == names(kegg_module2gene_lib)[i],]$kegg_module_name, " (", 
# #                         moduleId_2_moduleName[moduleId_2_moduleName$kegg_module_id == names(kegg_module2gene_lib)[i],]$kegg_module_id, ")")
# #   names(kegg_module2gene_lib)[i] <- moduleName
# # }
# if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
# saveRDS(kegg_module2gene_lib, "InputFiles/Enrichment_analysis_libraries/kegg_module2gene_lib.rds")





## Enrichment of hallmark of cancer related KEGG pathways -------------------------------------------------------------

# Download the HOC associated pathways
if(!dir.exists("Databases/CHG/")){dir.create("Databases/CHG/", recursive = TRUE)}
if(!file.exists("Databases/CHG/Table_1.xls")){
  download.file(url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7013921/bin/Table_1.xls",
                destfile = "Databases/CHG/Table_1.xls", method = "wget")
}

CHG_pathways <- read_excel("Databases/CHG/Table_1.xls", skip = 2)
CHG_pathways <- CHG_pathways$Pathway
CHG_pathways <- sort(unique(unlist(str_split(CHG_pathways, "\n"))))
CHG_pathways <- CHG_pathways[CHG_pathways != ""]

# Create mapping from Entrez gene ID to Ensembl gene ID
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)

# Create KEGG pathway ID to pathway name mapping
keggId_2_keggName <- GET("https://rest.kegg.jp/list/pathway/hsa")
keggId_2_keggName <- content(keggId_2_keggName)
keggId_2_keggName <- data.frame(str_split(keggId_2_keggName, pattern = "\n"))
colnames(keggId_2_keggName) <- "kegg_pathway_id"
keggId_2_keggName <- na.exclude(separate(keggId_2_keggName, col = "kegg_pathway_id", into = c("kegg_pathway_id", "kegg_pathway_name"), sep = "\t"))



# Retrieve genes for the pathways
CHG_keggPath2Gene_lib <- list()
for(pathway in CHG_pathways){
  print(paste0("-- retrieving for ", pathway))
  response <- GET(paste0("https://rest.kegg.jp/link/hsa/", pathway))
  response <- content(response)
  response <- data.frame(str_split(response, pattern = "\n"))
  colnames(response) <- "kegg_pathway_id"
  response <- na.exclude(separate(response, col = "kegg_pathway_id", into = c("kegg_pathway_id", "entrez_gene_id"), sep = "\t"))
  response$kegg_pathway_id <- gsub("^path:", "", response$kegg_pathway_id)
  response$entrez_gene_id <- gsub("^hsa:", "", response$entrez_gene_id)
  response$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(response$entrez_gene_id, entrezId_2_ensemblId$gene_id)]
  tmp <- paste0(keggId_2_keggName[keggId_2_keggName$kegg_pathway_id == pathway, "kegg_pathway_name"], " (", pathway, ")")
  tmp <- gsub("* - Homo sapiens \\(human\\)*", "", tmp)
  CHG_keggPath2Gene_lib[[tmp]] <- sort(na.exclude(unique(response$ensembl_gene_id)))
}
if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(CHG_keggPath2Gene_lib, "InputFiles/Enrichment_analysis_libraries/CHG_keggPath2Gene_lib.rds")





# ## Enrichment of REACTOME pathway -------------------------------------------------------------
# # Refer for categories and sub-categories: http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp
# msigdb_Gene_ReactomePath <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
# msigdb_ReactomePath2Gene_lib <- split(x = msigdb_Gene_ReactomePath$ensembl_gene, f = msigdb_Gene_ReactomePath$gs_exact_source)
# 
# pathwayName_2_pathwayId <- unique(as.data.frame(msigdb_Gene_ReactomePath[, c("gs_name", "gs_exact_source")]))
# for(i in 1:length(msigdb_ReactomePath2Gene_lib)){
#   pathwayName <- paste0(pathwayName_2_pathwayId[pathwayName_2_pathwayId$gs_exact_source == names(msigdb_ReactomePath2Gene_lib)[i],]$gs_name, " (", 
#                         pathwayName_2_pathwayId[pathwayName_2_pathwayId$gs_exact_source == names(msigdb_ReactomePath2Gene_lib)[i],]$gs_exact_source, ")")
#   names(msigdb_ReactomePath2Gene_lib)[i] <- pathwayName
# }
# if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
# saveRDS(msigdb_ReactomePath2Gene_lib, "InputFiles/Enrichment_analysis_libraries/msigdb_ReactomePath2Gene_lib.rds")





## Enrichment of SMPDB pathway -------------------------------------------------------------

if(!dir.exists("Databases/SMPDB")){dir.create("Databases/SMPDB", recursive = TRUE)}

if(!file.exists("Databases/SMPDB/smpdb_proteins.csv.zip")){
  download.file(url = "https://smpdb.ca/downloads/smpdb_proteins.csv.zip",
                destfile = "Databases/SMPDB/smpdb_proteins.csv.zip", method = "wget")
  unzip(zipfile = "Databases/SMPDB/smpdb_proteins.csv.zip", exdir = "Databases/SMPDB/smpdb_proteins")
}

if(!file.exists("Databases/SMPDB/smpdb_pathways.csv.zip")){
  download.file(url = "https://smpdb.ca/downloads/smpdb_pathways.csv.zip",
                destfile = "Databases/SMPDB/smpdb_pathways.csv.zip", method = "wget")
  unzip(zipfile = "Databases/SMPDB/smpdb_pathways.csv.zip", exdir = "Databases/SMPDB")
}

SMPDb_pathways <- read.csv("Databases/SMPDB/smpdb_pathways.csv", header = TRUE)
SMPDb_Pathway2Gene <- list()

for(pathway_subject in unique(SMPDb_pathways$Subject)){
  SMPDb_pathways_select <- SMPDb_pathways[SMPDb_pathways$Subject == pathway_subject, "SMPDB.ID"]
  file.list <- paste0("Databases/SMPDB/smpdb_proteins/", SMPDb_pathways_select, "_proteins.csv")
  
  for(path in file.list){
    if(file.exists(path)){
      file.content <- read.csv(file = path, header = TRUE)
      if(nrow(file.content) != 0){
        file.content$ensembl_gene_id <- geneSymbol_2_ensemblId$ensembl_id[match(file.content$`Gene.Name`, geneSymbol_2_ensemblId$symbol)]
        pathway <- paste0(unique(file.content$Pathway.Name), " (", unique(file.content$SMPDB.ID), ") (", pathway_subject, ")")
        SMPDb_Pathway2Gene[[pathway_subject]][[pathway]] <- as.vector(na.exclude(file.content$ensembl_gene_id))
        
      }else{
        print(paste0("Empty file: ", path))
      }
    }else{
      print(paste0("File not found: ", path))
    }
  }
}

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(SMPDb_Pathway2Gene, "InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds")



print(warnings())