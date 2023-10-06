set.seed(5081)



# Libraries for enrichment analysis of side-effects



# Load libraries
library(unixtools)
library(org.Hs.eg.db)
library(openxlsx)



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## Enrichment of ACReCS adverse drug reaction terms -------------------------------------------------------------
# Refer for ID info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383906/
# The ADR ID can be explained as hierchy. The data contains level 3 and 4 ADR IDs.



if(!dir.exists("Databases/ADReCS/")){dir.create("Databases/ADReCS/", recursive = TRUE)}
if(!file.exists("Databases/ADReCS/ADRAlert2GENE2ID.xlsx")){
  download.file(url = "http://bioinf.xmu.edu.cn/ADReCS-Target/files/download/ADRAlert2GENE2ID.xlsx",
                destfile = "Databases/ADReCS/ADRAlert2GENE2ID.xlsx", method = "wget")
}

ADReCS_Gene_ADR <- read.xlsx("Databases/ADReCS/ADRAlert2GENE2ID.xlsx")
ADRid_2_ADRterm <- unique(ADReCS_Gene_ADR[, c("ADR.ID", "ADR.Term")])
ADReCS_Gene_ADR$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(ADReCS_Gene_ADR$GeneID, entrezId_2_ensemblId$gene_id)]
ADReCS_Gene_ADR <- na.exclude(ADReCS_Gene_ADR)

### Level 4 ADR IDs library
ADReCS_Gene_ADR_level4 <- ADReCS_Gene_ADR[lengths(strsplit(ADReCS_Gene_ADR$ADR.ID, split = "\\."))==4,]
ADReCS_ADR2Gene_level4_lib <- list() # Create empty list
for(i in unique(ADReCS_Gene_ADR_level4$ADR.ID)) {             
  ADReCS_ADR2Gene_level4_lib[[i]] <- ADReCS_Gene_ADR_level4[ADReCS_Gene_ADR_level4$ADR.ID == i,]$ensembl_gene_id
}
for(i in 1:length(ADReCS_ADR2Gene_level4_lib)){
  adrName <- paste0(ADRid_2_ADRterm[ADRid_2_ADRterm$ADR.ID == names(ADReCS_ADR2Gene_level4_lib)[i],]$ADR.Term, " (", 
                    ADRid_2_ADRterm[ADRid_2_ADRterm$ADR.ID == names(ADReCS_ADR2Gene_level4_lib)[i],]$ADR.ID, ")")
  names(ADReCS_ADR2Gene_level4_lib)[i] <- adrName
}
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(ADReCS_ADR2Gene_level4_lib, "InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level4_lib.rds")


ADReCS_Gene_ADR_level3 <- ADReCS_Gene_ADR[lengths(strsplit(ADReCS_Gene_ADR$ADR.ID, split = "\\."))==3,]
ADReCS_ADR2Gene_level3_lib <- list() # Create empty list
for(i in unique(ADReCS_Gene_ADR_level3$ADR.ID)) {             
  ADReCS_ADR2Gene_level3_lib[[i]] <- ADReCS_Gene_ADR_level3[ADReCS_Gene_ADR_level3$ADR.ID == i,]$ensembl_gene_id
}
for(i in 1:length(ADReCS_ADR2Gene_level3_lib)){
  adrName <- paste0(ADRid_2_ADRterm[ADRid_2_ADRterm$ADR.ID == names(ADReCS_ADR2Gene_level3_lib)[i],]$ADR.Term, " (", 
                    ADRid_2_ADRterm[ADRid_2_ADRterm$ADR.ID == names(ADReCS_ADR2Gene_level3_lib)[i],]$ADR.ID, ")")
  names(ADReCS_ADR2Gene_level3_lib)[i] <- adrName
}
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(ADReCS_ADR2Gene_level3_lib, "InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level3_lib.rds")



print(warnings())