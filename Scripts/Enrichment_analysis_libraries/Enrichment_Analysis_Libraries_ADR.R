set.seed(5081)



# Libraries for enrichment analysis of side-effects



# Load libraries
library(unixtools)
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
  download.file(url = "http://www.bio-add.org/ADReCS-Target/files/download/ADRAlert2GENE2ID.xlsx",
                destfile = "Databases/ADReCS/ADRAlert2GENE2ID.xlsx", method = "wget")
}

ADReCS_Gene_ADR <- read.xlsx("Databases/ADReCS/ADRAlert2GENE2ID.xlsx")
ADReCS_Gene_ADR$ADR.Term <- paste0(ADReCS_Gene_ADR$ADR.Term, " (", ADReCS_Gene_ADR$ADR.ID, ")")
# ADRid_2_ADRterm <- unique(ADReCS_Gene_ADR[, c("ADR.ID", "ADR.Term")])
ADReCS_Gene_ADR$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(ADReCS_Gene_ADR$GeneID, entrezId_2_ensemblId$gene_id)]
ADReCS_Gene_ADR <- na.exclude(ADReCS_Gene_ADR)


### Level 4 ADR IDs library
ADReCS_Gene_ADR_level4 <- ADReCS_Gene_ADR[lengths(strsplit(ADReCS_Gene_ADR$ADR.ID, split = "\\."))==4,]

ADReCS_ADR2Gene_level4_lib <- split(x = ADReCS_Gene_ADR_level4$ensembl_gene_id, f = ADReCS_Gene_ADR_level4$ADR.Term)

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(ADReCS_ADR2Gene_level4_lib, "InputFiles/Enrichment_analysis_libraries/ADReCS_ADR2Gene_level4_lib.rds")

### Level 3 ADR IDs library
ADReCS_Gene_ADR_level3 <- ADReCS_Gene_ADR[lengths(strsplit(ADReCS_Gene_ADR$ADR.ID, split = "\\."))==3,]
ADReCS_ADR2Gene_level3_lib <- split(x = ADReCS_Gene_ADR_level3$ensembl_gene_id, f = ADReCS_Gene_ADR_level3$ADR.Term)

if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(ADReCS_ADR2Gene_level3_lib, "InputFiles/Enrichment_analysis_libraries/ADReCS_ADR2Gene_level3_lib.rds")



print(warnings())