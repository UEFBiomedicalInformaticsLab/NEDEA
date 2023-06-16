set.seed(5081)



# Enrichment analysis libraries for Lung Cancer


# Notes:
# Includes Lung Cancer and its two major types non-small cell lung cancer and small cell lung cancer [Ref: https://www.cancer.gov/types/stomach].
# NCBI MedGene Concept ID: C0684249, C0007131, C0149925
# Experimental Factor Ontology: EFO_0001071, EFO_0003060, EFO_0000702
# Mondo Disease Ontology: MONDO_0008903, MONDO_0005233, MONDO_0008433
# For Enricher, used Adenocarcinoma of lung. This is the closest available term





# Load libraries
library(org.Hs.eg.db)





# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")




## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_LungCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C0684249|C0007131|C0149925", x = names(DisGeNET_Disease2Gene_lib), ignore.case = TRUE)]
names(DisGeNET_LungCancer2Gene_lib) <- paste0(names(DisGeNET_LungCancer2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_LungCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "EFO_0001071|EFO_0003060|EFO_0000702|MONDO_0008903|MONDO_0005233|MONDO_0008433", x = names(OpenTargets_Disease2Gene_GA_lib), ignore.case = TRUE)]
names(OpenTargets_LungCancer2Gene_GA_lib) <- paste0(names(OpenTargets_LungCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_LungCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "EFO_0001071|EFO_0003060|EFO_0000702|MONDO_0008903|MONDO_0005233|MONDO_0008433", x = names(OpenTargets_Disease2Gene_RNA_lib), ignore.case = TRUE)]
names(OpenTargets_LungCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_LungCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_LungCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "EFO_0001071|EFO_0003060|EFO_0000702|MONDO_0008903|MONDO_0005233|MONDO_0008433", x = names(OpenTargets_Disease2Gene_lit_lib), ignore.case = TRUE)]
names(OpenTargets_LungCancer2Gene_lit_lib) <- paste0(names(OpenTargets_LungCancer2Gene_lit_lib), "[OpenTargets_literature]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoLungCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "Adenocarcinoma of lung", 
                                                                                     x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), ignore.case = TRUE)]
names(Enrichr_GeoLungCancerSignatures_Up) <- paste0(names(Enrichr_GeoLungCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoLungCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "Adenocarcinoma of lung", 
                                                                                         x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), ignore.case = TRUE)]
names(Enrichr_GeoLungCancerSignatures_Down) <- paste0(names(Enrichr_GeoLungCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")

## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_LungCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "Lung",
                                                             x = names(Intogen_Disease2Gene_lib), ignore.case = TRUE)]
names(Intogen_LungCancer2Gene_lib) <- paste0(names(Intogen_LungCancer2Gene_lib), "[Intogen]")


## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")

PharmGKB_LungCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "Carcinoma, Non-Small-Cell Lung", 
                                                               x = names(PharmGKB_Disease2Gene_lib), ignore.case = TRUE)]
names(PharmGKB_LungCancer2Gene_lib) <- paste0(names(PharmGKB_LungCancer2Gene_lib), "[PharmGKB_associated]")


## ThETA
library(ThETA)
source("Scripts/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
LungCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0001071")]]
LungCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0001071"),]
LungCancer_Tscores <- tissue.specific.scores(disease_genes = LungCancer_genes$entrez, 
                                             ppi_network = ppi_strdb_700, 
                                             directed_network = FALSE, 
                                             tissue_expr_data = gtexv7_zscore,
                                             dis_relevant_tissues = LungCancer_rel_tissue_scores, 
                                             W = centrality_score$borda.disc, selected_tissues = NULL, 
                                             cutoff = 4, verbose = TRUE)

ThETA_LungCancer2Gene_lib <- row.names(LungCancer_Tscores[order(LungCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_LungCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_LungCancer2Gene_lib, c("ensembl_id")])
names(ThETA_LungCancer2Gene_lib) <- paste0("lung carcinoma (EFO_0001071)", "[ThETA]")



# Merge and save into one

Enrichment_LungCancer2Gene_lib <- c(DisGeNET_LungCancer2Gene_lib, OpenTargets_LungCancer2Gene_GA_lib, 
                                    OpenTargets_LungCancer2Gene_RNA_lib, OpenTargets_LungCancer2Gene_lit_lib, 
                                    Enrichr_GeoLungCancerSignatures_Up, Enrichr_GeoLungCancerSignatures_Down,
                                    Intogen_LungCancer2Gene_lib,
                                    PharmGKB_LungCancer2Gene_lib, ThETA_LungCancer2Gene_lib)




# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_LungCancer2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_LungCancer_lib.rds")



print(warnings())