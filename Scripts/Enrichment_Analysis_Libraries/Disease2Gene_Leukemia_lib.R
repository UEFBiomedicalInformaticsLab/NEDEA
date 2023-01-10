set.seed(5081)
rm(list = ls())



# Enrichment analysis libraries for Leukemia


# Notes:
# NCBI MedGene Concept ID: C0023418
# Experimental Factor Ontology: EFO_0000565 
# Mondo Disease Ontology: MONDO_0005059 
# For Enricher, used Leukemia as search term and retrieved multiple types of Leukemia gene sets
# For TheTA used EFO:0000339 (chronic myelogenous leukemia) as leukemia not found

# Load libraries
library(org.Hs.eg.db)

# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_Leukemia2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C0023418", x = names(DisGeNET_Disease2Gene_lib), ignore.case = TRUE)]
names(DisGeNET_Leukemia2Gene_lib) <- paste0(names(DisGeNET_Leukemia2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_Leukemia2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "EFO_0000565|MONDO_0005059", x = names(OpenTargets_Disease2Gene_GA_lib), ignore.case = TRUE)]
names(OpenTargets_Leukemia2Gene_GA_lib) <- paste0(names(OpenTargets_Leukemia2Gene_GA_lib), "[OpenTargets_GA]")

# Not found
# OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
# OpenTargets_Leukemia2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "EFO_0000565|MONDO_0005059", x = names(OpenTargets_Disease2Gene_RNA_lib), ignore.case = TRUE)]
# names(OpenTargets_Leukemia2Gene_RNA_lib) <- paste0(names(OpenTargets_Leukemia2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_Leukemia2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "EFO_0000565|MONDO_0005059", x = names(OpenTargets_Disease2Gene_lit_lib), ignore.case = TRUE)]
names(OpenTargets_Leukemia2Gene_lit_lib) <- paste0(names(OpenTargets_Leukemia2Gene_lit_lib), "[OpenTargets_literature]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoLeukemiaSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "Leukemia", 
                                                                                     x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), ignore.case = TRUE)]
names(Enrichr_GeoLeukemiaSignatures_Up) <- paste0(names(Enrichr_GeoLeukemiaSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoLeukemiaSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "Leukemia", 
                                                                                         x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), ignore.case = TRUE)]
names(Enrichr_GeoLeukemiaSignatures_Down) <- paste0(names(Enrichr_GeoLeukemiaSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")

## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_Leukemia2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "Leukemia",
                                                             x = names(Intogen_Disease2Gene_lib), ignore.case = TRUE)]
names(Intogen_Leukemia2Gene_lib) <- paste0(names(Intogen_Leukemia2Gene_lib), "[Intogen]")


## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")

PharmGKB_Leukemia2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "^leukemia", 
                                                                  x = names(PharmGKB_Disease2Gene_lib), ignore.case = TRUE)]
names(PharmGKB_Leukemia2Gene_lib) <- paste0(names(PharmGKB_Leukemia2Gene_lib), "[PharmGKB_associated]")


## ThETA
library(ThETA)
source("ExternalTools/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
Leukemia_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0000339")]]
Leukemia_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0000339"),]
Leukemia_Tscores <- tissue.specific.scores(disease_genes = Leukemia_genes$entrez, 
                                             ppi_network = ppi_strdb_700, 
                                             directed_network = FALSE, 
                                             tissue_expr_data = gtexv7_zscore,
                                             dis_relevant_tissues = Leukemia_rel_tissue_scores, 
                                             W = centrality_score$borda.disc, selected_tissues = NULL, 
                                             cutoff = 4, verbose = TRUE)

ThETA_Leukemia2Gene_lib <- row.names(Leukemia_Tscores[order(Leukemia_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_Leukemia2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_Leukemia2Gene_lib, c("ensembl_id")])
names(ThETA_Leukemia2Gene_lib) <- paste0("leukemia (EFO_0000339)", "[ThETA]")


# Merge and save into one

Enrichment_Leukemia2Gene_lib <- c(DisGeNET_Leukemia2Gene_lib, OpenTargets_Leukemia2Gene_GA_lib, 
                                    OpenTargets_Leukemia2Gene_lit_lib, 
                                    Enrichr_GeoLeukemiaSignatures_Up, Enrichr_GeoLeukemiaSignatures_Down,
                                    Intogen_Leukemia2Gene_lib,
                                    PharmGKB_Leukemia2Gene_lib, ThETA_Leukemia2Gene_lib)



# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_Leukemia2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_Leukemia_lib.rds")

print(warnings())