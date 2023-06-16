set.seed(5081)
rm(list = ls())



# Enrichment analysis libraries for ColonCancer


# Notes:
# Includes: colon carcinoma, colon adenocarcinoma, 
# colon small cell neuroendocrine carcinoma, squamous cell carcinoma of colon, 
# epithelial tumor of colon

# NCBI MedGene Concept ID: C0007102, C0699790, C0338106, C4016849

# Experimental Factor Ontology: "EFO_1001950", "EFO_1001949", 
# Mondo Disease Ontology: "MONDO_0003978", "MONDO_0018513", "MONDO_0024479"
# For Enrichr, used term 



# Load libraries
library(org.Hs.eg.db)

# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_ColonCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C0007102|C0699790|C0338106|C4016849", 
                                                                   x = names(DisGeNET_Disease2Gene_lib), 
                                                                   ignore.case = TRUE)]
names(DisGeNET_ColonCancer2Gene_lib) <- paste0(names(DisGeNET_ColonCancer2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_ColonCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "EFO_1001950|EFO_1001949|MONDO_0003978|MONDO_0018513|MONDO_0024479", 
                                                                               x = names(OpenTargets_Disease2Gene_GA_lib), 
                                                                               ignore.case = TRUE)]
names(OpenTargets_ColonCancer2Gene_GA_lib) <- paste0(names(OpenTargets_ColonCancer2Gene_GA_lib), "[OpenTargets_GA]")


# OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
# OpenTargets_ColonCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "EFO_1001950|EFO_1001949|MONDO_0003978|MONDO_0018513|MONDO_0024479", 
#                                                                                  x = names(OpenTargets_Disease2Gene_RNA_lib), 
#                                                                                  ignore.case = TRUE)]
# names(OpenTargets_ColonCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_ColonCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_ColonCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "EFO_1001950|EFO_1001949|MONDO_0003978|MONDO_0018513|MONDO_0024479", 
                                                                                 x = names(OpenTargets_Disease2Gene_lit_lib), 
                                                                                 ignore.case = TRUE)]
names(OpenTargets_ColonCancer2Gene_lit_lib) <- paste0(names(OpenTargets_ColonCancer2Gene_lit_lib), "[OpenTargets_literature]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoColonCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "Cancer of Colon", 
                                                                                         x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                                                         ignore.case = TRUE)]
names(Enrichr_GeoColonCancerSignatures_Up) <- paste0(names(Enrichr_GeoColonCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoColonCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "Cancer of Colon", 
                                                                                             x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                                             ignore.case = TRUE)]
names(Enrichr_GeoColonCancerSignatures_Down) <- paste0(names(Enrichr_GeoColonCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")


## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_ColonCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "Colorectal adenocarcinoma",
                                                                 x = names(Intogen_Disease2Gene_lib), 
                                                                 ignore.case = TRUE)]
names(Intogen_ColonCancer2Gene_lib) <- paste0(names(Intogen_ColonCancer2Gene_lib), "[Intogen]")


## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")

PharmGKB_ColonCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "Colonic Neoplasms|Colorectal Neoplasms", 
                                                                   x = names(PharmGKB_Disease2Gene_lib), 
                                                                   ignore.case = TRUE)]
names(PharmGKB_ColonCancer2Gene_lib) <- paste0(names(PharmGKB_ColonCancer2Gene_lib), "[PharmGKB_associated]")


## ThETA
library(ThETA)
source("ExternalTools/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
ColonCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:1001950")]]


ColonCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:1001950"),]

ColonCancer_Tscores <- tissue.specific.scores(disease_genes = ColonCancer_genes$entrez, 
                                                 ppi_network = ppi_strdb_700, 
                                                 directed_network = FALSE, 
                                                 tissue_expr_data = gtexv7_zscore,
                                                 dis_relevant_tissues = ColonCancer_rel_tissue_scores, 
                                                 W = centrality_score$borda.disc, selected_tissues = NULL, 
                                                 cutoff = 4, verbose = TRUE)

ThETA_ColonCancer2Gene_lib <- row.names(ColonCancer_Tscores[order(ColonCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_ColonCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_ColonCancer2Gene_lib, c("ensembl_id")])
names(ThETA_ColonCancer2Gene_lib) <- paste0("Colon carcinoma (EFO_1001950)", "[ThETA]")



# Merge and save into one

Enrichment_ColonCancer2Gene_lib <- c(DisGeNET_ColonCancer2Gene_lib, OpenTargets_ColonCancer2Gene_GA_lib, 
                                        # OpenTargets_ColonCancer2Gene_RNA_lib, 
                                        OpenTargets_ColonCancer2Gene_lit_lib, 
                                        Enrichr_GeoColonCancerSignatures_Up, Enrichr_GeoColonCancerSignatures_Down,
                                        Intogen_ColonCancer2Gene_lib,
                                        PharmGKB_ColonCancer2Gene_lib, ThETA_ColonCancer2Gene_lib)



# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_ColonCancer2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_ColonCancer_lib.rds")

print(warnings())