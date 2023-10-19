set.seed(5081)



# Enrichment analysis libraries for SkinCancer
# Notes:
# (a) Includes melanoma (metastatic and cutaneous)


# Load libraries
library(org.Hs.eg.db)



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")



## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_SkinCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep("skin|C0025202|C0278883|C0151779", 
                                                               names(DisGeNET_Disease2Gene_lib), 
                                                               ignore.case = TRUE)]
DisGeNET_SkinCancer2Gene_lib <- DisGeNET_SkinCancer2Gene_lib[grep("cancer|carcinoma|sarcoma|melanoma", 
                                                                  names(DisGeNET_SkinCancer2Gene_lib), 
                                                                  ignore.case = TRUE)]
DisGeNET_SkinCancer2Gene_lib <- DisGeNET_SkinCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                  names(DisGeNET_SkinCancer2Gene_lib), 
                                                                  ignore.case = TRUE, 
                                                                  invert = TRUE)] 
names(DisGeNET_SkinCancer2Gene_lib) <- paste0(names(DisGeNET_SkinCancer2Gene_lib), "[DisGeNET_curated]")



## OpenTargets

OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_SkinCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep("skin|EFO_0000756|EFO_0002617|EFO_0000389", 
                                                                           names(OpenTargets_Disease2Gene_GA_lib), 
                                                                           ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_GA_lib <- OpenTargets_SkinCancer2Gene_GA_lib[grep("cancer|carcinoma|sarcoma|melanoma", 
                                                                              names(OpenTargets_SkinCancer2Gene_GA_lib), 
                                                                              ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_GA_lib <- OpenTargets_SkinCancer2Gene_GA_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                              names(OpenTargets_SkinCancer2Gene_GA_lib), 
                                                                              ignore.case = TRUE, 
                                                                              invert = TRUE)] 
names(OpenTargets_SkinCancer2Gene_GA_lib) <- paste0(names(OpenTargets_SkinCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_SkinCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep("skin|EFO_0000756|EFO_0002617|EFO_0000389",
                                                                             names(OpenTargets_Disease2Gene_RNA_lib),
                                                                             ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_RNA_lib <- OpenTargets_SkinCancer2Gene_RNA_lib[grep("cancer|carcinoma|sarcoma|melanoma",
                                                                                names(OpenTargets_SkinCancer2Gene_RNA_lib),
                                                                                ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_RNA_lib <- OpenTargets_SkinCancer2Gene_RNA_lib[grep("hereditary|familial|susceptibility|predisposition",
                                                                                names(OpenTargets_SkinCancer2Gene_RNA_lib),
                                                                                ignore.case = TRUE,
                                                                                invert = TRUE)]
names(OpenTargets_SkinCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_SkinCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_SkinCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep("skin|EFO_0000756|EFO_0002617|EFO_0000389",
                                                                             names(OpenTargets_Disease2Gene_lit_lib),
                                                                             ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_lit_lib <- OpenTargets_SkinCancer2Gene_lit_lib[grep("cancer|carcinoma|sarcoma|melanoma",
                                                                                names(OpenTargets_SkinCancer2Gene_lit_lib),
                                                                                ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_lit_lib <- OpenTargets_SkinCancer2Gene_lit_lib[grep("hereditary|familial|susceptibility|predisposition",
                                                                                names(OpenTargets_SkinCancer2Gene_lit_lib),
                                                                                ignore.case = TRUE,
                                                                                invert = TRUE)]
names(OpenTargets_SkinCancer2Gene_lit_lib) <- paste0(names(OpenTargets_SkinCancer2Gene_lit_lib), "[OpenTargets_literature]")


OpenTargets_Disease2Gene_SM_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_SM_lib.rds")
OpenTargets_SkinCancer2Gene_SM_lib <- OpenTargets_Disease2Gene_SM_lib[grep("skin|EFO_0000756|EFO_0002617|EFO_0000389",
                                                                           names(OpenTargets_Disease2Gene_SM_lib),
                                                                           ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_SM_lib <- OpenTargets_SkinCancer2Gene_SM_lib[grep("cancer|carcinoma|sarcoma|melanoma",
                                                                              names(OpenTargets_SkinCancer2Gene_SM_lib),
                                                                              ignore.case = TRUE)]
OpenTargets_SkinCancer2Gene_SM_lib <- OpenTargets_SkinCancer2Gene_SM_lib[grep("hereditary|familial|susceptibility|predisposition",
                                                                              names(OpenTargets_SkinCancer2Gene_SM_lib),
                                                                              ignore.case = TRUE,
                                                                              invert = TRUE)]
names(OpenTargets_SkinCancer2Gene_SM_lib) <- paste0(names(OpenTargets_SkinCancer2Gene_SM_lib), "[OpenTargets_SM]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoSkinCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep("skin", 
                                                                                     names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                                                     ignore.case = TRUE)]
Enrichr_GeoSkinCancerSignatures_Up <- Enrichr_GeoSkinCancerSignatures_Up[grep("cancer|carcinoma|sarcoma|melanoma", 
                                                                              names(Enrichr_GeoSkinCancerSignatures_Up), 
                                                                              ignore.case = TRUE)]
Enrichr_GeoSkinCancerSignatures_Up <- Enrichr_GeoSkinCancerSignatures_Up[grep("hereditary|familial|susceptibility|predisposition", 
                                                                              names(Enrichr_GeoSkinCancerSignatures_Up), 
                                                                              ignore.case = TRUE, 
                                                                              invert = TRUE)] 
names(Enrichr_GeoSkinCancerSignatures_Up) <- paste0(names(Enrichr_GeoSkinCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoSkinCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep("skin", 
                                                                                         names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                                         ignore.case = TRUE)]
Enrichr_GeoSkinCancerSignatures_Down <- Enrichr_GeoSkinCancerSignatures_Down[grep("cancer|carcinoma|sarcoma|melanoma", 
                                                                                  names(Enrichr_GeoSkinCancerSignatures_Down), 
                                                                                  ignore.case = TRUE)]
Enrichr_GeoSkinCancerSignatures_Down <- Enrichr_GeoSkinCancerSignatures_Down[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                  names(Enrichr_GeoSkinCancerSignatures_Down), 
                                                                                  ignore.case = TRUE, 
                                                                                  invert = TRUE)] 
names(Enrichr_GeoSkinCancerSignatures_Down) <- paste0(names(Enrichr_GeoSkinCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")



## Intogen

Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Intogen_Disease2Gene_lib.rds")
Intogen_SkinCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep("skin|melanoma", names(Intogen_Disease2Gene_lib), 
                                                             ignore.case = TRUE)]
Intogen_SkinCancer2Gene_lib <- Intogen_SkinCancer2Gene_lib[grep("Uveal", 
                                                                names(Intogen_SkinCancer2Gene_lib), 
                                                                ignore.case = TRUE, 
                                                                invert = TRUE)] 
Intogen_SkinCancer2Gene_lib <- Intogen_SkinCancer2Gene_lib[grep("cancer|carcinoma|sarcoma|melanoma", 
                                                                names(Intogen_SkinCancer2Gene_lib), 
                                                                ignore.case = TRUE)]
Intogen_SkinCancer2Gene_lib <- Intogen_SkinCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                names(Intogen_SkinCancer2Gene_lib), 
                                                                ignore.case = TRUE, 
                                                                invert = TRUE)] 
names(Intogen_SkinCancer2Gene_lib) <- paste0(names(Intogen_SkinCancer2Gene_lib), "[Intogen]")



## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/PharmGKB_Disease2Gene_lib.rds")
PharmGKB_SkinCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep("skin|melanoma",
                                                               names(PharmGKB_Disease2Gene_lib),
                                                               ignore.case = TRUE)]
PharmGKB_SkinCancer2Gene_lib <- PharmGKB_SkinCancer2Gene_lib[grep("cancer|carcinoma|sarcoma|melanoma",
                                                                  names(PharmGKB_SkinCancer2Gene_lib),
                                                                  ignore.case = TRUE)]
PharmGKB_SkinCancer2Gene_lib <- PharmGKB_SkinCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition",
                                                                  names(PharmGKB_SkinCancer2Gene_lib),
                                                                  ignore.case = TRUE,
                                                                  invert = TRUE)]
names(PharmGKB_SkinCancer2Gene_lib) <- paste0(names(PharmGKB_SkinCancer2Gene_lib), "[PharmGKB_associated]")



## ThETA [NONE FOUND]

library(ThETA)
source("Scripts/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
SkinCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0000756")]]


SkinCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0000756"),]

SkinCancer_Tscores <- tissue.specific.scores(disease_genes = SkinCancer_genes$entrez,
                                              ppi_network = ppi_strdb_700,
                                              directed_network = FALSE,
                                              tissue_expr_data = gtexv7_zscore,
                                              dis_relevant_tissues = SkinCancer_rel_tissue_scores,
                                              W = centrality_score$borda.disc, selected_tissues = NULL,
                                              cutoff = 4, verbose = TRUE)

ThETA_SkinCancer2Gene_lib <- row.names(SkinCancer_Tscores[order(SkinCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_SkinCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_SkinCancer2Gene_lib, c("ensembl_id")])
names(ThETA_SkinCancer2Gene_lib) <- paste0("melanoma (EFO_0000756)", "[ThETA]")



# Merge and save into one

Enrichment_SkinCancer2Gene_lib <- c(DisGeNET_SkinCancer2Gene_lib, OpenTargets_SkinCancer2Gene_GA_lib, 
                                    OpenTargets_SkinCancer2Gene_RNA_lib, OpenTargets_SkinCancer2Gene_lit_lib, 
                                    OpenTargets_SkinCancer2Gene_SM_lib, 
                                    Enrichr_GeoSkinCancerSignatures_Up, Enrichr_GeoSkinCancerSignatures_Down,
                                    PharmGKB_SkinCancer2Gene_lib, Intogen_SkinCancer2Gene_lib,
                                    ThETA_SkinCancer2Gene_lib)



# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(Enrichment_SkinCancer2Gene_lib, "InputFiles/Enrichment_analysis_libraries/Disease2Gene_SkinCancer_lib.rds")



print(warnings())