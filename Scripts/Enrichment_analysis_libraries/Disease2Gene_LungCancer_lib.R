set.seed(5081)



# Enrichment analysis libraries for Lung Cancer



# Load libraries
library(org.Hs.eg.db)



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")



## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_LungCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep("lung", 
                                                               names(DisGeNET_Disease2Gene_lib), 
                                                               ignore.case = TRUE)]
DisGeNET_LungCancer2Gene_lib <- DisGeNET_LungCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                                                                  names(DisGeNET_LungCancer2Gene_lib), 
                                                                  ignore.case = TRUE)]
DisGeNET_LungCancer2Gene_lib <- DisGeNET_LungCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                  names(DisGeNET_LungCancer2Gene_lib), 
                                                                  ignore.case = TRUE, 
                                                                  invert = TRUE)] 
names(DisGeNET_LungCancer2Gene_lib) <- paste0(names(DisGeNET_LungCancer2Gene_lib), "[DisGeNET_curated]")



## OpenTargets

OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_LungCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep("lung", 
                                                                           names(OpenTargets_Disease2Gene_GA_lib), 
                                                                           ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_GA_lib <- OpenTargets_LungCancer2Gene_GA_lib[grep("cancer|carcinoma|sarcoma", 
                                                                              names(OpenTargets_LungCancer2Gene_GA_lib), 
                                                                              ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_GA_lib <- OpenTargets_LungCancer2Gene_GA_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                              names(OpenTargets_LungCancer2Gene_GA_lib), 
                                                                              ignore.case = TRUE, 
                                                                              invert = TRUE)] 
names(OpenTargets_LungCancer2Gene_GA_lib) <- paste0(names(OpenTargets_LungCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_LungCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep("lung", 
                                                                             names(OpenTargets_Disease2Gene_RNA_lib), 
                                                                             ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_RNA_lib <- OpenTargets_LungCancer2Gene_RNA_lib[grep("cancer|carcinoma|sarcoma", 
                                                                                names(OpenTargets_LungCancer2Gene_RNA_lib), 
                                                                                ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_RNA_lib <- OpenTargets_LungCancer2Gene_RNA_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                names(OpenTargets_LungCancer2Gene_RNA_lib), 
                                                                                ignore.case = TRUE, 
                                                                                invert = TRUE)] 
names(OpenTargets_LungCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_LungCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_LungCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep("lung", 
                                                                             names(OpenTargets_Disease2Gene_lit_lib), 
                                                                             ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_lit_lib <- OpenTargets_LungCancer2Gene_lit_lib[grep("cancer|carcinoma|sarcoma", 
                                                                                names(OpenTargets_LungCancer2Gene_lit_lib), 
                                                                                ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_lit_lib <- OpenTargets_LungCancer2Gene_lit_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                names(OpenTargets_LungCancer2Gene_lit_lib), 
                                                                                ignore.case = TRUE, 
                                                                                invert = TRUE)] 
names(OpenTargets_LungCancer2Gene_lit_lib) <- paste0(names(OpenTargets_LungCancer2Gene_lit_lib), "[OpenTargets_literature]")


OpenTargets_Disease2Gene_SM_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_SM_lib.rds")
OpenTargets_LungCancer2Gene_SM_lib <- OpenTargets_Disease2Gene_SM_lib[grep("lung", 
                                                                           names(OpenTargets_Disease2Gene_SM_lib), 
                                                                           ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_SM_lib <- OpenTargets_LungCancer2Gene_SM_lib[grep("cancer|carcinoma|sarcoma", 
                                                                              names(OpenTargets_LungCancer2Gene_SM_lib), 
                                                                              ignore.case = TRUE)]
OpenTargets_LungCancer2Gene_SM_lib <- OpenTargets_LungCancer2Gene_SM_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                              names(OpenTargets_LungCancer2Gene_SM_lib), 
                                                                              ignore.case = TRUE, 
                                                                              invert = TRUE)] 
names(OpenTargets_LungCancer2Gene_SM_lib) <- paste0(names(OpenTargets_LungCancer2Gene_SM_lib), "[OpenTargets_SM]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoLungCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep("lung", 
                                                                                     names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                                                     ignore.case = TRUE)]
Enrichr_GeoLungCancerSignatures_Up <- Enrichr_GeoLungCancerSignatures_Up[grep("cancer|carcinoma|sarcoma", 
                                                                              names(Enrichr_GeoLungCancerSignatures_Up), 
                                                                              ignore.case = TRUE)]
Enrichr_GeoLungCancerSignatures_Up <- Enrichr_GeoLungCancerSignatures_Up[grep("hereditary|familial|susceptibility|predisposition", 
                                                                              names(Enrichr_GeoLungCancerSignatures_Up), 
                                                                              ignore.case = TRUE, 
                                                                              invert = TRUE)] 
names(Enrichr_GeoLungCancerSignatures_Up) <- paste0(names(Enrichr_GeoLungCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")


Enrichr_GeoLungCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep("lung", 
                                                                                         names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                                         ignore.case = TRUE)]
Enrichr_GeoLungCancerSignatures_Down <- Enrichr_GeoLungCancerSignatures_Down[grep("cancer|carcinoma|sarcoma", 
                                                                                  names(Enrichr_GeoLungCancerSignatures_Down), 
                                                                                  ignore.case = TRUE)]
Enrichr_GeoLungCancerSignatures_Down <- Enrichr_GeoLungCancerSignatures_Down[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                  names(Enrichr_GeoLungCancerSignatures_Down), 
                                                                                  ignore.case = TRUE, 
                                                                                  invert = TRUE)] 
names(Enrichr_GeoLungCancerSignatures_Down) <- paste0(names(Enrichr_GeoLungCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")



## Intogen

Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Intogen_Disease2Gene_lib.rds")
Intogen_LungCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep("lung", names(Intogen_Disease2Gene_lib), 
                                                             ignore.case = TRUE)]
Intogen_LungCancer2Gene_lib <- Intogen_LungCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                                                                names(Intogen_LungCancer2Gene_lib), 
                                                                ignore.case = TRUE)]
Intogen_LungCancer2Gene_lib <- Intogen_LungCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                names(Intogen_LungCancer2Gene_lib), 
                                                                ignore.case = TRUE, 
                                                                invert = TRUE)] 
names(Intogen_LungCancer2Gene_lib) <- paste0(names(Intogen_LungCancer2Gene_lib), "[Intogen]")



## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/PharmGKB_Disease2Gene_lib.rds")
PharmGKB_LungCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep("lung", 
                                                               names(PharmGKB_Disease2Gene_lib), 
                                                               ignore.case = TRUE)]
PharmGKB_LungCancer2Gene_lib <- PharmGKB_LungCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                                                                  names(PharmGKB_LungCancer2Gene_lib), 
                                                                  ignore.case = TRUE)]
PharmGKB_LungCancer2Gene_lib <- PharmGKB_LungCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                  names(PharmGKB_LungCancer2Gene_lib), 
                                                                  ignore.case = TRUE, 
                                                                  invert = TRUE)] 
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
                                    OpenTargets_LungCancer2Gene_SM_lib, 
                                    Enrichr_GeoLungCancerSignatures_Up, Enrichr_GeoLungCancerSignatures_Down,
                                    Intogen_LungCancer2Gene_lib,
                                    PharmGKB_LungCancer2Gene_lib, ThETA_LungCancer2Gene_lib)




# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(Enrichment_LungCancer2Gene_lib, "InputFiles/Enrichment_analysis_libraries/Disease2Gene_LungCancer_lib.rds")



print(warnings())