DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
select <- DisGeNET_Disease2Gene_lib[grep("toxic", names(DisGeNET_Disease2Gene_lib), ignore.case = TRUE)]
lengths(select)



Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")
select_up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep("toxic", names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), ignore.case = TRUE)]
select_down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep("toxic", names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), ignore.case = TRUE)]
lengths(select_up)
lengths(select_down)



Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
select <- Intogen_Disease2Gene_lib[grep("toxic", names(Intogen_Disease2Gene_lib), ignore.case = TRUE)]
lengths(select)


PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")
select <- PharmGKB_Disease2Gene_lib[grep("toxic", names(PharmGKB_Disease2Gene_lib), ignore.case = TRUE)]
lengths(select)



OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
select <- OpenTargets_Disease2Gene_GA_lib[grep("toxic", names(OpenTargets_Disease2Gene_GA_lib), ignore.case = TRUE)]
lengths(select)


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
select <- OpenTargets_Disease2Gene_RNA_lib[grep("toxic", names(OpenTargets_Disease2Gene_RNA_lib), ignore.case = TRUE)]
lengths(select)


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
select <- OpenTargets_Disease2Gene_lit_lib[grep("toxic", names(OpenTargets_Disease2Gene_lit_lib), ignore.case = TRUE)]
lengths(select)

select <- OpenTargets_Disease2Gene_lib[grep("toxic", names(OpenTargets_Disease2Gene_lib), ignore.case = TRUE)]
lengths(select)

ADReCS_ADR2Gene_level3_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level3_lib.rds")
ADReCS_ADR2Gene_level4_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level4_lib.rds")
ADReCS_ADR2Gene_lib <- c(ADReCS_ADR2Gene_level4_lib, ADReCS_ADR2Gene_level3_lib)
select <- ADReCS_ADR2Gene_lib[grep("toxic", names(ADReCS_ADR2Gene_lib), ignore.case = TRUE)]
lengths(select)


CTD_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/CTD_Disease2Gene_lib.rds")
select <- CTD_Disease2Gene_lib[grep("toxic", names(CTD_Disease2Gene_lib), ignore.case = TRUE)]
lengths(select)



msigdb_keggPath2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/msigdb_keggPath2Gene_lib.rds")
select <- msigdb_keggPath2Gene_lib[grep("drug", names(msigdb_keggPath2Gene_lib), ignore.case = TRUE)]
lengths(select)


rm(list = ls())
