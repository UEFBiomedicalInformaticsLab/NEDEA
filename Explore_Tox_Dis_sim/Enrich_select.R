DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
select <- DisGeNET_Disease2Gene_lib[grep("skin", 
                                         names(DisGeNET_Disease2Gene_lib), 
                                         ignore.case = TRUE)]
select <- select[grep("cancer|carcinoma|sarcoma|melanoma", 
                      names(select), 
                      ignore.case = TRUE)]
lengths(select)



Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

select_up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep("skin", 
                                                            names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                            ignore.case = TRUE)]
select_up <- select_up[grep("cancer|carcinoma|sarcoma|melanoma", 
                            names(select_up), 
                            ignore.case = TRUE)]

select_down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep("skin", 
                                                                names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                ignore.case = TRUE)]
select_down <- select_down[grep("cancer|carcinoma|sarcoma|melanoma", 
                                names(select_down), 
                                ignore.case = TRUE)]



Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
select <- Intogen_Disease2Gene_lib[grep("skin", names(Intogen_Disease2Gene_lib), 
                                        ignore.case = TRUE)]
select <- select[grep("cancer|carcinoma|sarcoma|melanoma", 
                      names(select), 
                      ignore.case = TRUE)]
lengths(select)


PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")
select <- PharmGKB_Disease2Gene_lib[grep("skin", 
                                         names(PharmGKB_Disease2Gene_lib), 
                                         ignore.case = TRUE)]
select <- select[grep("cancer|carcinoma|sarcoma|melanoma", 
                      names(select), 
                      ignore.case = TRUE)]
lengths(select)



OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
select <- OpenTargets_Disease2Gene_GA_lib[grep("skin", 
                                               names(OpenTargets_Disease2Gene_GA_lib), 
                                               ignore.case = TRUE)]
select <- select[grep("cancer|carcinoma|sarcoma|melanoma", 
                      names(select), 
                      ignore.case = TRUE)]
lengths(select)


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
select <- OpenTargets_Disease2Gene_RNA_lib[grep("skin", 
                                                names(OpenTargets_Disease2Gene_RNA_lib), 
                                                ignore.case = TRUE)]
select <- select[grep("cancer|carcinoma|sarcoma|melanoma", 
                      names(select), 
                      ignore.case = TRUE)]
lengths(select)


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
select <- OpenTargets_Disease2Gene_lit_lib[grep("skin", 
                                                names(OpenTargets_Disease2Gene_lit_lib), 
                                                ignore.case = TRUE)]
select <- select[grep("cancer|carcinoma|sarcoma|melanoma", 
                      names(select), 
                      ignore.case = TRUE)]
lengths(select)


rm(list = ls())