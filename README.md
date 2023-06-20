# DrugCombination_1

This repository contains the scripts used to train classifiers for prioritization of drug combinations for cancer based on efficacy and safety estimates.

## Preparing the environment

The following packages must be available for executing the scripts: biomaRt; dbparser; httr; igraph; jsonlite; msigdbr; openxlsx; org.Hs.eg.db; readxl; sparklyr; sparklyr.nested; tidyverse; unixtools

The packages can be installed as:

## Databases

The following databases have been used in the pipeline:

## Usage

The complete pipeline can be divided into five parts. The parts are described in sequential order as they must be executed. The following section also lists the primary input and output files for each scripts.

Several of the scripts require the following parameters as input:

**\--disease:** The disease for which the script to be executed. Possible values: BreastCancer, LungCancer, KidneyCancer, OvaryCancer, ProstateCancer, SkinCancer

### Part 1: Preparing the input files for training the classifiers

#### 1.1. Curating protein-protein interaction network

+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Network_build/STRING_PPI.R](Scripts/Network_build/STRING_PPI.R)                                                                                                                                                                                                                     |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | a.  9606.protein.links.detailed.v11.5.txt.gz (STRING)                                                                                                                                                                                                                                        |
|              | b.  9606.protein.aliases.v11.5.txt.gz (STRING)                                                                                                                                                                                                                                               |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | STRING_PPI_Net_database.rds                                                                                                                                                                                                                                                                  |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | Creates a protein-protein interaction network from interactions reported in the STRING database. The script retrieves the source interaction file for human from STRING directly. Only interactions supported by experiments or those curated in string from other databases are considered. |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

#### 1.2. Retrieving associations between drugs, targets and diseases

+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script**     | [Scripts/Associations/DrugBank_Drug_Target_Associations.R](Scripts/Associations/DrugBank_Drug_Target_Associations.R)                                                                               |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**        | drugbank_all_full_database.xml.zip                                                                                                                                                                 |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**       | DrugBank_Drug_Target_Net.rds                                                                                                                                                                       |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**      | The script retrieves the drug target association for the drugs in DrugBank. Only human targets are reported. The output is a dataframe containing drugs and its corresponding targets in each row. |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Troubleshoot** | Missing files: The DrugBank database as XML should be separately downloaded. Permission from DrugBank is required to download the file.                                                            |
|                  |                                                                                                                                                                                                    |
|                  | `curl -Lfv -o filename.zip -u EMAIL:PASSWORD https://go.drugbank.com/releases/5-1-8/downloads/all-full-database`                                                                                   |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+--------------+------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Associations/OpenTargets_Drug_Disease_Associations.R](Scripts/Associations/OpenTargets_Drug_Disease_Associations.R) |
+--------------+------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | Core annotation for drug molecules (OpenTarget parquest files)                                                               |
+--------------+------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | OpenTargets_Drug_Disease_Net.rds                                                                                             |
+--------------+------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The script retrieves the drug disease associations from OpenTargets database.                                                |
+--------------+------------------------------------------------------------------------------------------------------------------------------+

#### 1.3. Curating the drug combinations for training the classifiers

+--------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Reference_lists/FIMMDrugComb_drugCombinations.R](Scripts/Reference_lists/FIMMDrugComb_drugCombinations.R)                                                                                                                                                    |
+--------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | summary_v\_1_5.csv (DrugComb.org)                                                                                                                                                                                                                                     |
+--------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | a.  FimmDrugComb_drugCombinations.rds                                                                                                                                                                                                                                 |
|              | b.  FimmDrugComb\_[disease]\_drugCombinations.rds                                                                                                                                                                                                                     |
+--------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The scripts retrieves the drug pairs and the synergy scores from the FIMM DrugComb portal. The pairs are then categorized as synergistic or antagonistic based on the synergy scores. Drug pairs are filtered to keep only those tested on cancer related cell lines. |
+--------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+--------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Reference_lists/DrugBank_drugInteractions.R](Scripts/Reference_lists/DrugBank_drugInteractions.R)                                  |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | a.  drug_drug_interactions.csv                                                                                                              |
|              | b.  drug.csv                                                                                                                                |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | a.  DrugBank_drugInteractions.rds                                                                                                           |
|              | b.  DrugBank_drugInteractions_withRiskSeverity.rds                                                                                          |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The script retrieves drug-drug interactions from DrugBank that have been reported to be involved with risk or severity to any side effects. |
+--------------+---------------------------------------------------------------------------------------------------------------------------------------------+

+----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script**   | [Scripts/Drug_combinations/DrugCombinations.R](Scripts/Drug_combinations/DrugCombinations.R)                                                                                                                                                                                                                                                                                                                            |
+----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Parameters** | \--disease                                                                                                                                                                                                                                                                                                                                                                                                              |
+----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**      | a.  DrugBank_drugInteractions_withRiskSeverity.rds                                                                                                                                                                                                                                                                                                                                                                      |
|                | b.  DrugBank_Drug_Target_Net.rds                                                                                                                                                                                                                                                                                                                                                                                        |
|                | c.  FimmDrugComb\_[disease]\_drugCombinations.rds                                                                                                                                                                                                                                                                                                                                                                       |
+----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**     | DrugComb\_[disease].rds                                                                                                                                                                                                                                                                                                                                                                                                 |
+----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**    | The script prepares the effective and adverse drug combinations for each cancer type. The synergistic drug pairs retrieved from FIMM DrugComb portal are considered as effective. While the DrugBank drug-drug interaction pairs with risk/severity indications that contain both drugs involved in forming effective pairs are considered as adverse. Drug pairs are filtered to keep only those with reported targets |
+----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+--------------+------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Drug_combinations/DrugCombs_training_annotation.R](Scripts/Drug_combinations/DrugCombs_training_annotation.R) |
+--------------+------------------------------------------------------------------------------------------------------------------------+
| **Input**    | DrugComb\_[disease].rds                                                                                                |
+--------------+------------------------------------------------------------------------------------------------------------------------+
| **Output**   | DrugCombs_training_annotation.xlsx                                                                                     |
+--------------+------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The script summarises the drug pairs for the training set. Also includes DrugBank drug ID to drug name mapping.        |
+--------------+------------------------------------------------------------------------------------------------------------------------+

#### 1.4. Preparing the libraries for estimation of efficacy and safety

+--------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_ADR.R](Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_ADR.R) |
+--------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | ADRAlert2GENE2ID.xlsx                                                                                                                                  |
+--------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | a.  ADReCS_ADR2Gene_level4_lib.rds                                                                                                                     |
|              | b.  ADReCS_ADR2Gene_level3_lib.rds                                                                                                                     |
+--------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The script collates gene sets related to side effects from ADReCs. The gene sets for the two different levels of side effects are exported separately. |
+--------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+

+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Diseases.R](Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Diseases.R) |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | (a) curated_gene_disease_associations.tsv.gz (DisGeNET)                                                                                                          |
|              | (b) Data-type metrics for direct target-disease associations (OpenTarget parquest files)                                                                         |
|              | (c) Compendium_Cancer_Genes.tsv (Intogen)                                                                                                                        |
|              | (d) cohorts.tsv (Intogen)                                                                                                                                        |
|              | (e) Disease_Signatures_from_GEO_up_2014.txt (Enrichr)                                                                                                            |
|              | (f) Disease_Signatures_from_GEO_down_2014.txt (Enrichr)                                                                                                          |
|              | (g) CTD_genes_diseases.csv.gz (Comparative Toxicogenomics Database)                                                                                              |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | a.  DisGeNET_Disease2Gene_lib.rds                                                                                                                                |
|              | b.  DisGeNET_DiseaseGroup2Gene_lib.rds                                                                                                                           |
|              | c.  DisGeNET_Phenotype2Gene_lib.rds                                                                                                                              |
|              | d.  OpenTargets_Disease2Gene_GA_lib.rds                                                                                                                          |
|              | e.  OpenTargets_Disease2Gene_RNA_lib.rds                                                                                                                         |
|              | f.  OpenTargets_Disease2Gene_lit_lib.rds                                                                                                                         |
|              | g.  Intogen_Disease2Gene_lib.rds                                                                                                                                 |
|              | h.  Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds                                                                                                                   |
|              | i.  PharmGKB_Disease2Gene_lib.rds                                                                                                                                |
|              | j.  CTD_Disease2Gene_lib.rds                                                                                                                                     |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The script collates gene sets related to diseases from several databases like DisGeNET, OpenTargets, Intogen, Enrichr, PharmGKB.                                 |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_DrugWithdrawalADR.R](Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_DrugWithdrawalADR.R) |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | a.  ADReCS_ADR2Gene_level4_lib.rds                                                                                                                                                 |
|              | b.  ADReCS_ADR2Gene_level3_lib.rds                                                                                                                                                 |
|              | c.  CTD_Disease2Gene_lib.rds                                                                                                                                                       |
|              | d.  DisGeNET_Disease2Gene_lib.rds                                                                                                                                                  |
|              | e.  DisGeNET_DiseaseGroup2Gene_lib.rds                                                                                                                                             |
|              | f.  DisGeNET_Phenotype2Gene_lib.rds                                                                                                                                                |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | a.  drugWithdrawal_Adr2Gene_lib.rds                                                                                                                                                |
|              | b.  drugWithdrawal_Adr2Gene_libInfo.xlsx                                                                                                                                           |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The script prepares gene set library for adverse drug reactions related to drug withdrawal.                                                                                        |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script**     | [Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Miscellaneous.R](Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Miscellaneous.R)                       |
+------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**        | a.  categories.tsv (DGIdb)                                                                                                                                                                       |
|                  | b.  genes.tsv (PharmGKB)                                                                                                                                                                         |
|                  | c.  genage_human.csv (GenAge db)                                                                                                                                                                 |
|                  | d.  Halifax-curation.Table S2. Hallmark of Cancer Data - Gene and Pathway.xlsx                                                                                                                   |
+------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**       | miscellaneous_gene_lib.rds                                                                                                                                                                       |
+------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**      | The script compiles gene sets related to druggable geneome, very important pharmacogenes, ageing, and cancer hallmarks into a single library.                                                    |
+------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Troubleshoot** | Missing files: The Halifax-curation.Table S2. Hallmark of Cancer Data - Gene and Pathway.xlsx file must be downloaded manually from the journal website https://doi.org/10.1093/database/baaa045 |
+------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Pathways.R](Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Pathways.R) |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | a.  Cancer hallmark related pathways (Table_1.xls from [PMC7013921](https://doi.org/10.3389%2Ffgene.2020.00029))                                                 |
|              | b.  smpdb_proteins.csv.zip                                                                                                                                       |
|              | c.  smpdb_pathways.csv.zip                                                                                                                                       |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | a.  CHG_keggPath2Gene_lib.rds                                                                                                                                    |
|              | b.  SMPDb_Pathway2Gene_lib.rds                                                                                                                                   |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The script retrieves gene sets related to cancer hallmark linked cancer pathways and small molecular pathways from SMPDB.                                        |
+--------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------+

+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| **R Script** | [Scripts/Enrichment_Analysis_Libraries/Disease2Gene_BreastCancer_lib.R](Scripts/Enrichment_Analysis_Libraries/Disease2Gene_BreastCancer_lib.R)     |
|              |                                                                                                                                                    |
|              | [Scripts/Enrichment_Analysis_Libraries/Disease2Gene_KidneyCancer_lib.R](Scripts/Enrichment_Analysis_Libraries/Disease2Gene_KidneyCancer_lib.R)     |
|              |                                                                                                                                                    |
|              | [Scripts/Enrichment_Analysis_Libraries/Disease2Gene_LungCancer_lib.R](Scripts/Enrichment_Analysis_Libraries/Disease2Gene_LungCancer_lib.R)         |
|              |                                                                                                                                                    |
|              | [Scripts/Enrichment_Analysis_Libraries/Disease2Gene_OvaryCancer_lib.R](Scripts/Enrichment_Analysis_Libraries/Disease2Gene_OvaryCancer_lib.R)       |
|              |                                                                                                                                                    |
|              | [Scripts/Enrichment_Analysis_Libraries/Disease2Gene_ProstateCancer_lib.R](Scripts/Enrichment_Analysis_Libraries/Disease2Gene_ProstateCancer_lib.R) |
|              |                                                                                                                                                    |
|              | [Scripts/Enrichment_Analysis_Libraries/Disease2Gene_SkinCancer_lib.R](Scripts/Enrichment_Analysis_Libraries/Disease2Gene_SkinCancer_lib.R)         |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| **Input**    | a.  DisGeNET_Disease2Gene_lib.rds                                                                                                                  |
|              | b.  OpenTargets_Disease2Gene_GA_lib.rds                                                                                                            |
|              | c.  OpenTargets_Disease2Gene_RNA_lib.rds                                                                                                           |
|              | d.  OpenTargets_Disease2Gene_lit_lib.rds                                                                                                           |
|              | e.  Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds                                                                                                     |
|              | f.  Intogen_Disease2Gene_lib.rds                                                                                                                   |
|              | g.  PharmGKB_Disease2Gene_lib.rds                                                                                                                  |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| **Output**   | Disease2Gene\_[disease]\_lib.rds                                                                                                                   |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| **Summary**  | The scripts compiles disease specific gene set libraries from different databases and TheTA.                                                       |
+--------------+----------------------------------------------------------------------------------------------------------------------------------------------------+

### Part 2: Preparing the features to be used in training the classifiers

#### 2.1. Execute RWR on the drug combinations

+---------------+-------------------------------------------------+
| **R Script**  | Scripts/Feature_generation/DrugCombs_dNetRWR.R  |
+---------------+-------------------------------------------------+
| **Input**     |                                                 |
+---------------+-------------------------------------------------+
| **Output**    |                                                 |
+---------------+-------------------------------------------------+
| **Summary**   |                                                 |
+---------------+-------------------------------------------------+

#### 2.2. Execute FGSEA on the ranked gene list from RWR

+--------------+--------------------------------------------------------------+
| **R Script** | Scripts/Feature_generation/DrugCombs_fgsea_CommonLibraries.R |
|              |                                                              |
|              | Scripts/Feature_generation/DrugCombs_fgsea_EfficacySafety.R  |
+--------------+--------------------------------------------------------------+
| **Input**    |                                                              |
+--------------+--------------------------------------------------------------+
| **Output**   |                                                              |
+--------------+--------------------------------------------------------------+
| **Summary**  |                                                              |
+--------------+--------------------------------------------------------------+

#### 2.3. Calculate network proximity

+--------------+---------------------------------------------------------------------+
| **R Script** | Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugAdr.       |
|              |                                                                     |
|              | RScripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDisease.R |
|              |                                                                     |
|              | Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDrug.R     |
+--------------+---------------------------------------------------------------------+
| **Input**    |                                                                     |
+--------------+---------------------------------------------------------------------+
| **Output**   |                                                                     |
+--------------+---------------------------------------------------------------------+
| **Summary**  |                                                                     |
+--------------+---------------------------------------------------------------------+

### Part 3: Training and validation of the classifiers

#### 3.1. Prepare data for repeated cross validation

+---------------+---------------------------------------------------+
| **R Script**  | Scripts/Model_train/DrugCombs_trainML_rwrFgsea.R  |
+---------------+---------------------------------------------------+
| **Input**     |                                                   |
+---------------+---------------------------------------------------+
| **Output**    |                                                   |
+---------------+---------------------------------------------------+
| **Summary**   |                                                   |
+---------------+---------------------------------------------------+

#### 3.2. Train and validate the models

+--------------+------------------------------------------------------------------------+
| **R Script** | Scripts/Model_train/DrugCombs_trainML_rwrFgsea.R                       |
|              |                                                                        |
|              | Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugAdr.R        |
|              |                                                                        |
|              | Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDisease.R    |
|              |                                                                        |
|              | Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDrug.R       |
|              |                                                                        |
|              | Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDiseaseAdr.R |
+--------------+------------------------------------------------------------------------+
| **Input**    |                                                                        |
+--------------+------------------------------------------------------------------------+
| **Output**   |                                                                        |
+--------------+------------------------------------------------------------------------+
| **Summary**  |                                                                        |
+--------------+------------------------------------------------------------------------+

### Part 4: Summarizing the results
