
# DrugCombination_1

This repository contains the scripts used to train classifiers for
prioritization of drug combinations for cancer based on efficacy and
safety estimates.

## Preparing the environment

The following packages must be available for executing the scripts:
biomaRt; dbparser; httr; igraph; jsonlite; msigdbr; openxlsx;
org.Hs.eg.db; readxl; sparklyr; sparklyr.nested; tidyverse; unixtools

The packages can be installed as:

## Databases

The following databases have been used in the pipeline:

## Usage

The complete pipeline can be divided into five parts. The parts are
described in sequential order as they must be executed. The following
section also lists the primary input and output files for each scripts.

Several of the scripts require the following parameters as input:

**--disease:** The disease for which the script to be executed. Possible
values: BreastCancer, LungCancer, KidneyCancer, OvaryCancer,
ProstateCancer, SkinCancer

### Part 1: Preparing the input files for training the classifiers

#### 1.1. Curating protein-protein interaction network

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Network_build/STRING_PPI.R">Scripts/Network_build/STRING_PPI.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>9606.protein.links.detailed.v11.5.txt.gz (STRING)</li>
<li>9606.protein.aliases.v11.5.txt.gz (STRING)</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>STRING_PPI_Net_database.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>Creates a protein-protein interaction network from interactions
reported in the STRING database. The script retrieves the source
interaction file for human from STRING directly. Only interactions
supported by experiments or those curated in string from other databases
are considered.</td>
</tr>
</tbody>
</table>

#### 1.2. Retrieving associations between drugs, targets and diseases

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 91%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Associations/DrugBank_Drug_Target_Associations.R">Scripts/Associations/DrugBank_Drug_Target_Associations.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td>drugbank_all_full_database.xml.zip</td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>DrugBank_Drug_Target_Net.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script retrieves the drug target association for the drugs in
DrugBank. Only human targets are reported. The output is a dataframe
containing drugs and its corresponding targets in each row.</td>
</tr>
<tr class="odd">
<td><strong>Troubleshoot</strong></td>
<td><p>Missing files: The DrugBank database as XML should be separately
downloaded. Permission from DrugBank is required to download the
file.</p>
<p><code>curl -Lfv -o filename.zip -u EMAIL:PASSWORD https://go.drugbank.com/releases/5-1-8/downloads/all-full-database</code></p></td>
</tr>
</tbody>
</table>

|              |                                                                                                                              |
|--------------|------------------------------------------------------------------------------------------------------------------------------|
| **R Script** | [Scripts/Associations/OpenTargets_Drug_Disease_Associations.R](Scripts/Associations/OpenTargets_Drug_Disease_Associations.R) |
| **Input**    | Core annotation for drug molecules (OpenTarget parquest files)                                                               |
| **Output**   | OpenTargets_Drug_Disease_Net.rds                                                                                             |
| **Summary**  | The script retrieves the drug disease associations from OpenTargets database.                                                |

#### 1.3. Curating the drug combinations for training the classifiers

<table>
<colgroup>
<col style="width: 5%" />
<col style="width: 94%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Reference_lists/FIMMDrugComb_drugCombinations.R">Scripts/Reference_lists/FIMMDrugComb_drugCombinations.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td>summary_v_1_5.csv (DrugComb.org)</td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>FimmDrugComb_drugCombinations.rds</li>
<li>FimmDrugComb_[disease]_drugCombinations.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The scripts retrieves the drug pairs and the synergy scores from the
FIMM DrugComb portal. The pairs are then categorized as synergistic or
antagonistic based on the synergy scores. Drug pairs are filtered to
keep only those tested on cancer related cell lines.</td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col style="width: 9%" />
<col style="width: 90%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Reference_lists/DrugBank_drugInteractions.R">Scripts/Reference_lists/DrugBank_drugInteractions.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>drug_drug_interactions.csv</li>
<li>drug.csv</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>DrugBank_drugInteractions.rds</li>
<li>DrugBank_drugInteractions_withRiskSeverity.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script retrieves drug-drug interactions from DrugBank that have
been reported to be involved with risk or severity to any side
effects.</td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Drug_combinations/DrugCombinations.R">Scripts/Drug_combinations/DrugCombinations.R</a></td>
</tr>
<tr class="even">
<td><strong>Parameters</strong></td>
<td>--disease</td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>DrugBank_drugInteractions_withRiskSeverity.rds</li>
<li>DrugBank_Drug_Target_Net.rds</li>
<li>FimmDrugComb_[disease]_drugCombinations.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td>DrugComb_[disease].rds</td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script prepares the effective and adverse drug combinations for
each cancer type. The synergistic drug pairs retrieved from FIMM
DrugComb portal are considered as effective. While the DrugBank
drug-drug interaction pairs with risk/severity indications that contain
both drugs involved in forming effective pairs are considered as
adverse. Drug pairs are filtered to keep only those with reported
targets</td>
</tr>
</tbody>
</table>

|              |                                                                                                                        |
|--------------|------------------------------------------------------------------------------------------------------------------------|
| **R Script** | [Scripts/Drug_combinations/DrugCombs_training_annotation.R](Scripts/Drug_combinations/DrugCombs_training_annotation.R) |
| **Input**    | DrugComb\_\[disease\].rds                                                                                              |
| **Output**   | DrugCombs_training_annotation.xlsx                                                                                     |
| **Summary**  | The script summarises the drug pairs for the training set. Also includes DrugBank drug ID to drug name mapping.        |

#### 1.4. Preparing the libraries for estimation of efficacy and safety

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 91%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_ADR.R">Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_ADR.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td>ADRAlert2GENE2ID.xlsx</td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>ADReCS_ADR2Gene_level4_lib.rds</li>
<li>ADReCS_ADR2Gene_level3_lib.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script collates gene sets related to side effects from ADReCs.
The gene sets for the two different levels of side effects are exported
separately.</td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 91%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Diseases.R">Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Diseases.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>curated_gene_disease_associations.tsv.gz (DisGeNET)</li>
<li>Data-type metrics for direct target-disease associations (OpenTarget
parquest files)</li>
<li>Compendium_Cancer_Genes.tsv (Intogen)</li>
<li>cohorts.tsv (Intogen)</li>
<li>Disease_Signatures_from_GEO_up_2014.txt (Enrichr)</li>
<li>Disease_Signatures_from_GEO_down_2014.txt (Enrichr)</li>
<li>CTD_genes_diseases.csv.gz (Comparative Toxicogenomics Database)</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>DisGeNET_Disease2Gene_lib.rds</li>
<li>DisGeNET_DiseaseGroup2Gene_lib.rds</li>
<li>DisGeNET_Phenotype2Gene_lib.rds</li>
<li>OpenTargets_Disease2Gene_GA_lib.rds</li>
<li>OpenTargets_Disease2Gene_RNA_lib.rds</li>
<li>OpenTargets_Disease2Gene_lit_lib.rds</li>
<li>Intogen_Disease2Gene_lib.rds</li>
<li>Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds</li>
<li>PharmGKB_Disease2Gene_lib.rds</li>
<li>CTD_Disease2Gene_lib.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script collates gene sets related to diseases from several
databases like DisGeNET, OpenTargets, Intogen, Enrichr, PharmGKB.</td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 92%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_DrugWithdrawalADR.R">Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_DrugWithdrawalADR.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>ADReCS_ADR2Gene_level4_lib.rds</li>
<li>ADReCS_ADR2Gene_level3_lib.rds</li>
<li>CTD_Disease2Gene_lib.rds</li>
<li>DisGeNET_Disease2Gene_lib.rds</li>
<li>DisGeNET_DiseaseGroup2Gene_lib.rds</li>
<li>DisGeNET_Phenotype2Gene_lib.rds</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>drugWithdrawal_Adr2Gene_lib.rds</li>
<li>drugWithdrawal_Adr2Gene_libInfo.xlsx</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script prepares gene set library for adverse drug reactions
related to drug withdrawal.</td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 91%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Miscellaneous.R">Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Miscellaneous.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>categories.tsv (DGIdb)</li>
<li>genes.tsv (PharmGKB)</li>
<li>genage_human.csv (GenAge db)</li>
<li>Halifax-curation.Table S2. Hallmark of Cancer Data - Gene and
Pathway.xlsx</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>miscellaneous_gene_lib.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script compiles gene sets related to druggable geneome, very
important pharmacogenes, ageing, and cancer hallmarks into a single
library.</td>
</tr>
<tr class="odd">
<td><strong>Troubleshoot</strong></td>
<td>Missing files: The Halifax-curation.Table S2. Hallmark of Cancer
Data - Gene and Pathway.xlsx file must be downloaded manually from the
journal website <a href="https://doi.org/10.1093/database/baaa045"
class="uri">https://doi.org/10.1093/database/baaa045</a></td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 91%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Pathways.R">Scripts/Enrichment_Analysis_Libraries/Enrichment_Analysis_Libraries_Pathways.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>Cancer hallmark related pathways (Table_1.xls from <a
href="https://doi.org/10.3389%2Ffgene.2020.00029">PMC7013921</a>)</li>
<li>smpdb_proteins.csv.zip</li>
<li>smpdb_pathways.csv.zip</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>CHG_keggPath2Gene_lib.rds</li>
<li>SMPDb_Pathway2Gene_lib.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script retrieves gene sets related to cancer hallmark linked
cancer pathways and small molecular pathways from SMPDB.</td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col style="width: 9%" />
<col style="width: 90%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><p><a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_BreastCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_BreastCancer_lib.R</a></p>
<p><a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_KidneyCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_KidneyCancer_lib.R</a></p>
<p><a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_LungCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_LungCancer_lib.R</a></p>
<p><a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_OvaryCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_OvaryCancer_lib.R</a></p>
<p><a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_ProstateCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_ProstateCancer_lib.R</a></p>
<p><a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_SkinCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_SkinCancer_lib.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>DisGeNET_Disease2Gene_lib.rds</li>
<li>OpenTargets_Disease2Gene_GA_lib.rds</li>
<li>OpenTargets_Disease2Gene_RNA_lib.rds</li>
<li>OpenTargets_Disease2Gene_lit_lib.rds</li>
<li>Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds</li>
<li>Intogen_Disease2Gene_lib.rds</li>
<li>PharmGKB_Disease2Gene_lib.rds</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>Disease2Gene_[disease]_lib.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The scripts compiles disease specific gene set libraries from
different databases and TheTA.</td>
</tr>
</tbody>
</table>

### Part 2: Preparing the features to be used in training the classifiers

#### 2.1. Execute RWR on the drug combinations

|              |                                                |
|--------------|------------------------------------------------|
| **R Script** | Scripts/Feature_generation/DrugCombs_dNetRWR.R |
| **Input**    |                                                |
| **Output**   |                                                |
| **Summary**  |                                                |

#### 2.2. Execute FGSEA on the ranked gene list from RWR

<table>
<colgroup>
<col style="width: 19%" />
<col style="width: 80%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><p>Scripts/Feature_generation/DrugCombs_fgsea_CommonLibraries.R</p>
<p>Scripts/Feature_generation/DrugCombs_fgsea_EfficacySafety.R</p></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td></td>
</tr>
</tbody>
</table>

#### 2.3. Calculate network proximity

<table>
<colgroup>
<col style="width: 17%" />
<col style="width: 82%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><p>Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugAdr.</p>
<p>RScripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDisease.R</p>
<p>Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDrug.R</p></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td></td>
</tr>
</tbody>
</table>

### Part 3: Training and validation of the classifiers

#### 3.1. Prepare data for repeated cross validation

|              |                                                  |
|--------------|--------------------------------------------------|
| **R Script** | Scripts/Model_train/DrugCombs_trainML_rwrFgsea.R |
| **Input**    |                                                  |
| **Output**   |                                                  |
| **Summary**  |                                                  |

#### 3.2. Train and validate the models

<table>
<colgroup>
<col style="width: 17%" />
<col style="width: 82%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><p>Scripts/Model_train/DrugCombs_trainML_rwrFgsea.R</p>
<p>Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugAdr.R</p>
<p>Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDisease.R</p>
<p>Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDrug.R</p>
<p>Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDiseaseAdr.R</p></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td></td>
</tr>
</tbody>
</table>

### Part 4: Summarizing the results
