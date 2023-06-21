
# DrugCombination_1

This repository contains the scripts used to train classifiers for
prioritization of drug combinations for cancer based on efficacy and
safety estimates.

## Preparing the environment

The following packages must be available for executing the scripts:
biomaRt; dbparser; doParallel; foreach; httr; igraph; jsonlite; msigdbr;
openxlsx; optparse; org.Hs.eg.db; readxl; sparklyr; sparklyr.nested;
tidyverse; unixtools

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

**--folds:** Number of folds to be created. Default: 5

**--repeats:** Number of repeats of fold splitting to be performed.
Default: 5

**--data_balance_method:** The method to be used to balance imbalanced
data. Possible values: none. Default: none.

**--proximity:** The proximity type to use. Possible values: closest,
shortest, centre, kernel, separation. Default: separation

**--nproc:** Number of processes to use. Default: NULL

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
<col style="width: 4%" />
<col style="width: 95%" />
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

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Associations/OpenTargets_Drug_Disease_Associations.R">Scripts/Associations/OpenTargets_Drug_Disease_Associations.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td>Core annotation for drug molecules (OpenTarget parquet files)</td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>OpenTargets_Drug_Disease_Net.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script retrieves the drug disease associations from OpenTargets database.</td>
</tr>
</tbody>
</table>

#### 1.3. Curating the drug combinations for training the classifiers

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
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
<col style="width: 4%" />
<col style="width: 95%" />
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
<col style="width: 4%" />
<col style="width: 95%" />
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

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Drug_combinations/DrugCombs_training_annotation.R">Scripts/Drug_combinations/DrugCombs_training_annotation.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td>DrugComb_[disease].rds</td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>DrugCombs_training_annotation.xlsx </td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script summarises the drug pairs for the training set. Also includes DrugBank drug ID to drug name mapping.</td>
</tr>
</tbody>
</table>

#### 1.4. Preparing the libraries for estimation of efficacy and safety

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
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
<col style="width: 4%" />
<col style="width: 95%" />
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
<col style="width: 4%" />
<col style="width: 95%" />
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
<col style="width: 4%" />
<col style="width: 95%" />
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
<col style="width: 4%" />
<col style="width: 95%" />
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
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_BreastCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_BreastCancer_lib.R</a>
<a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_KidneyCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_KidneyCancer_lib.R</a>
<a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_LungCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_LungCancer_lib.R</a>
<a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_OvaryCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_OvaryCancer_lib.R</a>
<a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_ProstateCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_ProstateCancer_lib.R</a>
<a
href="Scripts/Enrichment_Analysis_Libraries/Disease2Gene_SkinCancer_lib.R">Scripts/Enrichment_Analysis_Libraries/Disease2Gene_SkinCancer_lib.R</a></td>
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

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Feature_generation/DrugCombs_dNetRWR.R">Scripts/Feature_generation/DrugCombs_dNetRWR.R</a></td>
</tr>
<tr class="even">
<td><strong>Parameters</strong></td>
<td>--disease, --nproc</td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>DrugComb_[disease].rds</li>
<li>DrugBank_Drug_Target_Net.rds</li>
<li>STRING_PPI_Net_database.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td>dNetRWR050_[disease].rda</td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script executes a random walk with restart on the input network
for the effective and adverse drug combinations. The known targets of
the drugs are used as seed for the RWR. The output is a probability
score for the genes/proteins in the network for each drug
combination.</td>
</tr>
</tbody>
</table>

#### 2.2. Execute FGSEA on the ranked gene list from RWR

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><p><a
href="Scripts/Feature_generation/DrugCombs_fgsea_CommonLibraries.R">Scripts/Feature_generation/DrugCombs_fgsea_CommonLibraries.R</a></p>
<p><a
href="Scripts/Feature_generation/DrugCombs_fgsea_EfficacySafety.R">Scripts/Feature_generation/DrugCombs_fgsea_EfficacySafety.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Parameters</strong></td>
<td>--disease, --nproc</td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>dNetRWR050_[disease].rda</li>
<li>CHG_keggPath2Gene_lib.rds</li>
<li>SMPDb_Pathway2Gene_lib.rds</li>
<li>SMPDb_Pathway2Gene_lib.rds</li>
<li>miscellaneous_gene_lib.rds</li>
<li>Disease2Gene_[disease]_lib.rds</li>
<li>drugWithdrawal_Adr2Gene_lib.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>fgseaProbCut_CommonLib_[disease].rds</li>
<li>fgseaProbCut_EfficacySafety_[disease].rds</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>For each drug combination, the ranked gene/protein list obtained
from RWR is used for fast gene set enrichment analysis. Only the genes
that occurs above the 90th quantile probability score is used for
FGSEA.</td>
</tr>
</tbody>
</table>

#### 2.3. Calculate network proximity

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><p><a
href="Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugAdr.R">Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugAdr.R</a></p>
<p><a
href="Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDisease.R">Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDisease.R</a></p>
<p><a
href="Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDrug.R">Scripts/Feature_generation/DrugCombs_BarabasiMetrics_DrugDrug.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Parameters</strong></td>
<td>--disease, --nproc</td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>DrugComb_[disease].rds</li>
<li>DrugBank_Drug_Target_Net.rds</li>
<li>drugWithdrawal_Adr2Gene_lib.rds</li>
<li>Disease2Gene_[disease]_lib.rds</li>
<li>STRING_PPI_Net_database.rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>BarabasiProx_DrugAdr_[disease].rds</li>
<li>BarabasiProx_DrugDisease_[disease].rds</li>
<li>BarabasiProx_DrugDrug_[disease].rds</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script calculates various types of network proximity between
drug targets, disease associated genes and adverse drug reaction related
genes.</td>
</tr>
</tbody>
</table>

### Part 3: Training and validation of the classifiers

#### 3.1. Prepare data for repeated cross validation

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><a
href="Scripts/Model_train/DrugCombs_ML_dataSplit.R">Scripts/Model_train/DrugCombs_ML_dataSplit.R</a></td>
</tr>
<tr class="even">
    <td><strong>Parameters</strong></td>
    <td>--disease, --folds, --repeats</td>
    </tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td>DrugComb_[disease].rds</td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td>ML_dataSplit_[disease].rds</td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script segregates the data into training and test set by splitting it into k folds and n repeats. For each turn, n-1 fold is used for training.</td>
</tr>
</tbody>
</table>

#### 3.2. Train and validate the models

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>R Script</strong></td>
<td><p><a
href="Scripts/Model_train/DrugCombs_trainML_rwrFgsea.R">Scripts/Model_train/DrugCombs_trainML_rwrFgsea.R</a></p>
<p><a
href="Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugAdr.R">Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugAdr.R</a></p>
<p><a
href="Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDisease.R">Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDisease.R</a></p>
<p><a
href="Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDrug.R">Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDrug.R</a></p>
<p><a
href="Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDiseaseAdr.R">Scripts/Model_train/DrugCombs_trainML_BarabasiMetrics_DrugDiseaseAdr.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Parameters</strong></td>
<td>--disease, --data_balance_method, --proximity, --nproc</td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ol type="a">
<li>ML_dataSplit_[disease].rds</li>
<li>fgseaProbCut_EfficacySafety_[disease].rds</li>
<li>fgseaProbCut_CommonLib_[disease].rds</li>
<li>BarabasiProx_DrugAdr_[disease].rds</li>
<li>BarabasiProx_DrugDisease_[disease].rds</li>
<li>BarabasiProx_DrugDrug_[disease].rds</li>
</ol></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ol type="a">
<li>models_[data_balance_method]_rwrFgsea_[disease].rds</li>
<li>models_[data_balance_method]_BarabasiProx_DrugAdr_[disease].rds</li>
<li>models_[data_balance_method]_BarabasiProx_DrugDisease_[disease].rds</li>
<li>models_[data_balance_method]_BarabasiProx_DrgDisAdr_[disease]_[proximity].rds</li>
</ol></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The features generated in part 2 are used here to train four
different types of classifiers.</td>
</tr>
</tbody>
</table>

### Part 4: Summarizing the results
