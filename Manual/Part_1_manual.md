Part 1: Preparing the data
================
2024-07-08

This part deals with the scripts that retrieved and process the data
that form the base of our framework. Here we provide details of the
scripts used for building the protein-protein interaction network,
parsing drug related data from DrugBank and curating the drug
combination data set from FIMM DrugComb portal.

### 1.1. Building the network

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 97%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Network_build/STRING_PPI.R">Scripts/Network_build/STRING_PPI.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td>9606.protein.links.detailed.v12.0.txt.gz</td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ul>
<li><p>STRING_PPI_Net.rds</p></li>
<li><p>STRING_PPI_Net_params.csv</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>Creates a protein-protein interaction network from interactions
reported in the <a href="https://string-db.org/">STRING</a> database.
The script retrieves the source interaction file for human from STRING
directly. The interactions with combined score in STRING greater than
400 are used to build the network. Only the largest connected component
is retained. The output is an igraph network stored in an rds file.
Additional, a csv file with the basic properties of the network is also
generated.</td>
</tr>
</tbody>
</table>

### 1.2. Retrieving the drug target associations

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 93%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Associations/DrugBank_Drug_Target_Associations.R">Scripts/Associations/DrugBank_Drug_Target_Associations.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td>drugbank_all_full_database.xml.zip</td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>DrugBank_Drug_Target_associations.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td><p>The script extracts the drug-target associations from <a
href="https://go.drugbank.com/">DrugBank</a> and stores it as a R data
frame with two columns mentioning the DrugBank drug IDs and the Ensembl
gene IDs for the targets.</p>
<p>Note: For this script to be successfully executed, the input file
must be manually downloaded from DrugBank (requires access permission)
and placed in <code>Databases/DrugBank/</code> directory.</p></td>
</tr>
</tbody>
</table>

### 1.3. Processing the FIMM drug combination data

<table>
<colgroup>
<col style="width: 5%" />
<col style="width: 94%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Reference_lists/FIMMDrugComb_drugCombinations.R">Scripts/Reference_lists/FIMMDrugComb_drugCombinations.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ul>
<li><p>summary_v_1_5.csv</p></li>
<li><p>Thesaurus.FLAT.zip</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>FimmDrugComb_drugCombinations.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script downloads the data from <a
href="https://drugcomb.fimm.fi/">FIMM DrugComb portal</a>, extracts the
drug combinations that have been tested on cancer cell lines and assigns
a cancer-specific synergy level by utilising the synergy scores reported
in the portal.</td>
</tr>
</tbody>
</table>

### 1.4. Retrieve the DrugBank drug-drug interactions

|             |                                                                                                                                                                                                                                |
|-------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Rscript** | [Scripts/Reference_lists/DrugBank_DDI.R](../Scripts/Reference_lists/DrugBank_DDI.R)                                                                                                                                            |
| **Input**   | drugbank_all_full_database.xml.zip                                                                                                                                                                                             |
| **Output**  | DrugBank_DDI_processed.rds                                                                                                                                                                                                     |
| **Summary** | The script extracts the drug-drug interactions reported in [DrugBank](https://go.drugbank.com/) and categorises drug pairs into various classes including metabolic effect, therapeutic effect and cancer-specific ADR status. |

### 1.5. Extract the targets of the drug combinations

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 97%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Drug_combinations/Extract_drugCombs_and_targets.R">Scripts/Drug_combinations/Extract_drugCombs_and_targets.R</a></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td>Mandatory: <code>--disease</code></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>FimmDrugComb_drugCombinations.rds</p></li>
<li><p>STRING_PPI_Net.rds</p></li>
<li><p>DrugBank_Drug_Target_associations.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>Drug_combination_data/drugCombs_data_[disease].rds</p></li>
<li><p>drugCombs_targets_extended_[disease].rds</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script first subsets the drug combinations for the cancer type
specified using <code>--disease</code> argument from the data retrieved
from <a href="https://drugcomb.fimm.fi/">FIMM DrugComb portal</a> (see
section 1.3). It then extracts all the known targets of the drugs
involved in the combination. Additionally, it also attempts to extend
the list of targets by considering the immediate neighbours of the known
targets in various kinds of functionally linked protein-protein
interaction networks retrieved from <a
href="https://omnipathdb.org/">OmniPath</a> database. Two output file
are generated. The first one contains a R data frame of the with the
list of drug combinations subsetted from processed DrugComb portal data.
The second file also contains a data frame with information about the
targets.</td>
</tr>
</tbody>
</table>
