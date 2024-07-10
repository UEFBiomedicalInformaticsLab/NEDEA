Part 2: Preparing the enrichment libraries
================
2024-07-08

This part deals with the scripts that are used to download associations
of genes to adverse effects, diseases, pathways and other relevant
phenotypes to build gene set libraries that are used for performing
enrichment analysis or for calculating proximity within our framework.
Here, library refer to a group of gene sets. The libraries are stored as
R named list with list element being the gene sets. The genes are
represented by their Ensembl IDs.

### 2.1. Retrieving the databases for curating the disease and adverse effect gene sets

<table>
<colgroup>
<col style="width: 5%" />
<col style="width: 94%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><p><a
href="../Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_ADR.R">Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_ADR.R</a></p>
<p><a
href="../Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_Diseases.R">Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_Diseases.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ul>
<li><p>ADRAlert2GENE2ID.xlsx (ADReCS)</p></li>
<li><p>curated_gene_disease_associations.tsv.gz (DisGeNet)</p></li>
<li><p>associationByDatatypeDirect (Open Targets)</p></li>
<li><p>IntOGen-Drivers-20230531.zip (Intogen)</p></li>
<li><p>Disease_Signatures_from_GEO_up_2014.txt (Enrichr)</p></li>
<li><p>relationships.zip (PharmGKB)</p></li>
<li><p>CTD_genes_diseases.csv.gz (Comparative Toxicogenomics
Database)</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ul>
<li><p>ADReCS_ADR2Gene_level4_lib.rds</p></li>
<li><p>ADReCS_ADR2Gene_level3_lib.rds</p></li>
<li><p>DisGeNET_Disease2Gene_lib.rds</p></li>
<li><p>DisGeNET_Phenotype2Gene_lib.rds</p></li>
<li><p>DisGeNET_DiseaseGroup2Gene_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_GA_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_RNA_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_lit_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_SM_lib.rds</p></li>
<li><p>Intogen_Disease2Gene_lib.rds</p></li>
<li><p>Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds</p></li>
<li><p>PharmGKB_Disease2Gene_lib.rds</p></li>
<li><p>CTD_Disease2Gene_lib.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The scripts retrieves associations of genes to adverse effect and
diseases as curated in different databases. Where ever possible, the
types or the sources of the associations have been segregated and are
stored separately.</td>
</tr>
</tbody>
</table>

### 2.2. Preparing disease specific libraries

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 91%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><p><a
href="../Scripts/Enrichment_analysis_libraries/Disease2Gene_BreastCancer_lib.R">Scripts/Enrichment_analysis_libraries/Disease2Gene_BreastCancer_lib.R</a></p>
<p><a
href="../Scripts/Enrichment_analysis_libraries/Disease2Gene_KidneyCancer_lib.R">Scripts/Enrichment_analysis_libraries/Disease2Gene_KidneyCancer_lib.R</a></p>
<p><a
href="../Scripts/Enrichment_analysis_libraries/Disease2Gene_LungCancer_lib.R">Scripts/Enrichment_analysis_libraries/Disease2Gene_LungCancer_lib.R</a></p>
<p><a
href="../Scripts/Enrichment_analysis_libraries/Disease2Gene_OvaryCancer_lib.R">Scripts/Enrichment_analysis_libraries/Disease2Gene_OvaryCancer_lib.R</a></p>
<p><a
href="../Scripts/Enrichment_analysis_libraries/Disease2Gene_ProstateCancer_lib.R">Scripts/Enrichment_analysis_libraries/Disease2Gene_ProstateCancer_lib.R</a></p>
<p><a
href="../Scripts/Enrichment_analysis_libraries/Disease2Gene_SkinCancer_lib.R">Scripts/Enrichment_analysis_libraries/Disease2Gene_SkinCancer_lib.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ul>
<li><p>DisGeNET_Disease2Gene_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_GA_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_RNA_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_lit_lib.rds</p></li>
<li><p>OpenTargets_Disease2Gene_SM_lib.rds</p></li>
<li><p>Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds</p></li>
<li><p>Intogen_Disease2Gene_lib.rds</p></li>
<li><p>PharmGKB_Disease2Gene_lib.rds</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ul>
<li><p>Disease2Gene_BreastCancer_lib.rds</p></li>
<li><p>Disease2Gene_KidneyCancer_lib.rds</p></li>
<li><p>Disease2Gene_LungCancer_lib.rds</p></li>
<li><p>Disease2Gene_OvaryCancer_lib.rds</p></li>
<li><p>Disease2Gene_ProstateCancer_lib.rds</p></li>
<li><p>Disease2Gene_SkinCancer_lib.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The scripts prepares six cancer type specific gene set libraries
using the gene disease associations retrieved in section 2.1.</td>
</tr>
</tbody>
</table>

### 2.3. Prepare the adverse effect linked gene set library

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 95%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_curatedADR.R">Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_curatedADR.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ul>
<li><p>ADReCS_ADR2Gene_level3_lib.rds</p></li>
<li><p>ADReCS_ADR2Gene_level4_lib.rds</p></li>
<li><p>CTD_Disease2Gene_lib.rds</p></li>
<li><p>DisGeNET_Disease2Gene_lib.rds</p></li>
<li><p>DisGeNET_DiseaseGroup2Gene_lib.rds</p></li>
<li><p>DisGeNET_Phenotype2Gene_lib.rds</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ul>
<li><p>curatedAdr2Gene_lib.rds</p></li>
<li><p>curatedAdr2Gene_libInfo.xlsx</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script curates gene set library with the genes linked to adverse
drug reactions that are frequently associated with drug withdrawal from
the market and other relevant organ toxicities. The script also
calculates the Jaccard similarity between the gene sets in this
library,</td>
</tr>
</tbody>
</table>

### 2.4. Retrieving and prepare pathway related gene set libraries

<table>
<colgroup>
<col style="width: 1%" />
<col style="width: 98%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_Pathways.R">Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_Pathways.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ul>
<li><p>PMC7013921/Table_1.xls (PubMed Central)</p></li>
<li><p>smpdb_proteins.csv.zip (SMPDB)</p></li>
<li><p>smpdb_pathways.csv.zip (SMPDB)</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td><ul>
<li><p>CHG_keggPath2Gene_lib.rds</p></li>
<li><p>SMPDb_Pathway2Gene_lib.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td>The script curates two library consisting of pathways from two
different sources. The first library consists of cancer hallmark
associated KEGG pathways. The list of characteristic KEGG pathways of
different cancer hallmarker were retrieved from the research article <a
href="https://doi.org/10.3389/fgene.2020.00029">PMC7013921</a> while the
genes in the KEGG pathway were sourced from the <a
href="https://www.genome.jp/kegg/">KEGG database</a>. The second library
is based on the pathways curated in the <a
href="https://www.smpdb.ca/">Small Molecule Pathway Database
(SMPDB)</a>. The library consists of two sub-libraries — one for the
pathways associated to drug action while the other for pathways linked
to drug metabolism.</td>
</tr>
</tbody>
</table>

### 2.5. Retrieve and prepare library consisting of miscellaneous gene sets

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_Miscellaneous.R">Scripts/Enrichment_analysis_libraries/Enrichment_Analysis_Libraries_Miscellaneous.R</a></td>
</tr>
<tr class="even">
<td><strong>Input</strong></td>
<td><ul>
<li><p>categories.tsv (DGIdb)</p></li>
<li><p>genes.zip (PharmGKB)</p></li>
<li><p>human_genes.zip (GenAge)</p></li>
<li><p>suppl_data_baaa045.zip (PubMed)</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Output</strong></td>
<td>miscellaneous_gene_lib.rds</td>
</tr>
<tr class="even">
<td><strong>Summary</strong></td>
<td><p>The script curates a library by combining four independent gene
sets — genes in drugable geneome from <a
href="https://www.dgidb.org/">Drug Gene Interaction Database
(DGIdb)</a>, very important pharmacogenes from <a
href="https://www.pharmgkb.org/">PharmGKB</a>, genes related to ageing
from <a
href="https://genomics.senescence.info/genes/index.html">GenAge</a>
database and hallmarkers of cancer (HOC) genes from HOC database.</p>
<p>Note: To compile the HOC genes, the supplementary file
‘Halifax-curation.Table S2. Hallmark of Cancer Data - Gene and
Pathway.xlsx’ must be manually downloaded from the research article <a
href="https://doi.org/10.1093/database/baaa045">PMC7294774</a> and
placed inside the <code>Databases/HOCdb/</code> directory.</p></td>
</tr>
</tbody>
</table>
