Part 3: Feature generation and exploration
================
2024-07-08

This part describes the scripts that are used to execute the RWR-FGSEA
framework to obtain the normalised enrichment scores and calculation of
network proximity between the drug targets to various gene sets.

### 3.1. Execute RWR

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 97%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Feature_generation/Execute_RWR_on_drugCombs.R">Scripts/Feature_generation/Execute_RWR_on_drugCombs.R</a></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td><p>Mandatory: <code>--disease</code></p>
<p>Optional: <code>--drug_target_type</code>,
<code>--nproc</code></p></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>STRING_PPI_Net.rds</p></li>
<li><p>drugCombs_targets_extended_[disease].rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td>rwrProbs_[disease]_[drug_target_type].rds</td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script executes random walk with restart (RWR) for each drug
combination using its targets as the seeds. By default, the known drug
targets are used but it can be modified to use the extended drug targets
using the <code>drug_target_type</code> argument. The output is a sparse
matrix containing the visiting probabilities of all nodes in the input
network. The number of rows in the matrix equals the number of nodes in
the input network while the number of columns equals the number of drug
combinations for which the RWR was executed.</td>
</tr>
</tbody>
</table>

### 3.2. Calculate NES

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><p><a
href="../Scripts/Feature_generation/Execute_FGSEA_EfficacySafety_on_drugCombs.R">Scripts/Feature_generation/Execute_FGSEA_EfficacySafety_on_drugCombs.R</a></p>
<p><a
href="../Scripts/Feature_generation/Execute_FGSEA_Pathways_on_drugCombs.R">Scripts/Feature_generation/Execute_FGSEA_Pathways_on_drugCombs.R</a></p>
<p><a
href="../Scripts/Feature_generation/Execute_FGSEA_Miscellaneous_on_drugCombs.R">Scripts/Feature_generation/Execute_FGSEA_Miscellaneous_on_drugCombs.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td><p>Mandatory: <code>--disease</code></p>
<p>Optional: <code>--drug_target_type</code>,
<code>--nproc</code></p></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>rwrProbs_[disease]_[drug_target_type].rds</p></li>
<li><p>Disease2Gene_[disease]_lib.rds</p></li>
<li><p>curatedAdr2Gene_lib.rds</p></li>
<li><p>CHG_keggPath2Gene_lib.rds</p></li>
<li><p>SMPDb_Pathway2Gene_lib.rds</p></li>
<li><p>miscellaneous_gene_lib.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>fgseaNES_EfficacySafety_[disease]_[drug_target_type].rds</p></li>
<li><p>fgseaNES_Pathway_[disease]_[drug_target_type].rds</p></li>
<li><p>fgseaNES_Miscellaneous_[disease]_[drug_target_type].rds</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The scripts performs fast gene set enrichment analysis (FGSEA) of
the genes selected from RWR against the mentioned libraries. The script
internally uses an elbow method on the visiting probability of nodes
obtained from RWR (section 3.1) to identify the genes to be used for
FGSEA. The output is a list of data frames for each feature type. The
data frames contains normalised enrichment scores for the drug
combinations in columns against the gene sets in rows.</td>
</tr>
</tbody>
</table>

### 3.3. Calculate network proximity

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><p><a
href="../Scripts/Feature_generation/Calculate_proximity_EfficacySafety_for_drugCombs.R">Scripts/Feature_generation/Calculate_proximity_EfficacySafety_for_drugCombs.R</a></p>
<p><a
href="../Scripts/Feature_generation/Calculate_proximity_EfficacySafety_for_drugCombs.R">Scripts/Feature_generation/Calculate_proximity_Pathways_for_drugCombs.R</a></p>
<p><a
href="../Scripts/Feature_generation/Calculate_proximity_Miscellaneous_for_drugCombs.R">Scripts/Feature_generation/Calculate_proximity_Miscellaneous_for_drugCombs.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td><p>Mandatory: <code>--disease</code></p>
<p>Optional: <code>--drug_target_type</code>,
<code>--nproc</code></p></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>STRING_PPI_Net.rds</p></li>
<li><p>Disease2Gene_[disease]_lib.rds</p></li>
<li><p>curatedAdr2Gene_lib.rds</p></li>
<li><p>CHG_keggPath2Gene_lib.rds</p></li>
<li><p>SMPDb_Pathway2Gene_lib.rds</p></li>
<li><p>miscellaneous_gene_lib.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>netProx_EfficacySafety_[disease]_[drug_target_type].rds</p></li>
<li><p>netProx_Pathway_[disease]_[drug_target_type].rds</p></li>
<li><p>netProx_Miscellaneous_[disease]_[drug_target_type].rds</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The scripts calculates network based proximities between the targets
of the drug combination and the genes sets of the mentioned feature. The
primary output is a data frame containing the proximity values of the
drug combinations in columns to the gene sets in rows. The data frame is
stored as list of list, with the first level being the feature type and
the second level being the proximity type.</td>
</tr>
</tbody>
</table>
