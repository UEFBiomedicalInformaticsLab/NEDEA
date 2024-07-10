Part 5: De novo predictions
================
2024-07-08

This part describes the scripts used to identify novel anti-cancer drug
combinations from existing licensed anti-cancer drugs using the
predictive system developed in the study.

### 5.1. Data for *de novo* prediction

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/DeNovo_predictions/DeNovo_data_1.R">Scripts/DeNovo_predictions/DeNovo_data_1.R</a></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td>Mandatory: <code>--disease</code></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>parsed_DrugBank_data.rds</p></li>
<li><p>cancerdrugsdb.txt (Cancer Drug Database)</p></li>
<li><p>drugCombs_cat_effVadv_[disease].rds</p></li>
<li><p>09.04.2024.zip (C-DCDB)</p></li>
<li><p>STRING_PPI_Net.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td>drugCombs_denovo1_[disease].rds</td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script generates an initial pool of drug pairs which will be
scouted to identify the novel potential effective and safe drug
combinations. The drug pairs are formed by combining all possible pairs
of licensed anti-cancer drugs followed by filtering our those that were
included in the training set, have been into clinical trials for the
specific cancer or have reported adverse drug-drug interactions.</td>
</tr>
</tbody>
</table>

### 5.2. Feature generation on *de novo* data

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 97%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/DeNovo_predictions/DeNovo_featureGeneration_1.R">Scripts/DeNovo_predictions/DeNovo_featureGeneration_1.R</a></td>
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
<li><p>drugCombs_denovo1_[disease].rds</p></li>
<li><p>NES_EfficacySafety_selectedFeatures_[disease]_[drug_target_type].csv</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>fgseaNES_combinedEfficacySafety_[disease]_[drug_target_type].rds
(in DeNovo_data_1 directory)</p></li>
<li><p>plot_DeNovo_1_AD_combinedEfficacySafety_[disease]_[drug_target_type].tiff</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The scripts takes in the drug pairs generated in section 5.1 and
implements the RWR-FGSEA (as described in part 3) framework to return a
data frame of the NES for the drug combinations in column against the
gene sets in rows. However, since our predictive system is based on the
combined efficacy-safety library, the NES for the <em>de novo</em> data
is generated only for the same. Additionally, the scripts plot the
applicibility domain of the <em>de novo</em> drug combinations w.r.t the
training set drug combinations by considering one efficacy and one
safety feature for the scatter plot.</td>
</tr>
</tbody>
</table>

### 5.3. Prediction on *de novo* data

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/DeNovo_predictions/DeNovo_prediction_1.R">Scripts/DeNovo_predictions/DeNovo_prediction_1.R</a></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td><p>Mandatory: <code>--disease</code></p>
<p>Optional: <code>--drug_target_type</code>,
<code>--select_top_combs</code></p></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>parsed_DrugBank_data.rds</p></li>
<li><p>model_NES_combinedEfficacySafety_[disease]_[drug_target_type].rds</p></li>
<li><p>fgseaNES_combinedEfficacySafety_[disease]_[drug_target_type].rds
(from DeNovo_data_1 directory)</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>predictions_NES_combinedEfficacySafety_[disease]_[drug_target_type].csv</p></li>
<li><p>priorityDrugCombs_NES_combinedEfficacySafety_[disease]_[drug_target_type].csv</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script uses the predictive system developed in section 4.2 and
the efficacy and safety estimates calculated in section 5.2 to generate
predictions for the <em>de novo</em> data set drug combinations. The
primary output is a data frame with the final score and the predicted
categories along with other drug information. The script also identifies
the top n (set using <code>select_top_combs</code> argument) candidates
to be prioritised for experimental validation.</td>
</tr>
</tbody>
</table>
