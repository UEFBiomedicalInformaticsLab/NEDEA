Part 4: Training and validation of predictive system
================
2024-07-08

This part describes the scripts used to develop and validate the
predictive systems for classification of effective and adverse
anti-cancer drug combinations using its efficacy and safety estimates.

### 4.1. Feature selection

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 97%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Predictive_model/Feature_selection_EfficacySafety.R">Scripts/Predictive_model/Feature_selection_EfficacySafety.R</a></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td><p>Mandatory: <code>--disease</code></p>
<p>Optional: <code>--drug_target_type</code></p></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>fgseaNES_EfficacySafety_[disease]_[drug_target_type].rds</p></li>
<li><p>drugCombs_cat_effVadv_[disease].rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td>NES_EfficacySafety_selectedFeatures_[disease]_[drug_target_type].csv</td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script intakes the NES generated from FGSEA (section 3.2) and
performs a one-sided Wilcoxon’s test to identify the features that have
statistically significant (p&lt;0.001) differences in the distribution
of NES between the effective and adverse drug combinations. The script
has been set to use the NES for the combined library consisting the
disease and adverse effect related gene sets. For disease related
features, the test checks if the effective drug combinations have a
higher NES while for the adverse effect related features, it checks if
the adverse drug combinations have higher NES.</td>
</tr>
</tbody>
</table>

### 4.2. Development of predictive system

<table>
<colgroup>
<col style="width: 1%" />
<col style="width: 98%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><a
href="../Scripts/Predictive_model/train_predictiveSystem_combinedEfficacySafety_NES.R">Scripts/Predictive_model/train_predictiveSystem_combinedEfficacySafety_NES.R</a></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td><p>Mandatory: <code>--disease</code></p>
<p>Optional: <code>--drug_target_type</code></p></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>fgseaNES_EfficacySafety_[disease]_[drug_target_type].rds</p></li>
<li><p>drugCombs_cat_effVadv_[disease].rds</p></li>
<li><p>NES_EfficacySafety_selectedFeatures_[disease]_[drug_target_type].csv</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>model_characterization_[disease]_[drug_target_type].tiff</p></li>
<li><p>training_data_partition_[disease]_[drug_target_type].tiff</p></li>
<li><p>model_NES_combinedEfficacySafety_[disease]_[drug_target_type].rds</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script trains a predictive system to distinguish effective from
adverse drug combinations using the efficacy and safety estimates.
Ideally, for each feature, the script identifies a threshold that
partitions the effective and adverse drug combinations. The threshold
forms the basis of the voting scheme that underlies the predictive
system. Additionally, the scripts implements a repeated cross-fold
validation framework to understand the variability of the threshold and
accuracy of the model on changing data. The final predictive system is
trained on the complete data. The output is an R list that contains the
logistic regression models for each features that was used for threshold
identification, the identified thresholds, and the results from the
repeated cross-fold validation framework.</td>
</tr>
</tbody>
</table>

### 4.3. Preparing the validation datasets

<table>
<colgroup>
<col style="width: 5%" />
<col style="width: 94%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><p><a
href="../Scripts/Predictive_model/Validation_data_1.R">Scripts/Predictive_model/Validation_data_1.R</a></p>
<p><a
href="../Scripts/Predictive_model/Validation_data_2.R">Scripts/Predictive_model/Validation_data_2.R</a></p>
<p><a
href="../Scripts/Predictive_model/Validation_data_2.R">Scripts/Predictive_model/Validation_data_3.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td>Mandatory: <code>--disease</code></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>drugCombs_data_[disease].rds</p></li>
<li><p>DrugBank_DDI_processed.rds</p></li>
<li><p>drugCombs_cat_effVadv_[disease].rds</p></li>
<li><p>STRING_PPI_Net.rds</p></li>
<li><p>DrugBank_Drug_Target_associations.rds</p></li>
<li><p>cancerdrugsdb.txt (Cancer Drug Database)</p></li>
<li><p>parsed_DrugBank_data.rds</p></li>
<li><p>09.04.2024.zip (C-DCDB)</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>drugCombs_validation1_[disease].rds</p></li>
<li><p>drugCombs_validation2_[disease].rds</p></li>
<li><p>drugCombs_validation3_[disease].rds</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The scripts generate labelled drug combination data that will be
used to validate the predictive systems. The output is a data frame with
each row representing a drug combination. The columns include the
DrugBank drug IDs for the combination, drug combination class and the
drug targets.</td>
</tr>
</tbody>
</table>

### 4.4. Feature generation for validation datasets

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 97%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><p><a
href="../Scripts/Predictive_model/Validation_featureGeneration_1.R">Scripts/Predictive_model/Validation_featureGeneration_1.R</a></p>
<p><a
href="../Scripts/Predictive_model/Validation_featureGeneration_2.R">Scripts/Predictive_model/Validation_featureGeneration_2.R</a></p>
<p><a
href="../Scripts/Predictive_model/Validation_featureGeneration_3.R">Scripts/Predictive_model/Validation_featureGeneration_3.R</a></p></td>
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
<li><p>drugCombs_validation[1/2/3]_[disease].rds</p></li>
<li><p>Disease2Gene_[disease]_lib.rds</p></li>
<li><p>curatedAdr2Gene_lib.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>fgseaNES_combinedEfficacySafety_[disease]_[drug_target_type].rds
(in Validation_data_[1/2/3] directory)</p></li>
<li><p>plot_validation[1/2/3]_AD_combinedEfficacySafety_[disease]_[drug_target_type].tiff</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The scripts takes in the drug combinations identified in section 4.3
and implements the RWR-FGSEA (as described in part 3) framework to
return a data frame of the NES for the drug combinations in column
against the gene sets in rows. However, since our predictive system is
based on the combined efficacy-safety library, the NES for the
validation data is generated only for the same. Additionally, the
scripts plot the applicibility domain of the validation drug
combinations w.r.t the training set drug combinations by considering one
efficacy and one safety feature for the scatter plot.</td>
</tr>
</tbody>
</table>

### 4.5. Prediction on validation data sets

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Rscript</strong></td>
<td><p><a
href="../Scripts/Predictive_model/Validation_prediction_1.R">Scripts/Predictive_model/Validation_prediction_1.R</a></p>
<p><a
href="../Scripts/Predictive_model/Validation_prediction_2.R">Scripts/Predictive_model/Validation_prediction_2.R</a></p>
<p><a
href="../Scripts/Predictive_model/Validation_prediction_3.R">Scripts/Predictive_model/Validation_prediction_3.R</a></p></td>
</tr>
<tr class="even">
<td><strong>Arguments</strong></td>
<td><p>Mandatory: <code>--disease</code></p>
<p>Optional: <code>--drug_target_type</code></p></td>
</tr>
<tr class="odd">
<td><strong>Input</strong></td>
<td><ul>
<li><p>model_NES_combinedEfficacySafety_[disease]_[drug_target_type].rds</p></li>
<li><p>drugCombs_validation[1/2/3]_[disease].rds</p></li>
<li><p>fgseaNES_combinedEfficacySafety_[disease]_[drug_target_type].rds</p></li>
<li><p>parsed_DrugBank_data.rds</p></li>
</ul></td>
</tr>
<tr class="even">
<td><strong>Output</strong></td>
<td><ul>
<li><p>predictions_NES_combinedEfficacySafety_[disease]_[drug_target_type].csv</p></li>
<li><p>predictionMetrics_NES_combinedEfficacySafety_[disease]_[drug_target_type].csv</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><strong>Summary</strong></td>
<td>The script uses the predictive system developed in section 4.2 and
the efficacy and safety estimates calculated in section 4.4 to generate
predictions for the validation data set drug combinations. The primary
output is a data frame with the final score and the predicted categories
along with other drug information. The script also calculates the
accuracy of the prediction and writes it as a separate file.</td>
</tr>
</tbody>
</table>
