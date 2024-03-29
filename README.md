# Analyze OAI database with unsupervised learning algorithm
## Weihua Guo, Ph.D.
## 2022.9.18

## Script descriptions
### Step 0: Learn that each file includes what variables
sum_data.py
"Outputs: file_sum_excel.xlsx, file_sum_dictionary.pkl, col_sum_excel.xlsx, col_sum_dictionary.pkl"

### Step 1: Clean the data
#### Step 1.1: Parse the input data, put the data into a matrix
concatenate_data.py
"Outputs: input_direct_merge_dataframe_12112020_data_clean.csv, input_direct_merge_dataframe_12112020_data_clean.pkl"

#### Step 1.2: Data cleaning, Check the number of missing items in each variables and each subject. Remove the variable and subjects with too many missing items
data_clean_input.py
"Outputs: input_direct_merge_dataframe_12112020_data_clean_sbj_clean_v25.csv, input_direct_merge_dataframe_12112020_data_clean_sbj_clean_v25.pkl"

#### Step 1.3: Self correlation analysis
Calculate correlation coefficients: self_corr_input.py
"Outputs: input_direct_merge_dataframe_12112020_data_clean_correlation_v25.csv, input_direct_merge_dataframe_12112020_data_clean_pvalue_v25.csv, input_direct_merge_dataframe_12112020_data_clean_na_num_v25.csv"
Visualize correlation matrices: self_corr_vis.R
"Outputs: input_direct_merge_dataframe_12112020_data_clean_correlation_vis_v25.png (Fig. Sx), input_direct_merge_dataframe_12112020_data_clean_correlation_order_v25.csv"

#### Step 1.4: Parse the outcome data
concatenate_data_outcome.py
"Outputs: outcome_all_dataframe_12112020_data_clean.csv, outcome_all_dataframe_12112020_data_clean.pkl"

#### Step 1.5: Parse the total knee replacement data
knee_replacement_extraction.py
"Outputs: total_lastfollowup_merge_patient_basic_outcome_information.csv"

#### Step 1.6: Align the outcome data into real time
Extract the real dates: oai_time_vdate.py
"Outputs: real_dates.csv"
Convert the real dates into time period and integrate with outcome data: convert_to_survival_data.py
"Outputs: xxx_used_clean_dataframe.csv, xxx_used_clean_dataframe_with_real_time.csv, xxx_v00_real_time_survival_dataframe.csv"

#### Step 1.7: Categorical variable one-hot spot encoding
Encode the categorical variables with one-hot spot method: categorical_onehot_encoding.py
"Outputs: onehotspot_merge_dataframe_12112020_data_clean_v25_without_na_cols.pkl"

### Step 2: Unsupervised clustering
#### Step 2.1: Imputation, clustering, UMAP
impute_cluster_umap.py
"Outputs: clean_v25_cluster4_kmean_pca_umap_res.xlsx"

#### Step 2.2: Cluster annotation and visualization
cluster_annotate_vis.R
"Outputs: Fig 2, Fig S3, Data S4 (cluster_annotation_vis_outputs), Demographic table"
cluster_outcome_visualization.R
"Outputs: Fig S3 (Outcome variables including TKR along all the visit times crossing clusters)"

#### Step 2.3: Survival analysis of outcomes crossing clusters
cluster_outcome_analysis.R
"Outputs: Fig 3, Fig S4, Fig 4, Fig S5"
knee_replacement_cross_cluster.R
"Outputs: Fig 4A, Fig S5"

### Step 3: Supervised learning for clustering results
predict_cluster_classic.py
"Outputs: All the accuracy scores for predicting the clusters"
predict_cluster_result_vis.R
"Outputs: Fig 4"

### Step 4: Supervised learning for WOMAC total score
#### Step 4.1: Correlation analysis of WOMAC total scores to all input variables (predictors)
womts_correlation_calculation.py
"Outputs: Data S11"
womts_correlation_visualization.R
"Outputs: Fig 5A&B"
#### Step 4.2: Selection of optimal machine learning model for WOMTS prediction
womts_super_learning_prediction.py
"Outputs: Fig 6B&C"
#### Step 4.3: Implementation of WOMTS prediction with CNN model and random selection of input variables
womts_cnn_prediction.py/womts_cnn_prediction_batch.sh
"Outputs: prediction&measurements for each test, accuracy for each test, selected variables"
womts_cnn_prediction_result_merge.py
"Outputs: Data S12"
womts_cnn_prediction_result_vis.R
"Outputs: Fig 6E"
