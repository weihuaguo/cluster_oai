# Analyze OAI database with unsupervised learning algorithm
## Weihua Guo, Ph.D.
## 2022.9.18

## Script descriptions
### Step 0: Learn that each file includes what variables
sum_data.py -- Done
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
oai_outcome_v99.py
"Outputs: xray_only_hr_between_knee_replacement_and_outcomes.csv, xray_visit_merge_patient_basic_outcome_information.csv"

#### Step 1.6: Align the outcome data into real time
Extract the real dates: oai_time_vdate.py
"Outputs: real_dates.csv"
Convert the real dates into time period and integrate with outcome data: convert_to_survival_data.py
"Outputs: Need to be cleaned!!! Select a way to quantify the survival. Should be V00"
TODO: left here
