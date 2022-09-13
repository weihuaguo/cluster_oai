#Analyze OAI database with unsupervised learning algorithm
##Weihua Guo, Ph.D.
##2022.9.18

##Script descriptions
###Step 0: Learn that each file includes what variables
sum_data.py -- Done
"Outputs: file_sum_excel.xlsx, file_sum_dictionary.pkl, col_sum_excel.xlsx, col_sum_dictionary.pkl"

###Step 1: Clean the data
####Step 1.1: Parse the input data, put the data into a matrix
concatenate_data.py
"Outputs: input_direct_merge_dataframe_12112020_data_clean.csv, input_direct_merge_dataframe_12112020_data_clean.pkl"

###Step 1.2: Data cleaning, Check the number of missing items in each variables and each subject. Remove the variable and subjects with too many missing items
data_clean_input.py
"Outputs: input_direct_merge_dataframe_12112020_data_clean_sbj_clean_v25.csv, input_direct_merge_dataframe_12112020_data_clean_sbj_clean_v25.pkl"

TODO: left here
###Step 1.3: Self correlation analysis
Calculate correlation coefficients: self_corr_input.py
Visualize correlation matrices: self_corr_vis.R

