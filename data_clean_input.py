#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 19:38:27 2020

@author: weihua
further clean up the input matrix
    Count #NA in variable and subject
"""

import os
os.system('clear')
import pandas as pd
import matplotlib.pyplot as plt

mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
exp_id = "12112020_data_clean"
clean_cutoff = 0.25
clean_co = "v25"
ann_xlsx = "AllClinical00_V5_column_annotation_WG_clean.xlsx"
enroll_txt = "Enrollees.txt"

data_df=pd.read_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+".pkl")
all_col_oi = pd.read_excel(mainDir+'/'+ann_xlsx, sheet_name='Baseline', index_col='Variables')
all_pat_info = pd.read_csv(dataDir+"/"+enroll_txt, sep='|', index_col=0)

# Check the variable side
print("Checking variables...")
vars_cols=['definition', 'na_num','type', 'na_rel_num']
vars_check_df=pd.DataFrame(index=all_col_oi.index, columns=vars_cols)
for iv, ir in all_col_oi.iterrows():
    vars_check_df.loc[iv, 'definition']=all_col_oi.loc[iv, 'Label']
    equal_check=sum(data_df.columns==iv.upper())
    if equal_check==1:
        if pd.api.types.is_string_dtype(data_df[iv.upper()]):
            tmp_na=data_df[iv.upper()].str.contains("Missing").sum()
            tmp_clean_col=data_df[iv.upper()].loc[~data_df[iv.upper()].str.contains("Missing")]
            colon_num=tmp_clean_col.str.contains(':').sum()
            if colon_num==0 and tmp_na<data_df.shape[0]:
                data_df[iv.upper()]=pd.to_numeric(data_df[iv.upper()], errors='coerce')
                print("WARNING: "+iv.upper()+'!!!')
        if pd.api.types.is_numeric_dtype(data_df[iv.upper()]):
            tmp_na=sum(data_df[iv.upper()].isna())
        vars_check_df.loc[iv, 'na_num']=tmp_na
        vars_check_df.loc[iv, 'type']=data_df[iv.upper()].dtype
    else:
        raise ValueError("No matching column!!!")
vars_check_df['na_rel_num']=vars_check_df['na_num']/all_pat_info.shape[0]
vars_check_df.to_excel(resDir+'/variable_check_dataframe_'+exp_id+'.xlsx')
ax=vars_check_df['na_rel_num'].hist(bins=100)
fig=ax.get_figure()
fig.savefig(resDir+'/variable_check_na_rel_hist_'+exp_id+'.png')
plt.close()

# Check the subject
print("Checking subjects...")
sbj_cols=['na_num', 'na_rel_num']
sbj_check_df=pd.DataFrame(index=data_df.index, columns=sbj_cols)
sub_data_df=data_df[[i.upper() for i in all_col_oi.index.tolist()]]
for isbj, ir in sub_data_df.iterrows():
    tmp_num_na=sum(ir.isna())
    tmp_obj_na=ir.str.contains("Missing").sum()
    sbj_check_df.loc[isbj, 'na_num']=tmp_num_na+tmp_obj_na
sbj_check_df['na_rel_num']=sbj_check_df['na_num']/sub_data_df.shape[1]
sbj_check_df.to_excel(resDir+'/subject_check_dataframe_'+exp_id+'.xlsx')
ax=sbj_check_df['na_rel_num'].hist(bins=100)
fig=ax.get_figure()
fig.savefig(resDir+'/subject_check_na_rel_hist_'+exp_id+'.png')
plt.close()

rm_sbj_id=sbj_check_df.index[sbj_check_df['na_rel_num']>0.5]
rm_var_id=vars_check_df.index[vars_check_df['na_rel_num']>clean_cutoff]
sub_clean_data_df=sub_data_df.drop(rm_sbj_id)
sub_clean_data_df=sub_clean_data_df.drop(rm_var_id,axis=1)
sub_clean_data_df.to_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+"_sbj_clean_"+clean_co+".pkl")
sub_clean_data_df.to_csv(resDir+"/input_direct_merge_dataframe_"+exp_id+"_sbj_clean_"+clean_co+".csv") # SD2
