#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:16:01 2020
@author: weihua
"""
import os
import pandas as pd
import seaborn as sn
import numpy as np
import gower
import umap
from datetime import datetime as dt
from scipy import stats
import matplotlib.pyplot as plt
st = dt.now()
# mainDir = "G:/My Drive/Weihua/OAI_Public_Data"
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
exp_id = "12112020_data_clean"
clean_co = "v25"
data_df=pd.read_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+"_sbj_clean_"+clean_co+".pkl")

ann_xlsx = "AllClinical00_V5_column_annotation_WG_clean.xlsx"
enroll_txt = "Enrollees.txt"

method = "onehot"

direct_merge_df = pd.read_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+"_sbj_clean_"+clean_co+".pkl")
print('Total time cost '+str(dt.now()-st))
use_merge_df = direct_merge_df

# Separate numeric columns and categorical columns and remove the categorical columns without variance...
droped_cols = []
num_cols = []
str_cols = []
for coln, col in use_merge_df.iteritems():
    if str(col.dtype) == "object":
        unique_items = use_merge_df[coln].value_counts().shape[0]
        if unique_items <= 1:
            droped_cols.append(coln)
        else:
            str_cols.append(coln)
    else:
        num_cols.append(coln)
clean_merge_df = use_merge_df.drop(droped_cols, axis=1)

if method == "onehot":
    print("Use ONEHOT Spot method for category encoding...")
    true_cate_cols = []
    na_item = '.: Not available'
    obj_merge_df = clean_merge_df.select_dtypes(include=['object']).copy()
    obj_merge_df.replace('.: Missing Form/Incomplete Workbook', na_item, inplace=True)
    num_merge_df = clean_merge_df.select_dtypes(exclude=['object']).copy()
    # Check whether some columns are mixture of numerical variable and categorical variable
    for ic in obj_merge_df.columns.values.tolist():
        num_mask = obj_merge_df[ic].str.isnumeric()
        if num_mask.sum() != 0:
            obj_merge_df[ic] =  pd.to_numeric(obj_merge_df[ic], errors='coerce')
        else:
            true_cate_cols.append(ic)
    cate_merge_df=obj_merge_df[true_cate_cols]
    if obj_merge_df.shape[1] > len(true_cate_cols):
        obj_num_df=obj_merge_df.drop(true_cate_cols, axis=1)
        encoded_merge_df = pd.concat([obj_num_df, num_merge_df], axis=1, sort=False)
    else:
        encoded_merge_df = num_merge_df
    cate_encoded_df = pd.get_dummies(cate_merge_df, columns=true_cate_cols, prefix=true_cate_cols) # THIS IS ONE-HOT ENOCODING
    encoded_merge_df = pd.concat([encoded_merge_df, cate_encoded_df], axis=1, sort=False)
    check_obj_df = encoded_merge_df.select_dtypes(include=['object']).copy()
    if check_obj_df.shape[1] != 0:
        raise ValueError("Still has non-numeric columns!!!")
    
    print("Remove the NA categorical columns")
    clean_encoded_merge_df = encoded_merge_df
    for iob in true_cate_cols:
        na_col = iob+'_'+na_item
        if na_col in encoded_merge_df.columns.values.tolist():
            # Replace the corresponding items to nan, then remove the whole na column!
            col_oi = [x for x in encoded_merge_df.columns.values.tolist() if iob in x]
            col_use = [x for x in col_oi if x != na_col]
            tmp_mask = encoded_merge_df[na_col] == 1
            for icu in col_use:
                clean_encoded_merge_df.loc[tmp_mask, icu] = np.nan
            clean_encoded_merge_df = clean_encoded_merge_df.drop(columns=na_col)
    clean_encoded_merge_df.to_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols.pkl')
    encoded_merge_df.to_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_with_na_cols.pkl')
