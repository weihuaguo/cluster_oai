#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 20:57:23 2023

OAI

Predictor -- WOMTS correlation

@author: weihua
"""
   
def cramers_v(x, y):
    confusion_matrix = pd.crosstab(x,y)
    chi2, p, dof, expctd = stats.chi2_contingency(confusion_matrix)
    n = confusion_matrix.sum().sum()
    phi2 = chi2/n
    r,k = confusion_matrix.shape
    phi2corr = max(0, phi2-((k-1)*(r-1))/(n-1))
    rcorr = r-((r-1)**2)/(n-1)
    kcorr = k-((k-1)**2)/(n-1)
    c = np.sqrt(phi2corr/min((kcorr-1),(rcorr-1)))
    return c,p;

import os
import pickle
import pandas as pd
import numpy as np
import scipy.stats as stats
from datetime import datetime
from statsmodels.formula.api import ols

data_dir="/mnt/sda1/OAI_Data/data_summary"

out_dir = data_dir+"/outcome_sum"
exp_id = "12112020_data_clean"
clean_co = "v25"

result_folder="predictor_select_230213"
result_dir=os.path.join(data_dir, result_folder)
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
spf=result_dir+"/"+result_folder+"_"

output_df=pd.read_csv(data_dir+"/outcome_all_dataframe_"+exp_id+"_real_date_year.csv", index_col=1).iloc[:,1:]
input_df=pd.read_pickle(data_dir+'/input_direct_merge_dataframe_'+exp_id+'_sbj_clean_'+clean_co+'.pkl')

outcome_patterns=['WOMTSL', 'WOMTSR']

for iop in outcome_patterns:
    print('\t'+iop)
    tmp_df=output_df.loc[:,output_df.columns.str.contains(iop)]
    
    cor_df=pd.DataFrame(columns=input_df.columns, index=tmp_df.columns)
    p_df=pd.DataFrame(columns=input_df.columns, index=tmp_df.columns)   
    na_df=pd.DataFrame(columns=input_df.columns, index=tmp_df.columns)
    cvr_df=pd.DataFrame(columns=input_df.columns, index=tmp_df.columns) # Record Cramer's V, take KL grade as categorical variable
    cvp_df=pd.DataFrame(columns=input_df.columns, index=tmp_df.columns)
    anr_df=pd.DataFrame(columns=input_df.columns, index=tmp_df.columns) # Record Point-Biserial == ANOVA, take KL grade as continuous variable
    anp_df=pd.DataFrame(columns=input_df.columns, index=tmp_df.columns)
    
    tmp_merge_df=pd.merge(tmp_df, input_df, how='inner', left_index=True, right_index=True, suffixes=('', '_y'))
    
    for iip in input_df.columns.values.tolist():
        print('\t\t'+iip)
        if 'KL' in iop:
            if input_df[iip].dtype != "object":
                for ioc in tmp_df.columns.values.tolist():
                    keep_mask=tmp_merge_df[[ioc, iip]].isnull().sum(axis=1)==0
                    tmp_use_df=tmp_merge_df.loc[keep_mask]
                    if tmp_use_df.shape[0] >=3:
                        tmp_use_df[ioc]=tmp_use_df[ioc].astype('category')
                        r,p=stats.kendalltau(tmp_use_df[iip], tmp_use_df[ioc])
                        cor_df.loc[ioc, iip]=r
                        p_df.loc[ioc, iip]=p
                        cvr_df.loc[ioc, iip]=r
                        cvp_df.loc[ioc, iip]=p
                        anr_df.loc[ioc, iip]=r
                        anp_df.loc[ioc, iip]=p
                    na_df.loc[ioc, iip]=keep_mask.sum()
            else:
                for ioc in tmp_df.columns.values.tolist():
                    na_mask=tmp_merge_df[ioc].isnull()==0
                    miss_mask=~tmp_merge_df[iip].str.contains(".: Missing")
                    keep_mask=na_mask&miss_mask
                    tmp_use_df=tmp_merge_df.loc[keep_mask, [iip,ioc]]
                    if tmp_use_df.shape[0] >=3:
                        anova_form=ioc+' ~ C('+iip+')'
                        tmp_res=ols(anova_form, data=tmp_use_df).fit()
                        anr=np.sqrt(tmp_res.rsquared)
                        anp=tmp_res.f_pvalue
                        anr_df.loc[ioc, iip]=anr
                        anp_df.loc[ioc, iip]=anp
                        
                        cvr,cvp=cramers_v(tmp_use_df[iip], tmp_use_df[ioc])
                        cvr_df.loc[ioc, iip]=cvr
                        cvp_df.loc[ioc, iip]=cvp
                    na_df.loc[ioc, iip]=keep_mask.sum()
                    
            cvr_df.to_csv(spf+iop+"_cvr.csv")
            cvp_df.to_csv(spf+iop+"_cvp.csv")
            anr_df.to_csv(spf+iop+"_anr.csv")
            anp_df.to_csv(spf+iop+"_anp.csv")
            cvr_df.to_csv(spf+iop+"_cor.csv")
            cvp_df.to_csv(spf+iop+"_pval.csv")
        else:
            if input_df[iip].dtype != "object":
                for ioc in tmp_df.columns.values.tolist():
                    tmp_merge_df[ioc]=pd.to_numeric(tmp_merge_df[ioc], errors='coerce')
                    keep_mask=tmp_merge_df[[ioc, iip]].isnull().sum(axis=1)==0
                    tmp_use_df=tmp_merge_df.loc[keep_mask]
                    if tmp_use_df.shape[0] >=3:
                        r,p=stats.spearmanr(tmp_use_df[iip], tmp_use_df[ioc])
                        cor_df.loc[ioc, iip]=r
                        p_df.loc[ioc, iip]=p
                    na_df.loc[ioc, iip]=keep_mask.sum()
            else:
                for ioc in tmp_df.columns.values.tolist():
                    tmp_merge_df[ioc]=pd.to_numeric(tmp_merge_df[ioc], errors='coerce')
                    na_mask=tmp_merge_df[ioc].isnull()==0
                    miss_mask=~tmp_merge_df[iip].str.contains(".: Missing")
                    keep_mask=na_mask&miss_mask
                    tmp_use_df=tmp_merge_df.loc[keep_mask, [iip,ioc]]
                    if tmp_use_df.shape[0] >=3:
                        anova_form=ioc+' ~ C('+iip+')'
                        tmp_res=ols(anova_form, data=tmp_use_df).fit()
                        cor=np.sqrt(tmp_res.rsquared)
                        pval=tmp_res.f_pvalue
                        cor_df.loc[ioc, iip]=cor
                        p_df.loc[ioc, iip]=pval
                    na_df.loc[ioc, iip]=keep_mask.sum()
            cor_df.to_csv(spf+iop+"_cor.csv")
            p_df.to_csv(spf+iop+"_pval.csv")
    na_df.to_csv(spf+iop+"_smp_num.csv")
