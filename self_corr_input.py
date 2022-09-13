#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 22:23:23 2020

@author: weihua
self correlation analysis for input variables
"""
def cramers_v(x, y):
    confusion_matrix = pd.crosstab(x,y)
    chi2, p, dof, expctd = chi2_contingency(confusion_matrix)
    n = confusion_matrix.sum().sum()
    phi2 = chi2/n
    r,k = confusion_matrix.shape
    phi2corr = max(0, phi2-((k-1)*(r-1))/(n-1))
    rcorr = r-((r-1)**2)/(n-1)
    kcorr = k-((k-1)**2)/(n-1)
    c = np.sqrt(phi2corr/min((kcorr-1),(rcorr-1)))
    return c,p;

import os
os.system('clear')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import datetime as dt
from statsmodels.formula.api import ols
from scipy.stats import chi2_contingency, pearsonr

mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
exp_id = "12112020_data_clean"
clean_co = "v25"
data_df=pd.read_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+"_sbj_clean_"+clean_co+".pkl")

cor_df=pd.DataFrame(columns=data_df.columns, index=data_df.columns)
p_df=pd.DataFrame(columns=data_df.columns, index=data_df.columns)
na_df=pd.DataFrame(columns=data_df.columns, index=data_df.columns)
tmp_c=2.
tmp_p=-1.
tmp_na=6000
st=dt.datetime.now()
for ic1, icol1 in data_df.iteritems():
    print(ic1)
    for ic2, icol2 in data_df.iteritems():
#        print('\t'+ic2)
        if ic1==ic2:
            tmp_c=1.
            tmp_p=0.
            if pd.api.types.is_string_dtype(icol2):
                tmp_na=icol1.shape[0]-icol1.str.contains("Missing").sum()
            if pd.api.types.is_numeric_dtype(icol1):
                tmp_na=icol1.shape[0]-sum(icol2.isna())
        else:
            if pd.api.types.is_numeric_dtype(icol1)&pd.api.types.is_numeric_dtype(icol2):
                tmp_df=data_df[[ic1,ic2]]
                tmp_df=tmp_df.dropna()
                tmp_na=tmp_df.shape[0]
                if tmp_na>=3:
                    tmp_c, tmp_p=pearsonr(tmp_df[ic1], tmp_df[ic2])
            elif pd.api.types.is_string_dtype(icol1)&pd.api.types.is_string_dtype(icol2):
                sbj1=icol1.index[~icol1.str.contains("Missing")].tolist()
                sbj2=icol2.index[~icol2.str.contains("Missing")].tolist()
                ol_sbj=list(set(sbj1) & set(sbj2))
                tmp_na=len(ol_sbj)
                if tmp_na>0:
                    tmp_c, tmp_p=cramers_v(icol1.loc[ol_sbj], icol2.loc[ol_sbj])             
            else:
                if pd.api.types.is_string_dtype(icol1)&pd.api.types.is_numeric_dtype(icol2):
                    cate_col=ic1
                    num_col=ic2
                if pd.api.types.is_string_dtype(icol2)&pd.api.types.is_numeric_dtype(icol1):
                    cate_col=ic2
                    num_col=ic1
                tmp_df=data_df[[ic1,ic2]]
                sbj_cate=tmp_df.index[~tmp_df[cate_col].str.contains("Missing")].tolist()
                sbj_num=tmp_df.index[~tmp_df[num_col].isna()]
                ol_sbj=list(set(sbj_cate)&set(sbj_num))
                tmp_na=len(ol_sbj)
                if tmp_na>3:
                    tmp_df=tmp_df.loc[ol_sbj]
                    anova_form=num_col+' ~ C('+cate_col+')'
                    tmp_res=ols(anova_form, data=tmp_df).fit()
                    tmp_c=np.sqrt(tmp_res.rsquared)
                    tmp_p=tmp_res.f_pvalue
        cor_df.loc[ic1,ic2]=tmp_c
        p_df.loc[ic1,ic2]=tmp_p
        na_df.loc[ic1, ic2]=tmp_na
    print("\tTime cost"+str(dt.datetime.now()-st))
    print('\tCurrent time '+str(dt.datetime.now()))
cor_df.to_csv(resDir+"/input_direct_merge_dataframe_"+exp_id+"_correlation_"+clean_co+".csv")
p_df.to_csv(resDir+"/input_direct_merge_dataframe_"+exp_id+"_pvalue_"+clean_co+".csv")
na_df.to_csv(resDir+"/input_direct_merge_dataframe_"+exp_id+"_na_num_"+clean_co+".csv")
cor_df.to_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+"_correlation_"+clean_co+".pkl")
p_df.to_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+"_pvalue_"+clean_co+".pkl")
na_df.to_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+"_na_num_"+clean_co+".pkl")