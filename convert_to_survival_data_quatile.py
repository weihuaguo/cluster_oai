#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 17:09:11 2024
enerate survival data based on selected outcome variables
Using holistic quantile approach
@author: weihua
"""

import os
os.system('clear')
import pandas as pd
import pickle as pkl
import glob
import math
import numpy as np
from dateutil.relativedelta import relativedelta

# mainDir = "G:/My Drive/Weihua/OAI_Public_Data"
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
exp_id = "12112020_data_clean"

surv_folder="survival_ready_220927"
surv_folder="survival_ready_qtl_240422"
surv_dir=os.path.join(resDir, surv_folder)
if not os.path.exists(surv_dir):
    os.makedirs(surv_dir)
    
outcome_patterns = ['XRKLL', 'XRKLR','WOMADLL', 'WOMADLR', 'WOMKPL', 'WOMKPR', 'MCMJSWL', 'MCMJSWR', 'WOMSTFL', 'WOMSTFR', 'WOMTSL', 'WOMTSR']

all_outcome_df = pd.read_pickle(resDir+"/outcome_all_dataframe_"+exp_id+".pkl")
real_dates = pd.read_csv(resDir+"/real_dates.csv", index_col=0)

real_dates = real_dates.loc[:,real_dates.columns.str.contains("DATE")].apply(pd.to_datetime, errors='coerce')
for iop in outcome_patterns:
    print(iop)
    i0 = 0
    ic = 0
    iy = 0
    tmp_df = all_outcome_df.loc[:,all_outcome_df.columns.str.contains(iop)]
    tmp_df = tmp_df.apply(pd.to_numeric, errors='coerce')
    tmp_df = tmp_df.dropna(axis=0, how='all')
    tmp_df = tmp_df[sorted(tmp_df.columns.values.tolist())]
    tmp_max = max(tmp_df.max())
    tmp_ydf = pd.DataFrame(columns = ['ID', 'vid', 'tstart', 'tstop', 'value'], index = range(tmp_df.shape[0]*tmp_df.shape[1]))
    for idx in tmp_df.index.values.tolist():
        for col in sorted(tmp_df.columns.values.tolist()):
            tmp_ydf.iloc[iy,0] = idx
            tmp_ydf.iloc[iy,1] = col
            tmp_ydf.iloc[iy,2] = real_dates.loc[idx, 'V00EVDATE']
            if "V00" in col:
                tmp_ydf.iloc[iy,3] = real_dates.loc[idx, col[0:3]+'EVDATE']
            else:
                tmp_ydf.iloc[iy,3] = real_dates.loc[idx, col[0:3]+'FVDATE']
            tmp_ydf.iloc[iy,4] = tmp_df.loc[idx, col]
            iy+=1
            
    tmp_df.to_csv(surv_dir+"/"+iop+"_used_clean_dataframe.csv")
    tmp_ydf.dropna(axis = 0, how = "any", subset = ['value','tstart', 'tstop'], inplace=True)
    tmp_ydf['tstart'] = pd.to_datetime(tmp_ydf['tstart'], errors='coerce')
    tmp_ydf['tstop'] = pd.to_datetime(tmp_ydf['tstop'], errors='coerce')
    tmp_ydf['td'] = (tmp_ydf['tstop'] - tmp_ydf['tstart']).dt.days
    tmp_ydf['tm'] = tmp_ydf['tstop'].dt.to_period('M').astype(int) - tmp_ydf['tstart'].dt.to_period('M').astype(int)
    tmp_ydf['ty'] = tmp_ydf['tstop'].dt.to_period('Y').astype(int) - tmp_ydf['tstart'].dt.to_period('Y').astype(int)
    tmp_ydf.to_csv(surv_dir+"/"+iop+"_used_clean_dataframe_with_real_time.csv")

for iop in outcome_patterns:
    print(iop)
    i0 = 0
    ic = 0
    tmp_df = all_outcome_df.loc[:,all_outcome_df.columns.str.contains(iop)]
    tmp_df = tmp_df.apply(pd.to_numeric, errors='coerce')
    tmp_df = tmp_df.dropna(axis=0, how='all')
    tmp_max = max(tmp_df.max())
    cut_off = tmp_max/4
    tmp_v00_surv = pd.DataFrame(columns=['ID', 'vtime', 'event', 'tstart', 'tstop', 'diff', 'vstart', 'vstop', 'cutoff'], index = tmp_df.index)
    tmp_v00_surv['cutoff'] = cut_off
    for idx in tmp_df.index.values.tolist():
        f0 = True
        for col in sorted(tmp_df.columns.values.tolist()):
            if not math.isnan(tmp_df.loc[idx, col]):
                max_col = col
            if "V00" in col:
                v0 = tmp_df.loc[idx, col]
            else:
                if "KL" in iop:
                    d0 = tmp_df.loc[idx, col] - v0
                    if d0 >= cut_off and f0:
                        tmp_v00_surv.loc[idx, 'ID'] = idx
                        tmp_v00_surv.loc[idx, 'vtime'] = col
                        tmp_v00_surv.loc[idx, 'event'] = 1
                        tmp_v00_surv.loc[idx, 'tstart'] = real_dates.loc[idx, 'V00EVDATE']
                        tmp_v00_surv.loc[idx, 'tstop'] = real_dates.loc[idx, col[0:3]+'FVDATE']
                        tmp_v00_surv.loc[idx, 'diff'] = d0
                        tmp_v00_surv.loc[idx, 'vstart'] = v0
                        tmp_v00_surv.loc[idx, 'vstop'] = tmp_df.loc[idx, col]
                        f0 = False
                elif "JSW" in iop:
                    d0 = tmp_df.loc[idx, col] - v0 # TODO: how to handle pre = 0
                    if d0 <= -cut_off and f0: # NOTE: Already quantile but per patient
                        tmp_v00_surv.loc[idx, 'ID'] = idx
                        tmp_v00_surv.loc[idx, 'vtime'] = col
                        tmp_v00_surv.loc[idx, 'event'] = 1
                        tmp_v00_surv.loc[idx, 'tstart'] = real_dates.loc[idx, 'V00EVDATE']
                        tmp_v00_surv.loc[idx, 'tstop'] = real_dates.loc[idx, col[0:3]+'FVDATE']
                        tmp_v00_surv.loc[idx, 'diff'] = d0
                        tmp_v00_surv.loc[idx, 'vstart'] = v0
                        tmp_v00_surv.loc[idx, 'vstop'] = tmp_df.loc[idx, col]
                        f0 = False
                        i0 += 1
                else:
                    d0 = tmp_df.loc[idx, col] - v0 # TODO: how to handle pre = 0
                    if d0 >= cut_off and f0:
                        tmp_v00_surv.loc[idx, 'ID'] = idx
                        tmp_v00_surv.loc[idx, 'vtime'] = col
                        tmp_v00_surv.loc[idx, 'event'] = 1
                        tmp_v00_surv.loc[idx, 'tstart'] = real_dates.loc[idx, 'V00EVDATE']
                        tmp_v00_surv.loc[idx, 'tstop'] = real_dates.loc[idx, col[0:3]+'FVDATE']
                        tmp_v00_surv.loc[idx, 'diff'] = d0
                        tmp_v00_surv.loc[idx, 'vstart'] = v0
                        tmp_v00_surv.loc[idx, 'vstop'] = tmp_df.loc[idx, col]
                        f0 = False
                        i0 += 1
        if tmp_v00_surv.loc[idx, 'event'] is np.nan:
            if "V00" not in max_col:
                tmp_v00_surv.loc[idx, 'ID'] = idx
                tmp_v00_surv.loc[idx, 'vtime'] = max_col
                tmp_v00_surv.loc[idx, 'event'] = 0
                tmp_v00_surv.loc[idx, 'tstart'] = real_dates.loc[idx, 'V00EVDATE']
                tmp_v00_surv.loc[idx, 'tstop'] = real_dates.loc[idx, max_col[0:3]+'FVDATE']
                tmp_v00_surv.loc[idx, 'vstart'] = v0
                tmp_v00_surv.loc[idx, 'vstop'] = tmp_df.loc[idx, max_col]
        else:
            tmp_v00_surv.loc[idx, 'vstop'] = tmp_df.loc[idx, max_col]
    tmp_v00_surv.dropna(axis = 0, how = "all", inplace=True)
    tmp_v00_surv['tstart'] = pd.to_datetime(tmp_v00_surv['tstart'], errors='coerce')
    tmp_v00_surv['tstop'] = pd.to_datetime(tmp_v00_surv['tstop'], errors='coerce')
    
    tmp_v00_surv = tmp_v00_surv.merge(real_dates['V00EVDATE'], left_index=True, right_index=True, how = "left")
    tmp_v00_surv['dstart'] = (tmp_v00_surv['tstart'] - tmp_v00_surv['V00EVDATE']).dt.days
    tmp_v00_surv['dstop'] = (tmp_v00_surv['tstop'] - tmp_v00_surv['V00EVDATE']).dt.days
    tmp_v00_surv['mstart'] = tmp_v00_surv['tstart'].dt.to_period('M').astype(int) - tmp_v00_surv['V00EVDATE'].dt.to_period('M').astype(int)
    tmp_v00_surv['mstop'] = tmp_v00_surv['tstop'].dt.to_period('M').astype(int) - tmp_v00_surv['V00EVDATE'].dt.to_period('M').astype(int)
#    tmp_v00_surv['dstart'] = (tmp_v00_surv['tstart'] - tmp_v00_surv['V00EVDATE'])/np.timedelta64(1, 'D')
#    tmp_v00_surv['dstart'] = tmp_v00_surv['dstart'].astype(int)
#    tmp_v00_surv['dstop'] = (tmp_v00_surv['tstop'] - tmp_v00_surv['V00EVDATE'])/np.timedelta64(1, 'D')
#    tmp_v00_surv['dstop'] = tmp_v00_surv['dstop'].astype(int)
#    tmp_v00_surv['mstart'] = (tmp_v00_surv['tstart'] - tmp_v00_surv['V00EVDATE'])/np.timedelta64(1, 'M')
#    tmp_v00_surv['mstart'] = tmp_v00_surv['mstart'].astype(int)
#    tmp_v00_surv['mstop'] = (tmp_v00_surv['tstop'] - tmp_v00_surv['V00EVDATE'])/np.timedelta64(1, 'M')
#    tmp_v00_surv['mstop'] = tmp_v00_surv['mstop'].astype(int)
    
    tmp_v00_surv['td'] = (tmp_v00_surv['tstop'] - tmp_v00_surv['tstart']).dt.days
#    tmp_v00_surv['td'] = (tmp_v00_surv['tstop'] - tmp_v00_surv['tstart'])/np.timedelta64(1, 'D')
#    tmp_v00_surv['td'] = tmp_v00_surv['td'].astype(int)
    tmp_v00_surv['tm'] = tmp_v00_surv['tstop'].dt.to_period('M').astype(int) - tmp_v00_surv['tstart'].dt.to_period('M').astype(int)
#    tmp_v00_surv['tm'] = (tmp_v00_surv['tstop'] - tmp_v00_surv['tstart'])/np.timedelta64(1, 'M')
#    tmp_v00_surv['tm'] = tmp_v00_surv['tm'].astype(int)
    tmp_v00_surv['ty'] = tmp_v00_surv['tstop'].dt.to_period('Y').astype(int) - tmp_v00_surv['tstart'].dt.to_period('Y').astype(int)
    #tmp_v00_surv['ty'] = (tmp_v00_surv['tstop'] - tmp_v00_surv['tstart'])/np.timedelta64(1, 'Y')
    #tmp_v00_surv['ty'] = tmp_v00_surv['ty'].astype(float)
    tmp_v00_surv.to_csv(surv_dir+"/"+iop+"_v00_real_time_survival_dataframe.csv") #SD6
