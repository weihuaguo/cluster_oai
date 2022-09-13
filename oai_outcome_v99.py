import os
os.system('clear')
import pandas as pd
import pickle as pkl
import numpy as np
import glob

mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"

enroll_txt = "Enrollees.txt"
real_date_txt = "real_dates.csv"
outcome_txt = "Outcomes99.txt"
kr_criteria = "xray"

all_pat_info = pd.read_csv(resDir+"/"+real_date_txt, sep=',', index_col=0)
outcome_info = pd.read_csv(dataDir+"/"+outcome_txt, sep='|', index_col=0)

print(all_pat_info)
print(outcome_info)
enr_cols = all_pat_info.columns.values.tolist()
print(all_pat_info.loc[:,all_pat_info.columns.str.contains("DATE")])
print(all_pat_info.loc[:,all_pat_info.columns.str.contains("COHORT")])

merge_df = pd.merge(all_pat_info, outcome_info, left_index=True, right_index=True)

day_df = merge_df.loc[:,merge_df.columns.str.contains("DATE")]
all_dates = merge_df.loc[:,merge_df.columns.str.contains("DATE")].columns.values.tolist()
day_df = day_df.apply(pd.to_datetime, errors = "coerce")
for idt in all_dates:
    day_df[idt.replace("DATE", "DAY")] = (day_df[idt] - day_df['V00EVDATE'])/np.timedelta64(1, 'D')
    
total_df = pd.merge(merge_df, day_df, left_index=True, right_index=True)
total_df['left_time'] = total_df['V99ELKDAYS']
total_df['left_kr'] = total_df['V99ELKTLPR']
total_df['lkr_event'] = 0
total_df['right_time'] = total_df['V99ERKDAYS']
total_df['right_kr'] = total_df['V99ERKTLPR']
total_df['rkr_event'] = 0

# NOTE: Use the lastest follow-up visit date for cencored time point!!!
for iid in total_df.index.values.tolist():
    print('\t'+str(iid))
    latest_time_day = sorted(total_df.loc[iid, total_df.columns.str.contains("VDAY")].dropna())[-1]
    latest_time = total_df.columns[(total_df.columns.str.contains("VDAY") & (total_df.loc[iid,:] == latest_time_day).tolist())].values.tolist()[0]
    if "Total" in total_df.loc[iid, 'left_kr']:
        total_df.loc[iid, 'lkr_event'] = 1
    else:
        if ".:" in total_df.loc[iid, 'left_time']:
            if total_df.loc[iid, "V99EDDDAY"] > latest_time_day:
                total_df.loc[iid, 'left_time'] = day_df.loc[iid, "V99EDDDAY"]
            else:
                total_df.loc[iid, 'left_time'] = day_df.loc[iid, latest_time]
        else:
            if int(total_df.loc[iid, "left_time"]) < latest_time_day:
                total_df.loc[iid, 'left_time'] = latest_time_day
            
    if "Total" in total_df.loc[iid, 'right_kr']:
        total_df.loc[iid, 'rkr_event'] = 1
    else:
        if ".:" in total_df.loc[iid, 'right_time']:
            if total_df.loc[iid, "V99EDDDAY"] > latest_time_day:
                total_df.loc[iid, 'right_time'] = day_df.loc[iid, "V99EDDDAY"]
            else:
                total_df.loc[iid, 'right_time'] = day_df.loc[iid, latest_time]
        else:
            if int(total_df.loc[iid, "right_time"]) < latest_time_day:
                total_df.loc[iid, 'right_time'] = latest_time_day

total_df.to_csv(resDir+"/total_only_merge_patient_basic_outcome_information.csv")


partial_df = pd.merge(merge_df, day_df, left_index=True, right_index=True)
partial_df['left_time'] = partial_df['V99ELKDAYS']
partial_df['left_kr'] = partial_df['V99ELKTLPR']
partial_df['lkr_event'] = 0
partial_df['right_time'] = partial_df['V99ERKDAYS']
partial_df['right_kr'] = partial_df['V99ERKTLPR']
partial_df['rkr_event'] = 0
# NOTE: Use the lastest follow-up visit date for cencored time point!!!
for iid in partial_df.index.values.tolist():
    print('\t'+str(iid))
    latest_time_day = sorted(partial_df.loc[iid, partial_df.columns.str.contains("VDAY")].dropna())[-1]
    latest_time = partial_df.columns[(partial_df.columns.str.contains("VDAY") & (partial_df.loc[iid,:] == latest_time_day).tolist())].values.tolist()[0]
    if "Total" in partial_df.loc[iid, 'left_kr'] or "Partial" in partial_df.loc[iid, 'left_kr']:
        partial_df.loc[iid, 'lkr_event'] = 1
    else:
        if ".:" in partial_df.loc[iid, 'left_time']:
            if partial_df.loc[iid, "V99EDDDAY"] > latest_time_day:
                partial_df.loc[iid, 'left_time'] = day_df.loc[iid, "V99EDDDAY"]
            else:
                partial_df.loc[iid, 'left_time'] = day_df.loc[iid, latest_time]
            
    if "Total" in partial_df.loc[iid, 'right_kr'] or "Partial" in partial_df.loc[iid, 'right_kr']:
        partial_df.loc[iid, 'rkr_event'] = 1
    else:
        if ".:" in partial_df.loc[iid, 'right_time']:
            if partial_df.loc[iid, "V99EDDDAY"] > latest_time_day:
                partial_df.loc[iid, 'right_time'] = day_df.loc[iid, "V99EDDDAY"]
            else:
                partial_df.loc[iid, 'right_time'] = day_df.loc[iid, latest_time]

partial_df.to_csv(resDir+"/total_partial_merge_patient_basic_outcome_information.csv")

xray_df = pd.merge(merge_df, day_df, left_index=True, right_index=True)
xray_df['left_time'] = xray_df['V99ELKDAYS']
xray_df['left_kr'] = xray_df['V99ELKRPSN']
xray_df['lkr_event'] = 0
xray_df['right_time'] = xray_df['V99ERKDAYS']
xray_df['right_kr'] = xray_df['V99ERKRPSN']
xray_df['rkr_event'] = 0

# NOTE: Use the LATEST X-RAY date for cencored time point!!!
for iid in xray_df.index.values.tolist():
    print('\t'+str(iid))
    xray_mask = xray_df.loc[iid, xray_df.columns.str.contains("IMAGES")].str.contains("ray")
    if any(xray_mask):
        latest_xray_time = sorted(xray_df.columns[xray_df.columns.str.contains("IMAGES")][xray_mask].values.tolist())[-1]
        if "V00" in latest_xray_time:
            latest_xray_time_day = xray_df.loc[iid, latest_xray_time[0:3]+"EVDAY"]
        else:
            latest_xray_time_day = xray_df.loc[iid, latest_xray_time[0:3]+"FVDAY"]
        if "1:" in xray_df.loc[iid, 'left_kr']:
            xray_df.loc[iid, 'lkr_event'] = 1
        else:
            if ".:" in xray_df.loc[iid, 'left_time']:
                if "V00" in latest_xray_time:
                    xray_df.loc[iid, 'left_time'] = day_df.loc[iid, latest_xray_time[0:3]+"EVDAY"]
                else:
                    xray_df.loc[iid, 'left_time'] = day_df.loc[iid, latest_xray_time[0:3]+"FVDAY"]
            else:
                if int(xray_df.loc[iid, "left_time"]) < latest_xray_time_day:
                    xray_df.loc[iid, 'left_time'] = latest_xray_time_day
        if "1:" in xray_df.loc[iid, 'right_kr']:
            xray_df.loc[iid, 'rkr_event'] = 1
        else:
            if ".:" in xray_df.loc[iid, 'right_time']:
                if "V00" in latest_xray_time:
                    xray_df.loc[iid, 'right_time'] = day_df.loc[iid, latest_xray_time[0:3]+"EVDAY"]
                else:
                    xray_df.loc[iid, 'right_time'] = day_df.loc[iid, latest_xray_time[0:3]+"FVDAY"]
            else:
                if int(xray_df.loc[iid, "right_time"]) < latest_xray_time_day:
                    xray_df.loc[iid, 'right_time'] = latest_xray_time_day
    else:
        xray_df = xray_df.drop(iid)
        
xray_df.to_csv(resDir+"/xray_only_merge_patient_basic_outcome_information.csv")

xray_df = pd.merge(merge_df, day_df, left_index=True, right_index=True)
xray_df['left_time'] = xray_df['V99ELKDAYS']
xray_df['left_kr'] = xray_df['V99ELKRPSN']
xray_df['lkr_event'] = 0
xray_df['right_time'] = xray_df['V99ERKDAYS']
xray_df['right_kr'] = xray_df['V99ERKRPSN']
xray_df['rkr_event'] = 0

# NOTE: Use the LATEST visit date for cencored time point!!!
for iid in xray_df.index.values.tolist():
    print('\t'+str(iid))
    latest_time_day = sorted(xray_df.loc[iid, xray_df.columns.str.contains("VDAY")].dropna())[-1]
    latest_time = xray_df.columns[(xray_df.columns.str.contains("VDAY") & (xray_df.loc[iid,:] == latest_time_day).tolist())].values.tolist()[0]

    if "1:" in xray_df.loc[iid, 'left_kr']:
        xray_df.loc[iid, 'lkr_event'] = 1
    else:
        if ".:" in xray_df.loc[iid, 'left_time']:
            if xray_df.loc[iid, "V99EDDDAY"] > latest_time_day:
                xray_df.loc[iid, 'left_time'] = day_df.loc[iid, "V99EDDDAY"]
            else:
                xray_df.loc[iid, 'left_time'] = day_df.loc[iid, latest_time]
        else:
            if int(xray_df.loc[iid, "left_time"]) < latest_time_day:
                xray_df.loc[iid, 'left_time'] = latest_time_day
    if "1:" in xray_df.loc[iid, 'right_kr']:
        xray_df.loc[iid, 'rkr_event'] = 1
    else:
        if ".:" in xray_df.loc[iid, 'right_time']:
            if xray_df.loc[iid, "V99EDDDAY"] > latest_time_day:
                xray_df.loc[iid, 'right_time'] = day_df.loc[iid, "V99EDDDAY"]
            else:
                xray_df.loc[iid, 'right_time'] = day_df.loc[iid, latest_time]
        else:
            if int(xray_df.loc[iid, "right_time"]) < latest_time_day:
                xray_df.loc[iid, 'right_time'] = latest_time_day

        
xray_df.to_csv(resDir+"/xray_visit_merge_patient_basic_outcome_information.csv")
