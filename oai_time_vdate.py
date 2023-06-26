import os
os.system('clear')
import pandas as pd
import pickle as pkl
import glob

mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"

enroll_txt = "Enrollees.txt"
outcome_txt = "Outcomes99.txt"


all_pat_info = pd.read_csv(dataDir+"/"+enroll_txt, sep='|', index_col=0)
print(all_pat_info)
print()
enr_cols = all_pat_info.columns.values.tolist()
print(all_pat_info.loc[:,all_pat_info.columns.str.contains("DATE")])
print(all_pat_info.loc[:,all_pat_info.columns.str.contains("COHORT")])

for itf in all_txt_files:
    print(itf)
    tmp_df = pd.read_csv(itf, sep = "|", index_col=0)
    if any(tmp_df.columns.str.contains("DATE")):
        all_pat_info=pd.merge(all_pat_info, tmp_df.loc[:,tmp_df.columns.str.contains("VDATE")], left_index=True, right_index=True, how='left')
all_pat_info.to_csv(resDir+"/real_dates.csv")
