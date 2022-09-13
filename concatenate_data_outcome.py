"""
This script is only run in python 2.7!
Weihua Guo, Ph.D.
07/17/2020
"""
import os
os.system('clear')
import pandas as pd
import pickle as pkl
import glob
import warnings
from pandas.api.types import is_numeric_dtype as inumdt
from datetime import datetime as dt
import matplotlib.pyplot as plt

# mainDir = "G:/My Drive/Weihua/OAI_Public_Data"
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
exp_id = "12112020_data_clean"

enroll_txt = "Enrollees.txt"

total_file_version_merge = False

all_pat_info = pd.read_csv(dataDir+"/"+enroll_txt, sep='|', index_col=0)
print(all_pat_info)
direct_merge_df = all_pat_info

outcome_patterns = ['XRKL', 'WOMADLL', 'WOMADLR', 'WOMKPL', 'WOMKPR', 'MCMJSW', 'WOMSTFL', 'WOMSTFR', 'WOMTSL', 'WOMTSR']

all_txt_files = glob.glob(dataDir+"/*[0-9].txt")
# print(all_txt_files)
print(len(all_txt_files))

# TODO: Check where these two files come from
with open(resDir+'/file_sum_dictionary.pkl', 'rb') as handle:
    file_dict = pkl.load(handle)
with open(resDir+'/col_sum_dictionary.pkl', 'rb') as handle:
    col_dict = pkl.load(handle)

all_files = col_dict['ID']['filename']
all_cols = col_dict.keys()
col_oi = []
for y in outcome_patterns:
    col_oi+=[x for x in all_cols if y in x]

ict = 0
side_cts = 0
st = dt.now()
for icol in col_oi:
    print(icol+'\t'+str(ict))
    cst = dt.now()
    if icol not in direct_merge_df.columns.values.tolist():
        tmp_files = col_dict[icol]['filename']
        if len(tmp_files) == 2:
            print("\tTwo files have this column!")
            tmp_pat_num = col_dict[icol]['patient_num']
            tmp_na_num = col_dict[icol]['na_num']
            # NOTE: assume equal patient num represents the same index
            ## TODO: NEED check!!!
            # Unequal patient num will be iterated added
            if tmp_pat_num[0]==tmp_pat_num[1]:
                print("\t\tEqual patient number")
                tmp_filename = dataDir+'/'+tmp_files[0]
                if tmp_na_num[0] > tmp_na_num[1]:
                    tmp_filename = dataDir+'/'+tmp_files[1]
                else:
                    tmp_filename = dataDir+'/'+tmp_files[0]
                tmp_data = pd.read_csv(tmp_filename, sep='|')
                tmp_data.columns = tmp_data.columns.str.upper()
                if 'SIDE' not in tmp_data.columns.values.tolist():
                    print('\t\tUnique')
                    tmp_data.set_index('ID', inplace=True)
                    direct_merge_df = pd.concat([direct_merge_df, tmp_data[icol]], axis=1)
                else:
                    # NOTE: assume SIDE, READPRJ, ID is always there
                    print('\t\tSIDE_two')
                    side_cols = ['ID', 'SIDE', 'READPRJ', icol]
                    basic_info_cols = ['ID', 'SIDE', 'READPRJ']
                    tmp_data[basic_info_cols] = tmp_data[basic_info_cols].astype('str')
                    tmp_data.loc[tmp_data['SIDE'].str.contains('1'), 'SIDE'] = '1: Right'
                    tmp_data.loc[tmp_data['SIDE'].str.contains('2'), 'SIDE'] = '2: Left'
                    if side_cts == 0:
                        concat_side_df = tmp_data[side_cols]
                    else:
                        concat_side_df = pd.merge(concat_side_df, tmp_data[side_cols], how='outer', on=['ID', 'SIDE', 'READPRJ'])
                    side_cts += 1
            else:
                print("\t\tUnequal patient number")
                tmp_data_1 = pd.read_csv(dataDir+'/'+tmp_files[0], sep='|')
                tmp_data_1.columns = tmp_data_1.columns.str.upper()
                tmp_data_2 = pd.read_csv(dataDir+'/'+tmp_files[1], sep='|')
                tmp_data_2.columns = tmp_data_2.columns.str.upper()
                ol_pid = list(set(tmp_data_1['ID'])&set(tmp_data_2['ID']))
                if 'SIDE' in tmp_data_1.columns.values.tolist():
                    # NOTE: assume SIDE, READPRJ, ID is always there
                    print('\t\tSIDE_two')
                    side_cols = ['ID', 'SIDE', 'READPRJ', icol]
                    basic_info_cols = ['ID', 'SIDE', 'READPRJ']
                    tmp_data_1[basic_info_cols] = tmp_data_1[basic_info_cols].astype('str')
                    tmp_data_2[basic_info_cols] = tmp_data_2[basic_info_cols].astype('str')
                    
                    tmp_data_1.loc[tmp_data_1['SIDE'].str.contains('1'), 'SIDE'] = '1: Right'
                    tmp_data_1.loc[tmp_data_1['SIDE'].str.contains('2'), 'SIDE'] = '2: Left'
                    tmp_data_2.loc[tmp_data_2['SIDE'].str.contains('1'), 'SIDE'] = '1: Right'
                    tmp_data_2.loc[tmp_data_2['SIDE'].str.contains('2'), 'SIDE'] = '2: Left'
                    
                    tmp_concat_data = pd.concat([tmp_data_1[side_cols], tmp_data_2[side_cols]])
                    tmp_concat_data.drop_duplicates(tmp_concat_data.columns, inplace=True)
                    # TODO: update merge function
                    if side_cts == 0:
                        concat_side_df = tmp_concat_data
                    else:
                        concat_side_df = pd.merge(concat_side_df, tmp_concat_data, how='outer', on=['ID', 'SIDE', 'READPRJ'])
                    side_cts += 1
                else:
                    print('\t\tUnique')
                    tmp_data_1.set_index('ID', inplace=True)
                    tmp_data_2.set_index('ID', inplace=True)
                    tmp_merge_data = pd.concat([tmp_data_1, tmp_data_2], axis=0, sort=False)
                    direct_merge_df = pd.concat([direct_merge_df, tmp_merge_data[icol]], axis=1)
        if len(tmp_files) == 1:
            print("\tOnly one file has this column!")
            tmp_filename = dataDir+'/'+tmp_files[0]
            tmp_data = pd.read_csv(tmp_filename, sep='|')
            tmp_data.columns = tmp_data.columns.str.upper()
            if 'SIDE' not in tmp_data.columns.values.tolist():
                print('\t\tUnique one')
                tmp_data.set_index('ID', inplace=True)
                direct_merge_df = pd.concat([direct_merge_df, tmp_data[icol]], axis=1)
            else:
                print('\t\tSIDE_one')
                side_cols = ['ID', 'SIDE', 'READPRJ', icol]
                basic_info_cols = ['ID', 'SIDE', 'READPRJ']
                tmp_data[basic_info_cols] = tmp_data[basic_info_cols].astype('str')
                tmp_data.loc[tmp_data['SIDE'].str.contains('1'), 'SIDE'] = '1: Right'
                tmp_data.loc[tmp_data['SIDE'].str.contains('2'), 'SIDE'] = '2: Left'
                if side_cts == 0:
                    concat_side_df = tmp_data[side_cols]
                else:
                    concat_side_df = pd.merge(concat_side_df, tmp_data[side_cols], how='outer', on=['ID', 'SIDE', 'READPRJ'])
                side_cts += 1
        if len(tmp_files) == 0:
            print(tmp_files)
            raise ValueError("CANNOT FIND THIS COLUMN!")
        if len(tmp_files) > 2:
            print(tmp_files)
            raise ValueError("Too many existed files!")
        ict += 1
    else:
        print(icol+' is already existed!')
    print("\tTime cost "+str(dt.now()-cst))
    print("\tTotal time cost "+str(dt.now()-st)+"\n")


print(direct_merge_df)
direct_merge_df.to_csv(resDir+"/outcome_direct_merge_dataframe_"+exp_id+".csv")
direct_merge_df.to_pickle(resDir+"/outcome_direct_merge_dataframe_"+exp_id+".pkl")

print(concat_side_df)
concat_side_df.to_csv(resDir+"/outcome_concat_side_dataframe_"+exp_id+".csv")
concat_side_df.to_pickle(resDir+"/outcome_concat_side_dataframe_"+exp_id+".pkl")
print('Total time cost '+str(dt.now()-st))

side_cols = concat_side_df.columns.tolist()[3:-1]
new_side_cols = []
all_outcome_df=direct_merge_df

for i in side_cols:
    new_side_cols = new_side_cols+[i+'L']
    new_side_cols = new_side_cols+[i+'R']
all_outcome_df = pd.concat([direct_merge_df, pd.DataFrame(columns=new_side_cols, index=direct_merge_df.index)], axis=1)

for i in side_cols:
    tmp_side_df = pd.concat([concat_side_df.iloc[:,[0,1,2]], concat_side_df[i]], axis=1)
    for isd in [0, 1, 2, 3, 4]:
        tmp_side_df.replace(to_replace=str(isd)+': '+str(isd), value=isd, inplace=True)
    tmp_side_df[i] = pd.to_numeric(tmp_side_df[i], errors='coerce')
    left_i=i+'L'
    right_i=i+'R'
    pid_group = tmp_side_df.groupby('ID')
    for pid, gp in pid_group:
        gp.dropna(inplace=True)
        left_mask = gp['SIDE'].str.contains('Left')
        if left_mask.sum() == 1:
            all_outcome_df.loc[int(pid), left_i] = gp.loc[gp['SIDE'].str.contains('Left'), i].tolist()[0]
        elif left_mask.sum() > 1:
            left_gp = gp.loc[left_mask]
            if left_gp[i].unique().shape[0] == 1:
                all_outcome_df.loc[int(pid), left_i] = left_gp[i].unique().tolist()[0]
            else:
                all_outcome_df.loc[int(pid), left_i] = max(left_gp[i].unique().tolist())
                warnings.warn("Conflict results from different project!!!")

        right_mask = gp['SIDE'].str.contains('Right')
        if right_mask.sum() == 1:
            all_outcome_df.loc[int(pid), right_i] = gp.loc[gp['SIDE'].str.contains('Right'), i].tolist()[0]
        elif right_mask.sum() > 1:
            right_gp = gp.loc[right_mask]
            if right_gp[i].unique().shape[0] == 1:
                all_outcome_df.loc[int(pid), right_i] = right_gp[i].unique().tolist()[0]
            else:
                all_outcome_df.loc[int(pid), right_i] = max(right_gp[i].unique().tolist())
                warnings.warn("Conflict results from different project!!!")
print(all_outcome_df)
all_outcome_df.to_csv(resDir+"/outcome_all_dataframe_"+exp_id+".csv")
all_outcome_df.to_pickle(resDir+"/outcome_all_dataframe_"+exp_id+".pkl")