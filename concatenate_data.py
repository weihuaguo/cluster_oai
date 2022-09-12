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
from pandas.api.types import is_numeric_dtype as inumdt
from datetime import datetime as dt
import matplotlib.pyplot as plt

# mainDir = "G:/My Drive/Weihua/OAI_Public_Data"
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
exp_id = "12112020_data_clean"

enroll_txt = "Enrollees.txt"
ann_xlsx = "AllClinical00_V5_column_annotation_WG_clean.xlsx"

total_file_version_merge = False

all_pat_info = pd.read_csv(dataDir+"/"+enroll_txt, sep='|', index_col=0)
print(all_pat_info)
direct_merge_df = all_pat_info

all_col_oi = pd.read_excel(mainDir+'/'+ann_xlsx, sheet_name='Baseline', index_col='Variables')
print(all_col_oi)

all_txt_files = glob.glob(dataDir+"/*[0-9].txt")
# print(all_txt_files)
print(len(all_txt_files))

with open(resDir+'/file_sum_dictionary.pkl', 'rb') as handle:
    file_dict = pkl.load(handle)
with open(resDir+'/col_sum_dictionary.pkl', 'rb') as handle:
    col_dict = pkl.load(handle)

all_files = col_dict['ID']['filename']
all_cols = col_dict.keys()

all_col_oi.index = all_col_oi.index.str.upper()
col_oi = all_col_oi.index.values.tolist()
ict = 0
side_cts = 0
st = dt.now()
for icol in col_oi:
    print(icol+'\t'+str(ict))
    cst = dt.now()
    if icol not in direct_merge_df.columns.values.tolist():
        if all_col_oi.loc[icol,'Usage'] == "Input":
            tmp_files = col_dict[icol]['filename']
            if len(tmp_files) == 2:
                print("\tTwo files have this column!")
                tmp_pat_num = col_dict[icol]['patient_num']
                tmp_na_num = col_dict[icol]['na_num']
                # NOTE: assume equal patient num represents the same index
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
        #                print(direct_merge_df)
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

#                    print(tmp_files)
                    ol_pid = list(set(tmp_data_1['ID'])&set(tmp_data_2['ID']))
#                    print(len(ol_pid))
#                    print(tmp_data_1.shape)
#                    print(tmp_data_2.shape)
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
direct_merge_df.to_csv(resDir+"/input_direct_merge_dataframe_"+exp_id+".csv")
direct_merge_df.to_pickle(resDir+"/input_direct_merge_dataframe_"+exp_id+".pkl")
print('Total time cost '+str(dt.now()-st))
