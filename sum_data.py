import os
os.system('clear')
import itertools
import pandas as pd
import pickle as pkl
import glob
from pandas.api.types import is_numeric_dtype as inumdt
from datetime import datetime as dt
from collections import defaultdict
import matplotlib.pyplot as plt

# mainDir = "G:/My Drive/Weihua/OAI_Public_Data"
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"

enroll_txt = "Enrollees.txt"

all_pat_info = pd.read_csv(dataDir+"/"+enroll_txt, sep='|', index_col=0)
print(all_pat_info)

all_txt_files = glob.glob(dataDir+"/*[0-9].txt")
# print(all_txt_files)
print(len(all_txt_files))

file_dict = {}
col_dict = {}
all_colnames = []
ast = dt.now()
for itf in all_txt_files:
    st = dt.now()
    tmp_data = pd.read_csv(itf, sep='|', encoding = "ISO-8859-1")
    tmp_data.columns = tmp_data.columns.str.upper()
    tmp_name = itf.split('/')[-1].split('.')[0]
    print(tmp_name)
    print(tmp_data.shape)
    file_dict[tmp_name] = {}
    file_dict[tmp_name]['filename'] = itf.split('/')[-1]
    file_dict[tmp_name]['row_num'] = tmp_data.shape[0]
    tmp_colname = tmp_data.columns.values.tolist()

    for icol in tmp_data.columns.values.tolist():
        print(icol)
        if icol in col_dict:
#            print(icol+' is existed!')
            col_dict[icol]['filename'] = col_dict[icol]['filename'] + [itf.split('/')[-1]]
            col_dict[icol]['type'] =col_dict[icol]['type'] + [str(tmp_data[icol].dtype)]
            col_dict[icol]['unique_item_num'] = col_dict[icol]['unique_item_num'] + [tmp_data[icol].value_counts().shape[0]]
            col_dict[icol]['row_num'] = col_dict[icol]['row_num'] + [tmp_data.shape[0]]
            col_dict[icol]['patient_num'] = col_dict[icol]['patient_num'] + [len(tmp_data['ID'].unique().tolist())]
           
            if inumdt(tmp_data[icol]):
                col_dict[icol]['na_num'] = col_dict[icol]['na_num'] + [tmp_data[icol].isna().sum()]
            else:
                col_dict[icol]['na_num'] = col_dict[icol]['na_num'] + [tmp_data[icol].str.contains('Missing').sum()]

        else:
#            print('New column '+icol+'\n')
            col_dict[icol] = {}
            col_dict[icol]['filename'] = [itf.split('/')[-1]]
            col_dict[icol]['type'] = [str(tmp_data[icol].dtype)]
            col_dict[icol]['unique_item_num'] = [tmp_data[icol].value_counts().shape[0]]
            col_dict[icol]['row_num'] = [tmp_data.shape[0]]
            col_dict[icol]['patient_num'] = [len(tmp_data['ID'].unique().tolist())]
            if inumdt(tmp_data[icol]):
                col_dict[icol]['na_num'] = [tmp_data[icol].isna().sum()]
            else:
                col_dict[icol]['na_num'] = [tmp_data[icol].str.contains('Missing').sum()]

    if "ID" in tmp_colname:
        tmp_multi_spl_num = sum(tmp_data.ID.value_counts()>1)
        tmp_pat_num = tmp_data.ID.value_counts().shape[0]
        file_dict[tmp_name]['patient_num'] = tmp_pat_num
        file_dict[tmp_name]['multi_spl_num'] = tmp_multi_spl_num
        tmp_colname.remove("ID")
    else:
        file_dict[tmp_name]['patient_num'] = 0
        file_dict[tmp_name]['multi_spl_num'] = 0

    if "VERSION" in tmp_colname:
        tmp_version = tmp_data.VERSION.unique().tolist()
        file_dict[tmp_name]['version_num'] = len(tmp_version)
        file_dict[tmp_name]['version'] = tmp_version
        tmp_colname.remove("VERSION")
    else:
        file_dict[tmp_name]['version'] = ['NA']
        file_dict[tmp_name]['version_num'] = 0

    if "SIDE" in tmp_colname:
        tmp_side_num = tmp_data.SIDE.value_counts().shape[0]
        file_dict[tmp_name]['side_num'] = tmp_side_num
        tmp_colname.remove('SIDE')
    else:
        file_dict[tmp_name]['side_num'] = 0

    if "READPRJ" in tmp_colname:
        tmp_prj_num = tmp_data.READPRJ.value_counts().shape[0]
        file_dict[tmp_name]['readprj_num'] = tmp_prj_num
        file_dict[tmp_name]['readprj'] = tmp_data.READPRJ.unique().tolist()
        tmp_colname.remove('READPRJ')
    else:
        file_dict[tmp_name]['readprj'] = ['NA']
        file_dict[tmp_name]['readprj_num'] = 0

    file_dict[tmp_name]['column_num'] = len(tmp_colname)
    all_colnames = all_colnames + tmp_colname
    print("\tTime cost: "+str(dt.now()-st))


file_df = pd.DataFrame.from_dict(file_dict).T
file_df.to_excel(resDir+'/file_sum_excel.xlsx')
with open(resDir+'/file_sum_dictionary.pkl', 'wb') as handle:
    pkl.dump(file_dict, handle, protocol=pkl.HIGHEST_PROTOCOL)

for key in col_dict:
    value=col_dict[key]
    value["num"] = len(value["filename"])

col_df = pd.DataFrame.from_dict(col_dict).T
col_df.to_excel(resDir+'/column_sum_excel.xlsx')
with open(resDir+'/col_sum_dictionary.pkl', 'wb') as handle:
    pkl.dump(col_dict, handle, protocol=pkl.HIGHEST_PROTOCOL)

print(len(all_colnames))
print(len(list(set(all_colnames))))
print("Total time cost "+str(dt.now()-ast))
## All the columns * filesa
