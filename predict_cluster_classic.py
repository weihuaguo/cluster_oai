#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 13:59:48 2023

predict the clusters

@author: weihua
"""
import os
import pandas as pd
import seaborn as sn
import numpy as np
from datetime import datetime
import math
from datetime import datetime as dt
from sklearn import linear_model, svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import cross_val_score
from keras.utils.np_utils import to_categorical
from sklearn.metrics import classification_report, roc_auc_score, accuracy_score

st = dt.now()
mainDir = "/mnt/sda1/OAI_Data"
#mainDir="/home/weihua/mnts/smb_plee/Group/weihua/data_summary"

dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
#resDir = mainDir

clst_dir = mainDir + "/kmean_cluster_12252020"
exp_id = "12112020_data_clean"
clean_co = "v25"
plot_prefix = clst_dir+'/clean_'+clean_co+'_'
result_folder = "direct_predict_cluster_230109"

importance_source = "rfimp"
importance_number = range(2, 51)

lr_flag=True
rf_flag=False
sv_flag=False
model_name="lr"

test_size = 0.1

oh_merge_df=pd.read_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols.pkl')
print("Start to read the imputed data...")
imp_oh_merge_df = pd.read_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols_imp.pkl')
input_df = imp_oh_merge_df

cls_rd_res = pd.read_excel(plot_prefix+'cluster4_kmean_pca_umap_res.xlsx')
oh_cls_df = pd.get_dummies(cls_rd_res[['ID', 'kmean_pca']], columns=['kmean_pca'])
output_df = cls_rd_res
output_df.index = output_df['ID']

print("Start to scale the data...")
scaled_oh_data = StandardScaler().fit_transform(imp_oh_merge_df)
mm_scaler = MinMaxScaler().fit(input_df)
scale_df = pd.DataFrame(mm_scaler.transform(input_df), index=input_df.index, columns=input_df.columns)
input_df=scale_df.loc[output_df['ID']]

rsub_dir=os.path.join(clst_dir, result_folder)
if not os.path.exists(rsub_dir):
    os.makedirs(rsub_dir)
rpf = rsub_dir+"/"+result_folder+"_"
ona_mask=output_df['ID']
tmp_X=input_df
x=tmp_X
y=output_df[output_df.columns[output_df.columns.str.contains("kmean_pca")]] #np.ravel(
y_cat = to_categorical(y)

print("Using importance from random forest...")
rfr = RandomForestClassifier(random_state=0)
rfr.fit(x, np.ravel(y))
importance = rfr.feature_importances_

imp_df = pd.DataFrame(importance, index = tmp_X.columns, columns = ['sig'])
imp_df['varname'] = tmp_X.columns
imp_df['importance_source'] = importance_source
imp_df.to_excel(rpf+"importance_dataframe.xlsx")
#rmse, cor, rele, cvn, impn, model, scores = list(), list(), list(), list(), list(), list(), list()
mrpf = rpf + model_name + "_"
for iimp in importance_number:
    select_x = x.iloc[:,sorted(range(len(importance)), key=lambda i: importance[i])[-iimp:]]
    ss = ShuffleSplit(n_splits=int(2*1/test_size), test_size=0.10, random_state=0)
    cvi = 0
    for train_index, test_index in ss.split(range(x.shape[0])):
        train_x = select_x.iloc[train_index]
        train_y = y.iloc[train_index]
        test_x = select_x.iloc[test_index]
        test_y = y.iloc[test_index]
        
        test_res = output_df[output_df.columns[output_df.columns.str.contains("kmean_pca")]]
        test_res = test_res.iloc[test_index]
        test_res = test_res.rename(columns={"kmean_pca":"measurement"})
        
        if lr_flag:
            print("\tLogistic regression\t" + model_name+"\t"+str(cvi))
            model = linear_model.LogisticRegression(solver='newton-cg', max_iter=1000, n_jobs=-1, random_state=66)
            
            record = model.fit(train_x, np.ravel(train_y))
            predict_y_class = model.predict(test_x)
            predict_y = model.predict_proba(test_x)
            predict_ty_class = model.predict(train_x)
            predict_ty = model.predict_proba(train_x)
            
            test_res['model'] = model_name
            
            test_res['prediction'] = predict_y_class
            test_res['cv_counter'] = cvi
            test_res['roc_auc_ovo'] = roc_auc_score(to_categorical(test_y), predict_y, multi_class="ovo")
            test_res['roc_auc_ovr'] = roc_auc_score(to_categorical(test_y), predict_y, multi_class="ovr")
            test_res['accuracy'] = accuracy_score(test_res['measurement'], test_res['prediction'], normalize = True)
            
            test_res['roc_auc_ovo_train'] = roc_auc_score(to_categorical(train_y), predict_ty, multi_class="ovo")
            test_res['roc_auc_ovr_train'] = roc_auc_score(to_categorical(train_y), predict_ty, multi_class="ovr")
            test_res['accuracy_train'] = accuracy_score(train_y, predict_ty_class, normalize = True)
            
            
        if rf_flag:
            print("\tRandom forest\t" + model_name+"\t"+str(cvi))
            model = RandomForestClassifier(n_estimators=100, criterion = "entropy", random_state = 66)
            record = model.fit(train_x, np.ravel(train_y))
            predict_y_class = model.predict(test_x)
            predict_y = model.predict_proba(test_x)
            predict_ty_class = model.predict(train_x)
            predict_ty = model.predict_proba(train_x)
            
            test_res['model'] = model_name
            
            test_res['prediction'] = predict_y_class
            test_res['cv_counter'] = cvi
            test_res['roc_auc_ovo'] = roc_auc_score(to_categorical(test_y), predict_y, multi_class="ovo")
            test_res['roc_auc_ovr'] = roc_auc_score(to_categorical(test_y), predict_y, multi_class="ovr")
            test_res['accuracy'] = accuracy_score(test_res['measurement'], test_res['prediction'], normalize = True)
            
            test_res['roc_auc_ovo_train'] = roc_auc_score(to_categorical(train_y), predict_ty, multi_class="ovo")
            test_res['roc_auc_ovr_train'] = roc_auc_score(to_categorical(train_y), predict_ty, multi_class="ovr")
            test_res['accuracy_train'] = accuracy_score(train_y, predict_ty_class, normalize = True)
            
        if sv_flag:
            print("\tSVM\t" + model_name+"\t"+str(cvi))
            model = svm.SVC(kernel='sigmoid', probability=True)
            record = model.fit(train_x, np.ravel(train_y))
            predict_y_class = model.predict(test_x)
            predict_y = model.predict_proba(test_x)
            predict_ty_class = model.predict(train_x)
            predict_ty = model.predict_proba(train_x)
            
            test_res['model'] = model_name
            
            test_res['prediction'] = predict_y_class
            test_res['cv_counter'] = cvi
            test_res['roc_auc_ovo'] = roc_auc_score(to_categorical(test_y), predict_y, multi_class="ovo")
            test_res['roc_auc_ovr'] = roc_auc_score(to_categorical(test_y), predict_y, multi_class="ovr")
            test_res['accuracy'] = accuracy_score(test_res['measurement'], test_res['prediction'], normalize = True)
            
            test_res['roc_auc_ovo_train'] = roc_auc_score(to_categorical(train_y), predict_ty, multi_class="ovo")
            test_res['roc_auc_ovr_train'] = roc_auc_score(to_categorical(train_y), predict_ty, multi_class="ovr")
            test_res['accuracy_train'] = accuracy_score(train_y, predict_ty_class, normalize = True)

        if cvi == 0:
            merge_res = test_res
        else:
            merge_res = pd.concat([merge_res, test_res], axis=0)
        cvi += 1
    merge_res['importance_source'] = importance_source
    merge_res['importance_number'] = iimp
    merge_res.to_csv(mrpf+"feat"+str(iimp)+"_manual_cv_predict_measure_df.csv")      
