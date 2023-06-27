#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  14 13:37:32 2023
Test different supervised learning models for predicting WOMAC total scores
@author: weihua
""" 
def mape(actual, pred): 
    actual, pred = np.array(actual), np.array(pred)
    return np.mean(np.abs((actual - pred) / actual)) * 100

import os
import re
import math
import pickle
import pandas as pd
import numpy as np
import tensorflow as tf
import scipy.stats
import matplotlib.pyplot as plt

from datetime import datetime
from sklearn.svm import SVR
from sklearn import linear_model
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from collections import Counter
from sklearn.metrics import mean_squared_error
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor


data_dir="/mnt/sda1/OAI_Data/data_summary"
exp_id = "12112020_data_clean"
clean_co = "v25"

importance_source = "randomforest"
importance_number = range(2, 51)
test_size = 0.10
lr_flag = False
rf_flag = True
sv_flag = False
cnn_flag = False

dnm = 3/4
sqr_nlayer = 0
ln_rt = 0.001
n_epoch = 100

input_act_func = None
hiden_act_func = None
output_act_func = None


result_folder="impt_"+importance_source+"_sl_imputed_230214"
result_dir=os.path.join(data_dir, result_folder)
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

output_df=pd.read_csv(data_dir+"/outcome_all_dataframe_"+exp_id+"_real_date_year.csv", index_col=1)
input_df=pd.read_pickle(data_dir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols_imp.pkl')
input_df.describe(include="all").to_csv(data_dir+'/input_one_hot_manual_merge_dataframe_'+exp_id+'_sbj_clean_'+clean_co+'descriptive_stat.csv')

print("Clean binary categorical variable into one...")
cate_cols=input_df.columns[input_df.columns.str.contains(":")].tolist()
var_df = pd.DataFrame([re.split("_[0-9]",i) for i in cate_cols], index = cate_cols)
var_df['var_cts'] = var_df.groupby([0])[0].transform('count')
two_df = var_df.loc[var_df['var_cts'] == 2]
rm_df = two_df.drop_duplicates(subset=[0])

raw_input_df = input_df
input_df = input_df.drop(columns=rm_df.index)

mm_scaler = MinMaxScaler().fit(input_df)
scale_df = pd.DataFrame(mm_scaler.transform(input_df), index=input_df.index, columns=input_df.columns)
print(output_df.shape)
print(scale_df.shape)
input_df=scale_df.loc[output_df.index] # NOTE: the output df includes all the samples!
outcome_patterns=['WOMTSL', 'WOMTSR']
year_of_interest=["Y04", "Y08"]

i = 0
for yoi in year_of_interest:
    for op in outcome_patterns:
        rsub_dir=os.path.join(result_dir, result_folder+"_"+yoi+op)
        if not os.path.exists(rsub_dir):
            os.makedirs(rsub_dir)
        rpf = rsub_dir+"/"+yoi+op+"_sl_"
        print(yoi+op)
        ona_mask=~output_df[yoi+op].isna()
        tmp_X=input_df.loc[ona_mask,:]
        x=tmp_X
        y=output_df.loc[x.index,yoi+op]
        
        if importance_source == "linear":
            print("Using importance from linear regression...")
            lnrgr = linear_model.LinearRegression()
            lnrgr.fit(x, y)
            importance = lnrgr.coef_
        if importance_source == "randomforest":
            print("Using importance from random forest...")
            rfr = RandomForestRegressor(random_state=66)
            rfr.fit(x, y)
            importance = rfr.feature_importances_
        imp_df = pd.DataFrame(importance, index = tmp_X.columns, columns = ['sig'])
        imp_df['varname'] = tmp_X.columns
        imp_df['importance_source'] = importance_source
        imp_df['outname'] = yoi+op
        imp_df.to_excel(rpf+"importance.xlsx")
        rmse, cor, rele, cvn, impn, mdl = list(), list(), list(), list(), list(), list()
        for iimp in importance_number:
            select_x = x.iloc[:,sorted(range(len(importance)), key=lambda i: importance[i])[-iimp:]]
            if cnn_flag:
                mrpf = rpf + "cnnopt_"
                s=select_x.shape[1]
                cnn_shape=[]
                print("Determing CNN shape")
                while math.floor(s) > 1:
                    s=math.floor(s*dnm)
                    if math.floor(s) > 1:
                        cnn_shape.append(s)
                ss = ShuffleSplit(n_splits=int(2*1/test_size), test_size=test_size, random_state=0)
                cvi = 0
                for train_index, test_index in ss.split(range(x.shape[0])):
                    train_x = select_x.iloc[train_index]
                    train_y = y.iloc[train_index]
                    test_x = select_x.iloc[test_index]
                    test_y = y.iloc[test_index]
                    print("Building the CNN model")
                    model = Sequential()
                    model.add(Dense(select_x.shape[1], input_dim = select_x.shape[1], activation = input_act_func))
                    if sqr_nlayer > 0:
                        for ii in range(sqr_nlayer):
                            model.add(Dense(math.floor(select_x.shape[1]*dnm), activation = hiden_act_func))
                    for cs in cnn_shape:
                        model.add(Dense(cs, activation = hiden_act_func))
                    model.add(Dense(1, activation = output_act_func))
                    opt = tf.keras.optimizers.Adam(learning_rate=ln_rt)
                    model.compile(loss = 'mean_squared_error', optimizer=opt)
                    print(model.summary())
                    fst = datetime.now()
                    record = model.fit(train_x, train_y, epochs = n_epoch, batch_size = 10, use_multiprocessing=True)
                    plt.figure()
                    plt.xlabel("Number of Epochs")
                    plt.ylabel("Loss")
                    plt.plot(range(1, n_epoch+1), record.history['loss'])
                    plt.savefig(mrpf+"feat"+str(iimp)+"_fold"+str(cvi)+"_loss.png", dpi=300, bbox_inches='tight')
                    plt.clf()
                    plt.close('all')
                    predict_y = model.predict(test_x, use_multiprocessing=True)
                    test_res = pd.DataFrame(test_y)
                    test_res = test_res.rename(columns={yoi+op: "measurement"})
                    test_res['prediction'] = predict_y
                    test_res['outcome'] = yoi+op
                    test_res['model'] = 'cnnopt_epch'+str(n_epoch)
                    test_res['cv_counter'] = cvi
                    correct_y = [i[0] for i in predict_y]
                    rmse.append(math.sqrt(mean_squared_error(test_y, predict_y)))
                    cor.append(scipy.stats.pearsonr(test_y, correct_y)[0])
                    cvn.append(cvi)
                    impn.append(iimp)
                    mdl.append('cnn_epch'+str(n_epoch))
                    if cvi == 0:
                        merge_res = test_res
                    else:
                        merge_res = pd.concat([merge_res, test_res], axis=0)
                    cvi += 1
                merge_res['importance_source'] = importance_source
                merge_res['importance_number'] = iimp
                merge_res.to_csv(mrpf+"feat"+str(iimp)+"_manual_cv_predict_measure_df.csv")

            if lr_flag:
                mrpf = rpf + "linear_"
                ss = ShuffleSplit(n_splits=int(2*1/test_size), test_size=0.10, random_state=0)
                cvi = 0
                for train_index, test_index in ss.split(range(x.shape[0])):
                    train_x = select_x.iloc[train_index]
                    train_y = y.iloc[train_index]
                    test_x = select_x.iloc[test_index]
                    test_y = y.iloc[test_index]
                    lnrgr = linear_model.LinearRegression()
                    lnrgr.fit(train_x, train_y)
                    predict_y = lnrgr.predict(test_x)
                    test_res = pd.DataFrame(test_y)
                    test_res = test_res.rename(columns={yoi+op: "measurement"})
                    test_res['prediction'] = predict_y
                    test_res['outcome'] = yoi+op
                    test_res['model'] = 'linear'
                    test_res['cv_counter'] = cvi
                    rmse.append(math.sqrt(mean_squared_error(test_y, predict_y)))
                    cor.append(scipy.stats.pearsonr(test_y, predict_y)[0])
                    rele.append(mape(test_y, predict_y))
                    cvn.append(cvi)
                    impn.append(iimp)
                    mdl.append("linear")
                    
                    if cvi == 0:
                        merge_res = test_res
                    else:
                        merge_res = pd.concat([merge_res, test_res], axis=0)
                    cvi += 1
                merge_res['importance_source'] = importance_source
                merge_res['importance_number'] = iimp
                merge_res.to_csv(mrpf+"feat"+str(iimp)+"_manual_cv_predict_measure_df.csv")
                
            if rf_flag:
                trees=[10, 20, 40, 60, 80, 100]
                for tree in trees:
                    mrpf = rpf + "rft" + str(tree) + "_"
                    ss = ShuffleSplit(n_splits=int(2*1/test_size), test_size=0.10, random_state=66)
                    cvi = 0
                    for train_index, test_index in ss.split(range(x.shape[0])):
                        train_x = select_x.iloc[train_index]
                        train_y = y.iloc[train_index]
                        test_x = select_x.iloc[test_index]
                        test_y = y.iloc[test_index]
                        lnrgr = RandomForestRegressor(n_estimators=tree, random_state=66)
                        lnrgr.fit(train_x, train_y)
                        predict_y = lnrgr.predict(test_x)
                        test_res = pd.DataFrame(test_y)
                        test_res = test_res.rename(columns={yoi+op: "measurement"})
                        test_res['prediction'] = predict_y
                        test_res['outcome'] = yoi+op
                        test_res['model'] = 'rft'+str(tree)
                        test_res['cv_counter'] = cvi
                        rmse.append(math.sqrt(mean_squared_error(test_y, predict_y)))
                        cor.append(scipy.stats.pearsonr(test_y, predict_y)[0])
                        rele.append(mape(test_y, predict_y))
                        cvn.append(cvi)
                        impn.append(iimp)
                        mdl.append("rft"+str(tree))
                        
                        if cvi == 0:
                            merge_res = test_res
                        else:
                            merge_res = pd.concat([merge_res, test_res], axis=0)
                        cvi += 1
                    merge_res['importance_source'] = importance_source
                    merge_res['importance_number'] = iimp
                    merge_res.to_csv(mrpf+"feat"+str(iimp)+"_manual_cv_predict_measure_df.csv")
                
            if sv_flag:
                kernels = ['linear', 'poly', 'rbf', 'sigmoid']
                for kernel in kernels:
                    mrpf = rpf + "svr_" + kernel +"_"
                    ss = ShuffleSplit(n_splits=int(2*1/test_size), test_size=0.10, random_state=66)
                    cvi = 0
                    for train_index, test_index in ss.split(range(x.shape[0])):
                        train_x = select_x.iloc[train_index]
                        train_y = y.iloc[train_index]
                        test_x = select_x.iloc[test_index]
                        test_y = y.iloc[test_index]
                        lnrgr = SVR(kernel=kernel, C=100, gamma="auto")
                        lnrgr.fit(train_x, train_y)
                        predict_y = lnrgr.predict(test_x)
                        test_res = pd.DataFrame(test_y)
                        test_res = test_res.rename(columns={yoi+op: "measurement"})
                        test_res['prediction'] = predict_y
                        test_res['outcome'] = yoi+op
                        test_res['model'] = 'svr'+kernel
                        test_res['cv_counter'] = cvi
                        rmse.append(math.sqrt(mean_squared_error(test_y, predict_y)))
                        cor.append(scipy.stats.pearsonr(test_y, predict_y)[0])
                        rele.append(mape(test_y, predict_y))
                        cvn.append(cvi)
                        impn.append(iimp)
                        mdl.append("svr"+kernel)
                        
                        if cvi == 0:
                            merge_res = test_res
                        else:
                            merge_res = pd.concat([merge_res, test_res], axis=0)
                        cvi += 1
                    merge_res['importance_source'] = importance_source
                    merge_res['importance_number'] = iimp
                    merge_res.to_csv(mrpf+"feat"+str(iimp)+"_manual_cv_predict_measure_df.csv")
        if any([lr_flag, rf_flag, sv_flag, cnn_flag]):
            score_df = pd.DataFrame([rmse, cor, rele, cvn, impn, mdl], index=['rmse', 'r', 'mape', 'cvn', 'nfeats', 'model']).T
            score_df.to_csv(mrpf+"manual_cv_score_df.csv")
