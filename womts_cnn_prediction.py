#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 13:16:50 2022

cnn regression

from direct_super_learn_cnn_v3.py
@author: weihua
"""

import tensorflow as tf
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor

import os
import sys
import re
import time
import pickle
import math
import scipy
import shutil
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from sklearn.svm import SVC
from sklearn import linear_model
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import ShuffleSplit
from sklearn.impute import KNNImputer
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

data_dir="/mnt/sda1/OAI_Data/data_summary"
data_dir="/home/weihua/mnts/smb_plee/Group/weihua/data_summary"

exp_id = "12112020_data_clean"
slc_id = "predictor_vis_220125"
clean_co = "v25"

importance_source = "random" # "randomforest"/"linear"/"random"
importance_filter = "cluster"
importance_number = range(25,26)
test_size = 0.5
dnm = 3/4 #1/3, 1/2, 2/3, 3/4, 0.9
ln_rt = 0.001 # 1, 0.1, 0.01, 0.001

input_act_funcs = [None]
hidden_act_funcs = [None]
output_act_funcs = ["relu"]
general_init = ["random_normal"]
n_add_layers = [0]
r1 = int(sys.argv[1])
r2 = int(sys.argv[2])
random_var = range(r1,r2)

c = 0
hyperparam = {}
hyperparam['expr_id'] = []
hyperparam['input_act_func'] = []
hyperparam['hidden_act_func'] = []
hyperparam['output_act_func'] = []
hyperparam['sqr_nlayer'] = []
hyperparam['iniz'] = []
hyperparam['random_var_cts'] = []

for iiaf in input_act_funcs:
    for ihaf in hidden_act_funcs:
        for ioaf in output_act_funcs:
            for inal in n_add_layers:
                for igi in general_init:
                    for ir in random_var:
                        hyperparam['expr_id'].append("v"+str(ir))
                        hyperparam['input_act_func'].append(iiaf)
                        hyperparam['hidden_act_func'].append(ihaf)
                        hyperparam['output_act_func'].append(ioaf)
                        hyperparam['sqr_nlayer'].append(inal)
                        hyperparam['iniz'].append(igi)
                        hyperparam['random_var_cts'].append(ir)
                        c += 1
hyperparam_df = pd.DataFrame.from_dict(hyperparam)
hyperparam_df.to_csv(data_dir+"/impt_"+importance_source+"_hyperparam_sl_cnn_imputed_220902.csv")

output_df=pd.read_csv(data_dir+"/outcome_all_dataframe_"+exp_id+"_real_date_year.csv", index_col=1)
cor_df=pd.read_csv(data_dir+"/"+slc_id+"/"+slc_id+"_cor_merge.csv", index_col=0)
pval_cor_df=pd.read_csv(data_dir+"/"+slc_id+"/"+slc_id+"_pval_merge.csv", index_col=0)
sn_cor_df=pd.read_csv(data_dir+"/"+slc_id+"/"+slc_id+"_sn_merge.csv", index_col=0)
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
input_df=scale_df.loc[output_df.index]
outcome_patterns=['WOMTSL', 'WOMTSR']
year_of_interest=["Y04", "Y08"]

for icc in range(0, c):
    st = datetime.now()
    input_act_func = hyperparam['input_act_func'][icc]
    hidden_act_func = hyperparam['hidden_act_func'][icc]
    output_act_func = hyperparam['output_act_func'][icc]
    sqr_nlayer = hyperparam['sqr_nlayer'][icc]
    gnrl_iniz = hyperparam['iniz'][icc]

    if importance_filter != "none":
        result_folder="impt_"+importance_source+"_"+importance_filter+"_sl_cnn_imputed_221022_v"+hyperparam['expr_id'][icc] ## NOTE: Change here for different experiments!
    else:
        result_folder="impt_"+importance_source+"_sl_cnn_imputed_221017_v"+hyperparam['expr_id'][icc] ## NOTE: Change here for different experiments!
    result_dir=os.path.join(data_dir, result_folder)

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    i = 0
    n_epoch = 100
    for yoi in year_of_interest:
        for op in outcome_patterns:
            rsub_dir=os.path.join(result_dir, result_folder+"_"+yoi+op+"_"+str(n_epoch))
            if not os.path.exists(rsub_dir):
                os.makedirs(rsub_dir)
            rpf = rsub_dir+"/"+yoi+op+"_sl_cnn_nepoch_"+str(n_epoch)+"_"
            print(yoi+op)
            ona_mask=~output_df[yoi+op].isna()
            tmp_X=input_df.loc[ona_mask,:]
            x=tmp_X
            y=output_df.loc[x.index,yoi+op]
            if importance_source == "random":
                print("Using importance from random selection...")
                if importance_filter == "rf":
                    print("\tFiltered by random forest importances...")
                    rfr = RandomForestRegressor(random_state=0)
                    rfr.fit(x, y)
                    filter_importance = rfr.feature_importances_
                    select_order = sorted(range(len(filter_importance)), key=lambda i: filter_importance[i])[-30:]
                    importance = filter_importance
                    importance[select_order] = random.sample(range(1000,10000), 30)

                elif importance_filter == "cluster":
                    print("\tFiltered by clustering top10 variables...")
                    marker_df = pd.read_csv(data_dir+"/clean_v25_final_cluster4_kmeans_direct_knn2imp_12112020_data_clean_marker_df_largeB.csv", index_col=0)
                    marker_df['abs_logfc'] = abs(marker_df['logfc'])
                    sig_df = marker_df[abs(marker_df['adjp'])<0.05]
                    sig_top_df = sig_df.sort_values('abs_logfc', ascending=False).groupby(['cluster', 'type']).head(5)
                    importance = [1]*x.shape[1]
                    for iippp in range(len(importance)):
                        for iv in sig_top_df['variable']:
                            if iv in x.columns.tolist()[iippp]:
                                importance[iippp] = random.sample(range(10,10000), 1)[0]
                else:
                    importance = random.sample(range(x.shape[1]), x.shape[1])
            imp_df = pd.DataFrame(importance, index = tmp_X.columns, columns = ['sig'])
            imp_df['varname'] = tmp_X.columns
            imp_df['importance_source'] = importance_source
            imp_df['outname'] = yoi+op
            rmse, cor, rele, cvn, impn, models = list(), list(), list(), list(), list(), list()
            for iimp in importance_number:
                select_x = x.iloc[:,sorted(range(len(importance)), key=lambda i: importance[i])[-iimp:]]
                s=select_x.shape[1]
                cnn_shape=[]
                print("Determing CNN shape")
                while math.floor(s) > 1:
                    s=math.floor(s*dnm)
                    if math.floor(s) > 1:
                        cnn_shape.append(s)
                ss = ShuffleSplit(n_splits=int(2*1/test_size), test_size=test_size, random_state=0)
                ss = ShuffleSplit(n_splits=10, test_size=test_size, random_state=0)
                cvi = 0
                for train_index, test_index in ss.split(range(x.shape[0])):
                    train_x = select_x.iloc[train_index]
                    train_y = y.iloc[train_index]
                    test_x = select_x.iloc[test_index]
                    test_y = y.iloc[test_index]
                    print("Building the CNN model")
                    model = Sequential()
                    model.add(Dense(select_x.shape[1], input_dim = select_x.shape[1], activation = input_act_func, kernel_initializer=gnrl_iniz))
                    if sqr_nlayer > 0:
                        for ii in range(sqr_nlayer):
                            model.add(Dense(math.floor(select_x.shape[1]*dnm), kernel_initializer=gnrl_iniz, activation = hidden_act_func))
                    for cs in cnn_shape:
                        model.add(Dense(cs, kernel_initializer=gnrl_iniz, activation = hidden_act_func)) # NOTE: check default activation function
                    model.add(Dense(1, kernel_initializer=gnrl_iniz, activation = input_act_func)) # NOTE: check default activation function is a(x) = x
                    opt = tf.keras.optimizers.Adam(learning_rate=ln_rt)
                    model.compile(loss = 'mean_squared_error', optimizer=opt)
                    fst = datetime.now()
                    record = model.fit(train_x, train_y, epochs = n_epoch, batch_size = 10, use_multiprocessing=True)
                    predict_y = model.predict(test_x, use_multiprocessing=True)
                    test_res = pd.DataFrame(test_y)
                    test_res = test_res.rename(columns={yoi+op: "measurement"})
                    test_res['prediction'] = predict_y
                    test_res['outcome'] = yoi+op
                    test_res['model'] = 'cnn_epch'+str(n_epoch)
                    test_res['cv_counter'] = cvi
                    if any(np.isnan(predict_y)):
                        rmse.append(np.nan)
                        cor.append(np.nan)
                    else:
                        correct_y = [i[0] for i in predict_y]
                        rmse.append(math.sqrt(mean_squared_error(test_y, predict_y)))
                        cor.append(scipy.stats.pearsonr(test_y, correct_y)[0])

                    cvn.append(cvi)
                    impn.append(iimp)
                    models.append('cnn_epch'+str(n_epoch))
                    if cvi == 0:
                        merge_res = test_res
                    else:
                        merge_res = pd.concat([merge_res, test_res], axis=0)
                    cvi += 1
                merge_res['importance_source'] = importance_source
                merge_res['importance_number'] = iimp
                merge_res.to_csv(rpf+"feat"+str(iimp)+"_manual_cv_predict_measure_df.csv")
            score_df = pd.DataFrame([rmse, cor, rele, cvn, impn, models], index=['rmse', 'r', 'mape', 'cvn', 'nfeats', 'model']).T
            score_df.to_csv(rpf+"_manual_cv_score_df.csv")
            if i == 0:
                merge_imp_df = imp_df
            else:
                merge_imp_df = pd.concat([merge_imp_df, imp_df], axis=0)
            i += 1
    merge_imp_df.to_csv(result_dir+"/"+result_folder+"_importance_df.csv")
    print("Time cost for "+str(icc)+": "+str(datetime.now()-st))
    print("+"*20)
    print("\n\n")
