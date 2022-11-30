#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
import seaborn as sn
import numpy as np
import umap
import math
from datetime import datetime as dt
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, ShuffleSplit
from sklearn.metrics import silhouette_score, davies_bouldin_score,v_measure_score
from scipy import stats
from scipy.stats import kruskal
import matplotlib.pyplot as plt

st = dt.now()
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
#clst_dir = mainDir + "/cv_kmean_cluster_06112021"
clst_dir = mainDir + "/cluster_kmean_cv"
exp_id = "12112020_data_clean"
clean_co = "v25"
pca_umap = True
cluster_num = 4
general_prefix = clst_dir+'/clean_'+clean_co+'_'

enroll_txt = "Enrollees.txt"
ann_xlsx = "AllClinical00_V5_column_annotation.xlsx"

all_col_oi = pd.read_excel(mainDir+'/'+ann_xlsx, sheet_name='Baseline', index_col='Variables')
print(all_col_oi)
all_col_oi.index = all_col_oi.index.str.upper()

all_oh_merge_df=pd.read_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols.pkl')

kf = KFold(n_splits=10)
ss = ShuffleSplit(n_splits=10, test_size=0.10, random_state=0)
ss_cts=0
for train, test in ss.split(all_oh_merge_df.index.values):
    oh_merge_df = all_oh_merge_df.iloc[train]
    plot_prefix=general_prefix+'ss'+str(ss_cts)+'_'
    oh_merge_df.to_csv(plot_prefix+'merge_dataframe.csv')
    ist = dt.now()
    print("Start to imputate the NaN...")
    imputer = KNNImputer(n_neighbors=2, weights="uniform")
    imp_oh_merge_data = imputer.fit_transform(oh_merge_df)
    print("Imputation cost: "+str(dt.now()-ist))
    
    print("Start to scale the data...")
    scaled_oh_data = StandardScaler().fit_transform(imp_oh_merge_data)
    
    ust = dt.now()
    print("Start to run UMAP...")
    pca = PCA(n_components=30)
    pca_pcs = pca.fit_transform(scaled_oh_data)
    pc_elbow_df = pd.DataFrame({'var':pca.explained_variance_ratio_, 'PC': range(30)})
    sn.scatterplot(data = pc_elbow_df, x = "PC", y = "var", linewidth=0)
    plt.savefig(plot_prefix+'pca_check_elbow_plot.png', dpi=300)
    plt.clf()
    plt.close()
    
    umap_red = umap.UMAP()
    umap_emb = umap_red.fit_transform(pca_pcs[:,:16])
    print("UMAP cost: "+str(dt.now()-ust))

    kmst = dt.now()
    cls_col = "kmean_pca"
    print("Start to clustering with KMeans with PCA...")
    pca = PCA(n_components=16).fit(scaled_oh_data)
    pca.fit(scaled_oh_data)
    pca_score = pca.transform(scaled_oh_data)
    estimator = KMeans(init='k-means++', n_clusters=cluster_num)
    cluster_est = estimator.fit(pca_score)
    cluster_res = cluster_est.predict(pca_score)
    print("KMeans with PCA initial cost: "+str(dt.now()-kmst))
    
    cls_rd_res = pd.DataFrame(umap_emb, columns = ["UMAP1", "UMAP2"], index = oh_merge_df.index)
    cls_rd_res['kmean_pca'] = cluster_res
    cls_rd_res.to_excel(plot_prefix+'cluster'+str(cluster_num)+'_kmean_pca_umap_res.xlsx')
    cls_rd_res['kmean_pca'] = cls_rd_res['kmean_pca'].astype("str")
    
    sn.scatterplot(data = cls_rd_res, x = "UMAP1", y = "UMAP2", hue = "kmean_pca", linewidth=0, s=0.2)
    plt.savefig(plot_prefix+'cluster'+str(cluster_num)+'_kmean_pca_umap.png', dpi=300)
    plt.clf()
    plt.close()
    ss_cts+=1
    
