#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 13:31:13 2024
Subcluster the lack-of-supplements cluster

@author: weihua
"""

def clst_metric(df, df_name, clst_nums = range(2,21)):
    cmst = dt.now()
    print("Screening for "+df_name+" dataframe...")
    score_dict = {}
    for i in clst_nums:
        tmp_index = df_name+'_C'+str(i)
        score_dict[tmp_index] = {}
        km = KMeans(n_clusters=i, init='k-means++').fit(df)
        preds = km.predict(df)
        print("Score for number of cluster(s) {}: {}".format(i,km.score(df)))
        score_dict[tmp_index]['km_scores'] = -km.score(df)
        
        silhouette = silhouette_score(df, preds)
        score_dict[tmp_index]['silhouette_score'] = silhouette
        print("Silhouette score for number of cluster(s) {}: {}".format(i,silhouette))
        
        db = davies_bouldin_score(df, preds)
        score_dict[tmp_index]['davies_bouldin_score'] = db
        print("Davies Bouldin score for number of cluster(s) {}: {}".format(i,db))

        print("-+"*45)
    score_df = pd.DataFrame.from_dict(score_dict).T
    print("Total time cost "+str(dt.now()-cmst))
    return score_df

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
from sklearn.metrics import silhouette_score, davies_bouldin_score,v_measure_score
from scipy import stats
from scipy.stats import kruskal
import matplotlib.pyplot as plt

st = dt.now()
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
clst_dir = mainDir + "/kmean_cluster_12252020"
sub_clst_dir = mainDir + "/kmean_subcluster_07072024"
exp_id = "12112020_data_clean"
clean_co = "v25"
pca_umap = True
cluster_num = 6
plot_prefix = clst_dir+'/clean_'+clean_co+'_'
sub_prefix = sub_clst_dir+'/clean_'+clean_co+'_los_subcluster_'

enroll_txt = "Enrollees.txt"
ann_xlsx = "AllClinical00_V5_column_annotation.xlsx"

method = "onehot_plot"

all_col_oi = pd.read_excel(mainDir+'/'+ann_xlsx, sheet_name='Baseline', index_col='Variables')
print(all_col_oi)
all_col_oi.index = all_col_oi.index.str.upper()

oh_merge_df=pd.read_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols.pkl')
overall_clst_rslt=pd.read_excel(plot_prefix+'cluster4_kmean_pca_umap_res.xlsx')
sub_clst_rslt=overall_clst_rslt.loc[overall_clst_rslt['kmean_pca'] == 0]

ist = dt.now()
print("Start to imputate the NaN...")
imputer = KNNImputer(n_neighbors=2, weights="uniform")
imp_oh_merge_data = imputer.fit_transform(oh_merge_df)

sub_imp_oh_merge_data = imp_oh_merge_data[overall_clst_rslt['kmean_pca'] == 0,:]
sub_imp_oh_merge_df = pd.DataFrame(sub_imp_oh_merge_data, columns = oh_merge_df.columns, index=sub_clst_rslt['ID'])
sub_imp_oh_merge_df.to_pickle(sub_prefix+'onehotspot_merge_dataframe_without_na_cols_imp.pkl') # SD3
print("Imputation cost: "+str(dt.now()-ist))

print("Start to scale the data...")
scaled_oh_data = StandardScaler().fit_transform(sub_imp_oh_merge_data)

ust = dt.now()
print("Start to run UMAP...")
pca = PCA(n_components=30)
pca_pcs = pca.fit_transform(scaled_oh_data)
pc_elbow_df = pd.DataFrame({'var':pca.explained_variance_ratio_, 'PC': range(30)})
sn.scatterplot(data = pc_elbow_df, x = "PC", y = "var", linewidth=0)
plt.savefig(sub_prefix+'pca_check_elbow_plot.png', dpi=300) # Fig S2A
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

kmean_score_df = clst_metric(pca_score, df_name="direct_knn2imp")
kmean_score_df.to_excel(clst_dir+"/kmeans_metric_pca_score_direct_knn2imp_"+exp_id+".xlsx")

kmean_score_df['c'] = range(2,21)
kmean_score_df.set_index('c', inplace=True, drop = False)
kmean_score_df['km_scores'].plot.line()
plt.savefig(sub_prefix+'cluster'+str(cluster_num)+'_pca_score_km_score_elbow.png', dpi=300)
plt.clf()
plt.close()

kmean_score_df['silhouette_score'].plot.line()
plt.savefig(sub_prefix+'cluster'+str(cluster_num)+'_pca_score_silhouette_score_elbow.png', dpi=300)
plt.clf()
plt.close()

kmean_score_df['davies_bouldin_score'].plot.line()
plt.savefig(sub_prefix+'cluster'+str(cluster_num)+'_pca_score_db_score_elbow.png', dpi=300)
plt.clf()
plt.close()

estimator = KMeans(init='k-means++', n_clusters=cluster_num)
cluster_est = estimator.fit(pca_score)
cluster_res = cluster_est.predict(pca_score)
print("KMeans with PCA initial cost: "+str(dt.now()-kmst))

cls_rd_res = pd.DataFrame(umap_emb, columns = ["UMAP1", "UMAP2"], index = sub_clst_rslt['ID'])
cls_rd_res['kmean_pca'] = cluster_res
cls_rd_res.to_excel(sub_prefix+'cluster'+str(cluster_num)+'_kmean_pca_umap_res.xlsx')
cls_rd_res['kmean_pca'] = cls_rd_res['kmean_pca'].astype("str")

sn.scatterplot(data = cls_rd_res, x = "UMAP1", y = "UMAP2", hue = "kmean_pca", linewidth=0, s=0.2)
plt.savefig(sub_prefix+'cluster'+str(cluster_num)+'_kmean_pca_umap.png', dpi=300)
plt.clf()
plt.close()
