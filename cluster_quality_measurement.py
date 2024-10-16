#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 10:44:18 2024
Calculate various scores to measure the cluster qualities

@author: weihua
"""
import pandas as pd
import seaborn as sn
import numpy as np
import umap
from datetime import datetime as dt
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from scipy.stats import pointbiserialr
import matplotlib.pyplot as plt

st = dt.now()
mainDir = "/mnt/sda1/OAI_Data"
dataDir = mainDir + "/ASCII"
resDir = mainDir + "/data_summary"
clst_dir = mainDir + "/kmean_cluster_12252020"
exp_id = "12112020_data_clean"
clean_co = "v25"
pca_umap = True
cluster_num = 4
plot_prefix = clst_dir+'/clean_'+clean_co+'_'

enroll_txt = "Enrollees.txt"
ann_xlsx = "AllClinical00_V5_column_annotation.xlsx"

method = "onehot_plot"

all_col_oi = pd.read_excel(mainDir+'/'+ann_xlsx, sheet_name='Baseline', index_col='Variables')
print(all_col_oi)
all_col_oi.index = all_col_oi.index.str.upper()

oh_merge_df=pd.read_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols.pkl')

ist = dt.now()
print("Start to imputate the NaN...")
imputer = KNNImputer(n_neighbors=2, weights="uniform")
imp_oh_merge_data = imputer.fit_transform(oh_merge_df)
imp_oh_merge_df = pd.DataFrame(imp_oh_merge_data, columns = oh_merge_df.columns, index=oh_merge_df.index)
imp_oh_merge_df.to_pickle(resDir+'/onehotspot_merge_dataframe_'+exp_id+'_'+clean_co+'_without_na_cols_imp.pkl') # SD3
print("Imputation cost: "+str(dt.now()-ist))

print("Start to scale the data...")
scaled_oh_data = StandardScaler().fit_transform(imp_oh_merge_data)

ust = dt.now()
print("Start to run UMAP...")
pca = PCA(n_components=30)
pca_pcs = pca.fit_transform(scaled_oh_data)
pc_elbow_df = pd.DataFrame({'var':pca.explained_variance_ratio_, 'PC': range(30)})
sn.scatterplot(data = pc_elbow_df, x = "PC", y = "var", linewidth=0)
plt.savefig(plot_prefix+'pca_check_elbow_plot.png', dpi=300) # Fig S2A
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

cmst = dt.now()
df_name="direct_knn2imp"
df = pca_score
clst_nums = range(2,21)
print("Screening for "+df_name+" dataframe...")
score_dict = {}
for i in clst_nums:
    tmp_index = df_name+'_C'+str(i)
    score_dict[tmp_index] = {}
    km = KMeans(n_clusters=i, init='k-means++').fit(df)

    preds = km.predict(df)
    print("Score for number of cluster(s) {}: {}".format(i,km.score(df)))
    score_dict[tmp_index]['km_scores'] = -km.score(df)
    score_dict[tmp_index]['km_inertia'] = km.inertia_
    
    # Get unique cluster labels
    unique_labels = np.unique(km.labels_)
    # Calculate inter-cluster distances
    inter_cluster_distances = []
    for m in range(len(unique_labels)):
        for n in range(m + 1, len(unique_labels)):
            cluster_i = df[km.labels_ == unique_labels[m]]
            cluster_j = df[km.labels_ == unique_labels[n]]
            inter_cluster_distances.append(np.min(np.linalg.norm(cluster_i[:, None] - cluster_j, axis=2)))

    intra_cluster_distances = []
    for label in unique_labels:
        cluster = df[km.labels_ == label]
        intra_cluster_distances.append(np.max(np.linalg.norm(cluster[:, None] - cluster, axis=2)))

    # Calculate Generalized Dunn index
    score_dict[tmp_index]['generalized_dunn'] = np.min(inter_cluster_distances) / np.max(intra_cluster_distances)
    score_dict[tmp_index]['baker_hubert'] = np.max(inter_cluster_distances) / np.min(intra_cluster_distances)
    score_dict[tmp_index]['intra_cluster_diff'] = np.mean(intra_cluster_distances).item()
    score_dict[tmp_index]['inter_cluster_diff'] = np.mean(inter_cluster_distances).item()

    km_centroid = km.cluster_centers_
    pbr_list = []
    for ik in range(i):
        tmp_dist_list = []
        for ir in range(len(pca_score)):
            tmp_dist = np.linalg.norm(pca_score[ir]-km_centroid[ik]).item() # Distance between each point to each cluster centroid
            tmp_dist_list.append(tmp_dist)
        tmp_pbr = pointbiserialr(km.labels_==ik, tmp_dist_list)
        pbr_list.append(tmp_pbr[0].item())
    score_dict[tmp_index]['pointbiserial'] = np.mean(pbr_list).item() # Note: Use the average PBr as the score
    
    ch_index = calinski_harabasz_score(df, preds).item()
    score_dict[tmp_index]['calinski_harabsz_score'] = ch_index
    print("Calinski-Harabasz index for number of cluster(s) {}: {}".format(i,ch_index))

    silhouette = silhouette_score(df, preds)
    score_dict[tmp_index]['silhouette_score'] = silhouette
    print("Silhouette score for number of cluster(s) {}: {}".format(i,silhouette))
    
    db = davies_bouldin_score(df, preds)
    score_dict[tmp_index]['davies_bouldin_score'] = db
    print("Davies Bouldin score for number of cluster(s) {}: {}".format(i,db))

    print("-+"*45)
score_df = pd.DataFrame.from_dict(score_dict).T
print("Total time cost "+str(dt.now()-cmst))

score_df.to_excel(clst_dir+"/kmeans_metric_pca_score_direct_knn2imp_"+exp_id+"241016_updated.xlsx")