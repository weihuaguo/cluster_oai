import glob
import pandas as pd

data_dir="/home/weihua/mnts/smb_plee/Group/weihua/data_summary"

exp_id = "12112020_data_clean"
slc_id = "predictor_vis_220125"
clean_co = "v25"

importance_source = "random" # "randomforest"/"linear"/"random"
importance_filter = "cluster"
importance_number = 25
test_id = [*range(0, 9389), *range(9725,10336)] # NOTE: total 10,000 tests!
result_pattern = "impt_" + importance_source +"_" + importance_filter + "_sl_cnn_imputed_221022_v"


all_folders = [result_pattern+"v"+str(i) for i in test_id]
#print(all_folders)

icc = 0
iss = 0

for iti in test_id:
    print(iti)
    iaf = result_pattern+"v"+str(iti)
    imp_df = pd.read_csv(data_dir+"/"+iaf+"/"+iaf+"_importance_df.csv", index_col=0)
    imp_df['test_id'] = "v"+str(iti)
    iy = 0
    for iuo in imp_df['outname'].unique().tolist():
        print(iuo)
        tmp_imp_df = imp_df.loc[imp_df['outname'] == iuo]
        importance = tmp_imp_df['sig'].to_numpy()
        tmp_imp_df = tmp_imp_df.iloc[sorted(range(len(importance)), key=lambda i: importance[i])[-25:],:]
        if iy == 0:
            select_imp_df = tmp_imp_df
        else:
            select_imp_df = pd.concat([select_imp_df, tmp_imp_df], axis=0)
        iy += 1
    if icc == 0:
        merge_imp_df = select_imp_df
    else:
        merge_imp_df = pd.concat([merge_imp_df, select_imp_df], axis=0)
    icc += 1
    all_sub_folder = glob.glob(data_dir+"/"+iaf+"/"+iaf+"_*/")
    for isf in all_sub_folder:
        print("\t"+isf)
        score_file = glob.glob(isf+"/*_score_df.csv")
        tmp_score_df = pd.read_csv(score_file[0], index_col=0)
        if iti < 8550:
            tmp_score_df['r'] = tmp_score_df['r'].str.replace("[", "")
            tmp_score_df['r'] = tmp_score_df['r'].str.replace("]", "")
        tmp_score_df['test_id'] = "v"+str(iti)
        tmp_score_df['outname'] = score_file[0].split("/")[-1].split("_")[0]

        if iss == 0:
            merge_score_df = tmp_score_df
        else:
            merge_score_df = pd.concat([merge_score_df, tmp_score_df], axis=0)
        iss += 1
merge_imp_df.to_csv(data_dir+"/merge_select_importance_df_test_max_" + str(max(test_id)) + "_230202.csv")
merge_score_df.to_csv(data_dir+"/merge_score_df_test_max_" + str(max(test_id)) + "_230202.csv") # SD12
