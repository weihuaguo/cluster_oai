# Visualize the self_corr_input.py
# Weihua Guo, Ph.D.
# 12/13/2020
rm(list=ls())
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(circlize))

cat("Reading matrices...\n")
data_dir<-'/mnt/sda1/OAI_Data/data_summary'
clean_co<-'v25'
cor_df<-read.csv(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_correlation_", clean_co, ".csv", sep=''),row.names=1)
p_df<-read.csv(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_pvalue_",clean_co,".csv", sep=''),row.names=1)
num_df<-read.csv(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_na_num_",clean_co,".csv", sep=''),row.names=1)
cat("Replace NA's with abnormal values for clustering...\n")

cor_df[is.na(cor_df)] <- -2.0
cor_df<-do.call(data.frame, lapply(cor_df, function(x) replace(x, is.infinite(x), 2.0)))
print(cor_df[1:9,1:6])
cor_color<-colorRamp2(c(-2, -1, 0, 1), c("dodgerblue4", "dodgerblue",  "white", "firebrick"))
p_color<-colorRamp2(c(0,0.01,0.05,0.10,1.0), c("gold4", "gold3", "gold", "white", "purple4"))
num_color<-colorRamp2(c(0,1500,3000,4500), c("white", "darkolivegreen1", "darkolivegreen3", "darkolivegreen"))
cor_hm <- Heatmap(as.matrix(cor_df), name = 'Correlation', 
		  show_column_names=FALSE, 
		  show_row_names=FALSE, 
#		  row_km=5, row_km_repeats=50,
#		  column_km=5, column_km_repeats=50,
		  row_dend_width = unit(2, "cm"),
		  column_dend_height = unit(2, "cm"),
		  column_title='Correlations',
		  col=cor_color)
p_hm <- Heatmap(as.matrix(p_df), name = 'P-value',
		col=p_color,
		cluster_rows=FALSE, 
		cluster_columns=FALSE, 
		show_column_names=FALSE, 
		show_row_names=FALSE,
		row_order=row_order(cor_hm), 
		column_order=column_order(cor_hm), 
		column_title='P-values',
		heatmap_legend_param = list(at = c(0.0, 0.01, 0.05, 0.10, 1.0),
					    legend_height = unit(3, "cm")))
num_hm <- Heatmap(as.matrix(num_df), name = '#Pairs', 
		  col=num_color,
		  cluster_rows=FALSE, 
		  cluster_columns=FALSE,
		  show_column_names=FALSE, 
		  show_row_names=FALSE,
		  column_title='Number of available pairs',
		  row_order=row_order(cor_hm), 
		  column_order=column_order(cor_hm))
hm_list<-cor_hm+p_hm+num_hm
png(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_correlation_vis_",clean_co,".png", sep = ""),
width=20, height=7, res=300, units='in')
print(draw(hm_list))
gar<-dev.off()
# order_df<-as.data.frame(matrix(ncol=1,nrow=nrow(cor_df)))
order_df<-as.data.frame(row_order(cor_hm))
rownames(order_df)<-colnames(cor_df)
colnames(order_df)<-c("hclust_order")
print(head(order_df))

write.csv(order_df, paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_correlation_order_",clean_co,".csv", sep = ""))
