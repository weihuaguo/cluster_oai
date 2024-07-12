# Visualize UMAP and cluster markers through R
# Weihua Guo, Ph.D.
# 07/02/2020

rm(list = ls())

cat("Loading packages...\n")
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(viridis))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(viridis))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(gridExtra))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(scales))
suppressMessages(library(rstatix))
suppressMessages(library(broom))
suppressMessages(library(boot))
suppressMessages(library(table1))
suppressMessages(library(colorspace))
suppressMessages(library(ggforce))


clst_dir <- "/mnt/sda1/OAI_Data/kmean_subcluster_07072024/"
data_dir <- "/mnt/sda1/OAI_Data/data_summary/"

exp_id <- "kmean_pca_umap"
input_id <- "12112020_data_clean"
clean_co <- 'v25'
input_prefix <- paste(clst_dir, 'clean_', clean_co, '_los_subcluster_', sep = "")
output_prefix <- paste(input_prefix, "final_", sep = "")
cluster_num <- 3
scale_top <- 4
kr_type <- "total_lastfollowup"

png_res <- 600
top_m <- 10

demographic_flag <- FALSE
outcome_table_flag <- FALSE
cluster_annotation_flag <- TRUE
one2one_annotation_flag <- FALSE
numeric_vis_marker_flag <- FALSE
categorical_vis_marker_flag <- FALSE
umap_vis_flag <- FALSE
violin_vis_flag <- FALSE
volcano_flag <- FALSE
supplement_flag <- FALSE

cat("Reading python output results...\n")
umap_df <- as.data.frame(read_excel(paste(input_prefix, "cluster", cluster_num, "_kmean_pca_umap_res.xlsx", sep = "")))
rownames(umap_df) <- umap_df$ID
umap_df$name <- str_c("C", umap_df$kmean_pca)

data_df <- as.data.frame(read.csv(paste(data_dir, "input_direct_merge_dataframe_",input_id,"_sbj_clean_",clean_co,".csv", sep = ""), 
				  sep = ",", header = TRUE, row.names=1, stringsAsFactors = FALSE))
data_df <- data_df[rownames(umap_df),]
var_types <- sapply(data_df, class)

prog_files <- list.files(data_dir, pattern = 'survival_ready_results.csv')
var_coding <- as.data.frame(read_excel(paste(data_dir, "Variable_coding_v2.xlsx", sep = "")))

kr_df <- read.csv(paste(data_dir, kr_type, "_merge_patient_basic_outcome_information.csv", sep = ""), header = T)
rownames(kr_df) <- kr_df$ID

##### UMAP overlay with clusters [Fig 2A]
umap_df$cluster <- str_c("C", umap_df$kmean_pca)
umap_gg <- ggplot(umap_df, aes_string(x = "UMAP1", y = "UMAP2", color = "name")) + 
	geom_point(size = 1) +
	scale_color_brewer(palette = "Spectral") +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	labs(color='Cluster') +
	theme_classic()
ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_umap_r.png', sep = ""), umap_gg, 
       dpi = png_res, width = 7.5, height = 5)

if (FALSE) {
#metric_df <- as.data.frame(read_excel(paste(input_prefix, 'kmeans_metric_result_direct_knn2imp_', input_id, '.xlsx', sep = "")))
##### UMAP overlay with clusters [Fig S2A]
cat("Visualize metrics for clustering...\n")
metric_df$cluster_num <- as.numeric(str_split_fixed(metric_df[,1], '_C', n = 2)[,2])
gath_metric_df <- metric_df %>% gather(metric, scores, colnames(metric_df)[2:4])
metric_gg <- ggplot(gath_metric_df, aes(x=cluster_num, y = scores)) +
	geom_point() +
	facet_wrap(~metric, ncol=1, scales='free_y') +
	theme_bw()
ggsave(paste(output_prefix, 'kmeans_metric_result_direct_knn2imp_', input_id, '_vis.png', sep = ""), metric_gg,
dpi = png_res, width = 6, height=6)
}

###### Generate the marker table [Table S2]
if (cluster_annotation_flag) {
	cat("Generate marker table...\n")
	marker_types <- sapply(data_df, class)
	marker_cols <- c('variable', 'cluster', 'type', 'test','avg1', 'avg2', 'diff', 'logfc', 'per1', 'per2', 'pval')
	marker_df <- as.data.frame(matrix(nrow=(cluster_num*ncol(data_df)), ncol=length(marker_cols)))
	colnames(marker_df) <- marker_cols
	ic <- 1
	for (i in 0:(cluster_num-1)) {
		cat("\tCluster", i, '\n')
		tmp_umap_df <- umap_df
		tmp_umap_df$comp <- ifelse(tmp_umap_df$kmean_pca == i, 'COI', 'Others')
		for (im in names(marker_types)) {
	#		print(marker_types[[im]])
			marker_df[ic, 'variable'] <- im
			marker_df[ic, 'cluster'] <- i
			marker_df[ic, 'type'] <- marker_types[[im]]
			tmp_df <- cbind(tmp_umap_df, data_df[,im])
			tmp_force_num <- as.numeric(as.character(tmp_df[,ncol(tmp_df)]))
			if (sum(!is.na(tmp_force_num))>0) {
				marker_types[[im]] <- 'numeric'
				tmp_df[,ncol(tmp_df)] <- tmp_force_num
			}

			if (marker_types[[im]]=='numeric'|marker_types[[im]]=='integer'){
	#			cat("\t\tKruskal-Wallis Test...\n")
				marker_df[ic, 'test'] <- 'Kruskal-Wallis test' #NOTE: 5 clusters doesn't work
	#			print(head(tmp_df))
				tmp_sum <- tmp_df %>%
					group_by(comp) %>%
					summarize(mean = mean(`data_df[, im]`, na.rm=T),
						  n=n(),
						  naper = sum(!is.na(`data_df[, im]`))/n())
				marker_df[ic, 'avg1'] <- tmp_sum[1,'mean']
				marker_df[ic, 'avg2'] <- tmp_sum[2, 'mean']
				marker_df[ic, 'diff'] <- tmp_sum[1, 'mean']-tmp_sum[2,'mean']
				marker_df[ic, 'logfc'] <- log2(tmp_sum[1, 'mean']/tmp_sum[2, 'mean'])
				marker_df[ic, 'per1'] <- tmp_sum[1, 'naper']
				marker_df[ic, 'per2'] <- tmp_sum[2, 'naper']
				tmp_test <- kruskal.test(`data_df[, im]`~comp, data = tmp_df)
				marker_df[ic, 'pval'] <- tmp_test$p.value
			} else {
	#			cat("\t\tFisher exact test...\n")
				tmp_table <- as.data.frame.matrix(table(tmp_df$comp, tmp_df[,ncol(tmp_df)]))
				tmp_prop_table <- prop.table(as.matrix(tmp_table),margin=1)
	#			print(im)
	#			print(tmp_prop_table)
	#			print(dim(tmp_table))
				print(tmp_table)
				if (sum(tmp_prop_table[1,]==max(tmp_prop_table[1,]))>1|sum(tmp_prop_table[2,]==max(tmp_prop_table[2,]))>1) {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])][1]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])][1]
				} else {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])]
				}
				cat(im, '\n')
				if (ncol(tmp_table) >= 2) {
				tmp_test <- fisher.test(tmp_table, simulate.p.value=TRUE, B=1e7)#, workspace=2e9)
				tmp_chitest <- chisq.test(tmp_table, correct=FALSE)
				marker_df[ic, 'logfc'] <- sqrt(tmp_chitest$statistic/sum(tmp_table))
				marker_df[ic, 'test'] <- 'Fisher\'s exact test'
				}
				na_mask <- str_detect(colnames(tmp_table), 'Missing')
				if (sum(na_mask) == 1) {
					marker_df[ic, 'per1'] <- 1-tmp_table['COI', na_mask]/sum(tmp_table['COI',])
					marker_df[ic, 'per2'] <- 1-tmp_table['Others', na_mask]/sum(tmp_table['Others',])
				}
				marker_df[ic, 'pval'] <- tmp_test$p.value
			}
			ic <- ic + 1
		}
	}
	marker_df$adjp <- p.adjust(marker_df$pval, method = "BH")
	write.csv(marker_df, paste(output_prefix, "cluster", cluster_num,'_kmeans_direct_knn2imp_', input_id, '_marker_df_largeB.csv', sep = "")) # SD4
}

