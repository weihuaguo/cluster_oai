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
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))


clst_dir <- "/mnt/sda1/OAI_Data/kmean_cluster_12252020/"
data_dir <- "/mnt/sda1/OAI_Data/data_summary/"

exp_id <- "kmean_pca_umap"
input_id <- "12112020_data_clean"
clean_co <- 'v25'
input_prefix <- paste(clst_dir, 'clean_', clean_co, '_', sep = "")
output_prefix <- paste(input_prefix, "final_", sep = "")
score_df <- as.data.frame(read_excel(paste(clst_dir, "kmeans_metric_pca_score_direct_knn2imp_", input_id, "241016_updated.xlsx", sep = "")))
colnames(score_df)[1] <- "cluster"
score_df$c <- as.numeric(str_split_fixed(score_df$cluster, "_C", n = 2)[,2])
print(head(score_df))
print(unique(score_df$cluster))
print(colnames(score_df))

gath_df <- gather(score_df, "metric", "value", colnames(score_df)[3:(ncol(score_df)-1)])
print(head(gath_df))

gg <- ggplot(gath_df, aes(x = c, y = value)) +
	geom_point() +
	geom_vline(xintercept = 4, linetype = "dotdash") +
	facet_wrap(~metric, nrow = 3, scales = "free") +
	labs(x = "Cluster", y = "Metric values", title = "Cluster quality") +
	theme_bw()
ggsave(paste(clst_dir, "kmeans_metric_pca_score_direct_knn2imp_", input_id, "241016_updated.png", sep = ""), dpi = 300, width = 9, height = 6)
