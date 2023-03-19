# Visualize the results of predict_cluster_classic.py
# Weihua Guo, Ph.D.
# 2023-01-09

rm(list = ls())
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(rstatix))

mainDir <- "/mnt/sda1/OAI_Data"
clst_dir <- paste(mainDir, "kmean_cluster_12252020", sep = "/")

result_folder <- "direct_predict_cluster_230109"
score_type <- c("roc_auc_ovo", "roc_auc_ovr", "accuracy", "roc_auc_ovo_train", "roc_auc_ovr_train", "accuracy_train")
model_name <- c("lr", "rf10tree", "rf20tree", "rf40tree", "rf60tree", "rf80tree", "rf100tree", "svclinear", "svcpoly", "svcrbf", "svcsigmoid")

merge_score_flag <- FALSE

rpf <- paste(clst_dir, "/", result_folder, "/", result_folder, "_withcnn_", sep = "")
rpf <- paste(clst_dir, "/", result_folder, "/", result_folder, "_", sep = "")


if (merge_score_flag) {
	c <- 0
	for (imn in model_name) {
		cat(imn, "\n")
		tmp_pattern <- paste(result_folder, imn, sep = "_")
		all_files <- list.files(paste(clst_dir, result_folder, sep = "/"), pattern = tmp_pattern)
		for (iaf in all_files) {
			cat("\t", iaf, "\n")
			tmp_df <- read.csv(paste(clst_dir, result_folder, iaf, sep = "/"))
			tmp_score_df <- tmp_df[!duplicated(tmp_df$cv_counter),]
			tmp_gath_score <- gather(tmp_score_df[,c(score_type, "model", "importance_number", "cv_counter")], "score_type", "score", score_type)
			if (c == 0) {
				merge_gath_score <- tmp_gath_score
			} else {
				merge_gath_score <- rbind(merge_gath_score, tmp_gath_score)
			}
			c <- c + 1
		}
	}
	write.csv(merge_gath_score, paste(rpf, "merge_score_dataframe.csv", sep = ""))
} else {
	merge_gath_score <- read.csv(paste(rpf, "merge_score_dataframe.csv", sep = ""), row.names = 1)
}
#print(head(merge_gath_score))

cat("Visualization...\n")
merge_stat_score <- merge_gath_score %>%
	group_by(score_type, importance_number, model) %>%
	summarise(n = n(),
		  mean = mean(score, na.rm = T),
		  median = median(score, na.rm = T),
		  sd = sd(score, na.rm = T),
		  min = min(score, na.rm = T),
		  max = max(score, na.rm = T),
		  se = sd/sqrt(n),
		  mean_se_upper = mean+se,
		  mean_se_lower = mean-se
	)

val_gath_score <- merge_gath_score[!str_detect(merge_gath_score$score_type, "_train"),]
zoomin_val_gath_score <- val_gath_score[val_gath_score$importance_number >= 5 & val_gath_score$importance_number <= 15,]

val_gath_score$importance_number <- as.factor(val_gath_score$importance_number)
val_gath_score$model <- factor(val_gath_score$model, levels = model_name)
val_box_gg <- ggplot(val_gath_score, aes(x = importance_number, y = score, color = model)) +
	geom_boxplot(position = "dodge") +
	geom_point(position = position_dodge(width = 0.75), size = 0.1) +
	facet_wrap(~score_type, ncol = 1) +
	theme_bw()
ggsave(paste(rpf, "validation_boxplot.png", sep = ""), dpi = 300, width = 48, height = 16)

zoomin_val_gath_score$importance_number <- as.factor(zoomin_val_gath_score$importance_number)
zoomin_val_gath_score$model <- factor(zoomin_val_gath_score$model, levels = model_name)
val_box_gg <- ggplot(zoomin_val_gath_score, aes(x = importance_number, y = score, color = model)) +
	geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
	geom_point(position = position_dodge(width = 0.75), size = 0.1) +
	facet_wrap(~score_type, ncol = 1, scales = "free") +
	labs(x = "Input variable number", y = "Prediction metric (Mean/SEM)", color = "Predictive model") +
	theme_bw()
ggsave(paste(rpf, "zoomin_validation_boxplot.png", sep = ""), dpi = 300, width = 12, height = 6)

val_stat_score <- merge_stat_score[!str_detect(merge_stat_score$score_type, "_train"),]
zoomin_val_stat_score <- val_stat_score[val_stat_score$importance_number >= 5 & val_stat_score$importance_number <= 15,]

val_stat_score$importance_number <- as.factor(val_stat_score$importance_number)
val_stat_score$model <- factor(val_stat_score$model, levels = model_name)
val_bar_gg <- ggplot(val_gath_score, aes(x = importance_number, y = score)) +
	geom_jitter(aes(color = model), position=position_dodge(width = 0.6), size = 0.5) +
	geom_bar(data = val_stat_score, position=position_dodge(width = 0.6, preserve = "single"), stat="identity", 
		 aes(color = model, fill = model, x = importance_number, y = mean), alpha = 0.1, width = 0.5) +
	geom_errorbar(data = val_stat_score, aes(ymin = mean_se_lower, ymax = mean_se_upper, color = model, x = importance_number, y = mean), 
		      width = 0.2, position = position_dodge(width = 0.6)) + 
	facet_wrap(~ score_type, ncol = 1) +
	labs(x = "Input variable number", y = "Prediction metric (Mean/SEM)", color = "Predictive model", fill = "Predictive model") +
	stat_compare_means(aes(group = model), label = "p.signif") +
	theme_bw()
ggsave(paste(rpf, "validation_barplot.png", sep = ""), dpi = 300, width = 48, height = 16)

val_dotline_gg <- ggplot(val_stat_score, aes(x = importance_number, y = mean, color = model, group = model)) +
	geom_errorbar(aes(ymin = mean_se_lower, ymax = mean_se_upper), width = .1) +
	geom_line() +
	geom_point(size = 1) +
	facet_wrap(~score_type, ncol = 1) +
	labs(x = "Input variable number", y = "Prediction metric (Mean/SEM)", color = "Predictive model", fill = "Predictive model") +
	theme_bw()
ggsave(paste(rpf, "validation_dotline.png", sep = ""), dpi = 300, width = 12, height = 9)

zoomin_val_stat_score$importance_number <- as.factor(zoomin_val_stat_score$importance_number)
zoomin_val_stat_score$model <- factor(zoomin_val_stat_score$model, levels = model_name)
val_bar_gg <- ggplot(zoomin_val_gath_score, aes(x = importance_number, y = score)) +
	geom_jitter(aes(color = model), position=position_dodge(width = 0.6), size = 0.5) +
	geom_bar(data = zoomin_val_stat_score, position=position_dodge(width = 0.6, preserve = "single"), stat="identity", 
		 aes(color = model, fill = model, x = importance_number, y = mean), alpha = 0.1, width = 0.3) +
	geom_errorbar(data = zoomin_val_stat_score, aes(ymin = mean_se_lower, ymax = mean_se_upper, color = model, x = importance_number, y = mean), 
		      width = 0.2, position = position_dodge(width = 0.6)) + 
	facet_wrap(~ score_type, ncol = 1, scales = "free") +
	labs(x = "Input variable number", y = "Prediction metric (Mean/SEM)", color = "Predictive model", fill = "Predictive model") +
	stat_compare_means(aes(group = model), label = "p.signif") +
	theme_bw()
ggsave(paste(rpf, "zoomin_validation_barplot.png", sep = ""), dpi = 300, width = 9, height = 6)

val_dotline_gg <- ggplot(zoomin_val_stat_score, aes(x = importance_number, y = mean, color = model, group = model)) +
	geom_errorbar(aes(ymin = mean_se_lower, ymax = mean_se_upper), width = .1) +
	geom_line() +
	geom_point(size = 1) +
	facet_wrap(~score_type, ncol = 1, scales = "free") +
	labs(x = "Input variable number", y = "Prediction metric (Mean/SEM)", color = "Predictive model", fill = "Predictive model") +
	theme_bw()
ggsave(paste(rpf, "zoomin_validation_dotline.png", sep = ""), dpi = 300, width = 6, height = 9)
