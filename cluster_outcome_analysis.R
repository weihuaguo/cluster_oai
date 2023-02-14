# Examine the differences of outcomes between different clusters
# Weihua Guo, Ph.D.
# 03/17/2022

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
suppressMessages(library(forestmodel))

#clst_dir <- "/mnt/sda1/OAI_Data/ASCII_Plot_05282020/cluster_07102020/"
#data_dir <- "/mnt/sda1/OAI_Data/ASCII_Plot_05282020/"
#clst_dir <- "/mnt/sda1/OAI_Data/kmean_cluster_01012021/"
clst_dir <- "/mnt/sda1/OAI_Data/kmean_cluster_12252020/"
data_dir <- "/mnt/sda1/OAI_Data/data_summary/"

exp_id <- "kmean_pca_umap"
input_id <- "12112020_data_clean"
surv_folder <- "survival_ready_220927" # DO NOT change
surv_type <-  "v00"#"true_continuous" # "v00", "semi_continuous"

clean_co <- 'v25'
input_prefix <- paste(clst_dir, 'clean_', clean_co, '_', sep = "")
cluster_num <- 4
scale_top <- 4

png_res <- 300
top_m <- 10
compare_flag <- FALSE
pairwise_flag <- FALSE
surv_flag <- TRUE

cat("Reading python output results...\n")
umap_df <- as.data.frame(read_excel(paste(input_prefix, "cluster", cluster_num, "_kmean_pca_umap_res.xlsx", sep = "")))
rownames(umap_df) <- umap_df$ID
umap_df$raw_cluster <- umap_df$kmean_pca
umap_df$Cluster <- umap_df$kmean_pca
umap_df$Cluster[umap_df$Cluster == 0] <- "Low supplemental vitamins"
umap_df$Cluster[umap_df$Cluster == 1] <- "Poor knee & general health"
umap_df$Cluster[umap_df$Cluster == 2] <- "Good knee & general health"
umap_df$Cluster[umap_df$Cluster == 3] <- "Intermediate knee & general health"

umap_df$Cluster <- factor(umap_df$Cluster, levels = c("Poor knee & general health", "Intermediate knee & general health", "Good knee & general health", "Low supplemental vitamins"))

data_df <- as.data.frame(read.csv(paste(data_dir, "input_direct_merge_dataframe_",input_id,"_sbj_clean_",clean_co,".csv", sep = ""), 
				  sep = ",", header = TRUE, row.names=1, stringsAsFactors = FALSE))
var_types <- sapply(data_df, class)
metric_df <- as.data.frame(read_excel(paste(input_prefix, 'kmeans_metric_result_direct_knn2imp_', input_id, '.xlsx', sep = "")))
# print(data_df[1:9, 1:6])
# print(head(umap_df))


data_files <- list.files(paste(data_dir, surv_folder, sep = ""), pattern = '_used_clean_dataframe_with_real_time.csv')
print(data_files)

if (compare_flag) {
	c <- 0
	for (idata in data_files) {
		cat(idata, "\n")
		out_name <- str_split_fixed(idata, "_", n = 2)[1]
		tmp_data <- read.csv(paste(data_dir, surv_folder, "/", idata, sep = ""), header = T)
		merge_data <- merge(tmp_data, umap_df, by = "ID", all.x = T)
	#	print(head(merge_data))
		tmp_pf <- paste(data_dir, surv_folder, "/", out_name, "_final_", sep = "")

		year_data <- merge_data[!is.na(merge_data$Cluster),]

		year_data$year_mod <- year_data$tm %% 12
		year_data <- year_data[year_data$year_mod <= 3 | year_data$year_mod >= 9,]
		year_data$tm_round <- ifelse(year_data$year_mod <= 3, year_data$tm-year_data$year_mod, 
					     year_data$tm+12-year_data$year_mod)
		year_data$tyr <- year_data$tm_round/12
		year_data$ytyr <- str_c("Y", year_data$tyr)
	#	print(head(year_data))

		
		year_sum <- year_data %>%
			group_by(Cluster, ytyr) %>%
			summarize(n = n()) %>%
			filter(n > 0) %>%
			group_by(ytyr) %>%
			summarize(n = n())
		keep_ys <- year_sum$ytyr[year_sum$n == length(unique(umap_df$Cluster))]
		year_sum$tyr <- as.numeric(str_sub(year_sum$ytyr, -1, -1))

		###### Fig S3x
		all_gg <- ggplot(year_data[year_data$ytyr %in% keep_ys,], aes(x = ytyr, y = value, color = Cluster)) +
			geom_boxplot(outlier.shape=NA) +
			geom_point(position = position_jitterdodge(jitter.width=0.3), size = 0.2, alpha = 0.6) +
			scale_color_brewer(palette = "Spectral") +
			stat_compare_means(aes(group = Cluster), label = "p.signif") +
			labs(x = "Time (year)", y = out_name, color = "Cluster") +
			theme_bw()
		ggsave(paste(tmp_pf, "boxplot_jitter.png", sep = ""), all_gg, dpi = png_res, width = 9, height = 4)

		stat_df <- year_data %>%
			filter(ytyr %in% keep_ys) %>%
			group_by(Cluster, ytyr) %>%
			summarize(n = n(),
				    mean = mean(value, na.rm = T),
				    median = median(value, na.rm = T),
				    sd = sd(value, na.rm = T),
				    min = min(value, na.rm = T),
				    max = max(value, na.rm = T),
				    se = sd/sqrt(n),
				    mean_se_upper = mean+se,
				    mean_se_lower = mean-se)
		
	#	print(head(stat_df))
		##### Fig S3x+1
		stat_df$tyr <- as.numeric(str_sub(stat_df$ytyr, -1, -1))
		y_name <- expression(paste(Mean%+-%SEM, sep = ""))
		time_gg <- ggplot(stat_df, aes(x = tyr, y = mean, color = Cluster, group = Cluster)) +
			geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
			geom_line() +
			geom_point(size = 1.8) +
			scale_color_brewer(palette = "Spectral") +
			labs(y = y_name, color = "Cluster", x = "Time (year)", title = out_name) +
			theme_bw()
		ggsave(paste(tmp_pf, "mean_se_dot_line.png", sep = ""), time_gg, dpi = png_res, width = 6, height = 4)


		if (str_detect(out_name, "KL")) {
			cts_df <- year_data[year_data$ytyr %in% keep_ys,] %>%
				group_by(value, Cluster, ytyr) %>%
				summarise(n = n())
			cts_df$value <- factor(cts_df$value)
			all_gg <- ggplot(cts_df, aes(x = Cluster, y = n, color = value, fill = value)) +
				geom_bar(position="stack", stat="identity") +
				scale_color_brewer(palette = "RdYlBu") +
				scale_fill_brewer(palette = "RdYlBu") +
				facet_grid(.~ytyr) +
				labs(x = "Cluster", y = "Number of knees", color = out_name, fill = out_name) +
				theme_bw()
			ggsave(paste(tmp_pf, "bar.png", sep = ""), all_gg, dpi = png_res, width = 9, height = 4)
			all_gg <- ggplot(cts_df, aes(x = Cluster, y = n, color = value, fill = value)) +
				geom_bar(position="fill", stat="identity") +
				scale_color_brewer(palette = "RdYlBu") +
				scale_fill_brewer(palette = "RdYlBu") +
				facet_grid(.~ytyr) +
				labs(x = "Cluster", y = "Proportion of knees", color = out_name, fill = out_name) +
				theme_bw()
			ggsave(paste(tmp_pf, "rel_bar.png", sep = ""), all_gg, dpi = png_res, width = 9, height = 4)

		}

		if (c == 0) {
			all_year_data <- year_data[year_data$ytyr %in% keep_ys,]
		} else {
			all_year_data <- rbind(all_year_data, year_data[year_data$ytyr %in% keep_ys,])
		}
		c <- c+1
	}
} 
# q(save = "no")

if (pairwise_flag) {#TOOOOO time consuming
	all_year_data$outcome_var <- str_sub(all_year_data$vid, 4, -1)
	c <- 0
	for (ic in unique(all_year_data$Cluster)) {
		cat(ic, "\n")
		all_year_data$group <- ifelse(all_year_data$Cluster == ic, ic, "others")
		all_year_stat <- all_year_data %>%
			group_by(ytyr, outcome_var) %>%
			wilcox_test(value ~ group, ref.group = ic, detailed=T) %>%
			add_significance("p") %>%
			adjust_pvalue() %>%
			add_significance("p.adj")
		if (c == 0) {
			merge_year_stat <- all_year_stat
		} else {
			merge_year_stat <- rbind(merge_year_stat, all_year_stat)
		}
		c <- c+1
	}

	print(head(all_year_data))
	print(head(merge_year_stat))
	write.csv(merge_year_stat, paste(data_dir, surv_folder, "/outcome_cluster_markers_wilcox_result.csv", sep = ""))


	padj_alpha <- c("ns" = 0.05, "*" = 0.5, "**" = 0.75, "***" = 0.9, "****" = 1.0)
	padj_size <- c("ns" = 1, "*" = 4, "**" = 5, "***" = 6, "****" = 7)
	# TODO change color scheme and color scale to make it clear
	dot_gg <- ggplot(merge_year_stat, aes(x = group1, y = ytyr)) +
		geom_point(aes(color = estimate, size = p.adj.signif, alpha = p.adj.signif)) +
		scale_color_distiller(palette = "RdYlBu") +
		facet_grid(.~outcome_var) +
		scale_size_manual(values = padj_size) +
		scale_alpha_manual(values = padj_alpha) +
		labs(x = "Cluster", y = "Time", color = "Difference", size = "Statistical significance\n(Adjusted p-value)", alpha = "Statistical significance\n(Adjusted p-value)") +
		theme_bw()
	ggsave(paste(data_dir, surv_folder,  "/outcome_cluster_marker_dotplot.png", sep = ""), dot_gg, dpi = png_res, width = 15, height = 5)
}# TOOOOOOOOOOOOO

if (surv_flag) {
	prog_files <- list.files(paste(data_dir, surv_folder, sep = ""), pattern = paste(surv_type, 'real_time_survival_dataframe.csv', sep = "_"))
	print(prog_files)

	for (iprog in prog_files) {
		cat(iprog, "\n")
		if (surv_type == "v00") {
			tmp_prog <- read.csv(paste(data_dir, surv_folder, "/", iprog, sep = ""), header = T)
		} else {
			tmp_prog <- read.csv(paste(data_dir, surv_folder, "/", iprog, sep = ""), header = T, row.names = 1)
		}
		out_name <- str_split_fixed(iprog, "_", n = 2)[1]
	#	print(dim(tmp_prog))
	#	print(out_name)
	#	q(save = "no")
		clean_prog <- tmp_prog
	#	print(dim(clean_prog))
		merge_prog <- merge(clean_prog, umap_df, by = "ID", all.x = T)
		merge_prog$dstart <- as.numeric(merge_prog$dstart)
		merge_prog$dstop <- as.numeric(merge_prog$dstop)

	#	print(dim(merge_prog))
	#	print(head(merge_prog))
		tmp_pf <- paste(data_dir, surv_folder, "/", out_name, "_", surv_type, "_", sep = "")
		if (str_detect(iprog, "JSW")) {
			ylim <- c(0.5, 1)
			p_day_pos <- c(10,0.6)
			p_mon_pos <- c(2,0.6)
		} else {
			ylim <- c(0.4, 1)
			p_day_pos <- c(10,0.5)
			p_mon_pos <- c(2,0.5)
		}


		tmp_lfit <- survfit(Surv((dstop - dstart), event) ~ Cluster, data = merge_prog, id = ID)
		tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = p_day_pos, ggtheme = theme_bw(), legend = "right", 
				      legend.labs = levels(merge_prog$Cluster), 
				      palette = "Spectral", ylim = ylim) + 
			labs(title = out_name, x = "Time (day)", color = "Cluster") 		
		png(paste(tmp_pf, 'day_kmplot.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(tmp_lgg)
		gar <- dev.off()
		tmp_cox <- coxph(Surv((dstop - dstart), event) ~ Cluster, data=merge_prog, id=ID, ties="breslow")
		write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'day_cox.csv', sep = ""))

		fprog <- merge_prog
		fprog$Cluster <- factor(fprog$Cluster, 
					levels = c("Good knee & general health", "Intermediate knee & general health", "Poor knee & general health", "Low supplemental vitamins"))

		png(paste(tmp_pf, 'day_forestmodel_hr.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(forest_model(coxph(Surv((dstop - dstart), event) ~ Cluster, data=fprog))+labs(title = out_name))
		gar <- dev.off()

		tmp_lfit <- survfit(Surv((mstop - mstart), event) ~ Cluster, data = merge_prog, id = ID)
		tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = p_mon_pos, ggtheme = theme_bw(), legend = "right", 
				      legend.labs = levels(merge_prog$Cluster), 
				      palette = "Spectral", ylim = ylim) + 
			labs(title = out_name, x = "Time (month)", color = "Cluster") 		
		png(paste(tmp_pf, 'month_kmplot.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(tmp_lgg)
		gar <- dev.off()
		tmp_cox <- coxph(Surv((mstop - mstart), event) ~ Cluster, data=merge_prog, id=ID, ties="breslow")
		write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'month_cox.csv', sep = ""))

		fprog <- merge_prog
		fprog$Cluster <- factor(fprog$Cluster, 
					levels = c("Good knee & general health", "Intermediate knee & general health", "Poor knee & general health", "Low supplemental vitamins"))

		png(paste(tmp_pf, 'month_forestmodel_hr.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(forest_model(coxph(Surv((mstop - mstart), event) ~ Cluster, data=fprog))+labs(title = out_name))
		gar <- dev.off()

		merge_prog$half_year_mod <- merge_prog$tm %% 6
		merge_prog$tm_half_round <- ifelse(merge_prog$half_year_mod <= 3, merge_prog$tm-merge_prog$half_year_mod, 
					      merge_prog$tm+6-merge_prog$half_year_mod)
		merge_prog$tyhr <- merge_prog$tm_half_round/12

		tmp_lfit <- survfit(Surv(tyhr, event) ~ Cluster, data = merge_prog, id = ID)
		tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = p_mon_pos, ggtheme = theme_bw(), legend = "right", 
				      legend.labs = levels(merge_prog$Cluster), 
				      palette = "Spectral", ylim = ylim) + 
			labs(title = out_name, x = "Time (1/2 year)", color = "Cluster") 		
		png(paste(tmp_pf, 'half_year_kmplot.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(tmp_lgg)
		gar <- dev.off()
		tmp_cox <- coxph(Surv(tyhr, event) ~ Cluster, data=merge_prog, id=ID, ties="breslow")
		write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'half_year_cox.csv', sep = ""))

		fprog <- merge_prog
		fprog$Cluster <- factor(fprog$Cluster, 
					levels = c("Good knee & general health", "Intermediate knee & general health", "Poor knee & general health", "Low supplemental vitamins"))
		png(paste(tmp_pf, 'half_year_forestmodel_hr.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(forest_model(coxph(Surv(tyhr, event) ~ Cluster, data=fprog))+labs(title = out_name))
		gar <- dev.off()



		year_prog <- merge_prog
		year_prog$year_mod <- merge_prog$tm %% 12
		year_prog <- year_prog[year_prog$year_mod <= 3 | year_prog$year_mod >= 9,]
		year_prog$tm_round <- ifelse(year_prog$half_year_mod <= 3, year_prog$tm-year_prog$year_mod, 
					     year_prog$tm+12-year_prog$year_mod)
		year_prog$tyr <- year_prog$tm_round/12
	#	print(head(year_prog))
		tmp_lfit <- survfit(Surv(tyr, event) ~ Cluster, data = year_prog, id = ID)
		tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = p_mon_pos, ggtheme = theme_bw(), legend = "right", 
				      legend.labs = levels(year_prog$Cluster), 
				      palette = "Spectral", ylim = ylim) + 
			labs(title = out_name, x = "Time (1/2 year)", color = "Cluster") 		
		png(paste(tmp_pf, 'year_kmplot.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(tmp_lgg)
		gar <- dev.off()
		tmp_cox <- coxph(Surv(tyr, event) ~ Cluster, data=year_prog, id=ID, ties="breslow")
		write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'year_cox.csv', sep = ""))

		fprog <- year_prog
		fprog$Cluster <- factor(fprog$Cluster, 
					levels = c("Good knee & general health", "Intermediate knee & general health", "Poor knee & general health", "Low supplemental vitamins"))

		png(paste(tmp_pf, 'year_forestmodel_hr.png', sep = ""), res = png_res, width = 9, height = 4, units = 'in')
		print(forest_model(coxph(Surv(tyr, event) ~ Cluster, data=fprog))+labs(title = out_name))
		gar <- dev.off()
	}
}
