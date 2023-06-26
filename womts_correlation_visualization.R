# Merge, analyze, visualize and select the predictor-outcomes
# Weihua Guo, Ph.D.
# 02/13/2023

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(viridis))

data_dir <- "/mnt/sda1/OAI_Data/data_summary"
result_folder <- "predictor_vis_230213"
select_dir <- paste(data_dir, "predictor_select_230213", sep = "/")
spf <- paste(select_dir, "/predictor_select_230213_", sep = "")
result_dir <- paste(data_dir, result_folder, sep = "/")
dir.create(result_dir, showWarnings=F)
rpf <- paste(result_dir, "/", result_folder, "_", sep = "")

mt <- Sys.time()
cor_files <- list.files(select_dir, pattern = "_cor.csv")
i <- 0
for (icf in cor_files) {
	cat(icf, "\n")
	pf <- str_split_fixed(icf, "_c", n = 2)[1]
	cor_df <- read.csv(paste(select_dir, "/", icf, sep = ""), row.names = 1)
	predictor <- colnames(cor_df)
	cor_df$outcome <- rownames(cor_df)
	gath_cor_df <- gather(cor_df, "predictor", "cor", predictor)
	gath_cor_df$op <- str_c(gath_cor_df$outcome, "__", gath_cor_df$predictor)
	gath_cor_df$predictor <- NULL
	gath_cor_df$outcome <- NULL
	cor_df$outcome <- NULL

	p_df <- read.csv(paste(select_dir, "/", pf, "_pval.csv", sep = ""), row.names = 1)
	predictor <- colnames(p_df)
	p_df$outcome <- rownames(p_df)
	gath_p_df <- gather(p_df, "predictor", "pval", predictor)
	gath_p_df$op <- str_c(gath_p_df$outcome, "__", gath_p_df$predictor)
	gath_p_df$predictor <- NULL
	gath_p_df$outcome <- NULL
	p_df$outcome <- NULL

	gath_df <- merge(gath_cor_df, gath_p_df, by = "op")

	sn_df <- read.csv(paste(select_dir, "/", pf, "_smp_num.csv", sep = ""), row.names = 1)
	predictor <- colnames(sn_df)
	sn_df$outcome <- rownames(sn_df)
	gath_sn_df <- gather(sn_df, "predictor", "sample_num", predictor)
	gath_sn_df$op <- str_c(gath_sn_df$outcome, "__", gath_sn_df$predictor)
	gath_sn_df$predictor <- NULL
	gath_sn_df$outcome <- NULL
	sn_df$outcome <- NULL

	gath_df <- merge(gath_df, gath_sn_df, by = "op")

	if (i == 0) {
		merge_cor_df <- cor_df
		merge_p_df <- p_df
		merge_sn_df <- sn_df
		merge_gath_df <- gath_df
	} else {
		merge_cor_df <- rbind(merge_cor_df, cor_df)
		merge_p_df <- rbind(merge_p_df, p_df)
		merge_sn_df <- rbind(merge_sn_df, sn_df)
		merge_gath_df <- rbind(merge_gath_df, gath_df)

	}
	i <- i+1
}
write.csv(merge_cor_df, paste(rpf, "cor_merge.csv", sep = ""))
write.csv(merge_p_df, paste(rpf, "pval_merge.csv", sep = ""))
write.csv(merge_sn_df, paste(rpf, "sn_merge.csv", sep = ""))
write.csv(merge_gath_df, paste(rpf, "gath_merge.csv", sep = "")) #SD10


cat("Merge ")
print(Sys.time()-mt)

cat("Summarize the correlation...\n")
merge_gath_df$outcome <- str_split_fixed(merge_gath_df$op, "__", n = 2)[,1]
merge_gath_df$predictor <- str_split_fixed(merge_gath_df$op, "__", n = 2)[,2]
merge_gath_df$outcome_year <- str_sub(merge_gath_df$outcome, 1,3)
merge_gath_df$outcome_side <- str_sub(merge_gath_df$outcome, -1,-1)
merge_gath_df$outcome_category <- str_sub(merge_gath_df$outcome, 4,-2)

outcome_df <- as.data.frame(matrix(ncol = 1, nrow = nrow(merge_cor_df)))
colnames(outcome_df) <- c("name")
rownames(outcome_df) <- rownames(merge_cor_df)
outcome_df$year <- str_sub(rownames(merge_cor_df), 1,3)
outcome_df$side <- str_sub(rownames(merge_cor_df), -1, -1)
outcome_df[,1] <- str_sub(rownames(merge_cor_df), 4,-2)

cat("\tPer outcome category...\n")
print(head(merge_gath_df))
ft_cor_df <- merge_cor_df[str_detect(rownames(merge_cor_df), 
				     paste(c("Y04", "Y08"), collapse= "|")),]
ft_p_df <- merge_p_df[str_detect(rownames(merge_p_df), 
				     paste(c("Y04", "Y08"), collapse= "|")),]

predictors <- colnames(ft_cor_df)

for (iotcm in unique(merge_gath_df$outcome_category)) {
	cat("\t", iotcm, "\n")

	if (FALSE) {
	cat("\t\tYear 4 vs Year 8...\n")
	tmp_ft_cor_df <- ft_cor_df[str_detect(rownames(ft_cor_df), iotcm),]
	tmp_ft_cor_df$outcome_side <- str_sub(rownames(tmp_ft_cor_df), -1,-1)
	tmp_ft_cor_df$outcome_year <- str_sub(rownames(tmp_ft_cor_df), 1,3)
	tmp_ft_gath_cor <- gather(tmp_ft_cor_df, "predictor", "cor", predictors)
	tmp_ft_cor_diff <- spread(tmp_ft_gath_cor, "outcome_year", "cor")

	tmp_ft_p_df <- ft_p_df[str_detect(rownames(ft_p_df), iotcm),]
	tmp_ft_p_df$outcome_side <- str_sub(rownames(tmp_ft_p_df), -1,-1)
	tmp_ft_p_df$outcome_year <- str_sub(rownames(tmp_ft_p_df), 1,3)
	tmp_ft_gath_p <- gather(tmp_ft_p_df, "predictor", "p", predictors)
	tmp_ft_p_diff <- spread(tmp_ft_gath_p, "outcome_year", "p")
	
	tmp_ft_p_diff$diff_p005_flag <- ifelse(tmp_ft_p_diff$Y04 <=0.05,
					  ifelse(tmp_ft_p_diff$Y08 <= 0.05, "BS", "Y04"),
					  ifelse(tmp_ft_p_diff$Y08 <= 0.05, "Y08", "BNS"))
	tmp_ft_diff <- merge(tmp_ft_cor_diff, tmp_ft_p_diff[,c("outcome_side", "predictor", "diff_p005_flag")], 
			     by = c("outcome_side", "predictor"))

	print(table(tmp_ft_diff$diff_p005_flag))
	tmp_ft_diff$ft_diff <- tmp_ft_diff$Y04-tmp_ft_diff$Y08
	write.csv(tmp_ft_diff, paste(rpf, iotcm, "_ft_diff.csv", sep = ""))
	}

	tmp_gath_df <- merge_gath_df[merge_gath_df$outcome_category == iotcm,]
	cat("\t\tPer side, per predictor...\n")
	tmp_sum_df <- tmp_gath_df %>%
		group_by(outcome_side, predictor) %>%
		summarise(cor_avg = mean(cor, na.rm=T),
			  cor_sd = sd(cor, na.rm=T),
			  cor_min = min(cor, na.rm = T),
			  cor_max = max(cor, na.rm = T),
			  cor_range = max(cor, na.rm = T)-min(cor, na.rm=T),
			  n_time = n(),
			  sn_avg = mean(sample_num, na.rm=T),
			  sn_sd = sd(sample_num, na.rm=T),
			  sn_min = min(sample_num, na.rm = T),
			  sn_max = max(sample_num, na.rm = T),
			  sn_range = max(sample_num, na.rm = T)-min(sample_num, na.rm=T),
		)
	tmp_sum_df <- tmp_sum_df[!is.na(tmp_sum_df$cor_avg),]
	tmp_sum_df <- as.data.frame(tmp_sum_df[order(tmp_sum_df$cor_avg, decreasing = T),])
	tmp_sum_df$cor_avg_abs <- abs(tmp_sum_df$cor_avg)
	tmp_sum_df <- tmp_sum_df[order(tmp_sum_df$cor_avg_abs, decreasing=T),]
	tmp_sum_df$avg_abs_order <- 1:nrow(tmp_sum_df)

	tmp_sum_df$avg_order <- -10
	tmp_sum_df$avg_order[tmp_sum_df$cor_avg >0.0] <- order(tmp_sum_df$cor_avg[tmp_sum_df$cor_avg>0.0], decreasing = TRUE)
	tmp_sum_df$avg_order[tmp_sum_df$cor_avg <=0.0] <- order(tmp_sum_df$cor_avg[tmp_sum_df$cor_avg<=0.0])

	tmp_sum_df$avg_abs_label <- "NS"
	tmp_sum_df$avg_abs_label[tmp_sum_df$avg_abs_order <= 20] <- tmp_sum_df$predictor[tmp_sum_df$avg_abs_order <= 20]
	tmp_sum_df$avg_label <- "NS"
	tmp_sum_df$avg_label[tmp_sum_df$avg_order <= 10] <- tmp_sum_df$predictor[tmp_sum_df$avg_order <= 10]
	write.csv(tmp_sum_df, paste(rpf, iotcm, "_sum_df.csv", sep = ""))

	tmp_merge_df <- merge(tmp_gath_df, tmp_sum_df, by = c("outcome_side", "predictor"), all.x = T)
	tmp_merge_df <- tmp_merge_df[order(abs(tmp_merge_df$cor), decreasing = T),]

	avg_color <- turbo(length(unique(tmp_merge_df$avg_label))-1)
	names(avg_color) <- unique(tmp_merge_df$avg_label)[unique(tmp_merge_df$avg_label) != "NS"]
	avg_color[["NS"]] <- "grey"

	avg_alpha <- rep(1,length(unique(tmp_merge_df$avg_label)))
	names(avg_alpha) <- unique(tmp_merge_df$avg_label)
	avg_alpha[["NS"]] <- 0.1

	gg <- ggplot(tmp_merge_df, aes(x = outcome_year, y = cor, 
				       color = avg_label, group = predictor,
				       alpha = avg_label)) +
		geom_point() +
		geom_line() +
		scale_color_manual(values = avg_color) +
		scale_alpha_manual(values = avg_alpha) +
		facet_wrap(~outcome_side, nrow = 1, scales = "free") +
		labs(x = "Follow up time (year)", 
		     y = "Association between predictors and outcomes",
		     title = iotcm,
		     color = "Top 10\npredictors\n(Both sides)"
		) +
		scale_y_continuous(breaks=0.1*(-10:10)) +
		guides(alpha = "none") +
		theme_classic() 	
	ggsave(paste(rpf, iotcm, "_avg_line_dot.png", sep = ""), dpi = 300, width = 16, height = 6)

	avg_abs_color <- turbo(length(unique(tmp_merge_df$avg_abs_label))-1)
	names(avg_abs_color) <- unique(tmp_merge_df$avg_abs_label)[unique(tmp_merge_df$avg_abs_label) != "NS"]
	avg_abs_color[["NS"]] <- "grey"

	avg_abs_alpha <- rep(1,length(unique(tmp_merge_df$avg_abs_label)))
	names(avg_abs_alpha) <- unique(tmp_merge_df$avg_abs_label)
	avg_abs_alpha[["NS"]] <- 0.1

	gg <- ggplot(tmp_merge_df, aes(x = outcome_year, y = cor, 
				       color = avg_abs_label, group = predictor,
				       alpha = avg_abs_label)) +
		geom_point() +
		geom_line() +
		scale_color_manual(values = avg_abs_color) +
		scale_alpha_manual(values = avg_abs_alpha) +
		facet_wrap(~outcome_side, nrow = 1, scales = "free") +
		labs(x = "Follow up time (year)", 
		     y = "Association between predictors and outcomes",
		     title = iotcm,
		     color = "Top 10\npredictors\n(Absolute value)"
		) +
		scale_y_continuous(breaks=0.1*(-10:10)) +
		guides(alpha = "none") +
		theme_classic() 	
	ggsave(paste(rpf, iotcm, "_avg_abs_line_dot.png", sep = ""), dpi = 300, width = 16, height = 6)
}

cat("Heatmap visualization...\n")

for (iotcm in unique(outcome_df$name)) {
	cat("\t", iotcm, "\n")

	tmp_otcm <- outcome_df[outcome_df$name == iotcm,]

	col_fun <- colorRamp2(c(-1, 0, 1), c("forestgreen", "cornsilk", "darkorange"))
	hm_df <- as.matrix(merge_cor_df[rownames(tmp_otcm),])
	hm_df <- hm_df[rowSums(is.na(hm_df))<0.7*ncol(hm_df),]
#	hm_df <- hm_df[,colSums(is.na(hm_df))==0]
	tmp_otcm <- tmp_otcm[rownames(hm_df),]

	pfilter_df <- merge_p_df[rownames(tmp_otcm),]
	psig_df <- as.data.frame(merge_p_df[rownames(tmp_otcm),])
	psig_df <- as.data.frame(ifelse(psig_df <= 0.05, 1, 0))
	all_vars <- colnames(psig_df)
	psig_df$outcome <- rownames(psig_df)
	gath_psig_df <- gather(psig_df, "vars", "pflag", all_vars)
	gath_psig_df$side <- str_sub(gath_psig_df$outcome, -1,-1)
	sum_psig_df <- gath_psig_df %>%
		group_by(vars) %>%
		summarise(sig_num = sum(pflag, na.rm=T))
	print(head(sum_psig_df))
	psig_df$outcome <- NULL
	psig_df[is.na(psig_df)] <- 0.0
	print(psig_df[,1:6])

	row_ann <- HeatmapAnnotation(Side = tmp_otcm$side, which = "row",
				     col = list(Side = c("L"="goldenrod", "R" = "darkorchid")))
	psig_cor <- Heatmap(psig_df,
			  name = "1: p<=0.05; 0: p>0.05",
			  col = colorRamp2(c(0, 1), c("dodgerblue", "firebrick")),
			  cluster_rows = FALSE,
			  cluster_columns = TRUE, 
			  column_names_gp = gpar(fontsize = 12),
			  row_split = tmp_otcm$side,
			  row_names_gp = gpar(fontsize = 12),
			  left_annotation = row_ann,
			  row_order = order(rownames(tmp_otcm)),
			  show_column_names = FALSE,
			  column_dend_height = unit(4, "cm"), 
			  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3.6, "cm")))
	png(paste(rpf, iotcm, "_psig_flag_005_heatmap.png", sep = ""), res = 300, width = 16, height = 6, units = 'in')
	draw(psig_cor, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()

	top_vars <- sum_psig_df %>% top_n(wt = sig_num, 20)
	sig_hm_df <- hm_df[,top_vars$vars]

	row_ann <- HeatmapAnnotation(Side = tmp_otcm$side, which = "row",
				     col = list(Side = c("L"="goldenrod", "R" = "darkorchid")))
	sig_hm_cor <- Heatmap(sig_hm_df,
			  name = "Correlation/Association",
			  col = col_fun,
			  cluster_rows = FALSE,
			  cluster_columns = TRUE, 
			  column_names_gp = gpar(fontsize = 9),
			  row_split = tmp_otcm$side,
			  row_names_gp = gpar(fontsize = 9),
			  left_annotation = row_ann,
			  row_order = order(rownames(tmp_otcm)),
			  show_column_names = TRUE,
			  column_dend_height = unit(4, "cm"), 
			  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3.6, "cm")))
	png(paste(rpf, iotcm, "_sgnft_cor_heatmap.png", sep = ""), res = 300, width = 16, height = 6, units = 'in')
	draw(sig_hm_cor, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()

	hm_df <- hm_df[,sum_psig_df$vars[order(sum_psig_df$sig_num)]]
	row_ann <- HeatmapAnnotation(Side = tmp_otcm$side, which = "row",
				     col = list(Side = c("L"="goldenrod", "R" = "darkorchid")))
	hm_cor <- Heatmap(hm_df,
			  name = "Correlation/Association",
			  col = col_fun,
			  cluster_rows = FALSE,
			  cluster_columns = TRUE, 
			  column_names_gp = gpar(fontsize = 12),
			  row_split = tmp_otcm$side,
			  row_names_gp = gpar(fontsize = 12),
			  left_annotation = row_ann,
			  row_order = order(rownames(tmp_otcm)),
			  show_column_names = FALSE,
			  column_dend_height = unit(4, "cm"), 
			  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3.6, "cm")))
	png(paste(rpf, iotcm, "_cor_heatmap.png", sep = ""), res = 300, width = 16, height = 6, units = 'in')
	draw(hm_cor, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()

	if (FALSE) {
	hm_df <- as.matrix(merge_sn_df[rownames(tmp_otcm), colnames(hm_df)])
	hm_df <- log2(hm_df)
	col_fun <- colorRamp2(c(min(hm_df, na.rm=T), median(hm_df, na.rm=T), max(hm_df, na.rm=T)), 
			      c("dodgerblue", "cornsilk", "firebrick"))
	row_ann <- HeatmapAnnotation(Side = tmp_otcm$side, which = "row",
				     col = list(Side = c("L"="goldenrod", "R" = "darkorchid")))

	hm_obj <- Heatmap(hm_df,
			  name = "Available sample number (log2)",
#			  col = col_fun,
			  cluster_rows = FALSE,
			  cluster_columns = FALSE, 
			  column_names_gp = gpar(fontsize = 12),
			  row_split = tmp_otcm$side,
			  row_names_gp = gpar(fontsize = 12),
			  left_annotation = row_ann,
			  row_order = order(rownames(tmp_otcm)),
			  column_order = order(sum_psig_df$sig_num, decreasing=T),
			  show_column_names = FALSE,
			  column_dend_height = unit(5, "cm"), 
			  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3.6, "cm")))
	png(paste(rpf, iotcm, "_sn_heatmap.png", sep = ""), res = 300, width = 16, height = 6, units = 'in')
	draw(hm_obj, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()

	hm_obj <- Heatmap(hm_df,
			  name = "Available sample number (log2)",
#			  col = col_fun,
			  cluster_rows = FALSE,
			  cluster_columns = FALSE, 
			  column_order = column_order(hm_cor),
			  column_names_gp = gpar(fontsize = 12),
			  row_split = tmp_otcm$side,
			  row_names_gp = gpar(fontsize = 12),
			  left_annotation = row_ann,
			  row_order = order(rownames(tmp_otcm)),
			  show_column_names = FALSE,
			  column_dend_height = unit(5, "cm"), 
			  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3.6, "cm")))
	png(paste(rpf, iotcm, "_sn_align_heatmap.png", sep = ""), res = 300, width = 16, height = 6, units = 'in')
	draw(hm_obj, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()
	}
}

