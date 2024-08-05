# Visualize the outcomes crossing the clusters
# Weihua Guo, Ph.D.
# 03/14/2023

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

clst_dir <- "/mnt/sda1/OAI_Data/kmean_cluster_12252020/"
data_dir <- "/mnt/sda1/OAI_Data/data_summary/"

exp_id <- "kmean_pca_umap"
input_id <- "12112020_data_clean"
clean_co <- 'v25'
input_prefix <- paste(clst_dir, 'clean_', clean_co, '_', sep = "")
cluster_num <- 4
scale_top <- 4

png_res <- 300
top_m <- 6

ppf <- paste(data_dir, "outcome_visualization_", sep = "")

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
oc_df<-read.csv(paste(data_dir, "outcome_all_dataframe_", input_id, ".csv", sep=''),row.names=1)
rtday_df <- as.data.frame(read_excel(paste(data_dir, "/real_dates_manual.xlsx", sep = ""), sheet = "days"))
rtd_df <- as.data.frame(read_excel(paste(data_dir, "/real_dates_manual.xlsx", sep = ""), sheet = "days"))
rtm_df <- as.data.frame(read_excel(paste(data_dir, "/real_dates_manual.xlsx", sep = ""), sheet = "months"))

kr_type <- "total_lastfollowup"
kr_df <- read.csv(paste(data_dir, kr_type, "_merge_patient_basic_outcome_information.csv", sep = ""), header = T)
rownames(kr_df) <- kr_df$ID

kr_df$left_year <- kr_df$left_time/365
kr_df$left_round_year <- str_c("LEFT_TKR_Y", str_pad(round(kr_df$left_time/365), 2, pad = "0"))
kr_df$right_year <- kr_df$right_time/365
kr_df$right_round_year <- str_c("RIGHT_TKR_Y", str_pad(round(kr_df$right_time/365), 2, pad = "0"))

lkr_spread <- spread(kr_df[,c("ID", "left_round_year", "left_kr")], "left_round_year", "left_kr")
rkr_spread <- spread(kr_df[,c("ID", "right_round_year", "right_kr")], "right_round_year", "right_kr")
kr_spread <- merge(lkr_spread, rkr_spread, by = "ID")
kr_spread[is.na(kr_spread)] <- "3: No"
kr_gath <- gather(kr_spread, "time_side", "kr", colnames(kr_spread)[str_detect(colnames(kr_spread), "TKR_Y")])
kr_gath$side <- str_split_fixed(kr_gath$time_side, "_TKR_", n = 2)[,1]
kr_gath$time <- str_split_fixed(kr_gath$time_side, "_TKR_", n = 2)[,2]
kr_gath$tkr <- kr_gath$kr
kr_gath$tkr[str_detect(kr_gath$tkr, "\\.\\:\\ Missing")] <- "3: No"

kr_umap_df <- merge(kr_gath, umap_df, by = "ID", all.x = T)

kr_sum_df <- kr_umap_df %>%
	filter(!is.na(Cluster)) %>%
	group_by(side, time, Cluster, tkr) %>%
	summarize(n = n())
bar_gg <- ggplot(kr_sum_df, aes(x = time, y = n, fill = tkr)) +
	geom_bar(position = "fill", stat = "identity") +
	scale_fill_brewer(palette = "Set1") +
	facet_grid(side~Cluster) +
	labs(x = "Time (year)", y = "Relative number of knees", fill = "Total knee replacement") +
	theme_bw() +
	theme(legend.position = "top")
ggsave(paste(ppf, "tkr_bar_stack.png", sep = ""), dpi = png_res, width = 16, height = 6)

outcome_patterns<-c('XRKLL', 'XRKLR','WOMADLL', 'WOMADLR', 'WOMKPL', 'WOMKPR', 'MCMJSWL', 'MCMJSWR', 
		    'WOMSTFL', 'WOMSTFR', 'WOMTSL', 'WOMTSR')

oc_cols <- c()
for (iop in outcome_patterns) {
	tmp_mask <- str_detect(colnames(oc_df), iop)
	oc_cols <- c(oc_cols, colnames(oc_df)[tmp_mask])
}
#print(oc_cols)
#print(head(umap_df))

cat("Converting days to years\n")
#print(head(rtd_df))
rownames(rtd_df) <- rtd_df$ID
rtd_df$P01 <- NULL
rtd_df$ID <- NULL

rt_ydf <- rtd_df
max_y <- round(max(rtd_df/365, na.rm=T))
#print(max_y)

rty_vconv_df <- as.data.frame(matrix(nrow = nrow(rtd_df), ncol = max_y+1))
rownames(rty_vconv_df) <- rownames(rtd_df)
colnames(rty_vconv_df) <- str_c("Y", str_pad(0:max_y, 2, pad = "0"))

for (i in 1:nrow(rtd_df)) {
	cat(rownames(rtd_df)[i], "\n")
	tmp_r <- rtd_df[i,]
	tmp_yr <- tmp_r/365
	tmp_ryr <- round(tmp_yr)
	tmp_dyr <- tmp_r%%365
	prv_v <- "V"
	for (iv in 0:max_y) {
		yiv <- str_c("Y", str_pad(iv, 2, pad = "0"))
		cat("\t", yiv, "\n")
		tmp_diff <- tmp_yr-iv
		tmp_round_diff <- tmp_ryr-iv
		tmp_diff_mask <- tmp_round_diff == 0
		if (sum(tmp_diff_mask, na.rm = T) > 1) {
			tmp_diff_mask <- abs(tmp_diff)[1,] == min(abs(tmp_diff), na.rm=T)
		}		
		tmp_v <- colnames(tmp_diff)[tmp_diff_mask]
		tmp_v <- tmp_v[!is.na(tmp_v)]
		if (length(tmp_v) == 1) {
			if (tmp_v != prv_v) {
				rty_vconv_df[rownames(rtd_df)[i], yiv] <- tmp_v[!is.na(tmp_v)]
				prv_v <- tmp_v
			}
		}
#		print(abs(tmp_diff)[1,] == min(abs(tmp_diff), na.rm=T))
#		print(tmp_diff_mask)
	}
#	print(rty_oc_df[1:9,1:6])
#	print(head(rty_vconv_df))
#	print(tmp_yr)
}
write.csv(rty_vconv_df, paste(data_dir, "outcome_real_date_conversion_year.csv", sep = "")) # SD6
write.csv(rtd_df/365, paste(data_dir, "outcome_real_date_year_decimal.csv", sep = ""))
rty_df <- read.csv(paste(data_dir, "outcome_real_date_conversion_year.csv", sep = ""), header = T, row.names = 1)

cat("Merge cluster results with outcomes...\n")
use_oc_df <- oc_df[,oc_cols]
print(dim(use_oc_df))
use_oc_df$ID <- rownames(use_oc_df)
gath_df <- gather(use_oc_df, "outcome", "value", oc_cols)
gath_df$v <- str_sub(gath_df$outcome, 1,3)
gath_df$o <- str_sub(gath_df$outcome, 4,-2)
gath_df$side <- str_sub(gath_df$outcome, -1, -1)
print(use_oc_df[1:9,1:6])
q(save = "no")
#print(head(gath_df))
#print(unique(gath_df$outcome))

use_rt_df <- rtday_df
use_rt_df$P01 <- NULL
gath_drt_df <- gather(use_rt_df, "v", "day", colnames(use_rt_df)[str_detect(colnames(use_rt_df), "V")])

use_rt_df <- rtm_df
use_rt_df$P01 <- NULL
gath_mrt_df <- gather(use_rt_df, "v", "month", colnames(use_rt_df)[str_detect(colnames(use_rt_df), "V")])

rty_df$ID <- rownames(rty_df)
gath_yrt_df <- gather(rty_df, "y", "v", colnames(rty_df)[str_detect(colnames(rty_df), "Y")])
gath_yrt_df$year <- as.numeric(str_replace_all(gath_yrt_df$y, "Y", ""))
print(head(gath_yrt_df))

gath_df <- merge(gath_df, gath_drt_df, by = c("ID", "v"), all.x=T)
gath_df <- merge(gath_df, gath_mrt_df, by = c("ID", "v"), all.x=T)
gath_df <- merge(gath_df, gath_yrt_df, by = c("ID", "v"), all.x=T)
gath_df <- merge(gath_df, umap_df, by = "ID", all.x=T)

gath_ydf <- gath_df[!is.na(gath_df$y),]
gath_ydf$yo <- str_c(gath_ydf$y, gath_ydf$o, gath_ydf$side)
rty_oc_df <- spread(gath_ydf[,c("ID", "yo", "value")], "yo", "value")
#print(dim(rty_oc_df))
#print(rty_oc_df[1:9,1:6])
write.csv(rty_oc_df, paste(data_dir, "outcome_all_dataframe_", input_id, "_real_date_year.csv", sep=''))
print(head(gath_ydf))
#print(unique(gath_ydf$yo))

for (iotcm in unique(gath_ydf$o)) {
	cat(iotcm, "\n")
	tmp_gath_df <- gath_ydf[gath_ydf$o == iotcm,]
	if (str_detect(iotcm, "KL")) {
		cat("\t", iotcm, "\n")
		cat("\t\tCalculate the differences...\n")
#		print(head(tmp_gath_df))
		tmp_diff_df <- tmp_gath_df[,c("ID", "value", "o", "side", "y", "Cluster")] %>%
			spread(y, value)
		tmp_diff_df <- gather(tmp_diff_df, "y", "value", colnames(tmp_diff_df)[str_detect(colnames(tmp_diff_df), "Y") & (colnames(tmp_diff_df) != "Y00")])
		tmp_diff_df$diff <- tmp_diff_df$value - tmp_diff_df$Y00
#		print(head(tmp_diff_df))
		tmp_diff_df$diff[is.na(tmp_diff_df$diff)] <- "N/A"
		print(head(tmp_diff_df))
		tmp_diff_df$diff <- factor(tmp_diff_df$diff)
		tmp_diff_sum_df <- tmp_diff_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y, Cluster, diff) %>%
			summarize(n = n())
		write.csv(tmp_diff_sum_df, paste(ppf, iotcm, "_diff_descriptive_stats.csv", sep = ""))

		bar_gg <- ggplot(tmp_diff_sum_df, aes(x = y, y = n, fill = diff)) +
			geom_bar(position = "fill", stat = "identity") +
			scale_fill_brewer(palette = "Dark2") +
			facet_grid(side~Cluster) +
			labs(x = "Time (year)", y = "Relative number of knees", fill = "KL grade difference") +
			theme_bw() +
			theme(legend.position = "top")
		ggsave(paste(ppf, "diff_klg_bar_stack_facet_cluster.png", sep = ""), dpi = png_res, width = 16, height = 6)

		bar_gg <- ggplot(tmp_diff_sum_df, aes(x = y, y = n, fill = Cluster)) +
			geom_bar(position = "fill", stat = "identity") +
			scale_fill_brewer(palette = "Spectral") +
			facet_grid(side~diff) +
			labs(x = "Time (year)", y = "Relative number of knees", fill = "Cluster") +
			theme_bw() +
			theme(legend.position = "top")
		ggsave(paste(ppf, "diff_klg_bar_stack_facet_klg.png", sep = ""), dpi = png_res, width = 16, height = 6)

		bar_gg <- ggplot(tmp_diff_sum_df, aes(x = Cluster, y = n, fill = diff)) +
			geom_bar(position = "fill", stat = "identity") +
			scale_fill_brewer(palette = "RdYlBu") +
			facet_grid(side~y) +
			labs(x = "Cluster", y = "Relative number of knees", fill = "KL grade difference") +
			theme_bw() +
			theme(legend.position = "top", 
			      axis.text.x = element_text(angle = 45, hjust = 1),
			      plot.margin = margin(0, 0.5, 0, 2, "cm")
			)
		ggsave(paste(ppf, "diff_klg_bar_stack_facet_year.png", sep = ""), dpi = png_res, width = 16, height = 9)


		cat("\t\tVisualize the outcome variables...\n")
#		print(unique(tmp_gath_df$value))
		tmp_gath_df$value[is.na(tmp_gath_df$value)] <- "N/A"
		print(head(tmp_gath_df))
		tmp_gath_df$value <- factor(tmp_gath_df$value, levels = c("N/A",0,1,2,3,4))
		tmp_sum_df <- tmp_gath_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y, Cluster, value) %>%
			summarize(n = n())
		write.csv(tmp_sum_df, paste(ppf, iotcm, "_descriptive_stats.csv", sep = ""))		
#		print(head(tmp_sum_df))

		bar_gg <- ggplot(tmp_sum_df, aes(x = y, y = n, fill = value)) +
			geom_bar(position = "fill", stat = "identity") +
			scale_fill_brewer(palette = "Dark2") +
			facet_grid(side~Cluster) +
			labs(x = "Time (year)", y = "Relative number of knees", fill = "KL grade") +
			theme_bw() +
			theme(legend.position = "top")
		ggsave(paste(ppf, "klg_bar_stack_facet_cluster.png", sep = ""), dpi = png_res, width = 16, height = 6)

		bar_gg <- ggplot(tmp_sum_df, aes(x = y, y = n, fill = Cluster)) +
			geom_bar(position = "fill", stat = "identity") +
			scale_fill_brewer(palette = "Spectral") +
			facet_grid(side~value) +
			labs(x = "Time (year)", y = "Relative number of knees", fill = "Cluster") +
			theme_bw() +
			theme(legend.position = "top")
		ggsave(paste(ppf, "klg_bar_stack_facet_klg.png", sep = ""), dpi = png_res, width = 16, height = 6)

		bar_gg <- ggplot(tmp_sum_df, aes(x = Cluster, y = n, fill = value)) +
			geom_bar(position = "fill", stat = "identity") +
			scale_fill_brewer(palette = "YlOrRd") +
			facet_grid(side~y) +
			labs(x = "Cluster", y = "Relative number of knees", fill = "KL grade") +
			theme_bw() +
			theme(legend.position = "top", 
			      axis.text.x = element_text(angle = 45, hjust = 1),
			      plot.margin = margin(0, 0.5, 0, 2, "cm")
			)
		ggsave(paste(ppf, "klg_bar_stack_facet_year.png", sep = ""), dpi = png_res, width = 16, height = 9)
	} else {
		if (str_detect(iotcm, "JSW")) {tmp_gath_df <- tmp_gath_df[tmp_gath_df$y %in% c("Y00", "Y01", "Y02", "Y03", "Y04", "Y08"),]}

		cat("\t\tCalculate the differences...\n")
		print(head(tmp_gath_df))
		tmp_diff_df <- tmp_gath_df[,c("ID", "value", "o", "side", "y", "Cluster")] %>%
			spread(y, value)
		tmp_diff_df <- gather(tmp_diff_df, "y", "value", colnames(tmp_diff_df)[str_detect(colnames(tmp_diff_df), "Y") & (colnames(tmp_diff_df) != "Y00")])
		tmp_diff_df$diff <- tmp_diff_df$value - tmp_diff_df$Y00
		tmp_diff_df$rel_diff <- tmp_diff_df$diff/(tmp_diff_df$Y00+0.01)
		tmp_diff_df$fc <- (tmp_diff_df$value+0.01)/(tmp_diff_df$Y00+0.01)

#		print(head(tmp_diff_df))

		tmp_diff_sum_df <- tmp_diff_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y, Cluster) %>%
			summarise(n = n(),
				  avg = mean(diff, na.rm = T),
				  sd = sd(diff, na.rm = T),
				  median = median(diff, na.rm=T),
				  se = sd/sqrt(n),
				  mean_se_upper = avg+se,
				  mean_se_lower = avg-se
			)
		write.csv(tmp_diff_sum_df, paste(ppf, iotcm, "_diff_descriptive_stats.csv", sep = ""))

#		print(head(tmp_sum_df))
		tmp_kw_df <- tmp_diff_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y) %>%
			kruskal_test(diff~Cluster) %>%
			add_significance('p') %>%
			adjust_pvalue('p') %>%
			add_significance('p.adj')
		write.csv(tmp_kw_df, paste(ppf, iotcm, "_diff_kruskal_wallis_res.csv", sep = ""))
		tmp_prwlx_df <- tmp_diff_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y) %>%
			pairwise_wilcox_test(diff~Cluster, detailed=T) %>%
			add_significance('p') %>%
			adjust_pvalue('p') %>%
			add_significance('p.adj')
		write.csv(tmp_prwlx_df, paste(ppf, iotcm, "_diff_pairwise_wilcox_res.csv", sep = ""))	

		ysum_gg <- ggplot(tmp_diff_sum_df, aes(x = y, y = avg, color = Cluster, group = Cluster)) +
			geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.1) +
			geom_line() +
			geom_point(size = 0.5) +
			scale_y_continuous(trans = 'log2') +
			scale_color_brewer(palette = 'Spectral') +
			labs(x = "Time (year)", y = paste(iotcm, "(fold change, Yxx/Y00)")) +
			facet_wrap(side~., scales = 'free_y', nrow = 2) +
			theme_bw()
		ggsave(paste(ppf, iotcm, "_diff_real_date_year_sum.png", sep=''), ysum_gg, 
		       dpi = 300, width = 9, height = 6)

		q(save = "no")


		tmp_fc_sum_df <- tmp_diff_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y, Cluster) %>%
			summarise(n = n(),
				  avg = mean(fc, na.rm = T),
				  sd = sd(fc, na.rm = T),
				  median = median(fc, na.rm=T),
				  se = sd/sqrt(n),
				  mean_se_upper = avg+se,
				  mean_se_lower = avg-se
			)
		write.csv(tmp_fc_sum_df, paste(ppf, iotcm, "_fc_descriptive_stats.csv", sep = ""))

#		print(head(tmp_sum_df))
		tmp_kw_df <- tmp_diff_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y) %>%
			kruskal_test(fc~Cluster) %>%
			add_significance('p') %>%
			adjust_pvalue('p') %>%
			add_significance('p.adj')
		write.csv(tmp_kw_df, paste(ppf, iotcm, "_fc_kruskal_wallis_res.csv", sep = ""))
		tmp_prwlx_df <- tmp_diff_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y) %>%
			pairwise_wilcox_test(fc~Cluster, detailed=T) %>%
			add_significance('p') %>%
			adjust_pvalue('p') %>%
			add_significance('p.adj')
		write.csv(tmp_prwlx_df, paste(ppf, iotcm, "_fc_pairwise_wilcox_res.csv", sep = ""))	

		ysum_gg <- ggplot(tmp_fc_sum_df, aes(x = y, y = avg, color = Cluster, group = Cluster)) +
			geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.1) +
			geom_line() +
			geom_point(size = 0.5) +
			scale_y_continuous(trans = 'log2') +
			scale_color_brewer(palette = 'Spectral') +
			labs(x = "Time (year)", y = paste(iotcm, "(fold change, Yxx/Y00)")) +
			facet_wrap(side~., scales = 'free_y', nrow = 2) +
			theme_bw()
		ggsave(paste(ppf, iotcm, "_fc_real_date_year_sum.png", sep=''), ysum_gg, 
		       dpi = 300, width = 9, height = 6)

		cat("\t\tVisualize the outcome variables...\n")
		tmp_sum_df <- tmp_gath_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y, Cluster) %>%
			summarise(n = n(),
				  avg = mean(value, na.rm = T),
				  sd = sd(value, na.rm = T),
				  median = median(value, na.rm=T),
				  se = sd/sqrt(n),
				  mean_se_upper = avg+se,
				  mean_se_lower = avg-se
			)
		write.csv(tmp_sum_df, paste(ppf, iotcm, "_descriptive_stats.csv", sep = ""))

#		print(head(tmp_sum_df))
		tmp_kw_df <- tmp_gath_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y) %>%
			kruskal_test(value~Cluster) %>%
			add_significance('p') %>%
			adjust_pvalue('p') %>%
			add_significance('p.adj')
		write.csv(tmp_kw_df, paste(ppf, iotcm, "_kruskal_wallis_res.csv", sep = ""))
		tmp_prwlx_df <- tmp_gath_df %>%
			filter(!is.na(Cluster)) %>%
			group_by(side, y) %>%
			pairwise_wilcox_test(value~Cluster, detailed=T) %>%
			add_significance('p') %>%
			adjust_pvalue('p') %>%
			add_significance('p.adj')
		write.csv(tmp_prwlx_df, paste(ppf, iotcm, "_pairwise_wilcox_res.csv", sep = ""))	

		ysum_gg <- ggplot(tmp_sum_df, aes(x = y, y = avg, color = Cluster, group = Cluster)) +
			geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.1) +
			geom_line() +
			geom_point(size = 0.5) +
			scale_color_brewer(palette = 'Spectral') +
			labs(x = "Time (year)", y = iotcm) +
			facet_wrap(side~., scales = 'free_y', nrow = 2) +
			theme_bw()
		ggsave(paste(ppf, iotcm, "_real_date_year_sum.png", sep=''), ysum_gg, 
		       dpi = 300, width = 9, height = 6)
	}
}
