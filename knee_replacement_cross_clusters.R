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
surv_folder <- "survival_ready_220927"

kr_type <- "total_lastfollowup"

clean_co <- 'v25'
input_prefix <- paste(clst_dir, 'clean_', clean_co, '_', sep = "")
cluster_num <- 4

png_res <- 300

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

kr_df <- read.csv(paste(data_dir, kr_type, "_merge_patient_basic_outcome_information.csv", sep = ""), header = T)
print(dim(kr_df))
print(colnames(kr_df))

kr_umap_df <- merge(kr_df, umap_df, by = "ID")
print(dim(kr_umap_df))
print(colnames(kr_umap_df))
cohort_cts <- kr_umap_df %>%
	group_by(V00COHORT, Cluster) %>%
	summarise(n = n())

all_gg <- ggplot(cohort_cts, aes(x = Cluster, y = n, color = V00COHORT, fill = V00COHORT)) +
	geom_bar(position="stack", stat="identity") +
	scale_color_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	labs(x = "Cluster", y = "Number of patients", color = "Cohort", fill = "Cohort") +
	coord_flip() +
	theme_bw()
ggsave(paste(data_dir, surv_folder, "/cluster_", cluster_num, "_kmean_pca_cohort_bar.png", sep = ""), all_gg, dpi = png_res, width = 9, height = 2.7)
all_gg <- ggplot(cohort_cts, aes(x = Cluster, y = n, color = V00COHORT, fill = V00COHORT)) +
	geom_bar(position="fill", stat="identity") +
	scale_color_brewer(palette = "Dark2") +
	scale_fill_brewer(palette = "Dark2") +
	labs(x = "Cluster", y = "Proportion of patients", color = "Cohort", fill = "Cohort") +
	coord_flip() +
	theme_bw()
ggsave(paste(data_dir, surv_folder, "/cluster_", cluster_num, "_kmean_pca_cohort_rel_bar.png", sep = ""), all_gg, dpi = png_res, width = 9, height = 2.7)

##### TKR crossing clusters ####
tmp_pf <- paste(data_dir, surv_folder, "/", "left_", kr_type, "_cluster_", sep = "")
tmp_lfit <- survfit(Surv(left_time, lkr_event) ~ Cluster, data = kr_umap_df)
tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = c(10, 0.82), ggtheme = theme_bw(), palette = "Spectral", censor.size=2, ylim=c(0.8,1), 
		      legend = "right", legend.title = "Cluster") + 
	labs(title = "Left knee, total knee replacement/last follow-up")
png(paste(tmp_pf, 'kmplot.png', sep = ""), res = png_res, width = 8, height = 4, units = 'in')
print(tmp_lgg)
gar <- dev.off()
tmp_cox <- coxph(Surv(left_time, lkr_event) ~ Cluster, data=kr_umap_df, ties="breslow")
write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'cox.csv', sep = ""))

tmp_pf <- paste(data_dir, surv_folder, "/", "right_", kr_type, "_cluster_", sep = "")
tmp_lfit <- survfit(Surv(right_time, rkr_event) ~ Cluster, data = kr_umap_df)
tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = c(10, 0.82), ggtheme = theme_bw(), palette = "Spectral", censor.size=2, ylim=c(0.8,1), 
		      legend = "right", legend.title = "Cluster") + 
	labs(title = "Right knee, total knee replacement/last follow-up")
png(paste(tmp_pf, 'kmplot.png', sep = ""), res = png_res, width = 8, height = 4, units = 'in')
print(tmp_lgg)
gar <- dev.off()
tmp_cox <- coxph(Surv(right_time, rkr_event) ~ Cluster, data=kr_umap_df, ties="breslow")
write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'cox.csv', sep = ""))

##### TKR crossing cohorts ####
tmp_pf <- paste(data_dir, surv_folder, "/", "left_", kr_type, "_cohort_", sep = "")
tmp_lfit <- survfit(Surv(left_time, lkr_event) ~ V00COHORT, data = kr_umap_df)
tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = c(10, 0.82), ggtheme = theme_bw(), palette = "Dark2", censor.size=2, ylim=c(0.8,1), 
		      legend = "right", legend.title = "Cohort") + 
	labs(title = "Left knee, total knee replacement/last follow-up")
png(paste(tmp_pf, 'kmplot.png', sep = ""), res = png_res, width = 8, height = 4, units = 'in')
print(tmp_lgg)
gar <- dev.off()
tmp_cox <- coxph(Surv(left_time, lkr_event) ~ V00COHORT, data=kr_umap_df, ties="breslow")
write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'cox.csv', sep = ""))

tmp_pf <- paste(data_dir, surv_folder, "/", "right_", kr_type, "_cohort_", sep = "")
tmp_lfit <- survfit(Surv(right_time, rkr_event) ~ V00COHORT, data = kr_umap_df)
tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, pval.coord = c(10, 0.82), ggtheme = theme_bw(), palette = "Dark2", censor.size=2, ylim=c(0.8,1), 
		      legend = "right", legend.title = "Cohort") + 
	labs(title = "Right knee, total knee replacement/last follow-up")
png(paste(tmp_pf, 'kmplot.png', sep = ""), res = png_res, width = 8, height = 4, units = 'in')
print(tmp_lgg)
gar <- dev.off()

##### Cox regression model ####
tmp_cox <- coxph(Surv(right_time, rkr_event) ~ V00COHORT, data=kr_umap_df, ties="breslow")
write.csv(summary(tmp_cox)$coefficients, paste(tmp_pf, 'cox.csv', sep = ""))


data_files <- list.files(paste(data_dir, surv_folder, sep = ""), pattern = '_used_clean_dataframe_with_real_time.csv')
print(data_files)
hr_cols <- c("ytyr", "outcome_name", "kr_side", "hr", "hci95", "lci95", "lrtest", "p")
hr_sum_df <- as.data.frame(matrix(ncol = length(hr_cols), nrow = 11*length(data_files)*2))
colnames(hr_sum_df) <- hr_cols
c <- 1
for (idata in data_files) {
	cat(idata, "\n")
	out_name <- str_split_fixed(idata, "_", n = 2)[1]
	tmp_data <- read.csv(paste(data_dir, surv_folder, "/", idata, sep = ""), header = T)
	merge_data <- merge(tmp_data, kr_umap_df, by = "ID", all.x = T)
#	print(head(merge_data))
	tmp_pf <- paste(data_dir, surv_folder, "/", out_name, "_", sep = "")

	year_data <- merge_data
	year_data$year_mod <- year_data$tm %% 12
	year_data <- year_data[year_data$year_mod <= 3 | year_data$year_mod >= 9,]
	year_data$tm_round <- ifelse(year_data$year_mod <= 3, year_data$tm-year_data$year_mod, 
				     year_data$tm+12-year_data$year_mod)
	year_data$tyr <- year_data$tm_round/12
	year_data$ytyr <- str_c("Y", year_data$tyr)
#	print(head(year_data))

	for (iy in unique(year_data$ytyr)) {
		cat("\t", iy, "\n")
		tmp_data <- year_data[year_data$ytyr == iy,]
		tmp_cox <- coxph(Surv(left_time, lkr_event) ~ value, data=tmp_data, ties="breslow")
		tmp_sum_cox <- summary(tmp_cox)
#		print(tmp_sum_cox)
		hr_sum_df[c, "ytyr"] <- iy
		hr_sum_df[c, "outcome_name"] <- out_name
		hr_sum_df[c, "kr_side"] <- "Left knee"
		hr_sum_df[c, "hr"] <- tmp_sum_cox$coef[2]
#		print(names(tmp_sum_cox))
		hr_sum_df[c, "hci95"] <- tmp_sum_cox$conf.int[,"upper .95"]
		hr_sum_df[c, "lci95"] <- tmp_sum_cox$conf.int[,"lower .95"]
		hr_sum_df[c, "lrtest"] <- tmp_sum_cox$logtest["test"]
		hr_sum_df[c, "p"] <- tmp_sum_cox$logtest["pvalue"]
		c <- c+1

		tmp_cox <- coxph(Surv(right_time, rkr_event) ~ value, data=tmp_data, ties="breslow")
		tmp_sum_cox <- summary(tmp_cox)
#		print(tmp_sum_cox)
		hr_sum_df[c, "ytyr"] <- iy
		hr_sum_df[c, "outcome_name"] <- out_name
		hr_sum_df[c, "kr_side"] <- "Right knee"
		hr_sum_df[c, "hr"] <- tmp_sum_cox$coef[2]
#		print(names(tmp_sum_cox))
		hr_sum_df[c, "hci95"] <- tmp_sum_cox$conf.int[,"upper .95"]
		hr_sum_df[c, "lci95"] <- tmp_sum_cox$conf.int[,"lower .95"]
		hr_sum_df[c, "lrtest"] <- tmp_sum_cox$logtest["test"]
		hr_sum_df[c, "p"] <- tmp_sum_cox$logtest["pvalue"]
		c <- c+1
	}
}
hr_sum_df <- hr_sum_df %>% add_significance("p")
hr_sum_df$hr <- as.numeric(hr_sum_df$hr)
hr_sum_df$hci95 <- as.numeric(hr_sum_df$hci95)
hr_sum_df$lci95 <- as.numeric(hr_sum_df$lci95)
hr_sum_df$cate_hr <- ifelse(hr_sum_df$hr > 1, ifelse(hr_sum_df$lci95 > 1, "Mean>1,LCI>1", "Mean>1,LCI<1"), ifelse(hr_sum_df$hci95<1, "Mean<1,HCI<1", "Mean<1,HCI>1"))
write.csv(hr_sum_df, paste(data_dir, surv_folder, "/", kr_type, "_hr_between_knee_replacement_and_outcomes.csv", sep = ""))
print(head(hr_sum_df))

plot_df <- hr_sum_df[!is.na(hr_sum_df$hr),]
plot_df <- plot_df[plot_df$ytyr != "Y10",]

padj_alpha <- c("ns" = 0.05, "*" = 0.5, "**" = 0.75, "***" = 0.9, "****" = 1.0)
padj_size <- c("ns" = 0.1, "*" = 4, "**" = 5, "***" = 6, "****" = 7)

dot_gg <- ggplot(plot_df, aes(x = ytyr, y = outcome_name)) +
	geom_point(aes(color = hr, size = p.signif, alpha = p.signif)) +
	scale_colour_gradient2(low="darkorange", high="forestgreen", mid = "cornsilk", midpoint = 1.0) +
	facet_grid(.~kr_side) +
	scale_size_manual(values = padj_size) +
	scale_alpha_manual(values = padj_alpha) +
	labs(x = "Time", y = "Outcomes", color = "Hazard ratios", size = "Statistical significance\n(p-value)", alpha = "Statistical significance\n(p-value)") +
	theme_bw()
ggsave(paste(data_dir, surv_folder, "/", kr_type, "_logrank_hr_dotplot.png", sep = ""), dot_gg, dpi = png_res, width = 12, height = 6)
dot_gg <- ggplot(plot_df, aes(x = ytyr, y = outcome_name)) +
	geom_point(aes(color = cate_hr, size = p.signif, alpha = p.signif)) +
	scale_colour_manual(values = c("Mean>1,LCI>1" = "deepskyblue4", "Mean>1,LCI<1" = "deepskyblue1", "Mean<1,HCI>1" = "darkorange1", "Mean<1,HCI<1"="darkorange4")) +
	facet_grid(.~kr_side) +
	scale_size_manual(values = padj_size) +
	scale_alpha_manual(values = padj_alpha) +
	labs(x = "Time", y = "Outcomes", color = "Hazard ratios", size = "Statistical significance\n(p-value)", alpha = "Statistical significance\n(p-value)") +
	theme_bw()
ggsave(paste(data_dir, surv_folder, "/", kr_type, "_logrank_cate_hr_dotplot.png", sep = ""), dot_gg, dpi = png_res, width = 12, height = 6)

