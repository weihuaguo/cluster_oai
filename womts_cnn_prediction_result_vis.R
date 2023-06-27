
rm(list = ls())
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))



data_dir<-"/home/weihua/mnts/smb_plee/Group/weihua/data_summary"
importance_number <- 25
all_folders <- list.dirs(data_dir, full.names = F, recursive = F)
test_id <- c(0:9999)
top_n <- seq(10, 1000, 10)
#print(all_folders)

importance_source <- "random"
importance_filter <- "cluster"

merge_imp_df <- read.csv(paste(data_dir, "/merge_select_importance_df_test_max_", max(test_id), "_230202.csv", sep = ""))
merge_score_df <- read.csv(paste(data_dir, "/merge_score_df_test_max_", max(test_id), "_230202.csv", sep = ""))
print(head(merge_imp_df))
print(dim(merge_imp_df))
print(head(merge_score_df))
print(dim(merge_score_df))
print(unique(merge_score_df$outname))

avg_score_df <- merge_score_df %>%
	group_by(outname, test_id) %>%
	summarize(avg_rmse = mean(rmse), avg_r = mean(r))
write.csv(avg_score_df, paste(data_dir, "/average_test_max_", max(test_id), "_vx_score_df.csv", sep = ""))

final_imp_df <- merge(merge_imp_df, avg_score_df, by = c("outname", "test_id"), all.x = T)
write.csv(final_imp_df, paste(data_dir, "/final_test_max_", max(test_id), "_vx_importance_score_df.csv", sep = ""))

print(head(final_imp_df))

final_imp_df <- final_imp_df %>%
	group_by(outname, varname) %>%
	mutate(var_total_cts = n())

ic <- 0
for (itn in top_n) {
	cat(itn, "\n")
	tmp_top_df <- final_imp_df %>%
		group_by(outname) %>%
		slice_min(order_by = avg_rmse, n = itn*25)
#	print(tmp_top_df)
	tmp_cts_df <- tmp_top_df %>%
		group_by(outname, varname) %>%
		mutate(var_cts = n(), var_cts2test = n()/itn*100)
	tmp_cts_df$var_cts2occ = tmp_cts_df$var_cts/tmp_cts_df$var_total_cts*100
	tmp_cts_df$topn <- itn
#	print(head(as.data.frame(tmp_cts_df)))
	if (ic == 0) {
		merge_cts_df <- tmp_cts_df
	} else {
		merge_cts_df <- rbind(merge_cts_df, tmp_cts_df)
	}
	ic <- ic + 1
}

cat("To TOPN...\n")
dlgg <- ggplot(merge_cts_df, aes(x = topn, y = var_cts2test, color = varname, group = varname)) +
	geom_line() +
	geom_point() +
	facet_wrap(~outname, ncol = 1) +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2topn_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 48, height = 16)

top_cts_df <- merge_cts_df %>%
	group_by(varname, outname) %>%
	summarize(avg = mean(var_cts2test), n = n())
write.csv(top_cts_df, paste(data_dir, "/top_test_vs_cts2topn_test_max_", max(test_id), "_average.csv", sep = ""))
top_cts_df <- top_cts_df %>%
	group_by(outname) %>%
	slice_max(order_by = avg, n = 5)
print(top_cts_df)

dlgg <- ggplot(merge_cts_df[merge_cts_df$varname %in% unique(top_cts_df$varname),], aes(x = topn, y = var_cts2test, color = varname, group = varname)) +
	geom_line() +
	geom_point(size = 1) +
	facet_wrap(~outname, ncol = 1) +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2topn_top5_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 9, height = 9)

merge_cts_df$general_var <- str_split_fixed(merge_cts_df$varname, "_", n = 2)[,1]
sum_cts_df <- merge_cts_df %>%
	group_by(general_var, outname, topn) %>%
	summarise(n = n(), avg_cts2test = mean(var_cts2test), sd_cts2test = sd(var_cts2test))
print(head(as.data.frame(merge_cts_df)))
print(sum_cts_df)

dlgg <- ggplot(sum_cts_df, aes(x = topn, y = avg_cts2test, color = general_var, group = general_var)) +
	geom_line() +
	geom_point() +
	facet_wrap(~outname, ncol = 1) +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2topn_general_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 32, height = 16)

top_cts_df <- sum_cts_df %>%
	group_by(general_var, outname) %>%
	summarize(avg = mean(avg_cts2test), n = n())
write.csv(top_cts_df, paste(data_dir, "/top_test_vs_cts2topn_general_test_max_", max(test_id), "_average.csv", sep = ""))

top_cts_df <- top_cts_df %>%
	group_by(outname) %>%
	slice_max(order_by = avg, n = 5)
print(top_cts_df)

dlgg <- ggplot(sum_cts_df[sum_cts_df$general_var %in% unique(top_cts_df$general_var),], aes(x = topn, y = avg_cts2test, color = general_var, group = general_var)) +
	geom_line() +
	geom_point(size = 1) +
	scale_color_brewer(palette = "Paired") +
	labs(x = "Top N tests based on RMSE", y = "Relative occurence to top N tests", color = "Top 5 variables") +
	facet_wrap(~outname, ncol = 1) +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2topn_general_top5_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 9, height = 9)

cat("To TOTAL OCC...\n")
dlgg <- ggplot(merge_cts_df, aes(x = topn, y = var_cts2occ, color = varname, group = varname)) +
	geom_line() +
	geom_point() +
	facet_wrap(~outname, ncol = 1) +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2occ_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 48, height = 16)

top_cts_df <- merge_cts_df %>%
	group_by(varname, outname) %>%
	summarize(avg = mean(var_cts2occ), n = n())
write.csv(top_cts_df, paste(data_dir, "/top_test_vs_cts2occ_test_max_", max(test_id), "_average.csv", sep = ""))

top_cts_df <- top_cts_df %>%
	group_by(outname) %>%
	slice_max(order_by = avg, n = 5)
print(top_cts_df)

dlgg <- ggplot(merge_cts_df[merge_cts_df$varname %in% unique(top_cts_df$varname),], aes(x = topn, y = var_cts2occ, color = varname, group = varname)) +
	geom_line() +
	geom_point(size = 1) +
	facet_wrap(~outname, ncol = 1) +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2occ_top5_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 9, height = 9)

merge_cts_df$general_var <- str_split_fixed(merge_cts_df$varname, "_", n = 2)[,1]
# NOTE: SUM not AVERAGE
sum_cts_df <- merge_cts_df %>%
	group_by(general_var, outname, topn) %>%
	summarise(var_sum = sum(var_cts), var_occ_sum = sum(var_total_cts))
sum_cts_df$avg_cts2occ <- sum_cts_df$var_sum/sum_cts_df$var_occ_sum*100
print(head(as.data.frame(merge_cts_df)))
print(sum_cts_df)

dlgg <- ggplot(sum_cts_df, aes(x = topn, y = avg_cts2occ, color = general_var, group = general_var)) +
	geom_line() +
	geom_point() +
	facet_wrap(~outname, ncol = 1) +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2occ_general_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 32, height = 16)

top_cts_df <- sum_cts_df %>%
	group_by(general_var, outname) %>%
	summarize(avg = mean(avg_cts2occ), n = n())
write.csv(top_cts_df, paste(data_dir, "/top_test_vs_cts2occ_general_test_max_", max(test_id), "_average.csv", sep = ""))

top_cts_df <- top_cts_df %>%
	group_by(outname) %>%
	slice_max(order_by = avg, n = 5)
print(top_cts_df)

dlgg <- ggplot(sum_cts_df[sum_cts_df$general_var %in% unique(top_cts_df$general_var),], aes(x = topn, y = avg_cts2occ, color = general_var, group = general_var)) +
	geom_line() +
	geom_point(size = 1) +
	scale_color_brewer(palette = "Paired") +
	facet_wrap(~outname, ncol = 1) +
	labs(x = "Top N tests based on RMSE", y = "Relative occurence to total test numbers", color = "Top 5 variables") +
	theme_bw()
ggsave(paste(data_dir, "/top_test_vs_cts2occ_general_top5_test_max_", max(test_id), "_dotline.png", sep = ""), dpi = 300, width = 9, height = 9)
