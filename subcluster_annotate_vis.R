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
cluster_num <- 4
scale_top <- 4
kr_type <- "total_lastfollowup"

png_res <- 600
top_m <- 10

names <- c("C0" = "High pain/depression", "C1" = "High strength", "C2" = "High depression", "C3" = "High nutrition")
name_order <- c("High nutrition", "High pain/depression", "High depression", "High strength")

general_names <- c("G0" = "Unhealthy diet", "G1" = "Poor knee & general health", "G2" = "Good knee & general health", "G3" = "Intermediate knee & general health")
general_name_order <- c("Poor knee & general health", "Intermediate knee & general health", "Good knee & general health", "Unhealthy diet")

demographic_flag <- FALSE
outcome_table_flag <- FALSE
cluster_annotation_flag <- FALSE
one2one_annotation_flag <- FALSE
numeric_vis_marker_flag <- TRUE
categorical_vis_marker_flag <- FALSE
umap_vis_flag <- FALSE
violin_vis_flag <- FALSE
volcano_flag <- FALSE
supplement_flag <- FALSE

cat("Reading python output results...\n")
general_umap_df <- as.data.frame(read_excel("/mnt/sda1/OAI_Data/kmean_cluster_12252020/clean_v25_cluster4_kmean_pca_umap_res.xlsx"))
general_umap_df$gcluster <- str_c("G", general_umap_df$kmean_pca)
general_umap_df$gname <- "NA"
print(head(general_umap_df))
for (ic in unique(general_umap_df$gcluster)) {
	general_umap_df$gname[general_umap_df$gcluster == ic] <- general_names[ic]
}
#general_umap_df$gname <- factor(general_umap_df$gname, levels = c("Poor knee & general health", "Intermediate knee & general health", "Good knee & general health", "Unhealthy diet"))

umap_df <- as.data.frame(read_excel(paste(input_prefix, "cluster", cluster_num, "_kmean_pca_umap_res.xlsx", sep = "")))
rownames(umap_df) <- umap_df$ID
umap_df$name <- str_c("C", umap_df$kmean_pca)
for (ic in unique(umap_df$name)) {
	umap_df$name[umap_df$name == ic] <- names[ic]
}
#umap_df$name <- factor(umap_df$name, levels = c("High nutrition", "High pain/depression", "High depression", "High strength"))

data_df <- as.data.frame(read.csv(paste(data_dir, "input_direct_merge_dataframe_",input_id,"_sbj_clean_",clean_co,".csv", sep = ""), 
				  sep = ",", header = TRUE, row.names=1, stringsAsFactors = FALSE))
data_df <- data_df[rownames(umap_df),]
var_types <- sapply(data_df, class)

prog_files <- list.files(data_dir, pattern = 'survival_ready_results.csv')
var_coding <- as.data.frame(read_excel(paste(data_dir, "Variable_coding_v2.xlsx", sep = "")))

kr_df <- read.csv(paste(data_dir, kr_type, "_merge_patient_basic_outcome_information.csv", sep = ""), header = T)
rownames(kr_df) <- kr_df$ID

##### UMAP overlay with clusters [Fig 2A]
umap_df$name <- factor(umap_df$name, levels = c("High nutrition", "High pain/depression", "High depression", "High strength"))
umap_gg <- ggplot(umap_df, aes_string(x = "UMAP1", y = "UMAP2", color = "name")) + 
	geom_point(size = 1) +
	scale_color_brewer(palette = "Spectral") +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	labs(color='Cluster') +
	theme_classic()
ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_umap_r.png', sep = ""), umap_gg, 
       dpi = png_res, width = 7, height = 3.5)
umap_df$name <- as.character(umap_df$name)
general_umap_df <- merge(general_umap_df, umap_df, by = "ID", all.x = T)
print(head(general_umap_df))

general_umap_df$sub_cluster <- ifelse(is.na(general_umap_df$name), "Others", general_umap_df$name)
general_umap_df$sub_cluster <- factor(general_umap_df$sub_cluster, levels = c("High nutrition", "High pain/depression", "High depression", "High strength", "Others"))

general_umap_df$all_cluster <- ifelse(is.na(general_umap_df$name), general_umap_df$gname, general_umap_df$name)
general_umap_df$all_cluster <- factor(general_umap_df$all_cluster, 
				      levels = c("High nutrition", "High pain/depression", "High depression", "High strength",
						 "Poor knee & general health", "Intermediate knee & general health", "Good knee & general health"))


print(head(general_umap_df$name[!is.na(general_umap_df$name)]))

umap_gg <- ggplot(general_umap_df, aes_string(x = "UMAP1.x", y = "UMAP2.x", color = "sub_cluster")) + 
	geom_point(size = 0.5) +
	scale_color_brewer(palette = "Spectral") +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	labs(color='Cluster', x = "UMAP1", y = "UMAP2") +
	theme_classic()
ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_umap_r_sub_cluster_on_big_umap.png', sep = ""), umap_gg, 
       dpi = png_res, width = 7, height = 3.5)

umap_gg <- ggplot(general_umap_df, aes_string(x = "UMAP1.x", y = "UMAP2.x", color = "all_cluster")) + 
	geom_point(size = 0.5) +
	scale_color_brewer(palette = "Set1") +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	labs(color='Cluster', x= "UMAP1", y = "UMAP2") +
	theme_classic()
ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_umap_r_all_cluster_on_big_umap.png', sep = ""), umap_gg, 
       dpi = png_res, width = 8, height = 3.5)

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

if (demographic_flag) {
	print(data_df[1:9, 1:6])
	print(colnames(data_df)[str_detect(colnames(data_df), "INCOME")])
	print(colnames(data_df)[str_detect(colnames(data_df), "AGE")])
	print(colnames(data_df)[str_detect(colnames(data_df), "RACE")])
	print(colnames(data_df)[str_detect(colnames(data_df), "EDC")])
	print(colnames(data_df)[str_detect(colnames(data_df), "EMPLOY")])
	print(colnames(data_df)[str_detect(colnames(data_df), "BMI")])
	print(colnames(data_df)[str_detect(colnames(data_df), "SEX")])
	print(colnames(data_df)[str_detect(colnames(data_df), "COHORT")])
	print(colnames(data_df)[str_detect(colnames(data_df), "COMORB")])


	demo_df <- data_df[,c("V00AGE",  "P02SEX", "P02RACE", "P01BMI", "V00EDCV", "V00INCOME", "V00CEMPLOY", "V00COMORB")]
	colnames(demo_df) <- c("Age", "Sex", "Race", "BMI", "Education", "Income", "Employment", "Comorbidity")
	demo_df$Cohort <- kr_df[rownames(demo_df), "V00COHORT"]
	table_df <- merge(umap_df, demo_df, by = "row.names")
	demo_table <- table1(~Age+Sex+Race+BMI+Education+Income+Employment+Comorbidity+Cohort | name, data = table_df)
	write.csv(as.data.frame(demo_table), paste(output_prefix, "cluster", cluster_num, "_demographic_format_table.csv", sep = ""))
}

###### Generate the outcome table [Table S2]
if (outcome_table_flag) {
	print(dim(kr_df))
	print(colnames(kr_df))
	print(head(kr_df[,200:ncol(kr_df)]))
	print(unique(kr_df[,"left_kr"]))
	cate_kr_df <- kr_df[,c("left_time", "left_kr", "right_time", "right_kr")]
	cate_kr_df$left_y04 <- ifelse(cate_kr_df$left_time <= 365*4+1, cate_kr_df$left_kr, 
				      ifelse(str_detect(cate_kr_df$left_kr, "Missing"), ".: Missing Form/Incomplete Workbook", "3: No"))
	cate_kr_df$left_y08 <- ifelse(cate_kr_df$left_time <= 365*8+2, cate_kr_df$left_kr, 
				      ifelse(str_detect(cate_kr_df$left_kr, "Missing"), ".: Missing Form/Incomplete Workbook", "3: No"))

	cate_kr_df$right_y04 <- ifelse(cate_kr_df$right_time <= 365*4+1, cate_kr_df$right_kr, 
				      ifelse(str_detect(cate_kr_df$right_kr, "Missing"), ".: Missing Form/Incomplete Workbook", "3: No"))
	cate_kr_df$right_y08 <- ifelse(cate_kr_df$right_time <= 365*8+2, cate_kr_df$right_kr, 
				      ifelse(str_detect(cate_kr_df$right_kr, "Missing"), ".: Missing Form/Incomplete Workbook", "3: No"))


	print(head(cate_kr_df))
	print(unique(cate_kr_df[,"left_y04"]))

	out_tbl_df <- cate_kr_df[,str_detect(colnames(cate_kr_df), "y0")]
	colnames(out_tbl_df) <- c("Year4_Left", "Year8_Left", "Year4_Right", "Year8_Right")
	print(head(out_tbl_df))

	table_df <- merge(umap_df, out_tbl_df, by = "row.names")
	kr_table <- table1(~ Year4_Left+Year8_Left+Year4_Right+Year8_Right | name, data = table_df)
	write.csv(as.data.frame(kr_table), paste(input_prefix, "tkr_format_table.csv", sep = ""))
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

##### Use Heatmap/dot plot to visualize the top feature markers [Fig. 2B&C]
marker_df <- read.csv(paste(output_prefix, "cluster", cluster_num,'_kmeans_direct_knn2imp_', input_id, '_marker_df_largeB.csv', sep = ""), row.names = 1)
marker_df$name <- str_c("C", marker_df$cluster)

if (categorical_vis_marker_flag) { # [Fig 2B]
	cat("Dot plot for categorical markers...\n")
	cat_marker_df <- marker_df[marker_df$type != 'numeric',]
#	cat_marker_df <- cat_marker_df[cat_marker_df$logfc >=0,]
	cat_marker_df$logpadj <- -log10(cat_marker_df$adjp)
	for (ic in unique(cat_marker_df$name)) {
		cat_marker_df$name[cat_marker_df$name == ic] <- names[ic]
	}

	cat_marker_df$cluster <- as.factor(cat_marker_df$cluster)
	cat_marker_df <- cat_marker_df[!is.na(cat_marker_df$cluster),]
	top_cat_marker_df <- cat_marker_df %>%
		group_by(cluster) %>%
		filter(adjp <= 0.10) %>%
		top_n(top_m,logfc)
	plot_top_cat_marker_df <- cat_marker_df[cat_marker_df$variable %in% unique(top_cat_marker_df$variable),]

	write.csv(cat_marker_df, paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_all_cat_markers.csv', sep = ""))
	co <- 1:6
	for (ico in co) {
		sig_cat_marker_df <- cat_marker_df %>% group_by(cluster) %>% filter(adjp <= 10^-ico)
		dengg <- ggplot(sig_cat_marker_df, aes(x = logfc)) +
			geom_density(aes(fill = name)) +
			scale_fill_brewer(palette = "Spectral") +
			facet_wrap(~name, ncol = 1) +
			labs(fill = "Cluster", x = "Cramer's V") +
			theme_bw()
		ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_adjp_co_10mns', ico, '_cat_markers_logfc_dens.png', sep = ""), dpi = png_res,
		       width = 9, height = 9)
	
		sig_cat_marker_df <- cat_marker_df %>% group_by(cluster) %>% filter(adjp <= 10^-ico, abs(logfc) >= 0.375)
		write.csv(sig_cat_marker_df, paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_adjp_co_10mns', ico, '_cat_markers.csv', sep = ""))
#		print(head(sig_cat_marker_df))
#		q(save = "no")
	}


	plot_top_cat_marker_df$var_name <- plot_top_cat_marker_df$variable
	for (ivn in 1:nrow(var_coding)) {
		if (var_coding[ivn, "VID"] %in% plot_top_cat_marker_df$variable) {
			tmp_mask <- plot_top_cat_marker_df$variable == var_coding[ivn, "VID"]
			plot_top_cat_marker_df$var_name[tmp_mask] <- var_coding[ivn, "Name"]
		}
	}

	var_levels <- c()
	for (inm in unique(plot_top_cat_marker_df$name)) {
		tmp_df <- plot_top_cat_marker_df[plot_top_cat_marker_df$name == inm,]
		tmp_levels <- tmp_df$var_name[order(abs(tmp_df$logfc), decreasing = T)][1:top_m]
		var_levels <- c(var_levels, tmp_levels[!(tmp_levels %in% var_levels)])
	}

	plot_top_cat_marker_df <- plot_top_cat_marker_df[plot_top_cat_marker_df$var_name %in% var_levels,]
	plot_top_cat_marker_df$var_name <- factor(plot_top_cat_marker_df$var_name, levels = var_levels)
	plot_top_cat_marker_df$name <- factor(plot_top_cat_marker_df$name)#, levels = name_order)
	dot_mgg <- ggplot(plot_top_cat_marker_df, aes(y = name, x = var_name)) +
		geom_point(aes(color = logfc, size = logpadj)) +
		scale_color_continuous_divergingx(palette = 'RdBu', mid = 0.2, p3 = 1, p4 = 1) + 
#		scale_color_continuous(palette = "Spectral") +
		labs(color = "Cramer\'s V", size = '-log10\nadjusted p-value', y = "Cluster", x = "Key categorical variables") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_top', top_m, '_cate_markers.png', sep = ""), 
	       dot_mgg, dpi = png_res, width = 10, height = 6)
}


if (numeric_vis_marker_flag) { # [Fig 2C]
	cat("Dot plot for numeric markers...\n")
	num_marker_df <- marker_df[marker_df$type == 'numeric',]
#	num_marker_df <- num_marker_df[!is.na(num_marker_df$logfc),]
#	num_marker_df <- num_marker_df[num_marker_df$logfc >=0,]
	num_marker_df$logpadj <- -log10(num_marker_df$adjp)
	print(head(num_marker_df))
	for (ic in unique(num_marker_df$name)) {
		num_marker_df$name[num_marker_df$name == ic] <- names[ic]
	}

	num_marker_df$cluster <- as.factor(num_marker_df$cluster)
	num_marker_df <- num_marker_df[!is.na(num_marker_df$cluster),]
	top_num_marker_df <- num_marker_df %>%
		group_by(cluster) %>%
		filter(adjp <= 0.10) %>%
		top_n(top_m,logfc)
	plot_top_num_marker_df <- num_marker_df[num_marker_df$variable %in% unique(top_num_marker_df$variable),]
	write.csv(num_marker_df, paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_all_num_markers.csv', sep = ""))
	print(unique(top_num_marker_df$variable))
	order_top_num_marker_df <- top_num_marker_df %>% 
		group_by(variable) %>% 
		mutate(logfc_order <- order(logfc))
	write.csv(order_top_num_marker_df, paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_top_used_num_markers.csv', sep = ""))
	top_num_marker_df$name <- factor(top_num_marker_df$name, levels = name_order)
	var_order <- top_num_marker_df$variable[order(top_num_marker_df$name, top_num_marker_df$logfc)]
	rownames(var_coding) <- var_coding[,1]
	order_var_coding <- var_coding[var_order,]
	print(order_var_coding)
	q(save = "no")

	sub_data_df <- data_df[,unique(top_num_marker_df$variable)]
	sub_data_df$ID <- rownames(sub_data_df)
	merge_sub_data_df <- merge(sub_data_df, umap_df, by = "ID", all = F)
	gath_plot_df <- gather(merge_sub_data_df, "marker", "value", unique(top_num_marker_df$variable))
	gath_plot_df <- merge(gath_plot_df, top_num_marker_df[,c("variable", "name")], by.x = "marker", by.y = "variable", all.x = T, suffix = c("", "_feature"))
	gath_plot_df$highlight <- ifelse(gath_plot_df$name == gath_plot_df$name_feature, "Y", "N")
	gath_plot_df$name <- factor(gath_plot_df$name, levels = name_order)
	gath_plot_df$marker <- factor(gath_plot_df$marker, levels = var_order)
	vln_gg <- ggplot(gath_plot_df, aes(x = name, y = value, fill = name)) +
		geom_violin(aes(alpha = highlight, color = highlight)) +
		scale_alpha_manual(values = c("Y" = 1.0, "N" = 0.25)) +
		scale_color_manual(values = c("Y" = "black", "N" = "gray")) +
		facet_grid(.~marker, scales = "free") +
		labs(x = "Subclusters", y = "Values") +
		coord_flip() +
		theme_bw() +
		theme(legend.position = "none", 
		      axis.text.x = element_blank(),
		      strip.text.x = element_text(angle = 90)
		)
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_top', top_m, '_num_violins.png', sep = ""), 
	       vln_gg, dpi = png_res, width = 24, height = 3)
	co <- 1:6
	for (ico in co) {
		sig_num_marker_df <- num_marker_df %>% group_by(cluster) %>% filter(adjp <= 10^-ico)
		dengg <- ggplot(sig_num_marker_df, aes(x = logfc)) +
			geom_density(aes(fill = name)) +
			scale_fill_brewer(palette = "Spectral") +
			facet_wrap(~name, ncol = 1) +
			labs(fill = "Cluster", x = "log2 Fold changes") +
			theme_bw()
		ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_adjp_co_10mns', ico, '_num_markers_logfc_dens.png', sep = ""), dpi = png_res,
		       width = 9, height = 9)

		sig_num_marker_df <- num_marker_df %>% group_by(cluster) %>% filter(adjp <= 10^-ico, abs(logfc) >= 0.5)
		write.csv(sig_num_marker_df, paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_adjp_co_10mns', ico, '_num_markers.csv', sep = ""))
#		print(head(sig_num_marker_df))
	}

	plot_top_num_marker_df$var_name <- plot_top_num_marker_df$variable
	for (ivn in 1:nrow(var_coding)) {
		if (var_coding[ivn, "VID"] %in% plot_top_num_marker_df$variable) {
			tmp_mask <- plot_top_num_marker_df$variable == var_coding[ivn, "VID"]
			plot_top_num_marker_df$var_name[tmp_mask] <- var_coding[ivn, "Name"]
		}
	}

	var_levels <- c()
	for (inm in unique(plot_top_num_marker_df$name)) {
		tmp_df <- plot_top_num_marker_df[plot_top_num_marker_df$name == inm,]
		tmp_levels <- tmp_df$var_name[order(abs(tmp_df$logfc), decreasing = T)][1:top_m]
		var_levels <- c(var_levels, tmp_levels[!(tmp_levels %in% var_levels)])
	}

	plot_top_num_marker_df <- plot_top_num_marker_df[plot_top_num_marker_df$var_name %in% var_levels,]
	plot_top_num_marker_df$var_name <- factor(plot_top_num_marker_df$var_name, levels = var_levels)
	plot_top_num_marker_df$name <- factor(plot_top_num_marker_df$name, levels = name_order)
	dot_mgg <- ggplot(plot_top_num_marker_df, aes(y = name, x = var_name)) +
		geom_point(aes(color = logfc, size = logpadj)) +
		scale_color_continuous_divergingx(palette = 'RdBu', mid = 0.0) + 
		scale_size(breaks = c(1.0, 2.0, 3.0, 4.0)) +
		labs(color = 'log2FC', size = '-log10\nadjusted p-value', y = "Cluster", x = "Key numeric variables") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_top', top_m, '_num_markers.png', sep = ""), 
	       dot_mgg, dpi = png_res, width = 10, height = 6)
}


