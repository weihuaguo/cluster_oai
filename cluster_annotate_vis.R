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


clst_dir <- "/mnt/sda1/OAI_Data/kmean_cluster_12252020/"
data_dir <- "/mnt/sda1/OAI_Data/data_summary/"

exp_id <- "kmean_pca_umap"
input_id <- "12112020_data_clean"
clean_co <- 'v25'
input_prefix <- paste(clst_dir, 'clean_', clean_co, '_', sep = "")
output_prefix <- paste(input_prefix, "final_", sep = "")
cluster_num <- 4
scale_top <- 4
kr_type <- "total_lastfollowup"

png_res <- 600
top_m <- 10

names <- c("C0" = "Low supplemental vitamins", "C1" = "Poor knee & general health", "C2" = "Good knee & general health", "C3" = "Intermediate knee & general health")
name_order <- c("Low supplemental vitamins", "Poor knee & general health", "Intermediate knee & general health", "Good knee & general health")
names <- c("C0" = "Unhealthy diet", "C1" = "Poor knee & general health", "C2" = "Good knee & general health", "C3" = "Intermediate knee & general health")
name_order <- c("Unhealthy diet", "Poor knee & general health", "Intermediate knee & general health", "Good knee & general health")


demographic_flag <- FALSE
outcome_table_flag <- FALSE
cluster_annotation_flag <- FALSE
one2one_annotation_flag <- FALSE
numeric_vis_marker_flag <- FALSE
categorical_vis_marker_flag <- FALSE
umap_vis_flag <- FALSE
violin_vis_flag <- FALSE
volcano_flag <- FALSE
supplement_flag <- TRUE
life_act_flag <- FALSE
norm_flag <- FALSE

cat("Reading python output results...\n")
umap_df <- as.data.frame(read_excel(paste(input_prefix, "cluster", cluster_num, "_kmean_pca_umap_res.xlsx", sep = "")))
rownames(umap_df) <- umap_df$ID
umap_df$name <- str_c("C", umap_df$kmean_pca)
for (ic in unique(umap_df$name)) {
	umap_df$name[umap_df$name == ic] <- names[ic]
}
umap_df$name <- factor(umap_df$name, levels = c("Poor knee & general health", "Intermediate knee & general health", "Good knee & general health", "Low supplemental vitamins"))

data_df <- as.data.frame(read.csv(paste(data_dir, "input_direct_merge_dataframe_",input_id,"_sbj_clean_",clean_co,".csv", sep = ""), 
				  sep = ",", header = TRUE, row.names=1, stringsAsFactors = FALSE))
var_types <- sapply(data_df, class)

metric_df <- as.data.frame(read_excel(paste(input_prefix, 'kmeans_metric_result_direct_knn2imp_', input_id, '.xlsx', sep = "")))
prog_files <- list.files(data_dir, pattern = 'survival_ready_results.csv')
var_coding <- as.data.frame(read_excel(paste(data_dir, "Variable_coding_v2.xlsx", sep = "")))

kr_df <- read.csv(paste(data_dir, kr_type, "_merge_patient_basic_outcome_information.csv", sep = ""), header = T)
rownames(kr_df) <- kr_df$ID

##### UMAP overlay with clusters [Fig 2A]
umap_df$cluster <- str_c("C", umap_df$kmean_pca)
umap_gg <- ggplot(umap_df, aes_string(x = "UMAP1", y = "UMAP2", color = "name")) + 
	geom_point(size = 0.1) +
	scale_color_brewer(palette = "Spectral") +
	guides(color = guide_legend(override.aes = list(size = 3))) +
	labs(color='Cluster') +
	theme_classic()
ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_umap_r.png', sep = ""), umap_gg, 
       dpi = png_res, width = 7.5, height = 5)

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

###### Generate the demographic table [Table S1]
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
	write.csv(as.data.frame(demo_table), paste(input_prefix, "demographic_format_table.csv", sep = ""))
}

###### Compare each 2 clusters [Table S3]
if (one2one_annotation_flag) {
	cat("Generate marker table...\n")
	marker_types <- sapply(data_df, class)
	marker_cols <- c('variable', 'cluster1', 'cluster2', 'type', 'test','avg1', 'avg2', 'diff', 'logfc', 'per1', 'per2', 'pval')
	marker_df <- as.data.frame(matrix(nrow=(cluster_num*ncol(data_df)), ncol=length(marker_cols)))
	colnames(marker_df) <- marker_cols

	print(unique(umap_df$cluster))
	cluster_comb <- combn(unique(umap_df$cluster), 2)
	print(cluster_comb)
#	q(save = "no")

	for (i in 1:ncol(cluster_comb)) {
		cat("\t", cluster_comb[1,i], "vs", cluster_comb[2,i], '\n')
		ic <- 1
		marker_cols <- c('variable', 'cluster1', 'cluster2', 'type', 'test','avg1', 'avg2', 'diff', 'logfc', 'per1', 'per2', 'pval')
		marker_df <- as.data.frame(matrix(nrow=(cluster_num*ncol(data_df)), ncol=length(marker_cols)))
		colnames(marker_df) <- marker_cols
		marker_df$cluster1 <- cluster_comb[1,i]
		marker_df$cluster2 <- cluster_comb[2,i]

		tmp_umap_df <- umap_df[umap_df$cluster %in% cluster_comb[,i],]
		print(unique(tmp_umap_df$cluster))
		tmp_data_df <- data_df[rownames(tmp_umap_df),]

		tmp_umap_df$comp <- ifelse(tmp_umap_df$cluster == cluster_comb[1,i], 'COI', 'Others')
		for (im in names(marker_types)) {
	#		print(marker_types[[im]])
			marker_df[ic, 'variable'] <- im
#			marker_df[ic, 'cluster'] <- i
			marker_df[ic, 'type'] <- marker_types[[im]]
			tmp_df <- cbind(tmp_umap_df, tmp_data_df[,im])
			tmp_force_num <- as.numeric(as.character(tmp_df[,ncol(tmp_df)]))
			if (sum(!is.na(tmp_force_num))>0) {
				marker_types[[im]] <- 'numeric'
				tmp_df[,ncol(tmp_df)] <- tmp_force_num
			}

			if (marker_types[[im]]=='numeric'|marker_types[[im]]=='integer'){
	#			cat("\t\tKruskal-Wallis Test...\n")
				marker_df[ic, 'test'] <- 'Kruskal-Wallis test'
	#			print(head(tmp_df))
				tmp_sum <- tmp_df %>%
					group_by(comp) %>%
					summarize(mean = mean(`tmp_data_df[, im]`, na.rm=T),
						  n=n(),
						  naper = sum(!is.na(`tmp_data_df[, im]`))/n())
				marker_df[ic, 'avg1'] <- tmp_sum[1,'mean']
				marker_df[ic, 'avg2'] <- tmp_sum[2, 'mean']
				marker_df[ic, 'diff'] <- tmp_sum[1, 'mean']-tmp_sum[2,'mean']
				marker_df[ic, 'logfc'] <- log2(tmp_sum[1, 'mean']/tmp_sum[2, 'mean'])
				marker_df[ic, 'per1'] <- tmp_sum[1, 'naper']
				marker_df[ic, 'per2'] <- tmp_sum[2, 'naper']
				tmp_test <- kruskal.test(`tmp_data_df[, im]`~comp, data = tmp_df)
				marker_df[ic, 'pval'] <- tmp_test$p.value
			} else {
	#			cat("\t\tFisher exact test...\n")
				tmp_table <- as.data.frame.matrix(table(tmp_df$comp, tmp_df[,ncol(tmp_df)]))
				tmp_prop_table <- prop.table(as.matrix(tmp_table),margin=1)
	#			print(im)
	#			print(tmp_prop_table)
	#			print(dim(tmp_table))
	#			print(tmp_table)
				if (sum(tmp_prop_table[1,]==max(tmp_prop_table[1,]))>1|sum(tmp_prop_table[2,]==max(tmp_prop_table[2,]))>1) {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])][1]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])][1]
				} else {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])]
				}
				cat(im, '\n')
				tmp_test <- fisher.test(tmp_table, simulate.p.value=TRUE, B=1e7)#, workspace=2e9)
				tmp_chitest <- chisq.test(tmp_table, correct=FALSE)
				marker_df[ic, 'logfc'] <- sqrt(tmp_chitest$statistic/sum(tmp_table))
				marker_df[ic, 'test'] <- 'Fisher\'s exact test'
				na_mask <- str_detect(colnames(tmp_table), 'Missing')
				if (sum(na_mask) == 1) {
					marker_df[ic, 'per1'] <- 1-tmp_table['COI', na_mask]/sum(tmp_table['COI',])
					marker_df[ic, 'per2'] <- 1-tmp_table['Others', na_mask]/sum(tmp_table['Others',])
				}
				marker_df[ic, 'pval'] <- tmp_test$p.value
			}
			ic <- ic + 1
		}
		marker_df$adjp <- p.adjust(marker_df$pval, method = "BH")
		write.csv(marker_df, paste(output_prefix, cluster_comb[1, i], "_vs_", cluster_comb[2,i],'_kmeans_direct_knn2imp_', input_id, '_marker_df_largeB.csv', sep = ""))
	}
#	marker_df$adjp <- p.adjust(marker_df$pval, method = "BH")
#	write.csv(marker_df, paste(output_prefix, "cluster", cluster_num,'_kmeans_direct_knn2imp_', input_id, '_marker_df_largeB.csv', sep = ""))
}

if (norm_flag) {
	var_ann <- as.data.frame(read_excel(paste(data_dir, "AllClinical00_V6_column_annotation.xlsx", sep = "/"), sheet = "Baseline"))
	print(head(var_ann))
	print(colnames(data_df))
	print(colnames(data_df)[str_detect(colnames(data_df), "WEIGHT")])
	var_ann <- var_ann[var_ann$Norm == "Yes",]
	use_data_df <- data_df[,c(var_ann$Variables, "P01BMI", "P01WEIGHT")]
	use_data_df <- use_data_df[rownames(umap_df),]
	norm_df <- as.data.frame(matrix(ncol = nrow(var_ann)*2, nrow = nrow(use_data_df)))
	rownames(norm_df) <- rownames(use_data_df)
	c <- 1
	for (i in var_ann$Variables) {
		cat(i, "\n")
		colnames(norm_df)[c] <- str_c(i,"_BMI")
		norm_df[,c] <- use_data_df[,i]/use_data_df[,"P01BMI"]
		c <- c+1
		colnames(norm_df)[c] <- str_c(i,"_WEIGHT")
		norm_df[,c] <- use_data_df[,i]/use_data_df[,"P01WEIGHT"]
		c <- c+1
	}
	use_data_df <- norm_df
	marker_types <- sapply(use_data_df, class)
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
			tmp_df <- cbind(tmp_umap_df, use_data_df[,im])
			print(head(tmp_df))
			tmp_force_num <- as.numeric(as.character(tmp_df[,ncol(tmp_df)]))
			if (sum(!is.na(tmp_force_num))>0) {
				marker_types[[im]] <- 'numeric'
				tmp_df[,ncol(tmp_df)] <- tmp_force_num
			}
			tmp_df <- tmp_df[!is.na(tmp_df[,ncol(tmp_df)]),]

			if (marker_types[[im]]=='numeric'|marker_types[[im]]=='integer'){
	#			cat("\t\tKruskal-Wallis Test...\n")
				marker_df[ic, 'test'] <- 'Kruskal-Wallis test'
	#			print(head(tmp_df))
				tmp_sum <- tmp_df %>%
					group_by(comp) %>%
					summarize(mean = mean(`use_data_df[, im]`, na.rm=T),
						  n=n(),
						  naper = sum(!is.na(`use_data_df[, im]`))/n())
				marker_df[ic, 'avg1'] <- tmp_sum[1,'mean']
				marker_df[ic, 'avg2'] <- tmp_sum[2, 'mean']
				marker_df[ic, 'diff'] <- tmp_sum[1, 'mean']-tmp_sum[2,'mean']
				marker_df[ic, 'logfc'] <- log2(tmp_sum[1, 'mean']/tmp_sum[2, 'mean'])
				marker_df[ic, 'per1'] <- tmp_sum[1, 'naper']
				marker_df[ic, 'per2'] <- tmp_sum[2, 'naper']
				tmp_test <- kruskal.test(`use_data_df[, im]`~comp, data = tmp_df)
				marker_df[ic, 'pval'] <- tmp_test$p.value
			} else {
	#			cat("\t\tFisher exact test...\n")
				tmp_table <- as.data.frame.matrix(table(tmp_df$comp, tmp_df[,ncol(tmp_df)]))
				tmp_prop_table <- prop.table(as.matrix(tmp_table),margin=1)
	#			print(im)
	#			print(tmp_prop_table)
	#			print(dim(tmp_table))
	#			print(tmp_table)
				if (sum(tmp_prop_table[1,]==max(tmp_prop_table[1,]))>1|sum(tmp_prop_table[2,]==max(tmp_prop_table[2,]))>1) {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])][1]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])][1]
				} else {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])]
				}
				cat(im, '\n')
				tmp_test <- fisher.test(tmp_table, simulate.p.value=TRUE, B=1e7)#, workspace=2e9)
				tmp_chitest <- chisq.test(tmp_table, correct=FALSE)
				marker_df[ic, 'logfc'] <- sqrt(tmp_chitest$statistic/sum(tmp_table))
				marker_df[ic, 'test'] <- 'Fisher\'s exact test'
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
	write.csv(marker_df, paste(output_prefix, "cluster", cluster_num,'_kmeans_calories_normed_', input_id, '_marker_df_largeB.csv', sep = "")) # SD4
}

if (life_act_flag) {
	life_act_df <- as.data.frame(read.csv(paste(data_dir, "life_act_direct_merge_dataframe_",input_id,".csv", sep = ""), 
					  sep = ",", header = TRUE, row.names=1, stringsAsFactors = FALSE))
	print(head(life_act_df))
	var_ann <- as.data.frame(read_excel(paste(data_dir, "AllClinical00_V6_column_annotation.xlsx", sep = "/"), sheet = "life_activity_accelo"))
	print(head(var_ann))
	use_data_df <- life_act_df[,var_ann$Variables]
	use_data_df <- use_data_df[rownames(umap_df),]
	print(head(use_data_df))
	
	marker_types <- sapply(use_data_df, class)
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
			tmp_df <- cbind(tmp_umap_df, use_data_df[,im])
			print(head(tmp_df))
			tmp_force_num <- as.numeric(as.character(tmp_df[,ncol(tmp_df)]))
			if (sum(!is.na(tmp_force_num))>0) {
				marker_types[[im]] <- 'numeric'
				tmp_df[,ncol(tmp_df)] <- tmp_force_num
			}
			tmp_df <- tmp_df[!is.na(tmp_df[,ncol(tmp_df)]),]

			if (marker_types[[im]]=='numeric'|marker_types[[im]]=='integer'){
	#			cat("\t\tKruskal-Wallis Test...\n")
				marker_df[ic, 'test'] <- 'Kruskal-Wallis test'
	#			print(head(tmp_df))
				tmp_sum <- tmp_df %>%
					group_by(comp) %>%
					summarize(mean = mean(`use_data_df[, im]`, na.rm=T),
						  n=n(),
						  naper = sum(!is.na(`use_data_df[, im]`))/n())
				marker_df[ic, 'avg1'] <- tmp_sum[1,'mean']
				marker_df[ic, 'avg2'] <- tmp_sum[2, 'mean']
				marker_df[ic, 'diff'] <- tmp_sum[1, 'mean']-tmp_sum[2,'mean']
				marker_df[ic, 'logfc'] <- log2(tmp_sum[1, 'mean']/tmp_sum[2, 'mean'])
				marker_df[ic, 'per1'] <- tmp_sum[1, 'naper']
				marker_df[ic, 'per2'] <- tmp_sum[2, 'naper']
				tmp_test <- kruskal.test(`use_data_df[, im]`~comp, data = tmp_df)
				marker_df[ic, 'pval'] <- tmp_test$p.value
			} else {
	#			cat("\t\tFisher exact test...\n")
				tmp_table <- as.data.frame.matrix(table(tmp_df$comp, tmp_df[,ncol(tmp_df)]))
				tmp_prop_table <- prop.table(as.matrix(tmp_table),margin=1)
	#			print(im)
	#			print(tmp_prop_table)
	#			print(dim(tmp_table))
	#			print(tmp_table)
				if (sum(tmp_prop_table[1,]==max(tmp_prop_table[1,]))>1|sum(tmp_prop_table[2,]==max(tmp_prop_table[2,]))>1) {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])][1]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])][1]
				} else {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])]
				}
				cat(im, '\n')
				tmp_test <- fisher.test(tmp_table, simulate.p.value=TRUE, B=1e7)#, workspace=2e9)
				tmp_chitest <- chisq.test(tmp_table, correct=FALSE)
				marker_df[ic, 'logfc'] <- sqrt(tmp_chitest$statistic/sum(tmp_table))
				marker_df[ic, 'test'] <- 'Fisher\'s exact test'
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
	write.csv(marker_df, paste(output_prefix, "cluster", cluster_num,'_kmeans_life_act_', input_id, '_marker_df_largeB.csv', sep = "")) # SD4

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
				marker_df[ic, 'test'] <- 'Kruskal-Wallis test'
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
	#			print(tmp_table)
				if (sum(tmp_prop_table[1,]==max(tmp_prop_table[1,]))>1|sum(tmp_prop_table[2,]==max(tmp_prop_table[2,]))>1) {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])][1]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])][1]
				} else {
					marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])]
					marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])]
				}
				cat(im, '\n')
				tmp_test <- fisher.test(tmp_table, simulate.p.value=TRUE, B=1e7)#, workspace=2e9)
				tmp_chitest <- chisq.test(tmp_table, correct=FALSE)
				marker_df[ic, 'logfc'] <- sqrt(tmp_chitest$statistic/sum(tmp_table))
				marker_df[ic, 'test'] <- 'Fisher\'s exact test'
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
for (ic in unique(marker_df$name)) {
	marker_df$name[marker_df$name == ic] <- names[ic]
}
marker_df$name <- factor(marker_df$name, levels = c("Poor knee & general health", "Intermediate knee & general health", "Good knee & general health", "Low supplemental vitamins"))

if (categorical_vis_marker_flag) { # [Fig 2B]
	cat("Dot plot for categorical markers...\n")
	cat_marker_df <- marker_df[marker_df$type != 'numeric',]
#	cat_marker_df <- cat_marker_df[cat_marker_df$logfc >=0,]
	cat_marker_df$logpadj <- -log10(cat_marker_df$adjp)
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
	plot_top_cat_marker_df$name <- factor(plot_top_cat_marker_df$name, levels = name_order)
	dot_mgg <- ggplot(plot_top_cat_marker_df, aes(y = name, x = var_name)) +
		geom_point(aes(color = logfc, size = logpadj)) +
		scale_color_continuous_divergingx(palette = 'RdBu', mid = 0.0, p3 = 1, p4 = 1) + 
		labs(color = 'log2FC', size = '-log10\nadjusted p-value', y = "Cluster", x = "Key categorical variables") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_top', top_m, '_cate_markers.png', sep = ""), 
	       dot_mgg, dpi = png_res, width = 10, height = 6)
}


if (numeric_vis_marker_flag) { # [Fig 2C]
	cat("Dot plot for numeric markers...\n")
	num_marker_df <- marker_df[marker_df$type == 'numeric',]
#	num_marker_df <- num_marker_df[num_marker_df$logfc >=0,]
	num_marker_df$logpadj <- -log10(num_marker_df$adjp)
	num_marker_df$cluster <- as.factor(num_marker_df$cluster)
	num_marker_df <- num_marker_df[!is.na(num_marker_df$cluster),]
	top_num_marker_df <- num_marker_df %>%
		group_by(cluster) %>%
		filter(adjp <= 0.10) %>%
		top_n(top_m,logfc)
	plot_top_num_marker_df <- num_marker_df[num_marker_df$variable %in% unique(top_num_marker_df$variable),]

	write.csv(num_marker_df, paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_all_num_markers.csv', sep = ""))
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
		labs(color = 'log2FC', size = '-log10\nadjusted p-value', y = "Cluster", x = "Key numeric variables") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_top', top_m, '_num_markers.png', sep = ""), 
	       dot_mgg, dpi = png_res, width = 10, height = 6)
}

##### Volcano plot to visualize the markers of each cluster [Fig. S2x]
if (volcano_flag) { # Fig S3
	marker_annot_df <- as.data.frame(read_excel(paste(input_prefix, "cluster4_kmeans_direct_knn2imp_12112020_data_clean_marker_df_JM_ZH_v2.xlsx", sep = "")))
#	marker_df <- read.csv(paste(input_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_marker_df.csv', sep = ""), row.names=1)
	marker_annot_df$unique_label <- str_c(marker_annot_df$variable, "_", marker_annot_df$cluster)
	marker_df$unique_label <- str_c(marker_df$variable, "_", marker_df$cluster)

	merge_marker_df <- merge(marker_df, marker_annot_df[,c('Volcano plot', 'Short name', 'Varible type', 'unique_label')], by = "unique_label")
	merge_marker_df$logpadj <- -log10(merge_marker_df$adjp)
	merge_marker_df$cluster <- as.factor(merge_marker_df$cluster)
	merge_marker_df$show_variable <- NA
	vol_var_mask <- merge_marker_df$`Volcano plot`=='Y'
	merge_marker_df$show_variable[vol_var_mask] <- merge_marker_df[vol_var_mask, 'variable']
	print(dim(marker_annot_df))
	print(dim(marker_df))
	print(head(marker_annot_df))
	print(head(merge_marker_df))

	## Numeric [Fig S2x]
	num_marker_df <- merge_marker_df[merge_marker_df$type=='numeric',]
	volgg <- ggplot(num_marker_df, aes(x = logfc, y = logpadj)) +
		geom_point(aes(color=name), alpha=0.6, size=3) +
		geom_label_repel(aes(label=`Short name`), segment.color = 'grey50',
				 max.overlaps = 15)+ # box.padding   = 0.35, point.padding = 0.5, 
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = 1) +
		scale_color_brewer(palette = "Spectral") +
		guides(color = guide_legend(override.aes = list(size = 5, alpha=1))) +
		facet_wrap(~name, ncol=2, scales='free') +
		labs(color="Clusters", x='log2FC', y = '-log10 adjusted p-value') +
		theme_bw() +
		theme(axis.title = element_text(size=21),
		axis.text=element_text(size=18),
		legend.position = "top",
		legend.text=element_text(size=21),
		legend.title=element_text(size=24),
		strip.text=element_text(size=21))
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_top', input_id, '_volcano_num_all.png', sep = ""), 
	       volgg, dpi = png_res, width = 27, height = 18)
	volgg <- ggplot(num_marker_df, aes(x = logfc, y = logpadj)) +
		geom_point(aes(color=name), alpha=0.6, size=3) +
		geom_label_repel(aes(label=show_variable), segment.color = 'grey50',
				 max.overlaps = 15)+ # box.padding   = 0.35, point.padding = 0.5, 
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = 1) +
		scale_color_brewer(palette = "Spectral") +
		guides(color = guide_legend(override.aes = list(size = 5, alpha=1))) +
		facet_wrap(~name, ncol=2, scales='free') +
		labs(color="Clusters", x='log2FC', y = '-log10 adjusted p-value') +
		theme_bw() +
		theme(axis.title = element_text(size=21),
		axis.text=element_text(size=18),
		legend.position = "top",
		legend.text=element_text(size=21),
		legend.title=element_text(size=24),
		strip.text=element_text(size=21))
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_top', input_id, '_volcano_var_num_all.png', sep = ""),
	       volgg, dpi = png_res, width = 27, height = 18)

	## Categorical [Fig S2x]
	num_marker_df <- merge_marker_df[merge_marker_df$type!='numeric',]
	volgg <- ggplot(num_marker_df, aes(x = logfc, y = logpadj)) +
		geom_point(aes(color=name), alpha=0.6, size=3) +
		geom_label_repel(aes(label=`Short name`), segment.color = 'grey50',
				 max.overlaps = 15)+ # box.padding   = 0.35, point.padding = 0.5, 
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = 1) +
		scale_color_brewer(palette = "Spectral") +
		guides(color = guide_legend(override.aes = list(size = 5, alpha=1))) +
		facet_wrap(~name, ncol=2, scales='free') +
		labs(color="Clusters", x='Cramer\'s V', y = '-log10 adjusted p-value') +
		theme_bw() +
		theme(axis.title = element_text(size=21),
		axis.text=element_text(size=18),
		legend.position = "top",
		legend.text=element_text(size=21),
		legend.title=element_text(size=24),
		strip.text=element_text(size=21))
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_top', input_id, '_volcano_cate_all.png', sep = ""), 
	       volgg, dpi = png_res, width = 27, height = 18)
	volgg <- ggplot(num_marker_df, aes(x = logfc, y = logpadj)) +
		geom_point(aes(color=name), alpha=0.6, size=3) +
		geom_label_repel(aes(label=show_variable), segment.color = 'grey50',
				 max.overlaps = 15)+ # box.padding   = 0.35, point.padding = 0.5, 
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = 1) +
		scale_color_brewer(palette = "Spectral") +
		guides(color = guide_legend(override.aes = list(size = 5, alpha=1))) +
		facet_wrap(~name, ncol=2, scales='free') +
		labs(color="Clusters", x='Cramer\'s V', y = '-log10 adjusted p-value') +
		theme_bw() +
		theme(axis.title = element_text(size=21),
		axis.text=element_text(size=18),
		legend.position = "top",
		legend.text=element_text(size=21),
		legend.title=element_text(size=24),
		strip.text=element_text(size=21))
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_top', input_id, '_volcano_var_cate_all.png',sep=""),
	       volgg, dpi = png_res, width = 27, height = 18)
}

##### UMAP for each variable [Fig. 2Ds]
if (umap_vis_flag) {
	umap_markers <- c('P01BMI','P01BPTOT','P01KPACDCV','V00CESD','V00COMORB','V00LFSFR','V00RFSFR','V00LKALNMT','V00RKALNMT','V00RKFHDEG',
			  'V00RKFHDEG','V00VITACV','V00VITCCV','V00WOMADLL','V00WOMADLR','V00WOMKPL','V00WOMKPR','V00WOMSTFL','V00WOMSTFR',
			  'V00WOMTSL','V00WOMTSR')
	umap_markers <- c("V00AGE",  "P02SEX", "P02RACE", "P01BMI",'V00WOMKPL','V00WOMKPR','V00VITCCV', 'V00SUPVITC',"V00COMORB","V00HSPSS", "V00CESD", "V00HSMSS", 
			  "P01LXRKOA", "P01RXRKOA", "V00EDCV", "V00INCOME")
#	print(colnames(data_df)[str_detect(colnames(data_df), "MSS")]) # NOTE: where is the baseline KL grade?
#	q(save = "no")

	umap_plot_list <- vector(mode='list', length=length(umap_markers))
	names(umap_plot_list) <- 1:length(umap_markers)
	umap_plot_map <- matrix(1:length(umap_markers), ncol=8, nrow=2)
	iu <- 1
	rescale_fun <- function(x, to=c(0,1), from=NULL, probs=c(0.1,0.9)) {
		xedges <- quantile(x, prob=probs, na.rm=T, names=FALSE)
		prep_x <- x
		prep_x[prep_x<xedges[1]] <- xedges[1]
		prep_x[prep_x>xedges[2]] <- xedges[2]
		x <- scales::rescale(prep_x, to=to, from=xedges)
		return(x)
	}
	for (ium in umap_markers) {
		tmp_umap_df <- cbind(umap_df, data_df[,ium])
		colnames(tmp_umap_df)[ncol(tmp_umap_df)] <- ium
		tmp_umap_df <- tmp_umap_df[order(tmp_umap_df[,ium]),]
		print(head(tmp_umap_df))
		umap_gg <- ggplot(tmp_umap_df, aes_string(x = "UMAP1", y = "UMAP2", color = ium)) + 
			geom_point(size = 0.1, alpha=0.6) +
			labs(color=ium) +
			theme_classic()
		marker_types <- sapply(tmp_umap_df, class)
		print(marker_types)

		if (marker_types[[ium]] == 'numeric' | marker_types[[ium]] == "integer") {
			umap_gg <- umap_gg + 
				scale_color_viridis_c(option = "plasma", rescaler=rescale_fun)
		ggsave(paste(output_prefix,"cluster",cluster_num,'_kmeans_direct_knn2imp_',input_id,'_',ium, '_umap_r.png', sep = ""), umap_gg, 
		       dpi = png_res, width = 3.6, height = 2)
		} else {
			umap_gg <- umap_gg + scale_color_brewer(palette = "Paired") +
				guides(colour = guide_legend(override.aes = list(size = 1, alpha=1))) +
				theme(legend.title=element_text(size=9),
				      legend.spacing.y = unit(0.01, 'mm'),
				legend.text=element_text(size=3))
		ggsave(paste(output_prefix,"cluster",cluster_num,'_kmeans_direct_knn2imp_',input_id,'_',ium, '_umap_r.png', sep = ""), umap_gg, 
		       dpi = png_res, width = 3.6, height = 2.4)
		}
		umap_plot_list[[iu]] <- umap_gg
		iu <- iu + 1
	}
	q(save = "no")
	comb_umap_gg <- grid.arrange(grobs = umap_plot_list, layout_matrix = umap_plot_map)
	ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_selected_marker.png', sep = ""), comb_umap_gg, 
	       dpi = png_res, width = 25.6, height = 6.4)
}

##### Violin plots directly compare the top markers [Fig S2x]
if (violin_vis_flag) {
	cat("Specific visualization for markers...\n")
	plot_list <- vector(mode = 'list', length = cluster_num*top_m)
	names(plot_list) <- 1:(cluster_num*top_m)
	plot_m <- matrix(1:(cluster_num*top_m), ncol = cluster_num, nrow = top_m)
	ic <- 1
	for (c in 0:(cluster_num-1)) {
		cat("\tCluster", c, "\n")
		clst_marker_df <- marker_df[marker_df$logfc>0|is.na(marker_df$logfc),]
		clst_marker_df <- clst_marker_df[clst_marker_df$cluster==c,]
		clst_marker_df <- clst_marker_df[clst_marker_df$adjp <= 0.10,]
		top_num_marker_df <- clst_marker_df %>% filter(type=='numeric') %>% top_n(top_m, logfc)
		top_cate_marker_df <- clst_marker_df %>% filter(type!='numeric') %>% top_n(top_m, logfc)

		cat("\t\tCategorical...\n")
		if (nrow(top_cate_marker_df) == top_m) {
			for (ir in 1:nrow(top_cate_marker_df)) {
				tmp_marker <- top_cate_marker_df[ir, 'variable']
				cat("\t\t\t", tmp_marker, "\n")
				tmp_plot_df <- cbind(umap_df, data_df[, tmp_marker])
				colnames(tmp_plot_df)[ncol(tmp_plot_df)] <- tmp_marker
				tmp_plot_df$cluster <- as.factor(tmp_plot_df$cluster)
				tmp_table <- as.data.frame(table(tmp_plot_df[,'name'], tmp_plot_df[, tmp_marker]))
				colnames(tmp_table) <- c("name", tmp_marker, "Freq")
				tmp_table$name <- factor(tmp_table$name, 
							 levels = c("Poor knee & general health", "Intermediate knee & general health", 
								    "Good knee & general health", "Low supplemental vitamins"))
				tmp_gg <- ggplot(tmp_table, aes_string(x = 'name', y = 'Freq', fill = tmp_marker)) +
					geom_bar(stat = "identity", position = position_stack(), width = 0.5) +
					scale_fill_brewer(palette = 'Paired') +
					labs(x = 'Cluster', y = 'Number of subjects', fill = tmp_marker) +
					coord_flip() +
					theme_bw()
				ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_', tmp_marker, '_stack_bar.png', sep = ""), 
				       tmp_gg, dpi = png_res, width = 9, height = 3)

				tmp_gg <- ggplot(tmp_table, aes_string(x = 'name', y = 'Freq', fill = tmp_marker)) +
					geom_bar(stat = "identity", position = "fill", width = 0.5) +
					scale_fill_brewer(palette = 'Paired') +
					labs(x = 'Cluster', y = 'Number of subjects', fill = tmp_marker) +
					coord_flip() +
					theme_bw()
				ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_', tmp_marker, '_fill_bar.png', sep = ""), 
				       tmp_gg, dpi = png_res, width = 9, height = 3)
			}


		cat("\t\tNumeric...\n")
		if (nrow(top_num_marker_df) == top_m) {
			for (ir in 1:nrow(top_num_marker_df)) {
				tmp_marker <- top_num_marker_df[ir, 'variable']
				cat("\t\t\t", tmp_marker, "\n")
				tmp_plot_df <- cbind(umap_df, data_df[, tmp_marker])
				colnames(tmp_plot_df)[ncol(tmp_plot_df)] <- tmp_marker
				tmp_gg <- ggplot(tmp_plot_df, aes_string(x = "name", y = tmp_marker, color = "name")) +
					geom_violin() +
					geom_jitter(position = position_jitterdodge(0.45), alpha = 0.24, size = 0.3) +
					scale_color_brewer(palette = 'Spectral') +
					labs(x = "Cluster") +
					stat_compare_means(aes_string(group = "name"), label = "p.format", label.y = max(tmp_plot_df[,tmp_marker], na.rm = T)*0.8, label.x = 1.5) +
					coord_flip() +
					theme_bw() +
					theme(legend.position = "None")
				ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_', tmp_marker, '_violin.png', sep = ""), tmp_gg, 
				       dpi = png_res, width = 6, height = 3)

				tmp_gg <- ggplot(tmp_plot_df, aes_string(x = "name", y = tmp_marker, color = "name")) +
					geom_boxplot(outlier.shape=NA) +
					geom_jitter(position = position_jitterdodge(0.45), alpha = 0.24, size = 0.3) +
					scale_color_brewer(palette = 'Spectral') +
					labs(x = "Cluster") +
					stat_compare_means(aes_string(group = "name"), label = "p.format", label.y = max(tmp_plot_df[,tmp_marker], na.rm = T)*0.8, label.x = 1.5) +
					coord_flip() +
					theme_bw() +
					theme(legend.position = "None")
				ggsave(paste(output_prefix, "cluster", cluster_num, '_kmeans_direct_knn2imp_', input_id, '_', tmp_marker, '_boxplot.png', sep = ""), tmp_gg, 
				       dpi = png_res, width = 6, height = 3)
			}
		}
		}
	}
}

if (supplement_flag) {
	var_ann <- as.data.frame(read_excel(paste(data_dir, "AllClinical00_V6_column_annotation.xlsx", sep = "/"), sheet = "Baseline"))
	print(head(var_ann))
	sup_var_ann <- var_ann[var_ann$Sup == "Yes",]
	suppf <- paste(data_dir, "/supp_analysis/supp_analysis_", sep= "")
	print(head(marker_df))

	non_supp_marker_df <- marker_df[!(marker_df$variable %in% sup_var_ann$Variables),]
	non_supp_marker_df$var_type <- ifelse(is.na(non_supp_marker_df$diff), "Categorical", "Numerical")
	non_supp_marker_df <- non_supp_marker_df %>% add_significance('adjp')
	select_marker_df <- non_supp_marker_df %>% filter(cluster == 0, adjp.signif != "ns") %>% group_by(var_type)%>% top_n(wt = abs(logfc), n = 25)
	print(head(select_marker_df))
	select_non_supp_marker_df <- non_supp_marker_df[non_supp_marker_df$variable %in% select_marker_df$variable,]

	dot_gg <- ggplot(select_non_supp_marker_df, aes(x = name, y = variable)) +
		geom_point(aes(color = logfc, size = adjp.signif)) +
		scale_color_continuous_divergingx(palette = 'RdBu', mid = 0.0, p3 = 1, p4 = 1) + 
		scale_size_manual(values = c("ns" = 0.1, "*" = 1, "**" = 2, "***" = 4, "****" = 6)) +
		labs(color = 'log2FC', size = 'Adjusted\nstatistic\nsignificance', y = "Cluster", x = "Nutrition related variables") +
		facet_row(~var_type, scales = "free", space = "free") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
	ggsave(paste(suppf, "nonsupp_marker_dotplot.png", sep = ""), dot_gg, dpi = 300, width = 12, height = 12)

	print(head(non_supp_marker_df))

	supp_marker_df <- marker_df[marker_df$variable %in% sup_var_ann$Variables,]
	supp_marker_df <- merge(supp_marker_df, sup_var_ann, by.x = "variable", by.y = "Variables", all.x = T)
	supp_marker_df$info_label <- str_split_fixed(supp_marker_df$Label, "\\\n", n = 2)[1,]
	supp_marker_df$var_type <- ifelse(is.na(supp_marker_df$diff), "Categorical", "Numerical")
	supp_marker_df$Resource <- ifelse(str_detect(supp_marker_df$Label, "from food"), "From food", 
				    ifelse(str_detect(supp_marker_df$Label, "from vitamin"), "From suppplements", "Frequency"))
	supp_marker_df <- supp_marker_df %>% add_significance('adjp')

	dot_gg <- ggplot(supp_marker_df, aes(x = name, y = variable)) +
		geom_point(aes(color = logfc, size = adjp.signif)) +
		scale_color_continuous_divergingx(palette = 'RdBu', mid = 0.0, p3 = 1, p4 = 1) + 
		scale_size_manual(values = c("ns" = 0.1, "*" = 1, "**" = 2, "***" = 4, "****" = 6)) +
		labs(color = 'log2FC', size = 'Adjusted\nstatistic\nsignificance', y = "Cluster", x = "Nutrition related variables") +
		facet_col(~Resource+var_type, scales = "free_y", space = "free") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
	ggsave(paste(suppf, "marker_dotplot.png", sep = ""), dot_gg, dpi = 300, width = 9, height = 21)

	print(head(supp_marker_df))
}

