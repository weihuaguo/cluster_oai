# Visualize cross-validation results
# Weihua Guo, Ph.D.
# 06/12/2021

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
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(gridExtra))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(scales))
suppressMessages(library(rstatix))
suppressMessages(library(broom))
suppressMessages(library(corrplot))
suppressMessages(library(tibble))
suppressMessages(library(circlize))
suppressMessages(library(RColorBrewer))

clst_dir <- "/mnt/sda1/OAI_Data/cluster_kmean_cv/"
data_dir <- "/mnt/sda1/OAI_Data/data_summary/"

exp_id <- "kmean_pca_umap"
input_id <- "12112020_data_clean"
clean_co <- 'v25'
input_prefix <- paste(clst_dir, 'clean_', clean_co, '_', sep = "")
cluster_num <- 4

png_res <- 300
top_m <- 6
sss <- 0:9
names <- c("C0" = "Low supplemental vitamins", "C1" = "Poor knee & general health", "C2" = "Good knee & general health", "C3" = "Intermediate knee & general health")
cluster_order <- c("Poor knee & general health", "Intermediate knee & general health", "Good knee & general health", "Low supplemental vitamins")

marker_flag <- TRUE
merge_flag <- TRUE
vis_flag <- TRUE

if (marker_flag) {
	all_data_df <- as.data.frame(read.csv(paste(data_dir, "input_direct_merge_dataframe_",input_id,"_sbj_clean_",clean_co,".csv", sep = ""), 
				  sep = ",", header = TRUE, row.names=1, stringsAsFactors = FALSE))
	print(head(all_data_df[1:9,1:6]))
	for (isss in sss) {
		umap_df <- as.data.frame(read_excel(paste(clst_dir, "clean_", clean_co, "_ss", isss, "_cluster", cluster_num, "_kmean_pca_umap_res.xlsx", sep = "")))
		umap_df$ID <- as.character(umap_df$ID)
		data_df <- all_data_df[umap_df$ID,]
		marker_types <- sapply(data_df, class)
		marker_cols <- c('variable', 'cluster', 'type', 'test','avg1', 'avg2', 'diff', 'logfc', 'per1', 'per2', 'pval')
		marker_df <- as.data.frame(matrix(nrow=(cluster_num*ncol(data_df)), ncol=length(marker_cols)))
		colnames(marker_df) <- marker_cols
		ic <- 1
		output_prefix <- paste(clst_dir, "clean_", clean_co, "_ss", isss, "_cluster", cluster_num, "_", sep = "")
		for (i in 0:(cluster_num-1)) {
			cat("\tCluster", i, '\n')
			tmp_umap_df <- umap_df
			tmp_umap_df$comp <- ifelse(tmp_umap_df$kmean_pca == i, 'COI', 'Others')
			for (im in names(marker_types)) {
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
					marker_df[ic, 'test'] <- 'Kruskal-Wallis test'
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
					tmp_table <- as.data.frame.matrix(table(tmp_df$comp, tmp_df[,ncol(tmp_df)]))
					tmp_prop_table <- prop.table(as.matrix(tmp_table),margin=1)
					if (sum(tmp_prop_table[1,]==max(tmp_prop_table[1,]))>1|sum(tmp_prop_table[2,]==max(tmp_prop_table[2,]))>1) {
						marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])][1]
						marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])][1]
					} else {
						marker_df[ic, 'avg1'] <- colnames(tmp_prop_table)[tmp_prop_table[1,]==max(tmp_prop_table[1,])]
						marker_df[ic, 'avg2'] <- colnames(tmp_prop_table)[tmp_prop_table[2,]==max(tmp_prop_table[2,])]
					}
					cat(im, '\n')
					tmp_test <- fisher.test(tmp_table, simulate.p.value=TRUE, B=1e7)
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
		marker_df$ss <- isss
		write.csv(marker_df, paste(output_prefix, "marker_df_largeB.csv", sep = ""))
	}
}

if (merge_flag) {
	ic <- 1
	for (isss in sss) {
		output_prefix <- paste(clst_dir, "clean_", clean_co, "_ss", isss, "_cluster", cluster_num, "_", sep = "")
		tmp_df <- read.csv(paste(output_prefix, "marker_df_largeB.csv", sep = ""))
		umap_df <- as.data.frame(read_excel(paste(clst_dir, "clean_", clean_co, "_ss", isss, "_cluster", cluster_num, "_kmean_pca_umap_res.xlsx", sep = "")))
		umap_df$ss <- isss

		if (ic == 1) {
			merge_df <- tmp_df
			merge_umap_df <- umap_df
		} else {
			merge_df <- rbind(merge_df, tmp_df)
			merge_umap_df <- rbind(merge_umap_df, umap_df)
		}
		ic <- ic + 1
	}
	write.csv(merge_df, paste(clst_dir, "clean_", clean_co, "_cluster", cluster_num, "_merged_marker_df.csv", sep = ""))
	write.csv(merge_umap_df, paste(clst_dir, "clean_", clean_co, "_cluster", cluster_num, "_merged_umap_df.csv", sep = ""))
}

if (vis_flag) {
	merge_umap_df <- read.csv(paste(clst_dir, "clean_", clean_co, "_cluster", cluster_num, "_merged_umap_df.csv", sep = ""), row.names = 1)
	true_umap_df <- as.data.frame(read_excel(paste(clst_dir, "clean_", clean_co, "_cluster", cluster_num, "_kmean_pca_umap_res.xlsx", sep = "")))

	umap_df <- merge_umap_df
	umap_df$test <- str_c("Test", umap_df$ss)

	all_marker_df <- read.csv(paste(clst_dir, "clean_", clean_co, "_final_cluster", cluster_num, "_kmeans_direct_knn2imp_12112020_data_clean_marker_df_largeB.csv", sep = ""), 
				  row.names = 1)
	all_marker_df$ss <- "N"
	marker_df <- read.csv(paste(clst_dir, "clean_", clean_co, "_cluster", cluster_num, "_merged_marker_df.csv", sep = ""), row.names = 1)
	marker_df$X <- NULL
	marker_df <- rbind(all_marker_df, marker_df)
	marker_df$unique_name <- str_c("Test", marker_df$ss, "_Cluster", marker_df$cluster)

	spr_marker_df <- marker_df[,c('variable', 'adjp', 'unique_name')] %>% spread(unique_name, adjp)

	cat("Generate the cluster alignment with the true clusters...\n")
	cor_marker_df <- spr_marker_df %>% cor_test(unique(marker_df$unique_name))

	dens_cor_df <- cor_marker_df[str_detect(cor_marker_df$var1, "TestN_Cluster"),]
	dens_cor_df <- dens_cor_df[dens_cor_df$cor < 1,]
	dens_cor_df <- dens_cor_df[!str_detect(dens_cor_df$var2, "TestN_Cluster"),]
	top_cor_df <- dens_cor_df %>% group_by(var2) %>% arrange(desc(cor)) %>% filter(row_number()==1) %>% arrange(var2, var1)
	write.csv(dens_cor_df, paste(clst_dir, "fold_cluster_match_table.csv", sep = ""))
	write.csv(top_cor_df, paste(clst_dir, "top_fold_cluster_match_table.csv", sep = ""))

	cat("Visualize the cluster alignment with the true clusters...\n")
	cor_marker_df$abs_cor <- abs(cor_marker_df$cor)
	cor_marker_df$logp <- -log10(cor_marker_df$p)
	cor_dot <- ggplot(cor_marker_df, aes(x=var1, y=var2, size=cor,color=cor)) +
		geom_point() +
		scale_color_distiller(palette='RdYlBu')+
		labs(color = "PCC", size = "PCC", x = "Test_Cluster", y = "Test_Cluster") +
		theme_bw() +
		theme(axis.text.x=element_text(angle=45, hjust=1))
	ggsave(paste(clst_dir, "dotplot_cv_cor_pearson_all.png", sep = ""), cor_dot, dpi = png_res, width = 12, height = 12)

	dens_cor_df$test <- str_split_fixed(dens_cor_df$var2, "_", n=2)[,1]
	cor_dot <- ggplot(dens_cor_df, aes(x=var1, y=var2, size=cor,color=cor)) +
		geom_point() +
		scale_color_distiller(palette='RdYlBu')+
		labs(color = "PCC", size = "PCC", x = "True cluster", y = "Cross validation cluster") +
		facet_wrap(~test, ncol=5, scales='free') +
		theme_bw() +
		theme(axis.text.x=element_text(angle=45, hjust=1))
	ggsave(paste(clst_dir, "dotplot_cv_cor_pearson_align.png", sep = ""), cor_dot, dpi = png_res, width = 11.2, height = 4.2)

	marker_df$align_name <- marker_df$unique_name
	for (ir in 1:nrow(top_cor_df)) {
		tmp_mask <- marker_df$align_name == top_cor_df$var2[ir]
		marker_df$align_name[tmp_mask] <- top_cor_df$var1[ir]
	}
	marker_df$align_name <- str_replace(marker_df$align_name, "TestN_", "")
	marker_df$var_name <- str_c(marker_df$variable, "-", marker_df$align_name)
	marker_df$cv <- str_c("Test",marker_df$ss)
	marker_df$sig <- ifelse(marker_df$adjp <= 0.05, "Sig", "Insig")
	marker_df$cv[marker_df$cv == "TestN"] <- "NonCV"

	var_cts_df <- marker_df %>%
		group_by(align_name, variable, sig) %>%
		summarise(n=n())
	print(head(marker_df))
	align_marker_df <- marker_df[,c("var_name", "cv", "sig")] %>% spread(cv, sig)
	align_marker_df <- gather(align_marker_df, "cv_test", "cv_sig", str_c("Test", sss))
	align_marker_df$comp <- align_marker_df$NonCV == align_marker_df$cv_sig

	test_sum <- align_marker_df %>%
		group_by(cv_test) %>%
		summarize(n = n(), aligned = sum(comp), rel = sum(comp)/n()*100) 
	var_sum <- align_marker_df %>%
		group_by(var_name) %>%
		summarize(n = n(), aligned = sum(comp), rel = sum(comp)/n()*100) 

	write.csv(test_sum, paste(clst_dir, "marker_sig_alignment_with_test.csv", sep = ""))
	write.csv(var_sum, paste(clst_dir, "marker_alignment_with_variable_cluster.csv", sep = ""))
	rownames(test_sum) <- test_sum$cv_test
	test_sum <- as.data.frame(t(test_sum))
	test_sum$NonCV <- 100.00
	test_sum <- as.data.frame(t(test_sum))
	rownames(var_sum) <- var_sum$var_name

	hm_df <- marker_df[,c('var_name', 'adjp', 'cv')] %>% spread(var_name, adjp)
	rownames(hm_df) <- hm_df[,1]
	hm_df[,1] <- NULL

	print(dim(hm_df))
	print(hm_df[,1:9])

	var_df <- as.data.frame(str_split_fixed(colnames(hm_df), "-", n=2))
	rownames(var_df) <- colnames(hm_df)
	colnames(var_df) <- c("Variable", "Cluster")
	print(head(var_df))

	test_sum <- test_sum[rownames(hm_df),]
	test_sum$rel <- as.numeric(test_sum$rel)
	var_sum <- var_sum[colnames(hm_df),]
	print(head(var_sum))
	print(head(test_sum))

	true_annot <- HeatmapAnnotation(Consistency = anno_barplot(var_sum$rel))
	row_annot <- rowAnnotation(Consistency = anno_barplot(test_sum$rel, width = unit(3, 'cm')))

	col_fun = colorRamp2(c(min(marker_df$adjp), 0.05, 1.0), c('dodgerblue','white','firebrick'))
	var_hm <- Heatmap(hm_df, 
			  col=col_fun,
			  top_annotation = true_annot,
			  left_annotation = row_annot,
			  name='Adjusted\np-value',
			  show_column_names=FALSE,
			  cluster_rows=FALSE,
			  column_split = var_df$Cluster,
			  heatmap_legend_param = list(at = c(0, 0.05, 1), 
						      labels = c("(0,0.05)", "0.05", "(0.05,1]")
			  )
	)

	png(paste(clst_dir, "var_align_heatmap.png", sep = ""), res=png_res, width=12, height=4, units='in')
	print(var_hm)
	gar <- dev.off()

	cat("Align the clusters to true clusters on samples...\n")
	umap_df$unique_name <- str_c(umap_df$test, "_Cluster", merge_umap_df$kmean_pca)
	umap_df$align_name <- umap_df$unique_name
	for (ir in 1:nrow(top_cor_df)) {
		tmp_mask <- umap_df$align_name == top_cor_df$var2[ir]
		umap_df$align_name[tmp_mask] <- top_cor_df$var1[ir]
	}

	umap_df$align_name <- str_replace(umap_df$align_name, "TestN_", "")
	umap_df$test <- str_c("Test", umap_df$ss)
#	print(head(umap_df))
	align_umap_df <- merge(umap_df, true_umap_df, by = "ID", all.x = TRUE)
	align_umap_df$true_cluster <- str_c("Cluster", align_umap_df$kmean_pca.y)
	align_umap_df$comp <- align_umap_df$align_name == align_umap_df$true_cluster
#	print(head(align_umap_df))

	test_sum <- align_umap_df %>%
		group_by(test) %>%
		summarize(n = n(), aligned = sum(comp), rel = sum(comp)/n()*100) 
	id_sum <- align_umap_df %>%
		group_by(ID) %>%
		summarize(n = n(), aligned = sum(comp), rel = sum(comp)/n()*100) 

	write.csv(test_sum, paste(clst_dir, "sample_alignment_with_test.csv", sep = ""))
	write.csv(id_sum, paste(clst_dir, "sampe_alignment_with_subject.csv", sep = ""))
	rownames(test_sum) <- test_sum$test
	rownames(id_sum) <- id_sum$ID

	umap_df$align_name <- as.numeric(str_replace(umap_df$align_name, "Cluster", ""))
	hm_smp_df <- umap_df[,c("test", "align_name", "ID")] %>% spread(ID, align_name)
	rownames(hm_smp_df) <- hm_smp_df$test
	hm_smp_df$test <- NULL

	annot_df <- true_umap_df
	annot_df$cluster <- str_c("C", annot_df$kmean_pca)
	annot_names <- c()
	
	for (ic in unique(annot_df$cluster)) {
		annot_df$name[annot_df$cluster == ic] <- names[ic]
		annot_names <- c(annot_names, names[ic])
	}
	annot_df$name <- factor(annot_df$name, levels = cluster_order)
	rownames(annot_df) <- annot_df$ID
	annot_df <- annot_df[colnames(hm_smp_df),]
	test_sum <- test_sum[rownames(hm_smp_df),]
	id_sum <- id_sum[colnames(hm_smp_df),]
	print(head(id_sum))


	col_smp_annot <- structure(brewer.pal(n=4, name='Spectral'), 
			     names = cluster_order)
	true_annot <- HeatmapAnnotation(Cluster = annot_df$name, Consistency = anno_barplot(id_sum$rel), col = list(Cluster = col_smp_annot))
	row_annot <- rowAnnotation(Consistency = anno_barplot(test_sum$rel))


	col_smp <- structure(brewer.pal(n=4, name='Spectral'), 
			     names = c("1", "3", "2", "0"))
	smp_hm <- Heatmap(hm_smp_df, 
			  col=col_smp,
			  top_annotation = true_annot,
			  left_annotation = row_annot,
			  name='Cluster',
			  column_title='Subjects',
			  row_title = 'Cross validation tests',
			  show_column_names=FALSE,
			  cluster_rows = FALSE,
			  show_column_dend = FALSE, 
			  show_heatmap_legend = FALSE)

	png(paste(clst_dir, "sample_align_heatmap.png", sep = ""), res=png_res, width=9, height=4, units='in')
	print(draw(smp_hm, heatmap_legend_side = NULL))
	gar <- dev.off()
}
