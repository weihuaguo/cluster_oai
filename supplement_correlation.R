# Merge, analyze, visualize and select the predictor-outcomes
# Weihua Guo, Ph.D.
# 02/13/2023

suppressMessages(library(dplyr))
suppressMessages(library(readxl))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(viridis))
suppressMessages(library(rstatix))
suppressMessages(library(ggforce))
suppressMessages(library(RColorBrewer))


data_dir <- "/mnt/sda1/OAI_Data/data_summary"
var_ann <- as.data.frame(read_excel(paste(data_dir, "AllClinical00_V6_column_annotation.xlsx", sep = "/"), sheet = "Baseline"))
print(head(var_ann))
sup_var_ann <- var_ann[var_ann$Sup == "Yes",]
suppf <- paste(data_dir, "/supp_analysis/supp_analysis_240802_", sep= "")

clean_co<-'v25'
cor_df<-read.csv(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_correlation_", clean_co, ".csv", sep=''),row.names=1)
p_df<-read.csv(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_pvalue_",clean_co,".csv", sep=''),row.names=1)
num_df<-read.csv(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_na_num_",clean_co,".csv", sep=''),row.names=1)
cat("Replace NA's with abnormal values for clustering...\n")

print(cor_df[1:9,1:6])
print(dim(cor_df))
cat("Correlation to other variables\n")
ol_var <- intersect(sup_var_ann$Variable, rownames(cor_df))
use_cor_df <- cor_df[ol_var, !(colnames(cor_df) %in% ol_var)]
use_p_df <- p_df[ol_var,!(colnames(cor_df) %in% ol_var)]
print(dim(use_cor_df))

gath_cor_df <- use_cor_df
all_other_col <- colnames(use_cor_df)
gath_cor_df$supp_var <- rownames(gath_cor_df)
gath_cor_df <- gather(gath_cor_df, "other_var", "cor", all_other_col)

gath_p_df <- use_p_df
all_other_col <- colnames(use_p_df)
gath_p_df$supp_var <- rownames(gath_p_df)
gath_p_df <- gather(gath_p_df, "other_var", "p", all_other_col)

gath_df <- merge(gath_cor_df, gath_p_df, by= c("supp_var", "other_var"))

gath_df <- merge(gath_df, sup_var_ann, by.x = "supp_var", by.y = "Variables", all.x = T)

gath_df$Resource <- ifelse(str_detect(gath_df$Label, "from food"), "From food", 
			    ifelse(str_detect(gath_df$Label, "from vitamin"), "From supplements", "Frequency"))
gath_df <- gath_df[gath_df$Resource != "Frequency",]
gath_df$p[is.na(gath_df$p)] <- 1.0
gath_df <- gath_df %>% add_significance('p')
gath_df <- gath_df %>%
	group_by(other_var, Resource) %>%
	mutate(n_sig_supp_var = sum(p.signif != "ns"))
gath_df$readable_name <- str_split_fixed(gath_df$Label, ", ", n = 2)[,2]
gath_df$readable_name <- str_split_fixed(gath_df$readable_name, "\\(cal", n = 2)[,1]
write.csv(gath_df, paste(suppf, "gathered_correlation_counts.csv", sep = ""))

gath_df <- gath_df[!str_detect(gath_df$other_var, "FFQ"),]
print(head(gath_df))
gath_df$cor <- as.numeric(gath_df$cor)
top_df <- gath_df %>%
	group_by(supp_var, Resource) %>%
	top_n(wt = cor, n = 5)
btm_df <- gath_df %>%
	group_by(supp_var, Resource) %>%
	top_n(wt = cor, n = -5)
cor_var <- unique(c(top_df$other_var, btm_df$other_var))
use_gath_df <- gath_df[gath_df$other_var %in% cor_var,]
write.csv(use_gath_df, paste(suppf, "gathered_correlation_top5.csv", sep = ""))


hm_df <- as.data.frame(spread(use_gath_df[,c("supp_var", "other_var", "cor")], "other_var", "cor"))
print(hm_df[1:9,1:6])

rownames(hm_df) <- hm_df$supp_var
print(hm_df[1:9,1:6])

hm_df$supp_var <- NULL
print(dim(hm_df))

supp_ann <- as.data.frame(use_gath_df[!duplicated(use_gath_df$supp_var),])
rownames(supp_ann) <- supp_ann$supp_var
print(head(supp_ann))
col_ann <- HeatmapAnnotation(Resource = supp_ann$Resource,
			     col = list(Resource = c("From supplements" = "goldenrod", "From food" = "purple"))
)

cor_hm <- Heatmap(t(hm_df), 
		  name = "Correlation",
		  top_annotation = col_ann,
		  column_split = supp_ann$Resource,
		  show_row_names = TRUE,
		  show_column_names = TRUE
) 
png(paste(suppf, "supplement_top5_other_cor_heatmap.png", sep = ""), res = 300, width = 16, height = 12, units = 'in')
draw(cor_hm, heatmap_legend_side = "bottom", merge_legend = T)
gar <- dev.off()

pos_gath_df <- gath_df %>%
	filter(cor > 0) %>%
	group_by(other_var, Resource) %>%
	mutate(n_pos_sig_supp_var = sum(p.signif != "ns"))

top_df <- pos_gath_df %>%
	group_by(supp_var) %>%
	top_n(wt = n_sig_supp_var, n=5)

col_num <- length(unique(top_df$other_var))
manual_color <- colorRampPalette(brewer.pal(8, "Paired"))(col_num)
gg <- ggplot(top_df, aes(x = readable_name, y = cor, color = other_var)) +
	geom_point(aes(shape = p.signif), size = 3) +
	scale_color_manual(values = manual_color) +
	labs(x = "Variable related to nutrition and supplements", y = "Correlation/association", shape = "Statistical\nsignificance", color = "Variables") +
	facet_row(~Resource, scales = "free_x", space = "free") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(suppf, "positive_sum_point.png", sep = ""), dpi = 300, width = 16, height =  9)
print(top_df)

neg_gath_df <- gath_df %>%
	filter(cor < 0) %>%
	group_by(other_var, Resource) %>%
	mutate(n_neg_sig_supp_var = sum(p.signif != "ns"))


top_df <- neg_gath_df %>%
	group_by(supp_var) %>%
	top_n(wt = n_sig_supp_var, n=5)

col_num <- length(unique(top_df$other_var))
manual_color <- colorRampPalette(brewer.pal(8, "Paired"))(col_num)
gg <- ggplot(top_df, aes(x = readable_name, y = cor, color = other_var)) +
	geom_point(aes(shape = p.signif), size = 3) +
	scale_color_manual(values = manual_color) +
	labs(x = "Variable related to nutrition and supplements", y = "Correlation/association", shape = "Statistical\nsignificance", color = "Variables") +
	facet_row(~Resource, scales = "free_x", space = "free") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(suppf, "negative_sum_point.png", sep = ""), dpi = 300, width = 16, height =  9)
print(top_df)

cat("Supplements self-correlation\n")
ol_var <- intersect(sup_var_ann$Variable, rownames(cor_df))
use_cor_df <- cor_df[ol_var, ol_var]
use_p_df <- p_df[ol_var, ol_var]
use_num_df <- num_df[ol_var, ol_var]


supp_ann <- sup_var_ann[sup_var_ann$Variables %in% ol_var,]
print(head(supp_ann))
rownames(supp_ann) <- supp_ann$Variables
supp_ann$Resource <- ifelse(str_detect(supp_ann$Label, "from food"), "From food", 
			    ifelse(str_detect(supp_ann$Label, "from vitamin"), "From suppplements", "Frequency")
)
supp_ann <- supp_ann[supp_ann$Resource != "Frequency",]

print(head(supp_ann))
use_cor_df <- use_cor_df[rownames(supp_ann), rownames(supp_ann)]
use_p_df <- use_p_df[rownames(supp_ann), rownames(supp_ann)]
use_num_df <- use_num_df[rownames(supp_ann), rownames(supp_ann)]

print(dim(supp_ann))
print(dim(use_cor_df))
row_ann <- HeatmapAnnotation(Resource = supp_ann$Resource,
			     col = list(Resource = c("From food" = "forestgreen", "From suppplements" = "dodgerblue"))
)

use_cor_df[is.na(use_cor_df)] <- -2.0
use_cor_df<-do.call(data.frame, lapply(use_cor_df, function(x) replace(x, is.infinite(x), 2.0)))
print(use_cor_df[1:9,1:6])
cor_color<-colorRamp2(c(-2, -1, 0, 1), c("dodgerblue4", "dodgerblue",  "white", "firebrick"))
p_color<-colorRamp2(c(0,0.01,0.05,0.10,1.0), c("gold4", "gold3", "gold", "white", "purple4"))
num_color<-colorRamp2(c(3600,3900,4200,4500), c("white", "darkolivegreen1", "darkolivegreen3", "darkolivegreen"))
cor_hm <- Heatmap(as.matrix(use_cor_df), name = 'Correlation', 
		  show_column_names=TRUE, 
		  show_row_names=TRUE, 
		  bottom_annotation = row_ann,
#		  row_km=5, row_km_repeats=50,
#		  column_km=5, column_km_repeats=50,
		  row_dend_width = unit(2, "cm"),
		  column_dend_height = unit(2, "cm"),
		  column_title='Correlations',
		  col=cor_color)
p_hm <- Heatmap(as.matrix(use_p_df), name = 'P-value',
		col=p_color,
		cluster_rows=FALSE, 
		cluster_columns=FALSE, 
		show_column_names=FALSE, 
		show_row_names=FALSE,
		row_order=row_order(cor_hm), 
		column_order=column_order(cor_hm), 
		column_title='P-values',
		heatmap_legend_param = list(at = c(0.0, 0.01, 0.05, 0.10, 1.0),
					    legend_height = unit(3, "cm")))
num_hm <- Heatmap(as.matrix(use_num_df), name = '#Pairs', 
		  col=num_color,
		  cluster_rows=FALSE, 
		  cluster_columns=FALSE,
		  show_column_names=FALSE, 
		  show_row_names=FALSE,
		  column_title='Number of available pairs',
		  row_order=row_order(cor_hm), 
		  column_order=column_order(cor_hm))
hm_list<-cor_hm+p_hm+num_hm
png(paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_correlation_vis_supp_",clean_co,".png", sep = ""),
width=42, height=16, res=300, units='in')
print(draw(hm_list))
gar<-dev.off()
# order_df<-as.data.frame(matrix(ncol=1,nrow=nrow(use_cor_df)))
order_df<-as.data.frame(row_order(cor_hm))
rownames(order_df)<-colnames(use_cor_df)
colnames(order_df)<-c("hclust_order")
print(head(order_df))

write.csv(order_df, paste(data_dir, "/input_direct_merge_dataframe_12112020_data_clean_correlation_order_supp_240802_",clean_co,".csv", sep = ""))

cat("Outcome association\n")
result_folder <- "predictor_vis_230213"
select_dir <- paste(data_dir, "predictor_select_230213", sep = "/")
spf <- paste(select_dir, "/predictor_select_230213_", sep = "")
result_dir <- paste(data_dir, result_folder, sep = "/")
dir.create(result_dir, showWarnings=F)
rpf <- paste(result_dir, "/", result_folder, "_", sep = "")

# write.csv(merge_cor_df, paste(rpf, "cor_merge.csv", sep = ""))
# write.csv(merge_p_df, paste(rpf, "pval_merge.csv", sep = ""))
# write.csv(merge_sn_df, paste(rpf, "sn_merge.csv", sep = ""))
merge_gath_df <- read.csv(paste(rpf, "gath_merge.csv", sep = ""), row.names = 1) #SD10
print(head(merge_gath_df))

merge_gath_df$outcome_var <- str_split_fixed(merge_gath_df$op, "__", n = 2)[,1]
merge_gath_df$sup_var <- str_split_fixed(merge_gath_df$op, "__", n = 2)[,2]

use_gath_df <- merge_gath_df[merge_gath_df$sup_var %in% sup_var_ann$Variables,]
use_gath_df <- merge(use_gath_df, sup_var_ann, by.x = "sup_var", by.y = "Variables", all.x=T)
use_gath_df$cor <- as.numeric(use_gath_df$cor)

print(head(use_gath_df))

hm_df <- spread(use_gath_df[,c("cor", "outcome_var", "Label")], "outcome_var", "cor")
rownames(hm_df) <- str_split_fixed(hm_df$Label, "\\\n", n = 3)[,1]
hm_df$Label <- NULL
print(head(hm_df))

supp_ann <- as.data.frame(matrix(nrow = nrow(hm_df), ncol = 2))
rownames(supp_ann) <- rownames(hm_df)
colnames(supp_ann) <- c("ID", "Resource")
supp_ann$ID <- rownames(supp_ann)
supp_ann$Resource <- ifelse(str_detect(supp_ann$ID, "from food"), "Daily nutrients from food (calculated)", 
			    ifelse(str_detect(supp_ann$ID, "from vitamin"), "Average daily nutrients\nfrom vitamin supplements (calculated)", "Frequency")
)
supp_ann <- supp_ann[supp_ann$Resource != "Frequency",]

print(head(supp_ann))
hm_df <- hm_df[rownames(supp_ann),]
rownames(supp_ann) <- str_replace_all(rownames(supp_ann), "\\ \\(calc\\)", "")
rownames(supp_ann) <- str_replace_all(rownames(supp_ann), "\\ \\(cal", "")
rownames(supp_ann) <- str_replace_all(rownames(supp_ann), "\\ \\(c", "")
rownames(supp_ann) <- ifelse(str_detect(rownames(supp_ann), "from food"), str_c(str_split_fixed(rownames(supp_ann), ",\\ ", n = 2)[,2]),
			     str_c(str_split_fixed(rownames(supp_ann), ",\\ ", n = 2)[,2], " ")
)
substr(rownames(supp_ann), 1,1) <- toupper(substr(rownames(supp_ann),1,1))
print(head(supp_ann))
rownames(hm_df) <- str_replace_all(rownames(hm_df), "\\ \\(calc\\)", "")
rownames(hm_df) <- str_replace_all(rownames(hm_df), "\\ \\(cal", "")
rownames(hm_df) <- str_replace_all(rownames(hm_df), "\\ \\(c", "")
rownames(hm_df) <- ifelse(str_detect(rownames(hm_df), "from food"), str_c(str_split_fixed(rownames(hm_df), ",\\ ", n = 2)[,2]),
			     str_c(str_split_fixed(rownames(hm_df), ",\\ ", n = 2)[,2], " ")
)
substr(rownames(hm_df), 1,1) <- toupper(substr(rownames(hm_df),1,1))

row_ann <- HeatmapAnnotation(Resource = supp_ann$Resource,
			     col = list(Resource = c("Daily nutrients from food (calculated)" = "forestgreen", 
						     "Average daily nutrients\nfrom vitamin supplements (calculated)" = "dodgerblue"))
)

outcome_ann <- as.data.frame(matrix(nrow = ncol(hm_df), ncol = 3))
rownames(outcome_ann) <- colnames(hm_df)
colnames(outcome_ann) <- c("ID", "Side", "Year")
outcome_ann$ID <- rownames(outcome_ann)
outcome_ann$Side <- str_sub(outcome_ann$ID, -1,-1)
outcome_ann$Year <- str_sub(outcome_ann$ID, 1, 3)
print(head(outcome_ann))
outcome_ann <- outcome_ann[outcome_ann$Year != "Y10",]
hm_df <- hm_df[,rownames(outcome_ann)]
col_ann <- rowAnnotation(Side = outcome_ann$Side,
			     col = list(Side = c("L" = "goldenrod", "R" = "purple"))
)


cor_hm <- Heatmap(t(hm_df), 
		  name = "Correlation\nbetween WOMAC-TS and nutrients",
		  heatmap_height = unit(18, 'cm'),
		  left_annotation = col_ann,
		  top_annotation = row_ann,
		  row_split = outcome_ann$Side,
		  column_split = supp_ann$Resource
) 
png(paste(rpf, "supplement_outcome_cor_heatmap_clusters.png", sep = ""), res = 300, width = 18, height = 9, units = 'in')
draw(cor_hm, heatmap_legend_side = "right", merge_legend = T)
gar <- dev.off()
