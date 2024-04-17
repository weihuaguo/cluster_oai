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

data_dir <- "/mnt/sda1/OAI_Data/data_summary"
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
var_ann <- as.data.frame(read_excel(paste(data_dir, "AllClinical00_V6_column_annotation.xlsx", sep = "/"), sheet = "Baseline"))
print(head(var_ann))
sup_var_ann <- var_ann[var_ann$Sup == "Yes",]

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
supp_ann$Resource <- ifelse(str_detect(supp_ann$ID, "from food"), "From food", 
			    ifelse(str_detect(supp_ann$ID, "from vitamin"), "From suppplements", "Frequency")
)
print(head(supp_ann))
hm_df <- hm_df[rownames(supp_ann),]
row_ann <- rowAnnotation(Resource = supp_ann$Resource,
			     col = list(Resource = c("From food" = "forestgreen", "From suppplements" = "dodgerblue", "Frequency" = "gray"))
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
col_ann <- HeatmapAnnotation(Side = outcome_ann$Side,
			     col = list(Side = c("L" = "goldenrod", "R" = "purple"))
)


cor_hm <- Heatmap(hm_df, 
		  name = "Correlation",
		  heatmap_width = unit(20, 'cm'),
		  top_annotation = col_ann,
		  left_annotation = row_ann,
		  column_split = outcome_ann$Side,
		  row_split = supp_ann$Resource
) 
png(paste(rpf, "supplement_outcome_cor_heatmap_clusters.png", sep = ""), res = 300, width = 24, height = 16, units = 'in')
draw(cor_hm, heatmap_legend_side = "left", merge_legend = T)
gar <- dev.off()
