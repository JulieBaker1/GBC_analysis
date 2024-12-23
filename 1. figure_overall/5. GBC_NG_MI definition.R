# External Dataset: 10.1038/s41467-020-17880-4 (Source: Fu Lab)
## Expression Matrix: 115 individuals
## Metadata: 119 individuals
## Matched individuals: 83

setwd("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked")
set.seed(925)

library(Seurat)
library(dplyr)
library(patchwork)

library(survminer)
library(survival)

library(reshape2)
library(ggplot2)
library(ggpubr)
library(pheatmap)


library(plotly)
# library(SingleR)
#library(clusterProfiler) #
#library(org.Hs.eg.db) #
#library(enrichplot) #
library(tidyverse)
# library(monocle) # Need R 4.0
#library(GSVA) #
library(fgsea)
library(msigdbr)
library(ggrepel)
library(corrplot)
#library(IOBR) #
library(circlize)
library(RColorBrewer)
library(GGally)
library(data.table)
#library(ComplexHeatmap) #
library(dendextend)
library(stringr)
library(psych)
library(reshape)
#library(jjPlot) #
library(gg.gap)
library(gmodels)

#library(ggradar) #
#library(DESeq2) #
library(caret)



# 1 - Molecular Signature ####
## 1.1 Update our metadata ####
# PatientInfo <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/Table S1_Clinical Information_v240507 copy.xlsx")
# PatientInfo <- PatientInfo[PatientInfo$`scRNA-seq` %in% "Yes",]
# PatientInfo <- PatientInfo[,-((ncol(PatientInfo)-1):ncol(PatientInfo))]
# colnames(PatientInfo)[1] <- "NewSample.ID"
# 
# Surv_New <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/GBC样本整理信息20240916预后更新 copy.xlsx")
# Surv_New <- Surv_New[,c("NewSample.ID","OS_month","event")]
# 
# AdenoP_MI <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/AdenoP_MI.xlsx")
# colnames(AdenoP_MI)[1] <- "NewSample.ID"
# 
# PatientInfo_New <- left_join(PatientInfo,Surv_New,by="NewSample.ID")
# PatientInfo_New <- left_join(PatientInfo_New,AdenoP_MI,by="NewSample.ID")
# 
# saveRDS(PatientInfo_New, file = "./1 - GBC_input data/PatientInfo_New/PatientInfo_New.RDS")

PatientInfo_New <- readRDS("./1 - GBC_input data/PatientInfo_New/PatientInfo_New.RDS")

## 1.2 Fetch AdenoP SampleID ####
# AP_ID = PatientInfo_New[PatientInfo_New$Histological.type %in% "adenocarcinoma",]
# AP_ID = AP_ID[AP_ID$Site %in% "Primary",]
# AP_ID = AP_ID$NewSample.ID

## 1.3 Filter and Merge ####
# Cell_type <- read.csv(file = "./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
# 
# Epithelial <- readRDS("./1 - GBC_input data/Epithelial.RDS") # 268836 cells
# Epithelial = subset(Epithelial, cells = Cell_type[Cell_type$celltype %in% c("normal", "tumor"),]$cellid) # 258121 cells
# Epithelial = subset(Epithelial, subset = orig.ident %in% AP_ID) # 154633 cells
# 
# Fibroblast <- readRDS("./1 - GBC_input data/Fibroblasts.RDS") # 98012 cells
# Fibroblast = subset(Fibroblast, cells = Cell_type[Cell_type$celltype %in% c("fibroblast"),]$cellid) # 78485 cells
# Fibroblast = subset(Fibroblast, subset = orig.ident %in% AP_ID) # 40610 cells
# 
# Endothelial <- readRDS("./1 - GBC_input data/Endothelial.RDS") # 22501 cells
# Endothelial = subset(Endothelial, cells = Cell_type[Cell_type$celltype %in% c("Endothelium"),]$cellid) # 16886 cells
# Endothelial = subset(Endothelial, subset = orig.ident %in% AP_ID) # 10260 cells
# 
# Myeloid_cell <- readRDS("./1 - GBC_input data/Myeloid_cell.RDS") # 323569 cells
# Myeloid_cell = subset(Myeloid_cell, cells = Cell_type[Cell_type$celltype %in% c("DC","Mast","Mono_Macro","Neu"),]$cellid) # 307285 cells
# Myeloid_cell = subset(Myeloid_cell, subset = orig.ident %in% AP_ID) # 174162 cells
# 
# B_cell <- readRDS("./1 - GBC_input data/B_cell.RDS") # 66389 cells
# B_cell = subset(B_cell, cells = Cell_type[Cell_type$celltype %in% c("B cell"),]$cellid) # 58754 cells
# B_cell = subset(B_cell, subset = orig.ident %in% AP_ID) # 27361 cells
# 
# Plasma <- readRDS("./1 - GBC_input data/Plasma_cell.RDS") # 42222 cells
# Plasma = subset(Plasma, cells = Cell_type[Cell_type$celltype %in% c("Plasma cell"),]$cellid) # 37237 cells
# Plasma = subset(Plasma, subset = orig.ident %in% AP_ID) # 22445 cells
# 
# T_NK <- readRDS("./1 - GBC_input data/T_NK.RDS") # 374954 cells
# T_NK = subset(T_NK, cells = Cell_type[Cell_type$celltype %in% c("CD4 T cell","CD8 T cell","NK cell"),]$cellid) # 360477 cells
# T_NK = subset(T_NK, subset = orig.ident %in% AP_ID) # 195304 cells
# 
# AdenoP_Merged = merge(Epithelial, y = c(Fibroblast, Endothelial, Myeloid_cell, B_cell, Plasma, T_NK))

## 1.4 Add Idents and MIs ####
# Cell_type <- Cell_type[Cell_type$cellid %in% colnames(AdenoP_Merged),]
# 
# AdenoP_Merged@meta.data$cellid <- rownames(AdenoP_Merged@meta.data)
# temp <- AdenoP_Merged@meta.data
# temp <- left_join(temp, Cell_type, by = "cellid")
# 
# identical(temp$cellid, rownames(AdenoP_Merged@meta.data))
# identical(temp$orig.ident.y, AdenoP_Merged@meta.data$orig.ident)
# 
# AdenoP_Merged@meta.data$Celltype <- temp$celltype.y
# AdenoP_Merged@meta.data$Subtype <- temp$subtype

# AdenoP_MI <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/AdenoP_MI.xlsx")
# colnames(AdenoP_MI)[1] <- "orig.ident"
# 
# temp <- AdenoP_Merged@meta.data
# temp <- left_join(temp, AdenoP_MI, by = "orig.ident")
# 
# identical(temp$cellid, rownames(AdenoP_Merged@meta.data))
# identical(temp$orig.ident, AdenoP_Merged@meta.data$orig.ident)
# 
# AdenoP_Merged@meta.data$MI.group <- temp$MI.group
# 
# saveRDS(AdenoP_Merged, file = "./1 - GBC_input data/AdenoP_Merged.RDS")

AdenoP_Merged <- readRDS("./1 - GBC_input data/AdenoP_Merged.RDS")

## 1.5 Find MI-specific Signatures ####
# Idents(AdenoP_Merged) <- "MI.group"
# MI_Sig_All <- FindAllMarkers(AdenoP_Merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# MI_Sig_All %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> MI_Sig_All_top20
# 
# saveRDS(MI_Sig_All, file = "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/MI_Sig_All.RDS")
# 
# MI_Sig_All_top20$probeid <- MI_Sig_All_top20$gene
# MI_Sig_All_top20$GeneName <- MI_Sig_All_top20$gene
# MI_Sig_All_top20$Class <- MI_Sig_All_top20$cluster
# MI_Sig_All_top20 <- MI_Sig_All_top20[,c("probeid","GeneName","Class")]
# MI_Sig_All_top20 = MI_Sig_All_top20[!(MI_Sig_All_top20$GeneName %in% MI_Sig_All_top20$GeneName[duplicated(MI_Sig_All_top20$GeneName)]),]
# table(MI_Sig_All_top20$Class)
# write.table(MI_Sig_All_top20, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/MI_Sig_All_top20.txt",row.names = F, col.names = T,quote = F,sep = "\t")
# 
# saveRDS(MI_Sig_All_top20, file = "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/MI_Sig_All_top20.RDS")

MI_Sig_All_top20 <- readRDS("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/MI_Sig_All_top20.RDS")

## 1.6 Transfer into mini-bulk ####
Idents(AdenoP_Merged) <- "orig.ident"
AdenoP_Merged <- NormalizeData(AdenoP_Merged)
AvgExprs_BySample <- AverageExpression(AdenoP_Merged)
normalized_matrix <- AvgExprs_BySample$RNA

metadata <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/AdenoP_MI.xlsx")
rownames(metadata) <- metadata$Patient
metadata <- metadata[colnames(normalized_matrix), ]

seurat_obj <- CreateSeuratObject(counts = normalized_matrix, meta.data = metadata)
seurat_obj[["RNA"]]@data <- normalized_matrix

## 1.7 Ident MIs by GenePattern ####
### GenePattern - gct obj ####
# 提取表达矩阵（这里使用的是归一化后的矩阵）
expr_matrix <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))

# 确保行名是基因，列名是细胞
genes <- rownames(expr_matrix)
cells <- colnames(expr_matrix)

# 创建 GCT 文件头信息
num_genes <- nrow(expr_matrix)
num_cells <- ncol(expr_matrix)

# 创建 GCT 文件的内容
gct_content <- c(
  "#1.3", 
  paste(num_genes, num_cells, sep = "\t"), 
  paste("probeid", "gene", paste(cells, collapse = "\t"), sep = "\t")
)

# 添加每个基因的表达值
for (i in 1:num_genes) {
  gene_name <- genes[i]
  gene_expr <- paste(expr_matrix[i, ], collapse = "\t")
  gct_content <- c(gct_content, paste(gene_name, gene_name, gene_expr, sep = "\t"))
}

# 将内容写入 GCT 文件
writeLines(gct_content, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/Chen_scRNA-seq_mini-bulk/sc_mini_bulk.gct")

### GenePattern - ROC ####
library(pROC)

NTP_prediction <- read.csv("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/Chen_scRNA-seq_mini-bulk/output/NTP_prediction_result.csv")
colnames(NTP_prediction)[1] <- "Patient"
NTP_prediction <- left_join(NTP_prediction, metadata, by = "Patient")

roc_curve <- roc(NTP_prediction$MI.group, NTP_prediction$predict.label) 
plot(roc_curve, col = "blue", main = "ROC Curve", lwd = 2)

# 计算AUC
auc_value <- auc(roc_curve)
cat("AUC:", auc_value, "\n")

# 在图上标注AUC值
text(0.4, 0.4, paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)

### GenePattern - Accuracy ####
equal_rows <- NTP_prediction[,"MI.group"] == NTP_prediction[,"predict.label"]
num_equal_rows <- sum(equal_rows)
accuracy_rate <- num_equal_rows/nrow(NTP_prediction)

## 1.8 Ident MIs by Score ####
gene_set_list <- list(
  gene_set1 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 1,]$GeneName,
  gene_set2 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 2,]$GeneName,
  gene_set3 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 3,]$GeneName,
  gene_set4 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 4,]$GeneName,
  gene_set5 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 5,]$GeneName
)

seurat_obj <- AddModuleScore(object = seurat_obj, features = gene_set_list)

Score_prediction <- seurat_obj@meta.data[,paste0("Cluster",1:5)]
Score_prediction <- as.data.frame(scale(Score_prediction))

Score_prediction$predicted.label <- max.col(Score_prediction)
Score_prediction$Patient <- rownames(Score_prediction)
Score_prediction <- left_join(Score_prediction, seurat_obj@meta.data, by = "Patient")

### Score - ROC ####
library(pROC)

roc_curve <- roc(Score_prediction$MI.group, Score_prediction$predicted.label) 
plot(roc_curve, col = "blue", main = "ROC Curve", lwd = 2)

# 计算AUC
auc_value <- auc(roc_curve)
cat("AUC:", auc_value, "\n")

# 在图上标注AUC值
text(0.4, 0.4, paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)

### Score - Accuracy ####
equal_rows <- Score_prediction[,"MI.group"] == Score_prediction[,"predicted.label"]
num_equal_rows <- sum(equal_rows)
accuracy_rate <- num_equal_rows/nrow(Score_prediction)

## 1.9 External Validation by Score ####
### Score ####
library(DESeq2)

# 读取count矩阵，并进行过滤
ExRNAseq <- read.table("./GBC_external RNA-seq/gene_counts.txt", header = T)
duplicated_genes <- ExRNAseq$Gene[duplicated(ExRNAseq$Gene)]
ExRNAseq <- ExRNAseq[!(ExRNAseq$Gene %in% duplicated_genes),]
rownames(ExRNAseq) <- ExRNAseq$Gene
ExRNAseq <- ExRNAseq[,-1]
ExRNAseq <- as.matrix(ExRNAseq)
ExRNAseq <- ExRNAseq[rowMeans(ExRNAseq)>1,]

# 创建metadata
col_data <- data.frame(condition = factor(rep("placeholder", ncol(ExRNAseq))), row.names = colnames(ExRNAseq))

# 构建dds对象
dds <- DESeqDataSetFromMatrix(countData = ExRNAseq, colData = col_data, design = ~ 1)

# 估计大小因子
dds <- estimateSizeFactors(dds)

# 提取标准化的计数数据
normalized_counts <- counts(dds, normalized = T)

seurat_obj <- CreateSeuratObject(counts = normalized_counts)
seurat_obj[["RNA"]]@data <- normalized_counts

MI_Sig_All_top20 <- readRDS("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/MI_Sig_All_top20.RDS")
gene_set_list <- list(
  gene_set1 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 1,]$GeneName,
  gene_set2 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 2,]$GeneName,
  gene_set3 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 3,]$GeneName,
  gene_set4 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 4,]$GeneName,
  gene_set5 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 5,]$GeneName
)
seurat_obj <- AddModuleScore(object = seurat_obj, features = gene_set_list)

Score_prediction <- seurat_obj@meta.data[,paste0("Cluster",1:5)]
Score_prediction <- as.data.frame(scale(Score_prediction))

Score_prediction$predicted.label <- max.col(Score_prediction)
Score_prediction$Patient <- rownames(Score_prediction)

### Prognosis ####
metadata <- openxlsx::read.xlsx("../2 - Analysis_V3_checked/GBC_external RNA-seq/NC_gbc_survival_GBCKorea 20210618.xlsx")
colnames(metadata)[2] <- "Patient"
metadata <- metadata[metadata$Patient %in% colnames(ExRNAseq),]
metadata <- left_join(metadata, Score_prediction, by = "Patient")

metadata$vital_status = ifelse(metadata$vital_status=="dead", 1, 0)
table(metadata$predicted.label)

metadata$group = ifelse(metadata$predicted.label == 2,2,"others")
fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=metadata)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = metadata,             # data used to fit survival curves.
  palette = c("#34ebc3","#bbbcbd"),
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

metadata$group = ifelse(metadata$Cluster2 > median(metadata$Cluster2), "High_MI2", "Low_MI2")
fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=metadata)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = metadata,             # data used to fit survival curves.
  palette = c("#34ebc3","#bbbcbd"),
  #risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in months",   # customize X axis label.
  ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

# ### CC 
# RNA_DeconResult <- readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
# selected_subtypes = c("CD4T_C2_CCR7","B_C0_IGHA1","γdT_C12_TRDC","CD8T_C0_CCR7_GZMK")
# 
# RNA_DeconResult <- RNA_DeconResult[,selected_subtypes]
# RNA_DeconResult$Patient <- rownames(RNA_DeconResult)
# RNA_DeconResult <- left_join(RNA_DeconResult, Score_prediction, by = "Patient")
# RNA_DeconResult$group_plot <- ifelse(RNA_DeconResult$predicted.label == 2, "MI2", "others")
# 
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD4T_C2_CCR7", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD4T_C2_CCR7")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "B_C0_IGHA1", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("B_C0_IGHA1")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "γdT_C12_TRDC", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("γdT_C12_TRDC")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD8T_C0_CCR7_GZMK", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD8T_C0_CCR7_GZMK")
# 
# ##
# RNA_DeconResult <- readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
# selected_subtypes = c("CD4T_C2_CCR7","B_C0_IGHA1","γdT_C12_TRDC","CD8T_C0_CCR7_GZMK")
# 
# RNA_DeconResult <- RNA_DeconResult[,selected_subtypes]
# RNA_DeconResult$Patient <- rownames(RNA_DeconResult)
# RNA_DeconResult <- left_join(RNA_DeconResult, Score_prediction, by = "Patient")
# RNA_DeconResult$group_plot <- ifelse(RNA_DeconResult$Cluster2 > median(RNA_DeconResult$Cluster2), "High_MI2", "Low_MI2")
# 
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD4T_C2_CCR7", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD4T_C2_CCR7")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "B_C0_IGHA1", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("B_C0_IGHA1")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "γdT_C12_TRDC", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("γdT_C12_TRDC")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD8T_C0_CCR7_GZMK", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD8T_C0_CCR7_GZMK")


## 1.10 External Validation by GenePattern ####
### GenePattern - gct obj ####
# 提取表达矩阵（这里使用的是归一化后的矩阵）
expr_matrix <- as.matrix(normalized_counts)

# 确保行名是基因，列名是细胞
genes <- rownames(expr_matrix)
cells <- colnames(expr_matrix)

# 创建 GCT 文件头信息
num_genes <- nrow(expr_matrix)
num_cells <- ncol(expr_matrix)

# 创建 GCT 文件的内容
gct_content <- c(
  "#1.3", 
  paste(num_genes, num_cells, sep = "\t"), 
  paste("probeid", "gene", paste(cells, collapse = "\t"), sep = "\t")
)

# 添加每个基因的表达值
for (i in 1:num_genes) {
  gene_name <- genes[i]
  gene_expr <- paste(expr_matrix[i, ], collapse = "\t")
  gct_content <- c(gct_content, paste(gene_name, gene_name, gene_expr, sep = "\t"))
}

# 将内容写入 GCT 文件
writeLines(gct_content, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/External_bulk-RNA-seq/bulk.gct")

### Prognosis ####
metadata <- openxlsx::read.xlsx("../2 - Analysis_V3_checked/GBC_external RNA-seq/NC_gbc_survival_GBCKorea 20210618.xlsx")
colnames(metadata)[2] <- "Patient"

NTP_prediction <- read.csv("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/External_bulk-RNA-seq/output/NTP_prediction_result.xls.csv") 
colnames(NTP_prediction)[1] <- "Patient"

NTP_prediction <- NTP_prediction[NTP_prediction$Patient %in% metadata$Patient,]
metadata <- left_join(metadata,NTP_prediction, metadata, by = "Patient")

metadata$vital_status = ifelse(metadata$vital_status=="dead", 1, 0)
table(metadata$predict.label)

metadata$group = ifelse(metadata$predict.label == 2, "Group2", "Others")
fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=metadata)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = metadata,             # data used to fit survival curves.
  palette = c("#34ebc3","#bbbcbd"),
  #risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in months",   # customize X axis label.
  ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

# ### CC
# RNA_DeconResult <- readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
# selected_subtypes = c("CD4T_C2_CCR7","B_C0_IGHA1","γdT_C12_TRDC","CD8T_C0_CCR7_GZMK")
# 
# RNA_DeconResult <- RNA_DeconResult[,selected_subtypes]
# RNA_DeconResult$Patient <- rownames(RNA_DeconResult)
# 
# NTP_prediction <- read.csv("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/External_bulk-RNA-seq/output/NTP_prediction_result.csv") 
# colnames(NTP_prediction)[1] <- "Patient"
# 
# RNA_DeconResult <- left_join(RNA_DeconResult, NTP_prediction, by = "Patient")
# RNA_DeconResult$group_plot <- ifelse(RNA_DeconResult$predict.label == 2, "MI2", "others")
# 
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD4T_C2_CCR7", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD4T_C2_CCR7")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "B_C0_IGHA1", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("B_C0_IGHA1")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "γdT_C12_TRDC", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("γdT_C12_TRDC")
# ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD8T_C0_CCR7_GZMK", color = "group_plot", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD8T_C0_CCR7_GZMK")




# * Fu scRNA-seq Validation - Mapping and annotating ####

# > dim(Endothelium)
# [1] 28306 21973
# > dim(`03 Fibroblast`)
# [1] 28306 25777
# > dim(Fibro)
# [1] 28306 25777
# > dim(`05 Tcell`)
# [1] 28306 85687
# > dim(`07 Myeloid`)
# [1] 28306 38159
# > dim(`06 Bcell`)
# [1] 28306 13003
# > dim(Myeloid_predict)
# [1] 28306 38159
# > (dim(B))
# [1] 28306 13003
# > dim(NK)
# [1] 28306  9002
# > dim(CD4T)
# [1] 28306 32043
# > dim(CD8T)
# [1] 28306 34062
# > 34062+32043+9002
# [1] 75107

## 1) Myeloid ####
# ## Query
# Myeloid <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/1 - Split_by_Celltype/07 Myeloid.rds")
# Myeloid <- NormalizeData(Myeloid)
# Myeloid <- FindVariableFeatures(Myeloid)
# 
# ## Ref
# Myeloid_cell <- readRDS("./1 - GBC_input data/Myeloid_cell.RDS") # 323569 cells
# Cell_type <- read.csv(file = "./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
# Myeloid_cell = subset(Myeloid_cell, cells = Cell_type[Cell_type$celltype %in% c("DC","Mast","Mono_Macro","Neu"),]$cellid) # 307285 cells
# 
# Cell_type <- Cell_type[Cell_type$cellid %in% colnames(Myeloid_cell),]
# 
# Myeloid_cell@meta.data$cellid <- rownames(Myeloid_cell@meta.data)
# temp <- Myeloid_cell@meta.data
# temp <- left_join(temp, Cell_type, by = "cellid")
# 
# identical(temp$cellid, rownames(Myeloid_cell@meta.data))
# identical(temp$orig.ident.y, Myeloid_cell@meta.data$orig.ident)
# 
# Myeloid_cell@meta.data$Celltype <- temp$celltype.y
# Myeloid_cell@meta.data$Subtype <- temp$subtype
# 
# Idents(Myeloid_cell) = "Subtype"
# 
# Myeloid_cell <- NormalizeData(Myeloid_cell)
# Myeloid_cell <- FindVariableFeatures(Myeloid_cell)
# 
# ## anchor
# anchors <- FindTransferAnchors(reference = Myeloid_cell, query = Myeloid, dims = 1:30) # Mac
# saveRDS(anchors, file = "../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/Myeloid/anchors.RDS")
# 
# ## Transfer
# predicted.labels <- TransferData(anchorset = anchors, refdata = Myeloid_cell$Subtype, dims = 1:30)
# predicted.labels$cellid <- rownames(predicted.labels)
# 
# temp <- Myeloid@meta.data
# temp$cellid <- rownames(temp)
# temp <- left_join(temp, predicted.labels, by = "cellid")
# 
# identical(temp$cellid, rownames(Myeloid@meta.data))
# Myeloid$predict.label <- temp$predicted.id
# Idents(Myeloid) <- "predict.label"
# 
# saveRDS(Myeloid, file = "../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/Myeloid/Myeloid_predict.RDS")

## 2) Endothelium ####
# ## Query 
# Endothelium <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/Split_by_Celltype/02 Endothelial.rds")
# Endothelium <- NormalizeData(Endothelium)
# Endothelium <- FindVariableFeatures(Endothelium)
# 
# ## Ref
# Endothelial_Ref <- readRDS("./1 - GBC_input data/Endothelial.RDS") # 22501 cells
# Cell_type <- read.csv(file = "./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
# Endothelial_Ref = subset(Endothelial_Ref, cells = Cell_type[Cell_type$celltype %in% c("Endothelium"),]$cellid) # 16886 cells
# 
# Cell_type <- Cell_type[Cell_type$cellid %in% colnames(Endothelial_Ref),]
# 
# Endothelial_Ref@meta.data$cellid <- rownames(Endothelial_Ref@meta.data)
# temp <- Endothelial_Ref@meta.data
# temp <- left_join(temp, Cell_type, by = "cellid")
# 
# identical(temp$cellid, rownames(Endothelial_Ref@meta.data))
# identical(temp$orig.ident.y, Endothelial_Ref@meta.data$orig.ident)
# 
# Endothelial_Ref@meta.data$Celltype <- temp$celltype.y
# Endothelial_Ref@meta.data$Subtype <- temp$subtype
# 
# Idents(Endothelial_Ref) = "Subtype"
# 
# Endothelial_Ref <- NormalizeData(Endothelial_Ref)
# Endothelial_Ref <- FindVariableFeatures(Endothelial_Ref)
# 
# ## anchor
# anchors <- FindTransferAnchors(reference = Endothelial_Ref, query = Endothelium, dims = 1:30)
# 
# ## Transfer
# predicted.labels <- TransferData(anchorset = anchors, refdata = Endothelial_Ref$Subtype, dims = 1:30)
# predicted.labels$cellid <- rownames(predicted.labels)
# 
# temp <- Endothelium@meta.data
# temp$cellid <- rownames(temp)
# temp <- left_join(temp, predicted.labels, by = "cellid")
# 
# identical(temp$cellid, rownames(Endothelium@meta.data))
# Endothelium$predict.label <- temp$predicted.id
# Idents(Endothelium) <- "predict.label"
# 
# saveRDS(Endothelium, file = "../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/Endothelium/Endothelium.RDS")

## 1.11 Fu scRNA-seq Validation by Score ####
### Transfer into mini-bulk ####
Endothelium <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/Endothelium/Endothelium.RDS")
Myeloid_predict <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/Myeloid/Myeloid_predict.RDS")

CD4T <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer CD4T.rds")
CD4T$predict.label <- CD4T$predicted.id
CD8T <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer CD8T.rds")
CD8T$predict.label <- CD8T$predicted.id
Fibro <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer fibro.rds")
Fibro$predict.label <- Fibro$predicted.id
NK <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer NK.rds")
NK$predict.label <- NK$predicted.id
B <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer B.rds")
B$predict.label <- B$predicted.id

Epithelial <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/1 - Split_by_Celltype/01 Epithelial.rds")

All <- merge(x = Myeloid_predict, y = list(CD4T, CD8T, NK, B, Endothelium, Fibro, Epithelial))
patient_select <- unique(All$orig.ident)
patient_select = patient_select[grep("T$", patient_select)]

All <- subset(All, subset = orig.ident %in% patient_select)
Idents(All) <- "orig.ident"
All <- NormalizeData(All)
All_mini <- AverageExpression(All)
normalized_matrix <- All_mini$RNA

### Score [No MI2] ####
seurat_obj <- CreateSeuratObject(counts = normalized_matrix)
seurat_obj[["RNA"]]@data <- normalized_matrix

gene_set_list <- list(
  gene_set1 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 1,]$GeneName,
  gene_set2 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 2,]$GeneName,
  gene_set3 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 3,]$GeneName,
  gene_set4 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 4,]$GeneName,
  gene_set5 = MI_Sig_All_top20[MI_Sig_All_top20$Class == 5,]$GeneName
)

seurat_obj <- AddModuleScore(object = seurat_obj, features = gene_set_list)

Score_prediction <- seurat_obj@meta.data[,paste0("Cluster",1:5)]
Score_prediction <- as.data.frame(scale(Score_prediction))

Score_prediction$predicted.label <- max.col(Score_prediction)
Score_prediction$orig.ident <- rownames(Score_prediction)

table(Score_prediction$predicted.label)

### GenePattern ####
# 提取表达矩阵（这里使用的是归一化后的矩阵）
expr_matrix <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))

# 确保行名是基因，列名是细胞
genes <- rownames(expr_matrix)
cells <- colnames(expr_matrix)

# 创建 GCT 文件头信息
num_genes <- nrow(expr_matrix)
num_cells <- ncol(expr_matrix)

# 创建 GCT 文件的内容
gct_content <- c(
  "#1.3", 
  paste(num_genes, num_cells, sep = "\t"), 
  paste("probeid", "gene", paste(cells, collapse = "\t"), sep = "\t")
)

# 添加每个基因的表达值
for (i in 1:num_genes) {
  gene_name <- genes[i]
  gene_expr <- paste(expr_matrix[i, ], collapse = "\t")
  gct_content <- c(gct_content, paste(gene_name, gene_name, gene_expr, sep = "\t"))
}

# 将内容写入 GCT 文件
writeLines(gct_content, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/Fu_scRNA-seq/sc_mini_bulk_fu.gct")

### CC ####
NTP_prediction <- read.csv("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_Sig/GenePattern/Fu_scRNA-seq/output/NTP_prediction_result_fu.csv")
colnames(NTP_prediction)[1] <- "Patient"

All <- merge(x = Myeloid_predict, y = list(CD4T, CD8T, NK, B, Endothelium, Fibro))
patient_select <- unique(All$orig.ident)
patient_select = patient_select[grep("T$", patient_select)]
All <- subset(All, subset = orig.ident %in% patient_select)

df.cell = data.frame(patient = All$orig.ident, celltype = All$predict.label)
df.cell = dcast(df.cell,patient~celltype)
rownames(df.cell) = df.cell$patient
df.cell = df.cell[,-1]

norm_func = function(x){
  x/sum(x)
}

df.cell.norm = apply(df.cell, 1, norm_func)
df.cell.norm = t(df.cell.norm) %>%as.data.frame()

selected_subtypes = c("CD4T_C2_CCR7","B_C0_IGHA1","γdT_C12_TRDC","CD8T_C0_CCR7_GZMK")
df.cell.norm <- df.cell.norm[,selected_subtypes]
df.cell.norm$Patient <- rownames(df.cell.norm)

df.cell.norm <- left_join(df.cell.norm, NTP_prediction, by = "Patient")
df.cell.norm$group_plot <- ifelse(df.cell.norm$predict.label == 2, "MI2", "others")

ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD4T_C2_CCR7", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD4T_C2_CCR7")
ggboxplot(RNA_DeconResult, x = "group_plot", y = "B_C0_IGHA1", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("B_C0_IGHA1")
ggboxplot(RNA_DeconResult, x = "group_plot", y = "γdT_C12_TRDC", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("γdT_C12_TRDC")
ggboxplot(RNA_DeconResult, x = "group_plot", y = "CD8T_C0_CCR7_GZMK", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD8T_C0_CCR7_GZMK")



# !!!!! 2 - Cellular Signature ####
## !!!!! 2.1 Find MI-specific Signatures ####
GBC_AllSubtype <- read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
GBC_TIMESubtype = GBC_AllSubtype[!(GBC_AllSubtype$celltype %in% c("tumor","normal")),]
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% unique(GBC_TIMESubtype$orig.ident),]
PatientInfo_metastasis = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c("adeno"),]
PatientInfo_metastasis = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% c("P"),]
PatientInfo_metastasis$Group = PatientInfo_metastasis$metastasis.type
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group == "P_LI", "P", PatientInfo_metastasis$Group)
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group %in% c("P_LN","P_LM"), "P_Mets", PatientInfo_metastasis$Group)
Select_Patient = PatientInfo_metastasis$NewSample.ID
cycle = Select_Patient
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(GBC_TIMESubtype, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$subtype)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

Freq_adj = scale(Freq)

Freq$Patient <- rownames(Freq)
metadata <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/AdenoP_MI.xlsx")
Freq <- left_join(Freq, metadata, by = "Patient")

subtype = c("CD4T_C1_FOXP3",
            "CD4T_C7_ISG15",
            #"CD4T_C8_MKI67",
            "CD8T_C1_CXCL13",
            #"CD8T_C8_MKI67",
            #"CD8T_C11_PCLAF",
            "γdT_C12_TRDC",
            "NKT_C3_CAPG",
            #"NK_C6_MKI67",
            "F_C0_MMP11",
            "Per_C0_RGS5",
            "DC_C0_cDC2_IL1B",
            "DC_C2_cDC3_FSCN1",
            "DC_C5_cDC1",
            "DC_C6_pDC",
            "M_C2_SPP1",
            "M_C7_AGR2",
            "N_C9_Nc",
            "EC_C2_CXCR4",
            "CD4T_C2_CCR7",
            #"CD4T_C4_Unassign",
            "CD8T_C0_CCR7_GZMK",
            "B_C0_IGHA1",
            "F_C1_CFD",
            "Per_C1_MYH11",
            "DC_C4_cDC2_PPP1R14A",
            #"M_C4_others",
            "N_C0_N0",
            "N_C1_N2",
            "N_C2_N1",
            "N_C7_N0",
            "EC_C3_GJA5",
            "EC_C5_PROX1",
            "Per_C3_STEAP4",
            "M_C1_S100A8",
            "CD8T_C9_MT1X_MT1E",
            "EC_C0_ACKR1")

# plot_function <- function(subtype) {
#   ggboxplot(Freq, x = "MI.group", y = subtype, color = "MI.group", 
#             add = "jitter") + 
#     stat_compare_means(method = 'wilcox.test',ref.group = ".all.",label = "p.signif") +
#     theme(axis.title.x = element_blank()) +
#     theme(axis.text.x = element_blank()) +
#     ggtitle(subtype)
# }
# 
# pdf(file = '../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/AdenoP_CC_Pvalue_PickFeature.pdf',
#     height = 4, width = 6)
# invisible(lapply(subtype, function(column) {
#   plot <- plot_function(column)
#   print(plot)
# }))
# dev.off()

##
Freq <- Freq[,c(subtype, "MI.group")]
Freq$MI.group <- factor(Freq$MI.group, levels = c(1,2,3,4,5))

Freq$MI1.group <- ifelse(Freq$MI.group == 1, "MI1", "Others")
Freq$MI2.group <- ifelse(Freq$MI.group == 2, "MI2", "Others")
Freq$MI3.group <- ifelse(Freq$MI.group == 3, "MI3", "Others")
Freq$MI4.group <- ifelse(Freq$MI.group == 4, "MI4", "Others")
Freq$MI5.group <- ifelse(Freq$MI.group == 5, "MI5", "Others")

p_values_1 <- data.frame()
for (i in 1:31) {
  col_name <- colnames(Freq)[i]
  
  test_result <- compare_means(as.formula(paste(col_name, "~ MI1.group")), data = Freq, method = 'wilcox.test', ref.group = "Others")
  
  p_values_1 <- rbind(p_values_1, test_result)
}
p_values_1$compare <- "MI1_vs_Others"

p_values_2 <- data.frame()
for (i in 1:31) {
  col_name <- colnames(Freq)[i]
  
  test_result <- compare_means(as.formula(paste(col_name, "~ MI2.group")), data = Freq, method = 'wilcox.test', ref.group = "Others")
  
  p_values_2 <- rbind(p_values_2, test_result)
}
p_values_2$compare <- "MI2_vs_Others"

p_values_3 <- data.frame()
for (i in 1:31) {
  col_name <- colnames(Freq)[i]
  
  test_result <- compare_means(as.formula(paste(col_name, "~ MI3.group")), data = Freq, method = 'wilcox.test', ref.group = "Others")
  
  p_values_3 <- rbind(p_values_3, test_result)
}
p_values_3$compare <- "MI3_vs_Others"

p_values_4 <- data.frame()
for (i in 1:31) {
  col_name <- colnames(Freq)[i]
  
  test_result <- compare_means(as.formula(paste(col_name, "~ MI4.group")), data = Freq, method = 'wilcox.test', ref.group = "Others")
  
  p_values_4 <- rbind(p_values_4, test_result)
}
p_values_4$compare <- "MI4_vs_Others"

p_values_5 <- data.frame()
for (i in 1:31) {
  col_name <- colnames(Freq)[i]
  
  test_result <- compare_means(as.formula(paste(col_name, "~ MI5.group")), data = Freq, method = 'wilcox.test', ref.group = "Others")
  
  p_values_5 <- rbind(p_values_5, test_result)
}
p_values_5$compare <- "MI5_vs_Others"

p_values <- rbind(p_values_1,
                  p_values_2,
                  p_values_3,
                  p_values_4,
                  p_values_5)

mean_values_1 <- data.frame()
for (i in 1:31) {
  df <- Freq[,c(i,33)]
  colnames(df)[1] <- "value"
  
  group_means <- df %>% group_by(MI1.group) %>% summarise(mean = mean(value, na.rm = TRUE))
  
  group_means$Subtype <- colnames(Freq)[i]
  
  mean_values_1 <- rbind(mean_values_1, group_means)
}
mean_values_1 <- dcast(mean_values_1, Subtype ~ MI1.group, value.var = "mean")
mean_values_1$Sig.group <- ifelse(mean_values_1$MI1 > mean_values_1$Others, "Up",
                                  ifelse(mean_values_1$MI1 < mean_values_1$Others, "Down", "Equal"))
mean_values_1$MI.Cat <- "MI1"

mean_values_2 <- data.frame()
for (i in 1:31) {
  df <- Freq[,c(i,34)]
  colnames(df)[1] <- "value"
  
  group_means <- df %>% group_by(MI2.group) %>% summarise(mean = mean(value, na.rm = TRUE))
  
  group_means$Subtype <- colnames(Freq)[i]
  
  mean_values_2 <- rbind(mean_values_2, group_means)
}
mean_values_2 <- dcast(mean_values_2, Subtype ~ MI2.group, value.var = "mean")
mean_values_2$Sig.group <- ifelse(mean_values_2$MI2 > mean_values_2$Others, "Up",
                                  ifelse(mean_values_2$MI2 < mean_values_2$Others, "Down", "Equal"))
mean_values_2$MI.Cat <- "MI2"

mean_values_3 <- data.frame()
for (i in 1:31) {
  df <- Freq[,c(i,35)]
  colnames(df)[1] <- "value"
  
  group_means <- df %>% group_by(MI3.group) %>% summarise(mean = mean(value, na.rm = TRUE))
  
  group_means$Subtype <- colnames(Freq)[i]
  
  mean_values_3 <- rbind(mean_values_3, group_means)
}
mean_values_3 <- dcast(mean_values_3, Subtype ~ MI3.group, value.var = "mean")
mean_values_3$Sig.group <- ifelse(mean_values_3$MI3 > mean_values_3$Others, "Up",
                                  ifelse(mean_values_3$MI3 < mean_values_3$Others, "Down", "Equal"))
mean_values_3$MI.Cat <- "MI3"

mean_values_4 <- data.frame()
for (i in 1:31) {
  df <- Freq[,c(i,36)]
  colnames(df)[1] <- "value"
  
  group_means <- df %>% group_by(MI4.group) %>% summarise(mean = mean(value, na.rm = TRUE))
  
  group_means$Subtype <- colnames(Freq)[i]
  
  mean_values_4 <- rbind(mean_values_4, group_means)
}
mean_values_4 <- dcast(mean_values_4, Subtype ~ MI4.group, value.var = "mean")
mean_values_4$Sig.group <- ifelse(mean_values_4$MI4 > mean_values_4$Others, "Up",
                                  ifelse(mean_values_4$MI4 < mean_values_4$Others, "Down", "Equal"))
mean_values_4$MI.Cat <- "MI4"

mean_values_5 <- data.frame()
for (i in 1:31) {
  df <- Freq[,c(i,37)]
  colnames(df)[1] <- "value"
  
  group_means <- df %>% group_by(MI5.group) %>% summarise(mean = mean(value, na.rm = TRUE))
  
  group_means$Subtype <- colnames(Freq)[i]
  
  mean_values_5 <- rbind(mean_values_5, group_means)
}
mean_values_5 <- dcast(mean_values_5, Subtype ~ MI5.group, value.var = "mean")
mean_values_5$Sig.group <- ifelse(mean_values_5$MI5 > mean_values_5$Others, "Up",
                                  ifelse(mean_values_5$MI5 < mean_values_5$Others, "Down", "Equal"))
mean_values_5$MI.Cat <- "MI5"

mean_values <- rbind(mean_values_1[,c(1,4,5)],
                     mean_values_2[,c(1,4,5)],
                     mean_values_3[,c(1,4,5)],
                     mean_values_4[,c(1,4,5)],
                     mean_values_5[,c(1,4,5)])

p_values$leftjoin <- paste0(p_values$.y., p_values$group2)
mean_values$leftjoin <- paste0(mean_values$Subtype, mean_values$MI.Cat)

p_values <- left_join(p_values, mean_values, by = "leftjoin")

p_values$p_value_log <- -log10(p_values$p)

subtype_order = c("CD4T_C1_FOXP3",
                  "CD8T_C1_CXCL13",
                  "CD4T_C2_CCR7",
                  "B_C0_IGHA1",
                  "γdT_C12_TRDC",
                  "CD8T_C0_CCR7_GZMK",
                  "EC_C3_GJA5",
                  "EC_C0_ACKR1",
                  "F_C1_CFD",
                  "Per_C1_MYH11",
                  "Per_C3_STEAP4",
                  "DC_C0_cDC2_IL1B",
                  "DC_C4_cDC2_PPP1R14A",
                  "DC_C5_cDC1",
                  "M_C1_S100A8",
                  "N_C1_N2",
                  "N_C2_N1",
                  "N_C0_N0",
                  "N_C7_N0",
                  "NKT_C3_CAPG",
                  "N_C9_Nc",
                  "DC_C2_cDC3_FSCN1", 
                  "CD4T_C7_ISG15",
                  "DC_C6_pDC",
                  "M_C7_AGR2",
                  "Per_C0_RGS5",
                  "EC_C2_CXCR4",
                  "M_C2_SPP1",
                  "EC_C5_PROX1",
                  "F_C0_MMP11",
                  "CD8T_C9_MT1X_MT1E"
            )
p_values$Subtype <- factor(p_values$Subtype, levels = rev(subtype_order))

p_values$Sig.group1 <- ifelse(p_values$p < 0.05, p_values$Sig.group, "NS")
p_values$Sig.group1 <- factor(p_values$Sig.group1, levels = c("Up","Down","NS"))

p <- ggplot(p_values, aes(x = compare, y = Subtype, size = p_value_log)) +
  geom_point(aes(color = Sig.group1), alpha = 0.7) +  # 添加气泡
  scale_size_continuous(range = c(1, 10)) +  # 设置气泡大小范围
  labs(title = "MI-specific subtypes", x = "", size = "-Log10(P)") +  # 添加标题和标签
  theme_minimal()  # 使用简约主题

pdf(file = '../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/AdenoP_CC_Pvalue_PickFeature_summary_241125.pdf',
    height = 10, 
    width = 5.5)
print(p)
dev.off()

write.csv(p_values, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/AdenoP_CC_Pvalue_PickFeature_summary_241125.csv")

## 2.2 Find MI-specific Signatures ####
MI5_Feature <- c("CD4T_C1_FOXP3","CD8T_C1_CXCL13")
MI2_Feature <- c("γdT_C12_TRDC","CD4T_C2_CCR7","CD8T_C0_CCR7_GZMK","B_C0_IGHA1")
MI1_Feature <- c("DC_C0_cDC2_IL1B","DC_C5_cDC1","F_C1_CFD","Per_C1_MYH11","DC_C4_cDC2_PPP1R14A","EC_C3_GJA5","Per_C3_STEAP4","EC_C0_ACKR1")
MI4_Feature <- c("M_C2_SPP1","M_C7_AGR2","EC_C2_CXCR4","Per_C0_RGS5")
MI3_Feature <- c("N_C9_Nc","N_C0_N0","N_C1_N2","N_C2_N1","N_C7_N0")

MI5_Feature <- c("CD4T-C1-FOXP3","CD8T-C1-CXCL13")
MI2_Feature <- c("γdT-C12-TRDC","CD4T-C2-CCR7","CD8T-C0-CCR7-GZMK","B-C0-IGHA1")
MI1_Feature <- c("DC-C0-cDC2-IL1B","DC-C5-cDC1","F-C1-CFD","Per-C1-MYH11","DC-C4-cDC2-PPP1R14A","EC-C3-GJA5","Per-C3-STEAP4","EC-C0-ACKR1")
MI4_Feature <- c("M-C2-SPP1","M-C7-AGR2","EC-C2-CXCR4","Per-C0-RGS5")
MI3_Feature <- c("N-C9-Nc","N-C0-N0","N-C1-N2","N-C2-N1","N-C7-N0")

# ## 2.3 Transfer into mini-bulk - CIBERSORT input
# PatientInfo_New <- readRDS("./1 - GBC_input data/PatientInfo_New/PatientInfo_New.RDS")
# AdenoP_Merged <- readRDS("./1 - GBC_input data/AdenoP_Merged.RDS")
# 
# Idents(AdenoP_Merged) <- "orig.ident"
# AdenoP_Merged <- NormalizeData(AdenoP_Merged)
# AvgExprs_BySample <- AverageExpression(AdenoP_Merged)
# normalized_matrix <- AvgExprs_BySample$RNA
# 
# normalized_matrix <- cbind(Gene = rownames(normalized_matrix), normalized_matrix)
# write.table(normalized_matrix, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/Chen_scRNA-seq_mini-bulk/gene_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# # CIBERSORTx

## 2.3 Features - PCA ####
## subtype
Freq_adj_pca <- Freq_adj[,subtype]
pca_result <- prcomp(Freq_adj_pca, center = TRUE, scale = TRUE)

# 转换主成分得分为数据框
pca_data <- as.data.frame(pca_result$x)
pca_data$Patient <- rownames(pca_data)
metadata <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/AdenoP_MI.xlsx")
pca_data <- left_join(pca_data,metadata,by="Patient")
pca_data$MI.group <- paste0("MI", pca_data$MI.group)

# 绘制 PCA 散点图（PC1 vs PC2）
ggplot(pca_data, aes(x = PC1, y = PC2, color = MI.group)) +
  geom_point(size = 3) +
  labs(title = "PCA Scatter Plot", x = "PC1", y = "PC2") +
  theme_minimal() +
  geom_text(aes(label = rownames(pca_data)), vjust = -1, hjust = 1)

## features
Freq_adj_pca <- Freq_adj[,c(MI1_Feature,
                            MI2_Feature,
                            MI3_Feature,
                            MI4_Feature,
                            MI5_Feature)]
pca_result <- prcomp(Freq_adj_pca, center = TRUE, scale = TRUE)

# 转换主成分得分为数据框
pca_data <- as.data.frame(pca_result$x)
pca_data$Patient <- rownames(pca_data)
pca_data <- left_join(pca_data,metadata,by="Patient")
pca_data$MI.group <- paste0("MI", pca_data$MI.group)

# 绘制 PCA 散点图（PC1 vs PC2）
p = ggplot(pca_data, aes(x = PC1, y = PC2, color = MI.group)) +
  geom_point() +
  labs(title = "PCA Scatter Plot", x = "PC1", y = "PC2") +
  theme_minimal() 
  #geom_text(aes(label = rownames(pca_data)), vjust = -1, hjust = 1)

pdf(file = '../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/AdenoP_CC_Pvalue_Feature_PCA_241125.pdf',
    height = 4, 
    width = 5)
print(p)
dev.off()

## 2.4 Ident MIs by Score ####
GBC_AllSubtype <- read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
GBC_TIMESubtype = GBC_AllSubtype[!(GBC_AllSubtype$celltype %in% c("tumor","normal")),]
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% unique(GBC_TIMESubtype$orig.ident),]
PatientInfo_metastasis = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c("adeno"),]
PatientInfo_metastasis = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% c("P"),]
PatientInfo_metastasis$Group = PatientInfo_metastasis$metastasis.type
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group == "P_LI", "P", PatientInfo_metastasis$Group)
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group %in% c("P_LN","P_LM"), "P_Mets", PatientInfo_metastasis$Group)
Select_Patient = PatientInfo_metastasis$NewSample.ID
cycle = Select_Patient
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(GBC_TIMESubtype, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$subtype)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

ExternalCohort_PctScore = CreateSeuratObject(t(Freq))

MI5_Feature <- c("CD4T-C1-FOXP3","CD8T-C1-CXCL13")
MI2_Feature <- c("γdT-C12-TRDC","CD4T-C2-CCR7","CD8T-C0-CCR7-GZMK","B-C0-IGHA1")
MI1_Feature <- c("DC-C0-cDC2-IL1B","DC-C5-cDC1","F-C1-CFD","Per-C1-MYH11","DC-C4-cDC2-PPP1R14A","EC-C3-GJA5","Per-C3-STEAP4","EC-C0-ACKR1")
MI4_Feature <- c("M-C2-SPP1","M-C7-AGR2","EC-C2-CXCR4","Per-C0-RGS5")
MI3_Feature <- c("N-C9-Nc","N-C0-N0","N-C1-N2","N-C2-N1","N-C7-N0")

gene_set_list <- list(
  gene_set1 = MI1_Feature,
  gene_set2 = MI2_Feature,
  gene_set3 = MI3_Feature,
  gene_set4 = MI4_Feature,
  gene_set5 = MI5_Feature
)

ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = gene_set_list, ctrl = 3)
Score_prediction <- ExternalCohort_PctScore@meta.data[,paste0("Cluster",1:5)]

Score_prediction <- as.data.frame(scale(Score_prediction))
Score_prediction$predicted.label <- max.col(Score_prediction)
Score_prediction$Patient <- rownames(Score_prediction)

metadata <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/AdenoP_MI.xlsx")
rownames(metadata) <- metadata$Patient

Score_prediction <- left_join(Score_prediction, metadata, by = "Patient")

### Score - ROC ####
library(pROC)

roc_curve <- roc(Score_prediction$MI.group, Score_prediction$predicted.label) 
plot(roc_curve, col = "blue", main = "ROC Curve", lwd = 2)

# 计算AUC
auc_value <- auc(roc_curve)
cat("AUC:", auc_value, "\n")

# 在图上标注AUC值
text(0.4, 0.4, paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)

### Score - Accuracy ####
equal_rows <- Score_prediction[,"MI.group"] == Score_prediction[,"predicted.label"]
num_equal_rows <- sum(equal_rows)
accuracy_rate <- num_equal_rows/nrow(Score_prediction)

### Pseudo - Processing ####
GBC_AllSubtype <- read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
GBC_TIMESubtype = GBC_AllSubtype[!(GBC_AllSubtype$celltype %in% c("tumor","normal")),]
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% unique(GBC_TIMESubtype$orig.ident),]
PatientInfo_metastasis = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c("adeno"),]
PatientInfo_metastasis = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% c("P"),]
PatientInfo_metastasis$Group = PatientInfo_metastasis$metastasis.type
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group == "P_LI", "P", PatientInfo_metastasis$Group)
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group %in% c("P_LN","P_LM"), "P_Mets", PatientInfo_metastasis$Group)
Select_Patient = PatientInfo_metastasis$NewSample.ID

MI5_Feature <- c("CD4T_C1_FOXP3","CD8T_C1_CXCL13")
MI2_Feature <- c("γdT_C12_TRDC","CD4T_C2_CCR7","CD8T_C0_CCR7_GZMK","B_C0_IGHA1")
MI1_Feature <- c("DC_C0_cDC2_IL1B","DC_C5_cDC1","F_C1_CFD","Per_C1_MYH11","DC_C4_cDC2_PPP1R14A","EC_C3_GJA5","Per_C3_STEAP4","EC_C0_ACKR1")
MI4_Feature <- c("M_C2_SPP1","M_C7_AGR2","EC_C2_CXCR4","Per_C0_RGS5")
MI3_Feature <- c("N_C9_Nc","N_C0_N0","N_C1_N2","N_C2_N1","N_C7_N0")
subtype = c(MI1_Feature,MI2_Feature,MI3_Feature,MI4_Feature,MI5_Feature)

set.seed(925)
cycle = Select_Patient
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(GBC_TIMESubtype, orig.ident == cycle[i])
  
  sub_fix <- sub[sub$subtype %in% subtype,]
  sub_sample <- sub[!(sub$subtype %in% subtype),]
  
  values <- sample(1:100,length(unique(sub_sample$subtype)))
  values <- values / sum(values)
  sub_sample$subtype <- sample(unique(sub_sample$subtype), nrow(sub_sample), replace = TRUE, prob = values)
  # sub_sample$subtype <- sample(unique(sub_sample$subtype), nrow(sub_sample), replace = TRUE)
  
  sub_pseudo <- rbind(sub_fix,sub_sample)
  
  Freq = bind_rows(Freq, round(table(sub_pseudo$subtype)/nrow(sub_pseudo),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0


ExternalCohort_PctScore = CreateSeuratObject(t(Freq))

MI5_Feature <- c("CD4T-C1-FOXP3","CD8T-C1-CXCL13")
MI2_Feature <- c("γdT-C12-TRDC","CD4T-C2-CCR7","CD8T-C0-CCR7-GZMK","B-C0-IGHA1")
MI1_Feature <- c("DC-C0-cDC2-IL1B","DC-C5-cDC1","F-C1-CFD","Per-C1-MYH11","DC-C4-cDC2-PPP1R14A","EC-C3-GJA5","Per-C3-STEAP4","EC-C0-ACKR1")
MI4_Feature <- c("M-C2-SPP1","M-C7-AGR2","EC-C2-CXCR4","Per-C0-RGS5")
MI3_Feature <- c("N-C9-Nc","N-C0-N0","N-C1-N2","N-C2-N1","N-C7-N0")

gene_set_list <- list(
  gene_set1 = MI1_Feature,
  gene_set2 = MI2_Feature,
  gene_set3 = MI3_Feature,
  gene_set4 = MI4_Feature,
  gene_set5 = MI5_Feature
)

ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = gene_set_list, ctrl = 3)
Score_prediction <- ExternalCohort_PctScore@meta.data[,paste0("Cluster",1:5)]

Score_prediction <- as.data.frame(scale(Score_prediction))
Score_prediction$predicted.label <- max.col(Score_prediction)
Score_prediction$Patient <- rownames(Score_prediction)

metadata <- openxlsx::read.xlsx("./1 - GBC_input data/PatientInfo_New/AdenoP_MI.xlsx")
rownames(metadata) <- metadata$Patient

Score_prediction <- left_join(Score_prediction, metadata, by = "Patient")

### Score - ROC ####
library(pROC)

roc_curve <- roc(Score_prediction$MI.group, Score_prediction$predicted.label) 
plot(roc_curve, col = "blue", main = "ROC Curve", lwd = 2)

# 计算AUC
auc_value <- auc(roc_curve)
cat("AUC:", auc_value, "\n")

# 在图上标注AUC值
text(0.4, 0.4, paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)

### Score - Accuracy ####
equal_rows <- Score_prediction[,"MI.group"] == Score_prediction[,"predicted.label"]
num_equal_rows <- sum(equal_rows)
accuracy_rate <- num_equal_rows/nrow(Score_prediction)

## plot
table(Score_prediction[Score_prediction$MI.group == 1,]$predicted.label)
table(Score_prediction[Score_prediction$MI.group == 2,]$predicted.label)
table(Score_prediction[Score_prediction$MI.group == 3,]$predicted.label)
table(Score_prediction[Score_prediction$MI.group == 4,]$predicted.label)
table(Score_prediction[Score_prediction$MI.group == 5,]$predicted.label)

plot_accuracy <- matrix(c(13,0,0,1,0,
                              2,14,1,0,1,
                              0,1,8,0,2,
                              0,0,4,9,3,
                              0,1,0,1,14), nrow=5, byrow = T)
formatted_numbers <- matrix(as.integer(plot_accuracy), nrow = nrow(plot_accuracy))
p = pheatmap(plot_accuracy, 
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c("lightgrey", "purple"))(50),
         display_numbers = formatted_numbers)
pdf(file = '../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/AdenoP_CC_Pvalue_Feature_Accuracy_ColPred_241125.pdf',
    height = 2.5, 
    width = 3)
print(p)
dev.off()

write.csv(plot_accuracy, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/AdenoP_CC_Pvalue_Feature_Accuracy_ColPred_241125.csv")

# ## 2.5 Ident MIs by GenePattern ####
# ### GenePattern - gct obj ####
# # 提取表达矩阵（这里使用的是归一化后的矩阵）
# expr_matrix <- as.matrix(GetAssayData(ExternalCohort_PctScore, assay = "RNA", slot = "data"))
# 
# # 确保行名是基因，列名是细胞
# genes <- rownames(expr_matrix)
# cells <- colnames(expr_matrix)
# 
# # 创建 GCT 文件头信息
# num_genes <- nrow(expr_matrix)
# num_cells <- ncol(expr_matrix)
# 
# # 创建 GCT 文件的内容
# gct_content <- c(
#   "#1.3", 
#   paste(num_genes, num_cells, sep = "\t"), 
#   paste("probeid", "gene", paste(cells, collapse = "\t"), sep = "\t")
# )
# 
# # 添加每个基因的表达值
# for (i in 1:num_genes) {
#   gene_name <- genes[i]
#   gene_expr <- paste(expr_matrix[i, ], collapse = "\t")
#   gct_content <- c(gct_content, paste(gene_name, gene_name, gene_expr, sep = "\t"))
# }
# 
# # 将内容写入 GCT 文件
# writeLines(gct_content, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/GenePattern/sc_mini_bulk_cc.gct")
# 
# ##
# MI_Sig_CC <- data.frame(probeid = c(MI1_Feature, MI2_Feature, MI3_Feature, MI4_Feature, MI5_Feature),
#                         GeneName = c(MI1_Feature, MI2_Feature, MI3_Feature, MI4_Feature, MI5_Feature),
#                         Class = c(rep("1", length(MI1_Feature)),
#                                   rep("2", length(MI2_Feature)),
#                                   rep("3", length(MI3_Feature)),
#                                   rep("4", length(MI4_Feature)),
#                                   rep("5", length(MI5_Feature))))
# write.table(MI_Sig_CC, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/GenePattern/MI_Sig_CC.txt",row.names = F, col.names = T,quote = F,sep = "\t")
# 
# ### GenePattern - ROC ####
# library(pROC)
# 
# NTP_prediction <- read.csv("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/GenePattern/output/NTP_prediction_result.xls.csv")
# colnames(NTP_prediction)[1] <- "Patient"
# NTP_prediction <- left_join(NTP_prediction, metadata, by = "Patient")
# 
# roc_curve <- roc(NTP_prediction$MI.group, NTP_prediction$predict.label) 
# plot(roc_curve, col = "blue", main = "ROC Curve", lwd = 2)
# 
# # 计算AUC
# auc_value <- auc(roc_curve)
# cat("AUC:", auc_value, "\n")
# 
# # 在图上标注AUC值
# text(0.4, 0.4, paste("AUC =", round(auc_value, 3)), col = "red", cex = 1.2)
# 
# ### GenePattern - Accuracy ####
# equal_rows <- NTP_prediction[,"MI.group"] == NTP_prediction[,"predict.label"]
# num_equal_rows <- sum(equal_rows)
# accuracy_rate <- num_equal_rows/nrow(NTP_prediction)

# ## Explain the low AUC values ####
# ## Counts
# ExRNAseq <- read.table("./GBC_external RNA-seq/gene_counts.txt", header = T)
# colSums(ExRNAseq[,2:ncol(ExRNAseq)])
# 
# AdenoP_Merged <- readRDS("./1 - GBC_input data/AdenoP_Merged.RDS")
# Idents(AdenoP_Merged) <- "orig.ident"
# AvgExprs_BySample <- AverageExpression(AdenoP_Merged, slot = "counts")
# sc_counts <- AvgExprs_BySample$RNA
# colSums(sc_counts)



# ## 2.7 External Validation by Score ####
# RNA_DeconResult = readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
# RNA_DeconResult = RNA_DeconResult[,1:119]
# RNA_DeconResult = as.data.frame(t(RNA_DeconResult))
# RNA_DeconResult = apply(RNA_DeconResult,2,function(x){x/sum(x)})
# 
# ExternalCohort_PctScore = CreateSeuratObject(RNA_DeconResult)
# 
# MI5_Feature <- c("CD4T-C1-FOXP3","CD8T-C1-CXCL13")
# MI2_Feature <- c("γdT-C12-TRDC","CD4T-C2-CCR7","CD8T-C0-CCR7-GZMK","B-C0-IGHA1")
# MI1_Feature <- c("DC-C0-cDC2-IL1B","DC-C5-cDC1","F-C1-CFD","Per-C1-MYH11","DC-C4-cDC2-PPP1R14A","EC-C3-GJA5","Per-C3-STEAP4","EC-C0-ACKR1")
# MI4_Feature <- c("M-C2-SPP1","M-C7-AGR2","EC-C2-CXCR4")
# MI3_Feature <- c("N-C9-Nc","N-C0-N0","N-C1-N2","N-C2-N1","N-C7-N0")
# 
# gene_set_list <- list(
#   gene_set1 = MI1_Feature,
#   gene_set2 = MI2_Feature,
#   gene_set3 = MI3_Feature,
#   gene_set4 = MI4_Feature,
#   gene_set5 = MI5_Feature
# )
# 
# ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = gene_set_list, ctrl = 3)
# Score_prediction <- ExternalCohort_PctScore@meta.data[,paste0("Cluster",1:5)]
# Score_prediction$predicted.label <- max.col(Score_prediction)
# Score_prediction$Patient <- rownames(Score_prediction)
# 
# metadata <- openxlsx::read.xlsx("../2 - Analysis_V3_checked/GBC_external RNA-seq/NC_gbc_survival_GBCKorea 20210618.xlsx")
# colnames(metadata)[2] <- "Patient"
# metadata <- left_join(metadata, Score_prediction, by = "Patient")
# 
# metadata$vital_status = ifelse(metadata$vital_status=="dead", 1, 0)
# table(metadata$predicted.label)
# 
# metadata$group = ifelse(metadata$predicted.label == 2,2,"others")
# fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=metadata)
# ggsurvplot(
#   fit,                     # survfit object with calculated statistics.
#   data = metadata,             # data used to fit survival curves.
#   palette = c("#34ebc3","#bbbcbd"),
#   risk.table = TRUE,       # show risk table.
#   pval = TRUE,             # show p-value of log-rank test.
#   #conf.int = TRUE,         # show confidence intervals for
#   # palette = "npg",
#   xlab = "Time in days",   # customize X axis label.
#   ggtheme = theme_classic(),
#   #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
#   #surv.median.line = "hv",  # add the median survival pointer.
#   # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
# )
# 
# RNA_DeconResult <- t(RNA_DeconResult)
# RNA_DeconResult <- as.data.frame(RNA_DeconResult)
# RNA_DeconResult$Patient <- rownames(RNA_DeconResult)
# Score_prediction <- left_join(Score_prediction, RNA_DeconResult, by="Patient")
# Score_prediction$group = ifelse(Score_prediction$predicted.label == 2,2,"others")
# 
# ggboxplot(Score_prediction, x = "group", y = "CD4T_C2_CCR7", color = "group", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD4T_C2_CCR7")
# ggboxplot(Score_prediction, x = "group", y = "B_C0_IGHA1", color = "group", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("B_C0_IGHA1")
# ggboxplot(Score_prediction, x = "group", y = "γdT_C12_TRDC", color = "group", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("γdT_C12_TRDC")
# ggboxplot(Score_prediction, x = "group", y = "CD8T_C0_CCR7_GZMK", color = "group", 
#           add = "jitter") + 
#   stat_compare_means(method = 't.test') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.text.x = element_blank()) +
#   ggtitle("CD8T_C0_CCR7_GZMK")

## 2.7 External Validation by Score - scale ####
RNA_DeconResult = readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
RNA_DeconResult = RNA_DeconResult[,1:119]
RNA_DeconResult = as.data.frame(t(RNA_DeconResult))
RNA_DeconResult = apply(RNA_DeconResult,2,function(x){x/sum(x)})

ExternalCohort_PctScore = CreateSeuratObject(RNA_DeconResult)

MI5_Feature <- c("CD4T-C1-FOXP3","CD8T-C1-CXCL13")
MI2_Feature <- c("γdT-C12-TRDC","CD4T-C2-CCR7","CD8T-C0-CCR7-GZMK","B-C0-IGHA1")
MI1_Feature <- c("DC-C0-cDC2-IL1B","DC-C5-cDC1","F-C1-CFD","Per-C1-MYH11","DC-C4-cDC2-PPP1R14A","EC-C3-GJA5","Per-C3-STEAP4","EC-C0-ACKR1")
MI4_Feature <- c("M-C2-SPP1","M-C7-AGR2","EC-C2-CXCR4","Per-C0-RGS5")
MI3_Feature <- c("N-C9-Nc","N-C0-N0","N-C1-N2","N-C2-N1","N-C7-N0")

gene_set_list <- list(
  gene_set1 = MI1_Feature,
  gene_set2 = MI2_Feature,
  gene_set3 = MI3_Feature,
  gene_set4 = MI4_Feature,
  gene_set5 = MI5_Feature
)

ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = gene_set_list, ctrl = 3)
Score_prediction <- ExternalCohort_PctScore@meta.data[,paste0("Cluster",1:5)]
Score_prediction <- as.data.frame(scale(Score_prediction))
Score_prediction$predicted.label <- max.col(Score_prediction)
Score_prediction$Patient <- rownames(Score_prediction)

metadata <- openxlsx::read.xlsx("../2 - Analysis_V3_checked/GBC_external RNA-seq/NC_gbc_survival_GBCKorea 20210618.xlsx")
colnames(metadata)[2] <- "Patient"
metadata <- left_join(metadata, Score_prediction, by = "Patient")

metadata$vital_status = ifelse(metadata$vital_status=="dead", 1, 0)
table(metadata$predicted.label)

c("MI1" = "#915F19",
  "MI2" = "#FF7F0E",
  "MI3" = "#D754DE",
  "MI4" = "#09A5ED",
  "MI5" = "#34EBC3")

metadata$group = ifelse(metadata$predicted.label == 2,"MI2","other MIs")
fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=metadata)
p = ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = metadata,             # data used to fit survival curves.
  palette = c("#FF7F0E","grey"),
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in months",   # customize X axis label.
  risk.table.pos = "in",
  risk.table.col = "strata",
  fontsize = 4,
  pval.size = 4,
  pval.coord = c(9, 0.9),
  legend = "top",
  title = "MI2"
  #ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

pdf(file = paste0("../Publication/240923_Nature Genetics_Revise/2 - Analysis/External_Prognosis_MI2_241218", ".pdf"), width =6, height = 6, onefile = F) #
print(p)
dev.off()

RNA_DeconResult <- t(RNA_DeconResult)
RNA_DeconResult <- as.data.frame(RNA_DeconResult)
RNA_DeconResult$Patient <- rownames(RNA_DeconResult)
Score_prediction <- left_join(Score_prediction, RNA_DeconResult, by="Patient")
Score_prediction$group = ifelse(Score_prediction$predicted.label == 2,2,"others")

p = ggboxplot(Score_prediction, x = "group", y = "CD4T_C2_CCR7", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 'wilcox') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD4T_C2_CCR7")

pdf(file = paste0("../Publication/240923_Nature Genetics_Revise/2 - Analysis/External_Prognosis_MI2_CD4T_C2_CCR7_241125", ".pdf"), width =3.5, height = 4, onefile = F) #
print(p)
dev.off()

p = ggboxplot(Score_prediction, x = "group", y = "B_C0_IGHA1", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 'wilcox') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("B_C0_IGHA1")

pdf(file = paste0("../Publication/240923_Nature Genetics_Revise/2 - Analysis/External_Prognosis_MI2_B_C0_IGHA1_241125", ".pdf"), width =3.5, height = 4, onefile = F) #
print(p)
dev.off()

p = ggboxplot(Score_prediction, x = "group", y = "γdT_C12_TRDC", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 'wilcox') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("γdT_C12_TRDC")

pdf(file = paste0("../Publication/240923_Nature Genetics_Revise/2 - Analysis/External_Prognosis_MI2_gdT_C12_TRDC_241125", ".pdf"), width =3.5, height = 4, onefile = F) #
print(p)
dev.off()

p = ggboxplot(Score_prediction, x = "group", y = "CD8T_C0_CCR7_GZMK", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 'wilcox') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD8T_C0_CCR7_GZMK")

pdf(file = paste0("../Publication/240923_Nature Genetics_Revise/2 - Analysis/External_Prognosis_MI2_CD8T_C0_CCR7_GZMK_241125", ".pdf"), width =3.5, height = 4, onefile = F) #
print(p)
dev.off()

## 2.7 External Validation by Score - Fixed scale ####
## scRNA-seq
Score_prediction <- ExternalCohort_PctScore@meta.data[,paste0("Cluster",1:5)]

sc_score_chen <- data.frame(mean = colMeans(Score_prediction),
                      sd = apply(Score_prediction, 2, sd))

## Bulk
RNA_DeconResult = readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
RNA_DeconResult = RNA_DeconResult[,1:119]
RNA_DeconResult = as.data.frame(t(RNA_DeconResult))
RNA_DeconResult = apply(RNA_DeconResult,2,function(x){x/sum(x)})

ExternalCohort_PctScore = CreateSeuratObject(RNA_DeconResult)

MI5_Feature <- c("CD4T-C1-FOXP3","CD8T-C1-CXCL13")
MI2_Feature <- c("γdT-C12-TRDC","CD4T-C2-CCR7","CD8T-C0-CCR7-GZMK","B-C0-IGHA1")
MI1_Feature <- c("DC-C0-cDC2-IL1B","DC-C5-cDC1","F-C1-CFD","Per-C1-MYH11","DC-C4-cDC2-PPP1R14A","EC-C3-GJA5","Per-C3-STEAP4","EC-C0-ACKR1")
MI4_Feature <- c("M-C2-SPP1","M-C7-AGR2","EC-C2-CXCR4")
MI3_Feature <- c("N-C9-Nc","N-C0-N0","N-C1-N2","N-C2-N1","N-C7-N0")

gene_set_list <- list(
  gene_set1 = MI1_Feature,
  gene_set2 = MI2_Feature,
  gene_set3 = MI3_Feature,
  gene_set4 = MI4_Feature,
  gene_set5 = MI5_Feature
)

ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = gene_set_list, ctrl = 3)
Score_prediction <- ExternalCohort_PctScore@meta.data[,paste0("Cluster",1:5)]

sc_mini_cibersort_z <- c()
for (i in 1:5) {
  mean_value <- sc_score_chen$mean[i]
  sd_value <- sc_score_chen$sd[i]
  sc_mini_cibersort_z <- cbind(sc_mini_cibersort_z, (Score_prediction[,i] - mean_value)/sd_value)
}

colnames(sc_mini_cibersort_z) <- paste0("Cluster",1:5)
rownames(sc_mini_cibersort_z) <- rownames(Score_prediction)

sc_mini_cibersort_z <- as.data.frame(sc_mini_cibersort_z)
sc_mini_cibersort_z$predicted.label <- max.col(sc_mini_cibersort_z)
sc_mini_cibersort_z$Patient <- rownames(sc_mini_cibersort_z)

RNA_DeconResult = as.data.frame(t(RNA_DeconResult))
RNA_DeconResult_p = RNA_DeconResult
RNA_DeconResult_p$Patient <- rownames(RNA_DeconResult_p)
sc_mini_cibersort_z <- left_join(sc_mini_cibersort_z, RNA_DeconResult_p, by = "Patient")

metadata <- openxlsx::read.xlsx("../2 - Analysis_V3_checked/GBC_external RNA-seq/NC_gbc_survival_GBCKorea 20210618.xlsx")
colnames(metadata)[2] <- "Patient"
metadata <- left_join(metadata, sc_mini_cibersort_z, by = "Patient")

metadata$vital_status = ifelse(metadata$vital_status=="dead", 1, 0)
table(metadata$predicted.label)

metadata$group = ifelse(metadata$predicted.label == 2,2,"others")
fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=metadata)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = metadata,             # data used to fit survival curves.
  palette = c("#34ebc3","#bbbcbd"),
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

sc_mini_cibersort_z$group = ifelse(sc_mini_cibersort_z$predicted.label == 2,2,"others")

ggboxplot(sc_mini_cibersort_z, x = "group", y = "CD4T_C2_CCR7", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD4T_C2_CCR7")
ggboxplot(sc_mini_cibersort_z, x = "group", y = "B_C0_IGHA1", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("B_C0_IGHA1")
ggboxplot(sc_mini_cibersort_z, x = "group", y = "γdT_C12_TRDC", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("γdT_C12_TRDC")
ggboxplot(sc_mini_cibersort_z, x = "group", y = "CD8T_C0_CCR7_GZMK", color = "group", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD8T_C0_CCR7_GZMK")

# ## 2.8 External Validation by GenePattern ####
# RNA_DeconResult = readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
# RNA_DeconResult = RNA_DeconResult[,1:119]
# RNA_DeconResult = as.data.frame(t(RNA_DeconResult))
# RNA_DeconResult = apply(RNA_DeconResult,2,function(x){x/sum(x)})
# 
# expr_matrix <- as.matrix(RNA_DeconResult)
# 
# # 确保行名是基因，列名是细胞
# genes <- rownames(expr_matrix)
# cells <- colnames(expr_matrix)
# 
# # 创建 GCT 文件头信息
# num_genes <- nrow(expr_matrix)
# num_cells <- ncol(expr_matrix)
# 
# # 创建 GCT 文件的内容
# gct_content <- c(
#   "#1.3", 
#   paste(num_genes, num_cells, sep = "\t"), 
#   paste("probeid", "gene", paste(cells, collapse = "\t"), sep = "\t")
# )
# 
# # 添加每个基因的表达值
# for (i in 1:num_genes) {
#   gene_name <- genes[i]
#   gene_expr <- paste(expr_matrix[i, ], collapse = "\t")
#   gct_content <- c(gct_content, paste(gene_name, gene_name, gene_expr, sep = "\t"))
# }
# 
# # 将内容写入 GCT 文件
# writeLines(gct_content, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/GenePattern/External_bulk-RNA-seq/bulk_cibersort_cc.gct")
# 
# ##
# MI5_Feature <- c("CD4T_C1_FOXP3","CD8T_C1_CXCL13")
# MI2_Feature <- c("γdT_C12_TRDC","CD4T_C2_CCR7","CD8T_C0_CCR7_GZMK","B_C0_IGHA1")
# MI1_Feature <- c("DC_C0_cDC2_IL1B","DC_C5_cDC1","F_C1_CFD","Per_C1_MYH11","DC_C4_cDC2_PPP1R14A","EC_C3_GJA5","Per_C3_STEAP4","EC_C0_ACKR1")
# MI4_Feature <- c("M_C2_SPP1","M_C7_AGR2","EC_C2_CXCR4")
# MI3_Feature <- c("N_C9_Nc","N_C0_N0","N_C1_N2","N_C2_N1","N_C7_N0")
# 
# MI_Sig_CC <- data.frame(probeid = c(MI1_Feature, MI2_Feature, MI3_Feature, MI4_Feature, MI5_Feature),
#                         GeneName = c(MI1_Feature, MI2_Feature, MI3_Feature, MI4_Feature, MI5_Feature),
#                         Class = c(rep("1", length(MI1_Feature)),
#                                   rep("2", length(MI2_Feature)),
#                                   rep("3", length(MI3_Feature)),
#                                   rep("4", length(MI4_Feature)),
#                                   rep("5", length(MI5_Feature))))
# write.table(MI_Sig_CC, "../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/GenePattern/External_bulk-RNA-seq/MI_Sig_CC.txt",row.names = F, col.names = T,quote = F,sep = "\t")
# 
# ##
# NTP_prediction <- read.csv("../Publication/240923_Nature Genetics_Revise/2 - Analysis/MI_Ident_CC/GenePattern/External_bulk-RNA-seq/output/NTP_prediction_result.xls.csv")
# colnames(NTP_prediction)[1] <- "Patient"
# 
# metadata <- openxlsx::read.xlsx("../2 - Analysis_V3_checked/GBC_external RNA-seq/NC_gbc_survival_GBCKorea 20210618.xlsx")
# colnames(metadata)[2] <- "Patient"
# metadata <- left_join(metadata, NTP_prediction, by = "Patient")
# 
# metadata$vital_status = ifelse(metadata$vital_status=="dead", 1, 0)
# table(metadata$predict.label)
# 
# metadata$group = ifelse(metadata$predict.label == 2,2,"others")
# fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=metadata)
# ggsurvplot(
#   fit,                     # survfit object with calculated statistics.
#   data = metadata,             # data used to fit survival curves.
#   palette = c("#34ebc3","#bbbcbd"),
#   risk.table = TRUE,       # show risk table.
#   pval = TRUE,             # show p-value of log-rank test.
#   #conf.int = TRUE,         # show confidence intervals for
#   # palette = "npg",
#   xlab = "Time in days",   # customize X axis label.
#   ggtheme = theme_classic(),
#   #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
#   #surv.median.line = "hv",  # add the median survival pointer.
#   # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
# )

## 2.9 Fu scRNA-seq Validation by Score ####
Endothelium <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/Endothelium/Endothelium.RDS")
Myeloid_predict <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/Myeloid/Myeloid_predict.RDS")

CD4T <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer CD4T.rds")
CD4T$predict.label <- CD4T$predicted.id
CD8T <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer CD8T.rds")
CD8T$predict.label <- CD8T$predicted.id
Fibro <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer fibro.rds")
Fibro$predict.label <- Fibro$predicted.id
NK <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer NK.rds")
NK$predict.label <- NK$predicted.id
B <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/2 - GBC_Trasnfer/LS/Label Transfer B.rds")
B$predict.label <- B$predicted.id

All <- merge(x = Myeloid_predict, y = list(CD4T, CD8T, NK, B, Endothelium, Fibro))
patient_select <- unique(All$orig.ident)
patient_select = patient_select[grep("T$", patient_select)]

All <- subset(All, subset = orig.ident %in% patient_select)

GBC_AllSubtype <- All@meta.data
cycle = unique(GBC_AllSubtype$orig.ident)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(GBC_AllSubtype, orig.ident == cycle[i])
  Freq = bind_rows(Freq, table(sub$predict.label)/nrow(sub))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

Freq_PctScore = CreateSeuratObject(t(Freq))

MI5_Feature <- c("CD4T-C1-FOXP3","CD8T-C1-CXCL13")
MI2_Feature <- c("γdT-C12-TRDC","CD4T-C2-CCR7","CD8T-C0-CCR7-GZMK","B-C0-IGHA1")
MI1_Feature <- c("DC-C0-cDC2-IL1B","DC-C5-cDC1","F-C1-CFD","Per-C1-MYH11","DC-C4-cDC2-PPP1R14A","EC-C3-GJA5","Per-C3-STEAP4","EC-C0-ACKR1")
MI4_Feature <- c("M-C2-SPP1","M-C7-AGR2","EC-C2-CXCR4")
MI3_Feature <- c("N-C9-Nc","N-C0-N0","N-C1-N2","N-C2-N1","N-C7-N0")

gene_set_list <- list(
  gene_set1 = MI1_Feature,
  gene_set2 = MI2_Feature,
  gene_set3 = MI3_Feature,
  gene_set4 = MI4_Feature,
  gene_set5 = MI5_Feature
)

Freq_PctScore = AddModuleScore(Freq_PctScore, features = gene_set_list, ctrl = 3)

Score_prediction <- Freq_PctScore@meta.data[,paste0("Cluster",1:5)]
Score_prediction$predicted.label <- max.col(Score_prediction)
Score_prediction$Patient <- rownames(Score_prediction)

selected_subtypes = c("CD4T_C2_CCR7","B_C0_IGHA1","γdT_C12_TRDC","CD8T_C0_CCR7_GZMK")
Freq <- Freq[,selected_subtypes]
Freq$Patient <- rownames(Freq)

Freq <- left_join(Freq, Score_prediction, by = "Patient")
Freq$group_plot <- ifelse(Freq$predicted.label == 2, "MI2", "others")

ggboxplot(Freq, x = "group_plot", y = "CD4T_C2_CCR7", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD4T_C2_CCR7")
ggboxplot(Freq, x = "group_plot", y = "B_C0_IGHA1", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("B_C0_IGHA1")
ggboxplot(Freq, x = "group_plot", y = "γdT_C12_TRDC", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("γdT_C12_TRDC")
ggboxplot(Freq, x = "group_plot", y = "CD8T_C0_CCR7_GZMK", color = "group_plot", 
          add = "jitter") + 
  stat_compare_means(method = 't.test') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  ggtitle("CD8T_C0_CCR7_GZMK")
















