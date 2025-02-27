setwd("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked")
set.seed(925)

library(Seurat)
library(dplyr)
library(patchwork)

# 1 - Mapping and annotating query datasets ####
## 1. Myeloid ####
### Query ####
Myeloid <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/Split_by_Celltype/07 Myeloid.rds")
Myeloid <- NormalizeData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid)

### Ref ####
Myeloid_cell <- readRDS("./1 - GBC_input data/Myeloid_cell.RDS") # 323569 cells
Cell_type <- read.csv(file = "./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
Myeloid_cell = subset(Myeloid_cell, cells = Cell_type[Cell_type$celltype %in% c("DC","Mast","Mono_Macro","Neu"),]$cellid) # 307285 cells

Cell_type <- Cell_type[Cell_type$cellid %in% colnames(Myeloid_cell),]

Myeloid_cell@meta.data$cellid <- rownames(Myeloid_cell@meta.data)
temp <- Myeloid_cell@meta.data
temp <- left_join(temp, Cell_type, by = "cellid")

identical(temp$cellid, rownames(Myeloid_cell@meta.data))
identical(temp$orig.ident.y, Myeloid_cell@meta.data$orig.ident)

Myeloid_cell@meta.data$Celltype <- temp$celltype.y
Myeloid_cell@meta.data$Subtype <- temp$subtype

Idents(Myeloid_cell) = "Subtype"

Myeloid_cell <- NormalizeData(Myeloid_cell)
Myeloid_cell <- FindVariableFeatures(Myeloid_cell)

### anchor ####
anchors <- FindTransferAnchors(reference = Myeloid_cell, query = Myeloid, dims = 1:30)

### Transfer ####
predicted.labels <- TransferData(anchorset = anchors, refdata = Myeloid_cell$Subtype, dims = 1:30)
predicted.labels$cellid <- rownames(predicted.labels)

temp <- Myeloid@meta.data
temp$cellid <- rownames(temp)
temp <- left_join(temp, predicted.labels, by = "cellid")

identical(temp$cellid, rownames(Myeloid@meta.data))
Myeloid$predict.label <- temp$predicted.id
Idents(Myeloid) <- "predict.label"

## 2. Endothelium ####
### Query ####
Endothelium <- readRDS("../Publication/240923_Nature Genetics_Revise/Paired_GBC_NAT_FJ/Split_by_Celltype/02 Endothelial.rds")
Endothelium <- NormalizeData(Endothelium)
Endothelium <- FindVariableFeatures(Endothelium)

### Ref ####
Endothelial_Ref <- readRDS("./1 - GBC_input data/Endothelial.RDS") # 22501 cells
Cell_type <- read.csv(file = "./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
Endothelial_Ref = subset(Endothelial_Ref, cells = Cell_type[Cell_type$celltype %in% c("Endothelium"),]$cellid) # 16886 cells

Cell_type <- Cell_type[Cell_type$cellid %in% colnames(Endothelial_Ref),]

Endothelial_Ref@meta.data$cellid <- rownames(Endothelial_Ref@meta.data)
temp <- Endothelial_Ref@meta.data
temp <- left_join(temp, Cell_type, by = "cellid")

identical(temp$cellid, rownames(Endothelial_Ref@meta.data))
identical(temp$orig.ident.y, Endothelial_Ref@meta.data$orig.ident)

Endothelial_Ref@meta.data$Celltype <- temp$celltype.y
Endothelial_Ref@meta.data$Subtype <- temp$subtype

Idents(Endothelial_Ref) = "Subtype"

Endothelial_Ref <- NormalizeData(Endothelial_Ref)
Endothelial_Ref <- FindVariableFeatures(Endothelial_Ref)

### anchor ####
anchors <- FindTransferAnchors(reference = Endothelial_Ref, query = Endothelium, dims = 1:30)

### Transfer ####
predicted.labels <- TransferData(anchorset = anchors, refdata = Endothelial_Ref$Subtype, dims = 1:30)
predicted.labels$cellid <- rownames(predicted.labels)

temp <- Endothelium@meta.data
temp$cellid <- rownames(temp)
temp <- left_join(temp, predicted.labels, by = "cellid")

identical(temp$cellid, rownames(Endothelium@meta.data))
Endothelium$predict.label <- temp$predicted.id
Idents(Endothelium) <- "predict.label"

# 2 - Calculate CC: N vs T ####
## 1. Myeloid ####
Myeloid_predict$group = substr(Myeloid_predict$orig.ident,nchar(Myeloid_predict$orig.ident),nchar(Myeloid_predict$orig.ident))
table(Myeloid_predict$group)

df.cell = data.frame(patient = Myeloid_predict$orig.ident, celltype = Myeloid_predict$predict.label)
df.cell = dcast(df.cell,patient~celltype)
rownames(df.cell) = df.cell$patient
df.cell = df.cell[,-1]

norm_func = function(x){
  (x/sum(x))*100
}

df.cell.norm = apply(df.cell, 1, norm_func)
df.cell.norm = t(df.cell.norm) %>%as.data.frame()
df.cell.norm$group = substr(rownames(df.cell.norm),nchar(rownames(df.cell.norm)),nchar(rownames(df.cell.norm)))


plot_function <- function(column_name) {
  ggboxplot(df.cell.norm, x = "group", y = column_name, color = "group", 
            palette = "jco", add = "jitter") + 
    stat_compare_means(method = 'wilcox.test', label = "p.format") +
    theme(legend.position = "top") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    ggtitle(column_name)
}

# 获取你感兴趣的列名列表
columns_to_plot <- colnames(df.cell.norm)[1:29]

# 打开 PDF 文件
pdf(file = '~/Downloads/NT compared Mye.pdf',
    height = 3, width = 2.5)
lapply(columns_to_plot, function(column) {
  plot <- plot_function(column)
  print(plot)
})
dev.off()
### 改输出路径 ####

## 2. Endothelium ####
Endothelium$group = substr(Endothelium$orig.ident,nchar(Endothelium$orig.ident),nchar(Endothelium$orig.ident))
table(Endothelium$group)

df.cell = data.frame(patient = Endothelium$orig.ident, celltype = Endothelium$predict.label)
df.cell = dcast(df.cell,patient~celltype)
rownames(df.cell) = df.cell$patient
df.cell = df.cell[,-1]

norm_func = function(x){
  (x/sum(x))*100
}

df.cell.norm = apply(df.cell, 1, norm_func)
df.cell.norm = t(df.cell.norm) %>%as.data.frame()
df.cell.norm$group = substr(rownames(df.cell.norm),nchar(rownames(df.cell.norm)),nchar(rownames(df.cell.norm)))


plot_function <- function(column_name) {
  ggboxplot(df.cell.norm, x = "group", y = column_name, color = "group", 
            palette = "jco", add = "jitter") + 
    stat_compare_means(method = 'wilcox.test', label = "p.format") +
    theme(legend.position = "top") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    ggtitle(column_name)
}

# 获取你感兴趣的列名列表
columns_to_plot <- colnames(df.cell.norm)[1:9]

# 打开 PDF 文件
pdf(file = '~/Downloads/NT compared Endo.pdf',
    height = 3, width = 2.5)
lapply(columns_to_plot, function(column) {
  plot <- plot_function(column)
  print(plot)
})
dev.off()

### 改输出路径 ####
