

library(Seurat)

t.data = readRDS('...\\05 Tcell.rds')


t.data <- NormalizeData(t.data, normalization.method = "LogNormalize", scale.factor = 10000)
t.data <- FindVariableFeatures(t.data, selection.method = "vst", nfeatures = 2000)

t.data <- ScaleData(t.data, features = VariableFeatures(t.data))
t.data <- RunPCA(t.data, features = VariableFeatures(object = t.data))

t.data <- FindNeighbors(t.data, dims = 1:20)
t.data <- FindClusters(t.data, resolution = 0.8)
t.data <- RunUMAP(t.data, dims = 1:20)
DimPlot(t.data, reduction = "umap", label = TRUE)
FeaturePlot(t.data,features = c('CD4','CD8A','NKG7'))+DimPlot(t.data, reduction = "umap", label = TRUE)

cd4 = subset(t.data,seurat_clusters %in%c(1,4,5,6,8,13,15))
cd8 = subset(t.data,seurat_clusters %in%c(0,3,7,11,16,17,20,21,22))
nk = subset(t.data, seurat_clusters %in%c(2,9))

unassign.cluster = subset(unassign.cluster,seurat_clusters %in%c(12,14,18))
unassign.cluster <- NormalizeData(unassign.cluster, normalization.method = "LogNormalize", scale.factor = 10000)
unassign.cluster <- FindVariableFeatures(unassign.cluster, selection.method = "vst", nfeatures = 2000)

unassign.cluster <- ScaleData(unassign.cluster, features = VariableFeatures(unassign.cluster))
unassign.cluster <- RunPCA(unassign.cluster, features = VariableFeatures(object = unassign.cluster))

unassign.cluster <- FindNeighbors(unassign.cluster, dims = 1:20)
unassign.cluster <- FindClusters(unassign.cluster, resolution = 0.1)
unassign.cluster <- RunUMAP(unassign.cluster, dims = 1:20)
DimPlot(unassign.cluster, reduction = "umap", label = TRUE)
FeaturePlot(unassign.cluster,features = c('CD4','CD8A','NKG7'))+DimPlot(unassign.cluster, reduction = "umap", label = TRUE)
FeaturePlot(unassign.cluster,features = 'MKI67')

cd8.u = subset(unassign.cluster,seurat_clusters %in%c(1,2,4))
cd4.u = subset(unassign.cluster,seurat_clusters %in%c(0))


cd8 = merge(cd8,cd8.u)
cd4 = merge(cd4,cd4.u)


cd4 <- NormalizeData(cd4)
cd4 <- FindVariableFeatures(cd4)
#### Reference  ####

cd4.ref = readRDS('...01 CD4.rds')
anchors <- FindTransferAnchors(reference = cd4.ref, query = cd4, dims = 1:30)
predicted.labels <- TransferData(anchorset = anchors, refdata = cd4.ref$celltype, dims = 1:30)
table(predicted.labels$predicted.id)
cd4 <- AddMetaData(cd4, metadata = predicted.labels)
saveRDS(cd4,'...\\Label Transfer CD4T.rds')
saveRDS(anchors,'...\\Anchors CD4T.rds')


cd8.ref = readRDS('...\\02 CD8.rds')
anchors <- FindTransferAnchors(reference = cd8.ref, query = cd8, dims = 1:30)
predicted.labels <- TransferData(anchorset = anchors, refdata = cd8.ref$celltype, dims = 1:30)
table(predicted.labels$predicted.id)
cd8 <- AddMetaData(cd8, metadata = predicted.labels)
saveRDS(cd8,'...\\Label Transfer CD8T.rds')
saveRDS(anchors,'...\\Anchors CD8T.rds')


cd8 = readRDS('...\\Label Transfer CD8T.rds')
table(cd8$predicted.id)


nk.ref = readRDS('...\\03 NK.rds')
anchors <- FindTransferAnchors(reference = nk.ref, query = nk, dims = 1:30)
predicted.labels <- TransferData(anchorset = anchors, refdata = nk.ref$celltype, dims = 1:30)
table(predicted.labels$predicted.id)
nk <- AddMetaData(nk, metadata = predicted.labels)
saveRDS(nk,'...\\Label Transfer NK.rds')
saveRDS(anchors,'...\\Anchors NK.rds')

#### Normal Tumor ####
library(Seurat)
library(dplyr)
library(reshape2)
cd8 = readRDS('...\\Label Transfer NK.rds')

cd8$group = substr(cd8$orig.ident,nchar(cd8$orig.ident),nchar(cd8$orig.ident))
table(cd8$group)

df.cell = data.frame(patient = cd8$orig.ident, celltype = cd8$predicted.id)
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

columns_to_plot <- colnames(df.cell.norm)[1:11]

pdf(file = '...\\NT Compared\\NT compared B.pdf',
    height = 3, width = 2.5)
lapply(columns_to_plot, function(column) {
  plot <- plot_function(column)
  print(plot)
})
dev.off()

cd4= readRDS('...\\Label Transfer CD4T.rds')

cd4$group = substr(cd4$orig.ident,nchar(cd4$orig.ident),nchar(cd4$orig.ident))
df.cell = data.frame(patient = cd4$orig.ident, celltype = cd4$predicted.id)
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
columns_to_plot <- colnames(df.cell.norm)[1:10]

pdf(file = '...\\NT cd4 compared.pdf',
    height = 3, width = 2.5)
lapply(columns_to_plot, function(column) {
  plot <- plot_function(column)
  print(plot)
})
dev.off()


cd4 <- NormalizeData(cd4, normalization.method = "LogNormalize", scale.factor = 10000)
cd4 <- FindVariableFeatures(cd4, selection.method = "vst", nfeatures = 2000)

cd4 <- ScaleData(cd4, features = VariableFeatures(cd4))
cd4 <- RunPCA(cd4, features = VariableFeatures(object = cd4))

cd4 <- FindNeighbors(cd4, dims = 1:20)
cd4 <- FindClusters(cd4, resolution = 0.8)
cd4 <- RunUMAP(cd4, dims = 1:20)
DimPlot(cd4, reduction = "umap", label = TRUE)

DimPlot(cd4, reduction = "umap", group.by ='predicted.id')
cd4$predicted.id

cd4.sub = subset(cd4,predicted.id == 'CD4T_C2_CCR7')
DimPlot(cd4.sub, reduction = "umap", group.by ='predicted.id')

FeaturePlot(cd4.sub,features = c('AREG','CCR7','LEF1','SELL'))
Idents(cd4.sub) = cd4.sub$group
DimPlot(cd4.sub, reduction = "umap", split.by  ='group')
cd4.sub$group = substr(cd4.sub$orig.ident,nchar(cd4.sub$orig.ident),nchar(cd4.sub$orig.ident))


cd4.fuji=readRDS('...\\Fig6.Tcell.rds')
cd4.fuji$celltype=cd4.fuji@active.ident
DimPlot(cd4.fuji)
cd4.naive= subset(cd4.fuji,celltype=='CD4 TN')


####
b = readRDS('...\\06 Bcell.rds')
b.ref = readRDS('...\\04 B cells.rds')
anchors <- FindTransferAnchors(reference = b.ref, query = b, dims = 1:30)
predicted.labels <- TransferData(anchorset = anchors, refdata = b.ref$celltype, dims = 1:30)
table(predicted.labels$predicted.id)
b <- AddMetaData(b, metadata = predicted.labels)
saveRDS(b,'...\\Label Transfer B.rds')
saveRDS(anchors,'...\\Anchors B.rds')




