####PBMC 99####
pbmc=Read10X('...\\GBC_099_PBMC\\outs\\filtered_feature_bc_matrix')
pbmc=CreateSeuratObject(pbmc, min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc<- ScaleData(object = pbmc,verbose = F)
PCS = 30
k.param = 20
pbmc= RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = FALSE)
pbmc=FindNeighbors(pbmc,dims = 1:PCS,verbose = FALSE,k.param = k.param)
pbmc=RunUMAP(pbmc,dims = 1:PCS,verbose = FALSE) 
pbmc=FindClusters(pbmc,resolution = 0.5,verbose = FALSE)
saveRDS(pbmc,'...\\01 PBMC 99 All.rds')
pbmc99 = readRDS('...\\01 PBMC 99 All.rds')


FeaturePlot(pbmc99,features = c('AREG','CCR7','CD4','CD3D')) +DimPlot(pbmc99,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')

maker=FindAllMarkers(pbmc,only.pos = TRUE)
saveRDS(maker,'...\\02 Marker PBMC 99.rds')

cd4.pbmc99 = subset(pbmc,seurat_clusters %in% c(3,6))
cd4.pbmc99$group = 'PBMC_99'
saveRDS(cd4.pbmc99,'...\\03 pbmc099 naive.rds')


####PBMC 98 ####
pbmc=Read10X('...\\GBC_098_PBMC\\filtered_feature_bc_matrix')
pbmc=CreateSeuratObject(pbmc, min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc<- ScaleData(object = pbmc,verbose = F)
PCS = 30
k.param = 20
pbmc= RunPCA(pbmc,features = VariableFeatures(pbmc),verbose = FALSE)
pbmc=FindNeighbors(pbmc,dims = 1:PCS,verbose = FALSE,k.param = k.param)
pbmc=RunUMAP(pbmc,dims = 1:PCS,verbose = FALSE) 
pbmc=FindClusters(pbmc,resolution = 0.5,verbose = FALSE)
saveRDS(pbmc,'...\\01 PBMC 98 All.rds')

pbmc98 = readRDS('...\\01 PBMC 98 All.rds')
DimPlot(pbmc98,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')

maker=FindAllMarkers(pbmc98,only.pos = TRUE)
saveRDS(maker,'...\\02 Marker PBMC 98.rds')
marker = readRDS('...\\02 Marker PBMC 98.rds')
table(pbmc$seurat_clusters)
dim(cd4)
FeaturePlot(pbmc98,features = c('CD3D','CCR7','CD4','SELL','AREG')) +DimPlot(pbmc98,label = TRUE) 
cd4.pbmc98=subset(pbmc98,seurat_clusters %in% c(0,8))
cd4.pbmc98$group = 'PBMC_98'
saveRDS(cd4.pbmc98,'...\\03 pbmc098 naive.rds')

pbmc.total = merge(cd4.pbmc99,cd4.pbmc98,add.cell.ids =c('PBMC_99','PBMC_98'))
pbmc.total <- NormalizeData(pbmc.total, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.total <- FindVariableFeatures(pbmc.total, selection.method = "vst", nfeatures = 2000)
pbmc.total<- ScaleData(object = pbmc.total,verbose = F)
PCS = 30
k.param = 20
pbmc.total= RunPCA(pbmc.total,features = VariableFeatures(pbmc.total),verbose = FALSE)
pbmc.total=FindNeighbors(pbmc.total,dims = 1:PCS,verbose = FALSE,k.param = k.param)
pbmc.total=RunUMAP(pbmc.total,dims = 1:PCS,verbose = FALSE) 
pbmc.total=FindClusters(pbmc.total,resolution = 0.5,verbose = FALSE)
DimPlot(pbmc.total,group.by = 'group')
FeaturePlot(pbmc.total,features = c('AREG','CCR7','CD4','SELL','LEF1')) 
pbmc.total$group = 'PBMC'

------------------------------------------------------------------------------------


####PBMC data from normal donor

setwd('...\\GSE126030_RAW(1)')
cd4= readRDS('...\\cd4_newname.rds')

data1 = read.table('PP017swap.filtered.matrix.txt',header = TRUE)
data2 = read.table('PP018swap.filtered.matrix.txt',header = TRUE)
data3 = read.table('PP019swap.filtered.matrix.txt',header = TRUE)
data4 = read.table('PP020swap.filtered.matrix.txt',header = TRUE)

gene = unique(data1$Gene)
gene = intersect(gene,rownames(cd4))

data1 = data1[data1$Gene %in% gene,]
data1<-data1[!duplicated(data1$Gene, fromLast=TRUE), ] 
rownames(data1)=data1$Gene
data1 = data1[,-c(1,2)]

data2 = data2[data2$Gene %in% gene,]
data2<-data2[!duplicated(data2$Gene, fromLast=TRUE), ] 
rownames(data2)=data2$Gene
data2 = data2[,-c(1,2)]

data3 = data3[data3$Gene %in% gene,]
data3<-data3[!duplicated(data3$Gene, fromLast=TRUE), ] 
rownames(data3)=data3$Gene
data3 = data3[,-c(1,2)]

data4= data4[data4$Gene %in% gene,]
data4<-data4[!duplicated(data4$Gene, fromLast=TRUE), ] 
rownames(data4)=data4$Gene
data4 = data4[,-c(1,2)]

obj.1 = CreateSeuratObject(counts = data1, project = 'Blood donor A1')
obj.2 = CreateSeuratObject(counts = data2, project = 'Blood donor A2')
obj.3 = CreateSeuratObject(counts = data3, project = 'Blood donor B1')
obj.4 = CreateSeuratObject(counts = data4, project = 'Blood donor B2')

obj.merge = merge(obj.1,obj.2)

obj.1 = CreateSeuratObject(counts = data1, project = 'Blood donor A1')
obj.2 = CreateSeuratObject(counts = data2, project = 'Blood donor A2')
obj.merge = merge(obj.1, obj.2)

obj.3 = CreateSeuratObject(counts = data3, project = 'Blood donor B1')
obj.merge = merge(x = obj.merge, y = obj.3)

obj.4 = CreateSeuratObject(counts = data4, project = 'Blood donor B2')
obj.merge = merge(x = obj.merge, y = obj.4)
dim(obj.merge)

obj.merge$group = 'Normal_PBMC'
obj.merge <- NormalizeData(obj.merge)
obj.merge <- FindVariableFeatures(obj.merge, selection.method = "vst", nfeatures = 2000)
obj.merge <- ScaleData(object = obj.merge)

PCS = 30
k.param = 20
obj.merge = RunPCA(obj.merge,features = VariableFeatures(obj.merge),verbose = FALSE)
obj.merge=FindNeighbors(obj.merge,dims = 1:PCS,verbose = FALSE,k.param = k.param)
obj.merge=RunUMAP(obj.merge,dims = 1:PCS,verbose = FALSE) 
obj.merge=FindClusters(obj.merge,resolution = 0.5,verbose = FALSE)
DimPlot(obj.merge,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
marker = FindAllMarkers(obj.merge,only.pos = TRUE)

saveRDS(obj.merge,'...\\Normal PBMC All.rds')

cd4.sub = subset(obj.merge,seurat_clusters %in% c(0,3,5,6,7,8))
cd4t.naive = subset(cd4.sub,seurat_clusters %in% c(5,6,7))

gbc.pbmc.98 = readRDS('...\\03 pbmc098 naive.rds')
gbc.pbmc.99 = readRDS('...\\03 pbmc099 naive.rds')
gbc.pbmc = merge(gbc.pbmc.98,gbc.pbmc.99)
gbc.pbmc$group = 'GBC_PBMC'

pbmc.total = merge(gbc.pbmc,cd4t.naive)
pbmc.total <- NormalizeData(pbmc.total)
pbmc.total <- FindVariableFeatures(pbmc.total, selection.method = "vst", nfeatures = 2000)
pbmc.total <- ScaleData(object = pbmc.total)

PCS = 30
k.param = 20
pbmc.total = RunPCA(pbmc.total,features = VariableFeatures(pbmc.total),verbose = FALSE)
pbmc.total=FindNeighbors(pbmc.total,dims = 1:PCS,verbose = FALSE,k.param = k.param)
pbmc.total=RunUMAP(pbmc.total,dims = 1:PCS,verbose = FALSE) 
pbmc.total=FindClusters(pbmc.total,resolution = 0.5,verbose = FALSE)
DimPlot(pbmc.total,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
DimPlot(pbmc.total,reduction = 'umap',label = TRUE,group.by = 'group')
saveRDS(pbmc.total,'...\\02 Normal PBMC add GBC PBMC cd4t naive.rds')

Idents(pbmc.total) = pbmc.total$group
marker = FindAllMarkers(pbmc.total,only.pos = TRUE)

areg.expr = FetchData(pbmc.total,vars = c('AREG','group','CCR7','LEF1','SELL','IL7R','TCF7'),slot = 'data')
ggplot(areg.expr,aes(group,AREG))+geom_violin()

melt = melt(areg.expr)
melt.1 = subset(melt,group =='GBC_PBMC')
ggplot(melt.1, aes(x = variable, y = value,fill=variable)) +
 geom_jitter(width=0.3,alpha=0.1,size=0.1)+
  geom_violin()  +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  stat_compare_means(aes(label="p.format"),
                     method = "anova",label = "p.format")+
  ggtitle('GBC_PBMC')+
  scale_fill_manual(values = c( '#B4505C','#D6A1A9','#A6C5B5','#FFE8CF','#A3937C','#4D3D30')) +
  theme_pubr()

melt.2 = subset(melt,group =='Normal_PBMC')
p=ggplot(melt.2, aes(x = variable, y = value,fill=variable)) +
  geom_jitter(width=0.3,alpha=0.2,size=0.1)+
  geom_violin(position = position_dodge(width = 0.6))  +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  stat_compare_means(aes(label="p.format"),
                     method = "anova",label = "p.format")+
  ggtitle('Normal_PBMC')+
  scale_fill_manual(values = c( '#B4505C','#D6A1A9','#A6C5B5','#FFE8CF','#A3937C','#4D3D30')) +
  theme_pubr()

data= readRDS('...\\02 Normal PBMC add GBC PBMC cd4t naive.rds')
df = FetchData(data,vars = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R','group'))
df = melt(df)
df2 = df[sample(c(1:dim(df)[1]),2000),]

df.mean = df2 %>% group_by(variable,group) %>% summarise_each(funs = mean)
a=ggplot(df2, aes(x = variable, y = value, fill = variable)) +
  geom_jitter(width=0.3,alpha=0.3,size=0.1)+
  geom_violin(width = 1, position = position_dodge(width = 3)) +
  geom_point(data = df.mean, aes(x = variable, y = value), size = 2, shape = 21, fill = "white") +
  labs(title = "Violin Plot", x = "Groups", y = "Value") +
  ggtitle('Naive CD4 T') +
  scale_fill_manual(values = c('#B4505C', '#D6A1A9', '#A6C5B5', '#FFE8CF', '#A3937C', '#4D3D30')) +
  theme_pubr() +
  facet_wrap(~group, ncol = 1)













