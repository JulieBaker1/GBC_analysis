#### Part 1 ICC Primary naive T ####

data1 = Read10X('...\\230093ICC1\\230093 primary focus')
data2 = Read10X('...\\357818ICC2\\357818 primary focus')

dim(data1)
dim(data2)

data1 = CreateSeuratObject(data1,project = 'ICC1_primary')
data2 = CreateSeuratObject(data2,project = 'ICC2_primary')
data = merge(data1,data2)

data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(object = data)

PCS = 30
k.param = 20
data = RunPCA(data,features = VariableFeatures(data),verbose = FALSE)
data=FindNeighbors(data,dims = 1:PCS,verbose = FALSE,k.param = k.param)
data=RunUMAP(data,dims = 1:PCS,verbose = FALSE) 
data=FindClusters(data,resolution = 0.8,verbose = FALSE)

DimPlot(data,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
saveRDS(data,'01 ICC Primary focus Total.rds')

FeaturePlot(data,features = c('CD3D','CD4','CD8A','IL7R','CCR7','TCF7','LEF1','AREG'))+
  DimPlot(data,label = TRUE)

cd4.sub = subset(data,seurat_clusters %in% c(4,12,6,9,5,10))
cd4.sub <- NormalizeData(cd4.sub)
cd4.sub <- FindVariableFeatures(cd4.sub, selection.method = "vst", nfeatures = 2000)
cd4.sub <- ScaleData(object = cd4.sub)
PCS = 30
k.param = 20
cd4.sub = RunPCA(cd4.sub,features = VariableFeatures(cd4.sub),verbose = FALSE)
cd4.sub=FindNeighbors(cd4.sub,dims = 1:PCS,verbose = FALSE,k.param = k.param)
cd4.sub=RunUMAP(cd4.sub,dims = 1:PCS,verbose = FALSE) 
cd4.sub=FindClusters(cd4.sub,resolution = 0.8,verbose = FALSE)
DimPlot(cd4.sub,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
FeaturePlot(cd4.sub,features = c('CD3D','CD4','CD8A','IL7R','CCR7','TCF7','LEF1','AREG'))+
  DimPlot(cd4.sub,label = TRUE)
cd4.naive = subset(cd4.sub,seurat_clusters ==3)
saveRDS(cd4.naive,'02 CD4T Naive ICC.rds')

VlnPlot(object = cd4.naive, 
        features = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R'))

df = FetchData(cd4.naive,vars = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R'))
df = melt(df)

ggplot(df, aes(x = variable, y = value,fill=variable)) +
  geom_violin() +
  #geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  ggtitle('ICC_P')+
  #scale_fill_manual(values = c( '#6A91B7','#D18796')) +
  theme_pubr()

ggplot(df, aes(x = variable, y = value,fill=variable)) +
  geom_jitter(width=0.3,alpha=0.3,size=0.1)+
  geom_violin()  +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  #stat_compare_means(aes(label="p.format"),
                     #method = "anova",label = "p.format")+
  ggtitle('ICC_P')+
  scale_fill_manual(values = c( '#B4505C','#D6A1A9','#A6C5B5','#FFE8CF','#A3937C','#4D3D30')) +
  theme_pubr()


#### Part 2 ICC Primary Naive T ####

data1 = Read10X('...\\230093ICC1\\230093 PBMC')
data2 = Read10X('...\\357818ICC2\\357818 PBMC')

dim(data1)
dim(data2)

data1 = CreateSeuratObject(data1,project = 'ICC1_PBMC')
data2 = CreateSeuratObject(data2,project = 'ICC2_PBMC')
data = merge(data1,data2)

data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(object = data)

PCS = 30
k.param = 20
data = RunPCA(data,features = VariableFeatures(data),verbose = FALSE)
data=FindNeighbors(data,dims = 1:PCS,verbose = FALSE,k.param = k.param)
data=RunUMAP(data,dims = 1:PCS,verbose = FALSE) 
data=FindClusters(data,resolution = 0.8,verbose = FALSE)
saveRDS(data,'03 ICC Pbmc TOTAL.rds')

cd4.sub = subset(data,seurat_clusters %in% c(0,1,7,10,12,13))
cd4.sub <- NormalizeData(cd4.sub)
cd4.sub <- FindVariableFeatures(cd4.sub, selection.method = "vst", nfeatures = 2000)
cd4.sub <- ScaleData(object = cd4.sub)
PCS = 30
k.param = 20
cd4.sub = RunPCA(cd4.sub,features = VariableFeatures(cd4.sub),verbose = FALSE)
cd4.sub=FindNeighbors(cd4.sub,dims = 1:PCS,verbose = FALSE,k.param = k.param)
cd4.sub=RunUMAP(cd4.sub,dims = 1:PCS,verbose = FALSE) 
cd4.sub=FindClusters(cd4.sub,resolution = 0.8,verbose = FALSE)
DimPlot(cd4.sub,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
FeaturePlot(cd4.sub,features = c('CD3D','CD4','CD8A','IL7R','CCR7','TCF7','LEF1','AREG'))+
  DimPlot(cd4.sub,label = TRUE)

cd4.naive = subset(cd4.sub,seurat_clusters %in% c(0,1,2,3,4,5,8,9,10))
saveRDS(cd4.naive,'04 CD4T Naive ICC PBMC.rds')

VlnPlot(object = cd4.naive, 
        features = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R'))

df = FetchData(cd4.naive,vars = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R'))
df = melt(df)

df2 = df[sample(c(1:22314),3000),]
ggplot(df2, aes(x = variable, y = value,fill=variable)) +
  geom_violin() +
  #geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  ggtitle('ICC_P')+
  #scale_fill_manual(values = c( '#6A91B7','#D18796')) +
  theme_pubr()

ggplot(df, aes(x = variable, y = value,fill=variable))+
  geom_jitter(width=0.3,alpha=0.3,size=0.1)+
  geom_violin(width=0.3,position=position_dodge(width=3))+
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  ggtitle('ICC_PBMC')+
  scale_fill_manual(values = c( '#B4505C','#D6A1A9','#A6C5B5','#FFE8CF','#A3937C','#4D3D30')) +
  theme_pubr()

#### Part 3 DCC Primary Naive T ####

data1 = Read10X('...\\230498DCC\\230498 primary focus')

data1 = CreateSeuratObject(data1,project = 'ICC1_Primary')
data = data1
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(object = data)

PCS = 30
k.param = 20
data = RunPCA(data,features = VariableFeatures(data),verbose = FALSE)
data=FindNeighbors(data,dims = 1:PCS,verbose = FALSE,k.param = k.param)
data=RunUMAP(data,dims = 1:PCS,verbose = FALSE) 
data=FindClusters(data,resolution = 0.8,verbose = FALSE)

DimPlot(data,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
saveRDS(data,'01 DCC Primary TOTAL.rds')
data = readRDS('...\\01 DCC Primary TOTAL.rds')
data@project.name = 'DCC_Primary'
saveRDS(data,'...\\01 DCC Primary TOTAL.rds')

cd4.sub = subset(data,seurat_clusters %in% c(3,7,13))
cd4.sub <- NormalizeData(cd4.sub)
cd4.sub <- FindVariableFeatures(cd4.sub, selection.method = "vst", nfeatures = 2000)
cd4.sub <- ScaleData(object = cd4.sub)
PCS = 30
k.param = 20
cd4.sub = RunPCA(cd4.sub,features = VariableFeatures(cd4.sub),verbose = FALSE)
cd4.sub=FindNeighbors(cd4.sub,dims = 1:PCS,verbose = FALSE,k.param = k.param)
cd4.sub=RunUMAP(cd4.sub,dims = 1:PCS,verbose = FALSE) 
cd4.sub=FindClusters(cd4.sub,resolution = 0.8,verbose = FALSE)
DimPlot(cd4.sub,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
FeaturePlot(cd4.sub,features = c('CD3D','CD4','CD8A','IL7R','CCR7','TCF7','LEF1','AREG'))+
  DimPlot(cd4.sub,label = TRUE)

cd4.naive = subset(cd4.sub,seurat_clusters %in% c(1))
saveRDS(cd4.naive,'02 CD4T Naive DCC Primary.rds')

VlnPlot(object = cd4.naive, 
        features = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R'))

df = FetchData(cd4.naive,vars = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R'))
df = melt(df)

df2 = df[sample(c(1:22314),3000),]
ggplot(df2, aes(x = variable, y = value,fill=variable)) +
  geom_violin() +
  #geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  ggtitle('ICC_P')+
  #scale_fill_manual(values = c( '#6A91B7','#D18796')) +
  theme_pubr()

df.mean = df %>% group_by(variable) %>% summarise_each(funs = mean)
a=ggplot(df, aes(x = variable, y = value,fill=variable))+
  geom_jitter(width=0.3,alpha=0.3,size=0.1)+
  geom_violin(width=1,position=position_dodge(width=3))+
  geom_point(data=df.mean,aes(x = variable, y = value),size=1.5,shape=21, fill = "white")+
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  ggtitle('DCC_Primary')+
  scale_fill_manual(values = c( '#EDA60E','#588b71','#DFD4C0','#679EB8','#E393A0','#d24f6b')) +
  theme_pubr()

pdf(file='DNN Primary.pdf',
    height=2.5,width=4)
print(a)
dev.off()

#### Part 4 DCC PBMC Naive T ####

data = Read10X('...\\230498DCC\\230498 PBMC')
data = CreateSeuratObject(data,project = 'DCC_PBMC')
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(object = data)
PCS = 30
k.param = 20
data = RunPCA(data,features = VariableFeatures(data),verbose = FALSE)
data=FindNeighbors(data,dims = 1:PCS,verbose = FALSE,k.param = k.param)
data=RunUMAP(data,dims = 1:PCS,verbose = FALSE) 
data=FindClusters(data,resolution = 0.8,verbose = FALSE)

DimPlot(data,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
saveRDS(data,'...\\05 DCC Pbmc TOTAL.rds')

FeaturePlot(data,features = c('CD3D','CD4','CD8A','IL7R','CCR7','TCF7','LEF1','AREG'))+
  DimPlot(data,label = TRUE)

cd4.sub = subset(data,seurat_clusters %in% c(4,6))
cd4.sub <- NormalizeData(cd4.sub)
cd4.sub <- FindVariableFeatures(cd4.sub, selection.method = "vst", nfeatures = 2000)
cd4.sub <- ScaleData(object = cd4.sub)
PCS = 30
k.param = 20
cd4.sub = RunPCA(cd4.sub,features = VariableFeatures(cd4.sub),verbose = FALSE)
cd4.sub=FindNeighbors(cd4.sub,dims = 1:PCS,verbose = FALSE,k.param = k.param)
cd4.sub=RunUMAP(cd4.sub,dims = 1:PCS,verbose = FALSE) 
cd4.sub=FindClusters(cd4.sub,resolution = 0.8,verbose = FALSE)
DimPlot(cd4.sub,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
FeaturePlot(cd4.sub,features = c('CD3D','CD4','CD8A','IL7R','CCR7','TCF7','LEF1','AREG'))+
  DimPlot(cd4.sub,label = TRUE)

cd4.naive = cd4.sub
saveRDS(cd4.naive,'...\\06 DCC Pbmc CD4 NAIVE.rds')

#### Part 5 ICC DCC Primary Naive T ####
dcc.p= readRDS('...\\02 CD4T Naive DCC Primary.rds')
icc.p= readRDS('...\\02 CD4T Naive ICC.rds')

dcc.p$group = 'DCC_Primary'
icc.p$group = 'ICC_Primary'

data = merge(icc.p,dcc.p)

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

pdf(file='...\\ICC DCC Primary.pdf',
    height=4,width=4)
print(a)
dev.off()

#### Part 6 ICC DCC PBMC Naive T ####

dnn.pbmc = readRDS('...\\06 DCC Pbmc CD4 NAIVE.rds')
dnn.pbmc$group = 'DCC_PBMC'

icc.pbmc = readRDS('...\\04 CD4T Naive ICC PBMC.rds')
icc.pbmc$group = 'ICC_PBMC'

data = merge(dnn.pbmc,icc.pbmc)
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

pdf(file='...\\DCC ICC PBMC.pdf',
    height=4,width=4)
print(a)
dev.off()

#### Part 7 Split Violin Plot ####


library(Seurat)
dcc.p= readRDS('...\\02 CD4T Naive DCC Primary.rds')
icc.p= readRDS('...\\02 CD4T Naive ICC.rds')

icc.pbmc = readRDS('...\\04 CD4T Naive ICC PBMC.rds')
dcc.pbmc = readRDS('...\\06 DCC Pbmc CD4 NAIVE.rds')

dcc.p$group = 'Primary'
icc.p$group = 'Primary'
dcc.p$class = 'DCC'
icc.p$class = 'ICC'

icc.pbmc$group = 'PBMC'
dcc.pbmc$group = 'PBMC'
dcc.pbmc$class = 'DCC'
icc.pbmc$class = 'ICC'

data1 = merge(icc.p,dcc.p)
data2 = merge(icc.pbmc,dcc.pbmc)
data = merge(data1,data2)

df = FetchData(data,vars = c("AREG", "TCF7", "CCR7",'LEF1','SELL','IL7R','group','class'))
df = melt(df)
df2 = df[sample(c(1:dim(df)[1]),2000),]
df.mean = df2 %>% group_by(variable,group,class) %>% summarise_each(funs = mean)
library(ggpubr)
source("...\\split_violin_ggplot.R")
a=ggplot(df2, aes(x = variable, y = value, fill = group))+
  geom_split_violin(width=0.5)+
  geom_point(data = df.mean, aes(x = variable, y = value), size = 1, shape = 21, fill = "white",
             position=position_dodge(0.4))+
  theme_pubr()+
  stat_compare_means(aes(group = group),
                     label = "p.format",
                     method = "wilcox",
                     label.y = 4,
                     hide.ns = T)+
  facet_wrap(~class, ncol = 1)


pdf(file='...\\Split_violin.pdf',
    height=5,width=6)
print(a)
dev.off()

















