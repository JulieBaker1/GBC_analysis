library(Seurat)
library(dplyr)
library(reshape2)

expr=readRDS('...\\T_cell.RDS')
nms=rownames(expr)
ig_genes = c(grep("^IGJ", nms, v=T),  
             grep("^IGH",nms,v=T), 
             grep("^IGK", nms, v=T),  
             grep("^IGL", nms, v=T)) 
bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),
                     "NEAT1","TMSB4X", "TMSB10", ig_genes)) 

genes_nature= read.table('...\\Neutrophils\\genes.txt') # nature article gene list (removed)
bad_genes=unique(c(genes_nature$V1,bad_genes))

expr=expr[-which(rownames(expr) %in% bad_genes),]
expr <- NormalizeData(expr)
expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
expr <- ScaleData(object = expr, verbose = T)

PCS = 30
k.param = 20
expr = RunPCA(expr,features = VariableFeatures(expr),verbose = FALSE)
expr=FindNeighbors(expr,dims = 1:PCS,verbose = FALSE,k.param = k.param)
expr=RunUMAP(expr,dims = 1:PCS,verbose = FALSE) 
expr=FindClusters(expr,resolution = 0.5
                  ,verbose = FALSE)
DimPlot(expr)
FeaturePlot(expr,features = c('CD3D','CD4','CCR7','SELL','AREG','LEF1','TCF7'))
saveRDS(expr,'...\\nature_t_aftercluster.rds')

expr=readRDS('...\\nature_t_aftercluster.rds')
DimPlot(expr,label = TRUE)+FeaturePlot(expr,features = c('CD4','CD8A'))
marker = FindAllMarkers(expr,only.pos = TRUE)

plot.list = list()
j=1
for (i in c(0:27)){
  sub = subset(expr,seurat_clusters ==i)
  a= FeaturePlot(sub,features = c('CD4','CD8A'))
  plot.list[[j]] = a 
  j=j+1
}

expr.sub = subset(expr,seurat_clusters == 12)
FeaturePlot(expr.sub ,features = c('CD4','AREG','SELL','TCF7'),cols = c("lightgrey" ,"#DE1F1F"))
FeaturePlot(expr,features = c('AREG','FOXP3'))

cd4 = subset(expr, seurat_clusters %in% c(0,3,8,12,14,20))
cd4 <- FindVariableFeatures(cd4, selection.method = "vst", nfeatures = 2000)
cd4 <- ScaleData(object = cd4, verbose = T)
PCS = 30
k.param = 20
cd4 = RunPCA(cd4,features = VariableFeatures(cd4),verbose = FALSE)
cd4=FindNeighbors(cd4,dims = 1:PCS,verbose = FALSE,k.param = k.param)
cd4=RunUMAP(cd4,dims = 1:PCS,verbose = FALSE) 
cd4=FindClusters(cd4,resolution = 0.5
                 ,verbose = FALSE)

cd4.2 = subset(cd4,seurat_clusters %in% c(1:4,6:11,13:16,19:24))
cd4.2 <- FindVariableFeatures(cd4.2, selection.method = "vst", nfeatures = 2000)
cd4.2 <- ScaleData(object = cd4.2, verbose = T)
PCS = 30
k.param = 20
cd4.2 = RunPCA(cd4.2,features = VariableFeatures(cd4.2),verbose = FALSE)
cd4.2=FindNeighbors(cd4.2,dims = 1:PCS,verbose = FALSE,k.param = k.param)
cd4.2=RunUMAP(cd4.2,dims = 1:PCS,verbose = FALSE) 
cd4.2=FindClusters(cd4.2,resolution = 0.1
                   ,verbose = FALSE)
FeaturePlot(cd4.2,features = c('CD4','SELL','CCR7','LEF1','AREG'))+DimPlot(cd4.2,label = TRUE)

marker = FindAllMarkers(cd4.2,only.pos = TRUE)
saveRDS(marker,'D:\\GBC\\18 Final Final code file\\29 HCC\\01 CD4 marker.rds')
new.cluster.ids <- c('Tcm', #0
                     "Treg",#1
                     "Tem",#2 
                     "Th1",#3
                     "Tn_PBMC",#4
                     "Tn_HCC",#5
                     "Treg",#6
                     "Unknow",#7
                     "Treg"#8 
)

names(new.cluster.ids) <- levels(cd4.2)
cd4.2<- RenameIdents(cd4.2, new.cluster.ids)
cd4.2$celltype = cd4.2@active.ident

DimPlot(cd4.2,label = TRUE)
#cd4.2 = subset(cd4.2,Tissue == 'liver')

naive = c("CCR7","TCF7","SELL","IL7R" ,"IL2RG",'LEF1')
func.gene = list()
func.gene[['naive']]=naive
sign.name = names(func.gene)
for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(cd4.2))
  a=list(a)
  cd4.2<- AddModuleScore(cd4.2,features = a,
                         ctrl = 100,
                         name = i)
}

df = FetchData(cd4.2,vars = c('naive1','celltype','AREG'))
df = df %>% group_by(celltype) %>% summarise_each(funs = mean)
df = as.data.frame(df)
rownames(df) = df$celltype
df = df[,-1]
library(pheatmap)
df = scale(df)

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
df = t(df)
a=pheatmap(df,cluster_rows = FALSE,cluster_cols = FALSE,
           #gaps_col = c(5),
           #gaps_row = c(9,17),
           cellwidth = 15, cellheight = 15,
           color = c(colorRampPalette(colors = c("#467897","white"))(length(bk)/2),
                     colorRampPalette(colors = c("white","#BD3438"))(length(bk)/2)),
           #legend_breaks=seq(-1,0,1),
           breaks=bk)


saveRDS(cd4.2,'...\\01 CD4 Data.rds')
cd4.2  = readRDS('...\\01 CD4 Data.rds')

df.hcc = data.frame(patient= cd4.2$Sample.Name,celltype=cd4.2$celltype)
df.hcc$value = 1
df.hcc=dcast(df.hcc,patient~celltype,fun.aggregate=sum)
rownames(df.hcc)=df.hcc$patient
df.hcc = df.hcc[,-1]

fn = function(x){
  return(-x*log(x))
}
softmax=function(x){
  (x/sum(x))
} 
fn2 = function(x){
  return(sum(x)/log(length(x)))
}

a=apply(df.hcc,2,softmax)
a=apply(a, 2, fn)
a[is.na(a)]=0
a=apply(a, 2, fn2)
shang.value=as.data.frame(a)
shang.value$celltype=rownames(shang.value)

df = FetchData(cd4.2,vars = c('naive1','celltype','AREG'))
df = df %>% group_by(celltype) %>% summarise_each(funs = mean)
df = as.data.frame(df)
rownames(df) = df$celltype
df = df[,-1]
df = scale(df)
df = as.data.frame(df)
df$celltype = rownames(df)

df = df%>% left_join(shang.value)
colnames(df) = c('Naive','AREG','Celltype','Shanno.Value','Cluster')
df2= melt(df,id.vars = c('Shanno.Value','Celltype','Cluster'))
df2 = df2[-c(7,14),]
df2$Cluster = factor(df2$Cluster,levels = c('Tcm','Treg','Tem','Th1','Tn_PBMC','Tn_HCC'))
a=ggplot(df2, aes(Cluster,variable,fill=value,size=Shanno.Value,color=value)) +
  geom_point(alpha=1)+
  scale_color_gradient2(low = "#467897", mid = "gray", high = "#BD3438")+
  theme(plot.background = element_blank())+
  theme(panel.grid = element_blank(), panel.background = element_blank())+
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray", linetype = "dashed"))+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  scale_size(range = c(1, 10))








