
setwd('...\\GSE126030_RAW(1)')
data1 = read.table('PP005swap.filtered.matrix.txt',header = TRUE)
data2 = read.table('PP006swap.filtered.matrix.txt',header = TRUE)
cd4= readRDS('...\\cd4_newname.rds')

gene = unique(data1$Gene)
gene = intersect(gene,rownames(cd4))

data1 = data1[data1$Gene %in% gene,]
data1<-data1[!duplicated(data1$Gene, fromLast=TRUE), ] 
rownames(data1)=data1$Gene
data1 = data1[,-c(1,2)]

obj.1 = CreateSeuratObject(counts = data1, project = 'LN Donor 1')
obj.1 <- NormalizeData(obj.1)

obj.1 <- FindVariableFeatures(obj.1, selection.method = "vst", nfeatures = 2000)
obj.1 <- ScaleData(object = obj.1)

PCS = 30
k.param = 20
obj.1 = RunPCA(obj.1,features = VariableFeatures(obj.1),verbose = FALSE)
obj.1=FindNeighbors(obj.1,dims = 1:PCS,verbose = FALSE,k.param = k.param)
obj.1=RunUMAP(obj.1,dims = 1:PCS,verbose = FALSE) 
obj.1=FindClusters(obj.1,resolution = 0.5,verbose = FALSE)
DimPlot(obj.1,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
marker = FindAllMarkers(obj.1,only.pos = TRUE)

FeaturePlot(obj.1,features = c('CD4','CCR7','LEF1','AREG','SELL','CD8A'))


data2 = data2[data2$Gene %in% gene,]
data2<-data2[!duplicated(data2$Gene, fromLast=TRUE), ] 
rownames(data2)=data2$Gene
data2 = data2[,-c(1,2)]

obj.2 = CreateSeuratObject(counts = data2, project = 'LN Donor 2')
obj.2 <- NormalizeData(obj.2)
obj.2 <- FindVariableFeatures(obj.2, selection.method = "vst", nfeatures = 2000)
obj.2 <- ScaleData(object = obj.2)

PCS = 30
k.param = 20
obj.2 = RunPCA(obj.2,features = VariableFeatures(obj.2),verbose = FALSE)
obj.2=FindNeighbors(obj.2,dims = 1:PCS,verbose = FALSE,k.param = k.param)
obj.2=RunUMAP(obj.2,dims = 1:PCS,verbose = FALSE) 
obj.2=FindClusters(obj.2,resolution = 0.5,verbose = FALSE)

FeaturePlot(obj.2,features = c('CD4','CCR7','LEF1','AREG','SELL','CD8A'))
DimPlot(obj.2,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
marker = FindAllMarkers(obj.2,only.pos = TRUE)

donor = merge(obj.1,obj.2,add.cell.ids =c('LN_1','LN_2'))
donor <- FindVariableFeatures(donor, selection.method = "vst", nfeatures = 2000)
donor <- ScaleData(object = donor)
PCS = 30
k.param = 20
donor = RunPCA(donor,features = VariableFeatures(donor),verbose = FALSE)
donor=FindNeighbors(donor,dims = 1:PCS,verbose = FALSE,k.param = k.param)
donor=RunUMAP(donor,dims = 1:PCS,verbose = FALSE) 
donor=FindClusters(donor,resolution = 0.5,verbose = FALSE)

marker = FindAllMarkers(donor,only.pos = TRUE)
FeaturePlot(donor,features = c('CD4','CCR7','LEF1','AREG','SELL','CD8A'))+DimPlot(donor,label = TRUE)

saveRDS(donor,'Lymphoid Node Donor1.rds')
donor = readRDS('Lymphoid Node Donor1.rds')

naive = c("CCR7","TCF7","SELL","IL7R" ,"IL2RG",'LEF1')
func.gene = list()
func.gene[['naive']]=naive
sign.name = names(func.gene)
for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(donor))
  a=list(a)
  donor<- AddModuleScore(donor,features = a,
                         ctrl = 100,
                         name = i)
}

df = FetchData(donor,vars = c('CD4','naive1','seurat_clusters','AREG'))
df = df %>% group_by(seurat_clusters) %>% summarise_each(funs = mean)
df = as.data.frame(df)
rownames(df) = df$seurat_clusters
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
pdf(file = '03 Pheatmap Donor Three Sig.pdf',
    width =6,height = 4)
print(a)
dev.off()

naive.ln.t = subset(donor,seurat_clusters %in% c(0,4))
saveRDS(naive.ln.t,'01 LN Naive Donor1.rds')
naive.ln.t =readRDS('01 LN Naive Donor1.rds')
naive.ln.t$group = 'Normal_LN'
cd4= readRDS('...\\cd4_newname.rds')
cd4.naive = subset(cd4,seurat_clusters == 2)
DimPlot(cd4.naive)

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
df=data.frame(NewSample.ID =  cd4.naive$orig.ident)
df$cellid = rownames(df)
df = df %>% left_join(info.patient)
df = filter(df,Tumors.for.scRNA.seq.short=='LN',histological.type.short=='adeno')
cd4.naive.ln = cd4.naive[,df$cellid]
cd4.naive.ln$group = 'GBC'
table(cd4.naive.ln$orig.ident)

total = merge(cd4.naive.ln,naive.ln.t)
total <- NormalizeData(total)
total <- FindVariableFeatures(total, selection.method = "vst", nfeatures = 2000)
total <- ScaleData(object = total)
PCS = 30
k.param = 20
total = RunPCA(total,features = VariableFeatures(total),verbose = FALSE)
total=FindNeighbors(total,dims = 1:PCS,verbose = FALSE,k.param = k.param)
total=RunUMAP(total,dims = 1:PCS,verbose = FALSE) 
total=FindClusters(total,resolution = 0.5,verbose = FALSE)
DimPlot(total,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
saveRDS(total,'02 Total Merge GBC LN Naive.rds')
Idents(total) = total$group
marker = FindAllMarkers(total,only.pos = TRUE)

marker.1= subset(marker,cluster == 'GBC')
marker.2= subset(marker,cluster == 'Normal_LN')
colnames(marker.2)
marker.2 = marker.2[,c(1,2,4,3,5,6,7)]
colnames(marker.2)=colnames(marker.1)
marker.plot = rbind(marker.1,marker.2)
marker.plot = subset(marker.plot,p_val_adj<0.05)

ggplot(marker.plot,aes(pct.1,pct.2,size=avg_log2FC,color = cluster,shape = cluster))+
  geom_point(alpha=0.8)+
  scale_shape_manual(values = c(15, 16))+
  scale_size(range = c(0, 4))+
  scale_color_manual(values = c('#2b6a99','#1b7c3d'))

marker.plot$group = case_when(marker.plot$pct.1+0.25<marker.plot$pct.2~'Color',
                              marker.plot$pct.1-0.25>marker.plot$pct.2~'Color')
marker.plot[is.na(marker.plot$group),'group'] = 'gray'

marker.plot1 = subset(marker.plot,group == 'Color')
marker.plot2 = subset(marker.plot,group == 'gray')

a=ggplot(marker.plot, aes(pct.1, pct.2, size = avg_log2FC, color = cluster, shape = cluster)) +
  geom_point(alpha = 0.8) +
  geom_point(data = subset(marker.plot, !(pct.1 + 0.25 <= pct.2  | pct.1 - 0.25>= pct.2)), color = "gray") +
  geom_abline(intercept = 0.25, slope = 1, linetype = "dashed", color = "black") +
  geom_abline(intercept = -0.25, slope = 1, linetype = "dashed", color = "black") +
  scale_size(range = c(1, 5)) +
  scale_shape_manual(values = c(15, 16)) +
  scale_color_manual(values = c('#2b6a99', '#1b7c3d'))+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black"),
    axis.line.y = element_line(color = "black",
                               size = 0.5),
    axis.line.x =element_line(color = "black",
                              size =0.5),
    panel.grid = element_blank()
  )



FeaturePlot(total,features = 'CCR7',pt.size = 0.3,max.cutoff = 4)
FeaturePlot(total,features = 'LEF1',pt.size = 0.3,max.cutoff = 4)
FeaturePlot(total,features = 'SELL',pt.size = 0.3,max.cutoff = 4)
FeaturePlot(total,features = 'AREG',pt.size = 0.3,max.cutoff = 4)

library(org.Hs.eg.db)
library(clusterProfiler)
marker$enterz <- mapIds(org.Hs.eg.db, keys = marker$gene, keytype = "SYMBOL", column="ENTREZID")
enrich <- compareCluster(enterz ~cluster, data=marker,
                         fun="enrichGO",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP" ,
                         pAdjustMethod = "BH",
                         readabl = TRUE)
dotplot(enrich,showCategory=10)
result=enrich@compareClusterResult

interest.pathway=c('regulation of epidermal growth factor-activated receptor activity',
                   'regulation of epithelial cell proliferation',
                   'lymphocyte proliferation',
                   'aerobic respiration')

enrich2=filter(result,Description %in% interest.pathway)

enrich2.df <- data.frame(enrich2) 
a=enrich2.df$GeneRatio
a <- data.frame(name=a)
b <- apply(a,1,function(x) eval(parse(text=x)))
enrich2.df$GeneRatio=b

a=ggplot(enrich2.df,aes(x=cluster,y=Description,size=GeneRatio,color=pvalue))+
  geom_point()+
  scale_colour_gradient(high= "#408591",low="#C6736E")+
  theme(axis.text.x = element_text(size=10,colour = 'black',face='bold',angle = 30),
        axis.text.y = element_text(size=10))





