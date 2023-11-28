library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

expr=readRDS('...\\T_NK.RDS')
nms=rownames(expr)
ig_genes = c(grep("^IGJ", nms, v=T),  
             grep("^IGH",nms,v=T), 
             grep("^IGK", nms, v=T),  
             grep("^IGL", nms, v=T)) 
bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes)) 
expr=expr[-which(rownames(expr) %in% bad_genes),]

expr <- NormalizeData(expr)
vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent")
expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
expr <- ScaleData(object = expr, vars.to.regress = vars.to.regress, verbose = F)

PCS = 30
k.param = 20
expr = RunPCA(expr,features = VariableFeatures(expr),verbose = FALSE)
expr=FindNeighbors(expr,dims = 1:PCS,verbose = FALSE,k.param = k.param)
expr=RunUMAP(expr,dims = 1:PCS,verbose = FALSE) 
expr=FindClusters(expr,resolution = 0.5,verbose = FALSE)

# analysis cycling cluster
cycling=subset(expr,seurat_clusters=='8')
cycling <- FindVariableFeatures(cycling, selection.method = "vst", nfeatures = 2000)
cycling <- ScaleData(object = cycling, vars.to.regress = vars.to.regress, verbose = F)

PCS = 30
k.param = 20
cycling = RunPCA(cycling,features = VariableFeatures(cycling),verbose = FALSE)
cycling=FindNeighbors(cycling,dims = 1:PCS,verbose = FALSE,k.param = k.param)
cycling=RunUMAP(cycling,dims = 1:PCS,verbose = FALSE) 
cycling=FindClusters(cycling,resolution = 0.3,verbose = FALSE)
DimPlot(cycling,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')

## extract NK cell
cycling_nk=subset(cycling,seurat_clusters=='4')
DimPlot(cycling_nk,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
nk=subset(expr,seurat_clusters %in% c('7','13'))
DimPlot(nk,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
nk_total=merge(nk,cycling_nk)
nk_total <- FindVariableFeatures(nk_total, selection.method = "vst", nfeatures = 2000)
nk_total <- ScaleData(object = nk_total, vars.to.regress = vars.to.regress, verbose = F)
PCS = 30
k.param = 20
nk_total = RunPCA(nk_total,features = VariableFeatures(nk_total),verbose = FALSE)
nk_total=FindNeighbors(nk_total,dims = 1:PCS,verbose = FALSE,k.param = k.param)
nk_total=RunUMAP(nk_total,dims = 1:PCS,verbose = FALSE) 
nk_total=FindClusters(nk_total,resolution = 0.5,verbose = FALSE)
DimPlot(nk_total,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
saveRDS(nk_total,'...\\NK_0.5resolution.RDS')

##extract T cell 
nk_id=colnames(nk_total)
expr=readRDS('...\\T_NK.RDS')
nms=rownames(expr)
ig_genes = c(grep("^IGJ", nms, v=T),  
             grep("^IGH",nms,v=T), 
             grep("^IGK", nms, v=T),  
             grep("^IGL", nms, v=T)) 
bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes)) 
expr=expr[-which(rownames(expr) %in% bad_genes),]
expr=expr[,-which(colnames(expr) %in% nk_id)]
expr <- NormalizeData(expr)

vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent")
expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
expr <- ScaleData(object = expr, vars.to.regress = vars.to.regress, verbose = F)

PCS = 30
k.param = 20
expr = RunPCA(expr,features = VariableFeatures(expr),verbose = FALSE)
expr=FindNeighbors(expr,dims = 1:PCS,verbose = FALSE,k.param = k.param)
expr=RunUMAP(expr,dims = 1:PCS,verbose = FALSE) 
expr=FindClusters(expr,resolution = 0.5,verbose = FALSE)
DimPlot(expr,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
saveRDS(expr,'...\\T_0.5resolution.RDS')

###Logistic Regression###

setwd('...')
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
expr=readRDS('T_0.5resolution.RDS')

save_path='...'
marker = c('CD8A')
combName=marker
cluster=c('0','1','2','3','4','5','6','7','8','9','10',
          '11','12','13','14','15','16','17','18','19','20','21')
pdf(paste0(save_path,"/",combName,".pdf"))
for(i in cluster){
  cluster_sub=subset(expr,seurat_clusters==i)
  print(FeaturePlot(cluster_sub,features = 'CD8A',reduction = "umap"))
}
dev.off()

marker=FindMarkers(expr,ident.1 = c(1,4,7,8),ident.2 = c(2,3,5),
                   only.pos =TRUE,,min.pct=0.25,logfc.threshold=0.25,test.use = "MAST")

marker_2=FindMarkers(expr,ident.1 = c(2,3,5),ident.2 =c(1,4,7,8),
                     only.pos =TRUE,,min.pct=0.25,logfc.threshold=0.25,test.use = "MAST")

top50_1=marker%>%top_n(n=50,wt=avg_log2FC)
top50_2=marker_2%>%top_n(n=50,wt=avg_log2FC)


common_marker=union(rownames(top50_2),rownames(top50_1))
common_marker=c(common_marker,'CD4')
## Training dataset and validation dataset
training_data1=expr[,expr@meta.data$seurat_clusters %in% c(1,4,7,8)]
training_data2=expr[,expr@meta.data$seurat_clusters %in% c(2,3,5)]

set.seed(777)
training_data1=training_data1[,sample(c(1:100714),size=85381)]

index_training=sample(c(1:85381),size=59767)
index_validation=c(1:85381)[-index_training]


train_data1=training_data1[common_marker,index_training]
train_data2=training_data2[common_marker,index_training]

train_data1_fram=as.data.frame(t(as.matrix(train_data1@assays[["RNA"]]@data)))
train_data1_fram$celltype=1
train_data2_fram=as.data.frame(t(as.matrix(train_data2@assays[["RNA"]]@data)))
train_data2_fram$celltype=0

train_total=rbind(train_data1_fram,train_data2_fram)
model=glm(celltype~.,family = binomial(link = 'logit'),data=train_total)
summary(model)

# evaluate the classifier performance
pred<-predict(model)
prob<-exp(pred)/(1+exp(pred))
yhat<-1*(prob>0.5)
table(train_total$celltype,yhat)

validation_data1=training_data1[common_marker,index_validation]
validation_data2=training_data2[common_marker,index_validation]

vali_data1_fram=as.data.frame(t(as.matrix(validation_data1@assays[["RNA"]]@data)))
vali_data1_fram$celltype=1
vali_data2_fram=as.data.frame(t(as.matrix(validation_data2@assays[["RNA"]]@data)))
vali_data2_fram$celltype=0

vali_totol=rbind(vali_data1_fram,vali_data2_fram)

pred2<-predict(model,newdata = vali_totol)
prob2<-exp(pred2)/(1+exp(pred2))
yhat<-1*(prob2>0.5)
table(vali_totol$celltype,yhat)

library(ROCR)
pred2<-prediction(pred2,vali_totol$celltype)
performance(pred2,'auc')@y.values

perf<-performance(pred2,'tpr','fpr')
plot(perf,col=2,type="l",lwd=2)
f=function(x){ y=x;return(y)}
curve(f(x),0,1,col=4,lwd=2,lty=2,add=T)

saveRDS(model,'logistic_regression.rds')

##perform the classifier with unclassical clusters

save_path='...'
combName='0.5_classification'
pdf(paste0(save_path,"/",combName,".pdf"))
for (i in c(0,6,9,10,11,12,13,14,15,16,17,19,20,21)){
  cluster=subset(expr,seurat_clusters==i)
  DimPlot(cluster)
  FeaturePlot(cluster,features = 'CD4')
  FeaturePlot(cluster,features = 'CD8A')
  cluster=cluster[common_marker,]
  testing_data=as.data.frame(t(as.matrix(cluster@assays[["RNA"]]@data)))
  result<-predict(model,newdata = testing_data)
  prob<-exp(result)/(1+exp(result))
  prob=as.data.frame(prob)
  prob$cellname=colnames(cluster)
  prob$predict_define <- ifelse(prob$prob<0.5,"CD4","CD8")
  cluster@meta.data$predict=prob$prob
  cluster@meta.data$predict_define=prob$predict_define 
  FeaturePlot(cluster,features = 'predict')
  print(DimPlot(cluster,group.by = 'predict_define')+ggtitle(i) +
          theme(plot.title = element_text(hjust = 0.5)))
}
dev.off()


##extract CD4 T and CD8 T 
cluster_cd4_name=c()
cluster_cd8_name=c()
for (i in c(0,6,9,10,11,12,13,14,15,16,17,19,20,21)){
  cluster=subset(expr,seurat_clusters==i)
  DimPlot(cluster)
  FeaturePlot(cluster,features = 'CD4')
  FeaturePlot(cluster,features = 'CD8A')
  cluster=cluster[common_marker,]
  testing_data=as.data.frame(t(as.matrix(cluster@assays[["RNA"]]@data)))
  result<-predict(model,newdata = testing_data)
  prob<-exp(result)/(1+exp(result))
  prob=as.data.frame(prob)
  prob$cellname=colnames(cluster)
  prob$predict_define <- ifelse(prob$prob<0.5,"CD4","CD8")
  cluster@meta.data$predict=prob$prob
  cluster@meta.data$predict_define=prob$predict_define 
  cluster_cd4_name=union(cluster_cd4_name,colnames(subset(cluster,predict_define=='CD4'))) 
  cluster_cd8_name=union(cluster_cd8_name,colnames(subset(cluster,predict_define=='CD8')))
}


training_data1=expr[common_marker,expr@meta.data$seurat_clusters %in% c(1,4,7,8)]
testing_data=as.data.frame(t(as.matrix(training_data1@assays[["RNA"]]@data)))
result<-predict(model,newdata = testing_data)
prob<-exp(result)/(1+exp(result))
prob=as.data.frame(prob)
prob$cellname=colnames(training_data1)
prob$predict_define <- ifelse(prob$prob<0.5,"CD4","CD8")
training_data1@meta.data$predict=prob$prob
training_data1@meta.data$predict_define=prob$predict_define 
cluster_cd4_name2=union(cluster_cd4_name,colnames(subset(training_data1,predict_define=='CD4'))) 
cluster_cd8_name2=union(cluster_cd8_name,colnames(subset(training_data1,predict_define=='CD8')))

training_data2=expr[common_marker,expr@meta.data$seurat_clusters %in% c(2,3,5)]
testing_data=as.data.frame(t(as.matrix(training_data2@assays[["RNA"]]@data)))
result<-predict(model,newdata = testing_data)
prob<-exp(result)/(1+exp(result))
prob=as.data.frame(prob)
prob$cellname=colnames(training_data2)
prob$predict_define <- ifelse(prob$prob<0.5,"CD4","CD8")
training_data2@meta.data$predict=prob$prob
training_data2@meta.data$predict_define=prob$predict_define 
cluster_cd4_name2=union(cluster_cd4_name2,colnames(subset(training_data2,predict_define=='CD4'))) 
cluster_cd8_name2=union(cluster_cd8_name2,colnames(subset(training_data2,predict_define=='CD8')))

saveRDS(cluster_cd4_name2,'cd4_name.rds')
saveRDS(cluster_cd8_name2,'cd8_name.rds')
