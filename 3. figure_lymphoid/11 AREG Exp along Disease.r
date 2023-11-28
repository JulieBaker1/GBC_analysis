
library(Seurat)
library(reshape2)
library(dplyr)
cd4=readRDS('...\\cd4_newname.rds')
cd4 = subset(cd4,seurat_clusters %in% c(2))

meta.data = data.frame(NewSample.ID = cd4$orig.ident)
info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
meta.data = meta.data %>% left_join(info.patient)

cd4$clinical.stage = meta.data$Clinical.stage
cd4$metastasis.type =meta.data$metastasis.type
cd4$Tumors.for.scRNA.seq.short =meta.data$Tumors.for.scRNA.seq.short
cd4$histological.type.short =meta.data$histological.type.short

norm.df= FetchData(cd4,vars = c('AREG','Tumors.for.scRNA.seq.short','histological.type.short'))
HG=filter(norm.df,histological.type.short=='HG')
HG$group='Non-Tumor'

LG=filter(norm.df,histological.type.short=='LG')
LG$group='Non-Tumor'

XGC=filter(norm.df,Tumors.for.scRNA.seq.short=='XGC')
XGC$group='Non-Tumor'

CC=filter(norm.df,Tumors.for.scRNA.seq.short=='CC')
CC$group='Non-Tumor'

normal = rbind(HG,LG,XGC,CC)

tumor.df = FetchData(cd4,vars = c('AREG','Tumors.for.scRNA.seq.short','histological.type.short','clinical.stage',
                                  'metastasis.type'))
tumor = filter(tumor.df,Tumors.for.scRNA.seq.short == 'P',histological.type.short == 'adeno')
tumor$group = 'Adeno'
tumor[tumor$clinical.stage=='IIA','clinical.stage']='II'
tumor[tumor$clinical.stage=='IIB','clinical.stage']='II'

mean.df= tumor %>% group_by(clinical.stage) %>% summarise_each(fun=mean)
mean.df$clinical.stage=factor(mean.df$clinical.stage,
                              levels = c('I','II','IIIA','IIIB','IVB','IVA'))
mean.df = mean.df %>% select('AREG','clinical.stage')
mean.df$group = 'AREG Expression'
colnames(mean.df)[1]='Value'

total.count.info=readRDS('...\\01 Lymphoid Celltype Count.rds')
celltype = colnames(total.count.info)[1:10]
total.count.info = total.count.info[,celltype]

normalize_to1 <- function(x){
  (x/sum(x))*100
} 

total.count.info = apply(total.count.info,1,normalize_to1)
total.count.info = as.data.frame(t(total.count.info))
total.count.info$NewSample.ID = rownames(total.count.info)

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
total.count.info = total.count.info %>% left_join(info.patient)

tumor=total.count.info%>% filter(histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
colnames(tumor)
tumor = tumor %>% select('CD4T_C2_CCR7','Clinical.stage')
tumor[tumor$Clinical.stage=='IIA','Clinical.stage']='II'
tumor[tumor$Clinical.stage=='IIB','Clinical.stage']='II'
mean.df2= tumor %>% group_by(Clinical.stage) %>% summarise_each(fun=mean)
mean.df2 = mean.df2[,c(2,1)]
colnames(mean.df2) = c('Value','clinical.stage')
mean.df2$group = 'Naive T Number'

radar = cbind(mean.df,mean.df2)
rownames(radar) = radar$clinical.stage
radar = radar[,c(1,4)]
colnames(radar) = c('Expression','Number')
radar = data.frame(t(radar))
radar = apply(radar,1,scale)
radar = as.data.frame(t(radar))
colnames(radar) = c('I','II','IIIA','IIIB','IVA','IVB')
radar$group = rownames(radar)
radar = radar[,c(7,1,2,3,4,5,6)]

library(ggradar)
a=ggradar(
  radar, 
  values.radar = c("-2", "0", "2"),
  grid.min = -2, grid.mid =0, grid.max =2,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#3d3b4f", "#ff7500"),
  # Background and grid lines
  #background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
)

###CD4 Naive 

expr = cd4
vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent")

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
df=data.frame(NewSample.ID =  expr$orig.ident)
df$cellid = rownames(df)
df = df %>% left_join(info.patient)
df = filter(df,Tumors.for.scRNA.seq.short=='P',histological.type.short=='adeno')
expr= expr[,df$cellid]

expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
expr <- ScaleData(object = expr, vars.to.regress = vars.to.regress, verbose = F)

PCS = 30
k.param = 20
expr = RunPCA(expr,features = VariableFeatures(expr),verbose = FALSE)
expr=FindNeighbors(expr,dims = 1:PCS,verbose = FALSE,k.param = k.param)
expr=RunUMAP(expr,dims = 1:PCS,verbose = FALSE) 
expr=FindClusters(expr,resolution = 0.5,verbose = FALSE)
DimPlot(expr,reduction = 'umap',label = TRUE,group.by = 'seurat_clusters')
saveRDS(expr,'...\\GBC Tn.rds')
expr = expr[,sample(c(1:dim(expr)[2]),3000)]
