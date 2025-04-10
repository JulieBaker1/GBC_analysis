#### 20240627 ROE ####

library(pheatmap)
library(reshape2)
library(dplyr)

data = readRDS('01 Lymphoid Disease Count.rds')
data = data[,c(1:9)] # cd8t
data = data[,c(10:17)]#cd4t
data = data[,c(18:25)] # nk
data = data[,c(26:30)] #B 

data = readRDS('ZTSubtypes_Counts_Roe_p_230605.RDS')
row.name = rownames(data)
data = as.data.frame(data)
rownames(data) = row.name

data = data[,c(21:26)] #Macrophage
data = data[,c(1,2,3,8,10)]#Neutrophils
data = data[,c(13:15,17:20)] #DC
data = data[,c(30:36)] #Endothelial

fibro.data = readRDS('Fibro_disease_count.RDS')
rownames(fibro.data)
data= as.data.frame(fibro.data)
rownames(data)=rownames(fibro.data)
rownames(data) = c('Adeno','Adeno Squa','CC','HG','LG','Neuro','Squa','Undiff','XGC')
#data = data[,c(1,2,3,13,14,16)]
data = data[c(3,9,5,4,1,2,6,7,8),]

data[is.na(data)]=0
a=data
roe=data
for (i in 1:dim(a)[1]){
  for (j in 1:dim(a)[2]){
    onecelltype_onetissue=a[i,j]
    onetissue_othercelltype=sum(a[i,])-onecelltype_onetissue
    onecelltype_othertissue=sum(a[,j])-onecelltype_onetissue
    othercelltype_ohertissue=sum(a) - onecelltype_onetissue - onetissue_othercelltype - onecelltype_othertissue
    x=matrix(c(onecelltype_onetissue,onecelltype_othertissue,onetissue_othercelltype,othercelltype_ohertissue),
             nrow = 2,ncol = 2)
    b=chisq.test(x)$expected[1,1]
    m=onecelltype_onetissue/b
    roe[i,j]=m
  }
} 

bk <- c(seq(0,1.0,by=0.01),seq(1.1,2,by=0.01))
roe = t(roe)

order = c(3,4,9,5,1,6,8,7,2) #cd8
order = c(5,6,1,2,3,4,7,8) #cd4
order = c(2,5,8,4,6,7,1,3) #nk
order = c(5,2,3,4,1) #b
order = c(3,1,0,4,5,2)+1 #Macrophage
order = c(1,4,5,3,2) #Neutrophil
order = c(1,5,2,3,6,4,7) #dc
order = c(1,2,5,7,3,6,4) #endo
order = c(2,16,3,13,1,14) #Fibro

roe = roe[order,]
a=pheatmap(roe,cluster_rows = FALSE,cluster_cols = FALSE,
           #gaps_col = c(5),
           #gaps_row = c(9,17),
           cellwidth = 25, cellheight = 20,
           color = c(colorRampPalette(colors = c("#404398","white"))(length(bk)/2),
                     colorRampPalette(colors = c("white","#F18932"))(length(bk)/2)),
           #legend_breaks=seq(-1,0,1),
           breaks=bk)



#### site roe analysis ####

setwd('...\\18 ROE Disease and Site All Celltype')
temp.data = readRDS('ZTSubtypes_Counts_Roe_p_230605.RDS')
temp.data = readRDS('ZTSubtypes_Counts_Roe_site_230605.RDS')
temp.data = readRDS('Celltypes_Counts_Roe_p_230117.RDS')
temp.data = readRDS('Celltypes_Counts_Roe_site_230117.RDS')
temp.data = readRDS('01 Lymphoid Site Count.rds')


data= as.data.frame(temp.data)
rownames(data)=rownames(temp.data)
data = data[,c(10:17)] #cd4t
data = data[,c(21:26)] #Macrophage
data = data[,c(1,2,3,8,10)]#Neutrophils
data = data[,c(13:15,17:20)] #DC
data = data[,c(30:36)] #Endothelial

data[is.na(data)]=0
a=data
roe=data
for (i in 1:dim(a)[1]){
  for (j in 1:dim(a)[2]){
    onecelltype_onetissue=a[i,j]
    onetissue_othercelltype=sum(a[i,])-onecelltype_onetissue
    onecelltype_othertissue=sum(a[,j])-onecelltype_onetissue
    othercelltype_ohertissue=sum(a) - onecelltype_onetissue - onetissue_othercelltype - onecelltype_othertissue
    x=matrix(c(onecelltype_onetissue,onecelltype_othertissue,onetissue_othercelltype,othercelltype_ohertissue),
             nrow = 2,ncol = 2)
    b=chisq.test(x)$expected[1,1]
    m=onecelltype_onetissue/b
    roe[i,j]=m
  }
} 


bk <- c(seq(0,1.5,by=0.01),seq(1.6,3,by=0.01))
roe = t(roe)

#order = c(3,1,0,4,5,2)+1 #Macrophage
order = c(1,4,5,3,2) #Neutrophage
order = c(1,5,2,3,6,4,7)
order = c(1,2,5,7,3,6,4)
order = c(1,2,5,9,3,8,7,6,4)
order = c(1,2,6,10,3,9,8,7,5) # celltype Primary
roe = roe[-c(4,11),] # celltype site

roe = roe[order,]
a=pheatmap(roe,cluster_rows = FALSE,cluster_cols = FALSE,
           #gaps_col = c(5),
           #gaps_row = c(9,17),
           cellwidth =15, cellheight = 15,
           color = c(colorRampPalette(colors = c("#3373B2","white"))(length(bk)/2),
                     colorRampPalette(colors = c("white","#BD3438"))(length(bk)/2)),
           #legend_breaks=seq(-1,0,1),
           breaks=bk)









