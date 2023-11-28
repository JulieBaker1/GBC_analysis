
setwd('...\\ ROE Disease and Site All Celltype')
temp.data = readRDS('ZTSubtypes_Counts_Roe_p_230605.RDS')
temp.data = readRDS('ZTSubtypes_Counts_Roe_site_230605.RDS')
temp.data = readRDS('Celltypes_Counts_Roe_p_230117.RDS')
temp.data = readRDS('Celltypes_Counts_Roe_site_230117.RDS')
temp.data = readRDS('01 Lymphoid Disease Count.rds')
temp.data = readRDS('01 Lymphoid Site Count.rds')
fibro.data = readRDS('Fibro_disease_count.RDS')
fibro.data = readRDS('Fibro_Site_count.RDS')


data= as.data.frame(fibro.data)
rownames(data)=rownames(fibro.data)
rownames(data) = c('Adeno','Adeno Squa','CC','HG','LG','Neuro','Squa','Undiff','XGC')
data = data[rownames(temp.data),]
colnames(data)
data = data[,c(1,2,3,13,14,16)]
data = data[,-c(4,11)] #all celltype
data = data[,c(21:26)] #Macrophage
data = data[,c(1,2,3,8,10)]#Neutrophils
data = data[,c(13:15,17:20)] #DC
data = data[,c(30:36)] #Endothelial
data = data[,c(1:9)] # cd8t
data = data[,c(10:17)]#cd4t
data = data[,c(18:25)] # nk
data = data[,c(26:30)] #B 

data[is.na(data)]=0
a=data
roe=data
for (i in 1:dim(a)[1]){
  for (j in 1:dim(a)[2]){
    onecelltype_onetissue=a[i,j]
    onetissue_othercelltype=sum(a[i,])-onecelltype_onetissue
    onecelltype_othertissue=sum(a[,j])-onecelltype_onetissue
    othercelltype_ohertissue=sum(a)-onecelltype_onetissue
    x=matrix(c(onecelltype_onetissue,onecelltype_othertissue,onetissue_othercelltype,othercelltype_ohertissue),
             nrow = 2,ncol = 2)
    b=chisq.test(x)$expected[1,1]
    m=onecelltype_onetissue/b
    roe[i,j]=m
  }
} 


bk <- c(seq(0,1.4,by=0.01),seq(1.5,3,by=0.01))
roe = t(roe)
a=pheatmap(roe,cluster_rows = FALSE,cluster_cols = FALSE,
           #gaps_col = c(5),
           #gaps_row = c(9,17),
           cellwidth = 15, cellheight = 15,
           color = c(colorRampPalette(colors = c("#3373B2","white"))(length(bk)/2),
                     colorRampPalette(colors = c("white","#BD3438"))(length(bk)/2)),
           #legend_breaks=seq(-1,0,1),
           breaks=bk)

