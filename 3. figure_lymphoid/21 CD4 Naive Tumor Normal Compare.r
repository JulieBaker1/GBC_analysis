#### Tumor Normal Compared

library(clusterProfiler)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
cd4= readRDS('...\\cd4_newname.rds')

df.cd4 = data.frame(patient=cd4$orig.ident,celltype=cd4$celltype)
df.cd4$value = 1
df.cd4=dcast(df.cd4,patient~celltype,fun.aggregate=sum)
rownames(df.cd4)=df.cd4$patient
df.cd4 = df.cd4[,-1]

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
df.cd4$NewSample.ID = rownames(df.cd4)
total.count.info=df.cd4 %>% full_join(info.patient,by='NewSample.ID')

HG=filter(total.count.info,histological.type.short=='HG')
HG=HG[,c(11,3)]
HG$group='HG'

LG=filter(total.count.info,histological.type.short=='LG')
LG=LG[,c(11,3)]
LG$group='LG'

XGC=filter(total.count.info,Tumors.for.scRNA.seq.short=='XGC')
XGC=XGC[,c(11,3)]
XGC$group='XGC'

CC=filter(total.count.info,Tumors.for.scRNA.seq.short=='CC')
CC=CC[,c(11,3)]
CC$group='CC'

adeno=filter(total.count.info,histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
adeno=adeno[,c(11,3)]
adeno$group='adeno_p'

cd4.sub =subset(cd4,celltype == 'CD4T_C2_CCR7')
df= data.frame(patient = cd4.sub$orig.ident)
df$group =  case_when(df$patient %in%c(XGC$NewSample.ID,CC$NewSample.ID,HG$NewSample.ID,LG$NewSample.ID)~ 'Normal',
                      df$patient %in%adeno$NewSample.ID~ 'Tumor')

cd4.sub$group = df$group

Idents(cd4.sub) = cd4.sub$group
table(cd4.sub$group)
marker.naive = FindAllMarkers(cd4.sub,only.pos = TRUE)  
marker.naive = filter(marker.naive,p_val_adj<0.05)
marker.naive$enterz <- mapIds(org.Hs.eg.db, keys = marker.naive$gene, keytype = "SYMBOL", column="ENTREZID")
enrich <- compareCluster(enterz ~cluster, data=marker.naive,
                         fun="enrichGO",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP" ,
                         pAdjustMethod = "BH",
                         readabl = TRUE)
dotplot(enrich,showCategory=10)

result=enrich@compareClusterResult
interest.pathway=c('cellular response to heat',
                   'response to temperature stimulus',
                   'positive regulation of interleukin-10 production',
                   'response to heat',
                   'interleukin-10 production',
                   'cellular divalent inorganic cation homeostasis',
                   'regulation of sterol transport',
                   'regulation of cholesterol transport',
                   'cholesterol efflux')

enrich2=filter(result,Description %in% interest.pathway)
enrich2.df <- data.frame(enrich2) 
a=enrich2.df$GeneRatio
a <- data.frame(name=a)
b <- apply(a,1,function(x) eval(parse(text=x)))
enrich2.df$GeneRatio=b
enrich2.df$Description = factor(enrich2.df$Description,levels = enrich2.df$Description)

a=ggplot(enrich2.df,aes(x=cluster,y=Description,size=GeneRatio,color=pvalue))+
  geom_point()+
  scale_colour_gradient(high= "#408591",low="#C6736E")+
  theme(axis.text.x = element_text(size=10,colour = 'black',face='bold',angle = 30),
        axis.text.y = element_text(size=10))

df.expression = FetchData(cd4.sub,vars = c('AREG','group'))
df.expression=df.expression[!is.na(df.expression$group),]
ggplot(df.expression, aes(x = group, y = AREG,fill=group))+
  geom_jitter(width=0.3,alpha=0.3,size=0.1)+
  geom_violin(width=1,position=position_dodge(width=3))




#### CC XGC HG LG Adeno
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
cd4= readRDS('...\\cd4_newname.rds')

df.cd4 = data.frame(patient=cd4$orig.ident,celltype=cd4$celltype)
df.cd4$value = 1
df.cd4=dcast(df.cd4,patient~celltype,fun.aggregate=sum)
rownames(df.cd4)=df.cd4$patient
df.cd4 = df.cd4[,-1]

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
df.cd4$NewSample.ID = rownames(df.cd4)
total.count.info=df.cd4 %>% full_join(info.patient,by='NewSample.ID')

table(total.count.info$histological.type.short)

HG=filter(total.count.info,histological.type.short=='HG')
HG=HG[,c(11,2)]
HG$group='HG'

LG=filter(total.count.info,histological.type.short=='LG')
LG=LG[,c(11,2)]
LG$group='LG'

XGC=filter(total.count.info,Tumors.for.scRNA.seq.short=='XGC')
XGC=XGC[,c(11,2)]
XGC$group='XGC'

CC=filter(total.count.info,Tumors.for.scRNA.seq.short=='CC')
CC=CC[,c(11,2)]
CC$group='CC'

adeno=filter(total.count.info,histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
adeno=adeno[,c(11,2)]
adeno$group='adeno_p'

cd4.sub =subset(cd4,celltype == 'CD4T_C2_CCR7')
df= data.frame(patient = cd4.sub$orig.ident)
df$group =  case_when(df$patient %in%c(XGC$NewSample.ID,CC$NewSample.ID,HG$NewSample.ID,LG$NewSample.ID)~ 'Normal',
                      df$patient %in%adeno$NewSample.ID~ 'Tumor')

cd4.sub$group = case_when(cd4.sub$orig.ident %in% XGC$NewSample.ID~'XGC',
                          cd4.sub$orig.ident %in% HG$NewSample.ID~'HG',
                          cd4.sub$orig.ident %in% LG$NewSample.ID~'LG',
                          cd4.sub$orig.ident %in% adeno$NewSample.ID~'adeno')

Idents(cd4.sub) = cd4.sub$group
table(cd4.sub$group)
marker = FindAllMarkers(cd4.sub,only.pos = TRUE) 

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
interest.pathway=c('positive regulation of interleukin-10 production',
                   'interleukin-10 production',
                   'positive regulation of apoptotic signaling pathway',
                   'protein refolding',
                   'protein folding',
                   'response to reactive oxygen species',
                   'positive regulation of I-kappaB kinase/NF-kappaB signaling',
                   'I-kappaB kinase/NF-kappaB signaling',
                   'humoral immune response',
                   'antimicrobial humoral response'
                   )
enrich2=filter(result,Description %in% interest.pathway)

enrich2.df <- data.frame(enrich2) 
a=enrich2.df$GeneRatio
a <- data.frame(name=a)
b <- apply(a,1,function(x) eval(parse(text=x)))
enrich2.df$GeneRatio=b

enrich2.df$Description = factor(enrich2.df$Description,levels = interest.pathway)
enrich2.df$cluster = factor(enrich2.df$cluster,levels = c('LG','HG','XGC','adeno'))

table(enrich2.df$Description)
a=ggplot(enrich2.df,aes(x=cluster,y=Description,size=GeneRatio,color=pvalue))+
  geom_point()+
  scale_colour_gradient(high= "#408591",low="#C6736E")+
  theme(axis.text.x = element_text(size=10,colour = 'black',face='bold',angle = 30),
        axis.text.y = element_text(size=10))+
  theme_bw()




