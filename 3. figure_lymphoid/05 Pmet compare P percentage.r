
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

norm.total = readRDS('...\\02 Lymphoid Celltype Normalized.rds')
colnames(norm.total)[42]='NK_C10_HLA_DRA'
celltype = colnames(norm.total)[1:64]

XGC=filter(norm.total,Tumors.for.scRNA.seq.short=='XGC')
XGC=XGC[,celltype]
XGC$group='XGC'

CC=filter(norm.total,Tumors.for.scRNA.seq.short=='CC')
CC=CC[,celltype]
CC$group='CC'

HG=filter(norm.total,histological.type.short=='HG')
HG=HG[,celltype]
HG$group='HG'

LG=filter(norm.total,histological.type.short=='LG')
LG=LG[,celltype]
LG$group='LG'

normal = rbind(XGC,CC,HG,LG)
normal$group = 'Non-Tumor'

P_TOTAL= norm.total %>% filter(histological.type.short=='adeno',metastasis.type=='P' |metastasis.type=='P_LI')
P_TOTAL=P_TOTAL[,celltype]
P_TOTAL$group = 'P'
Metastasis= norm.total %>% filter(histological.type.short=='adeno',metastasis.type=='P_LM'|metastasis.type=='P_LN')
Metastasis=Metastasis[,celltype]
Metastasis$group = 'P(Met)'

box.plot = rbind(P_TOTAL,Metastasis)

a=ggplot(box.plot,aes(x=group,y=B_C0_IGHA1,fill=group))+geom_boxplot(outlier.shape = NA,width=0.8)+
  geom_jitter(aes(fill=group,shape = group),width =0.2,size=1.5)+ 
  scale_shape_manual(values = c(15, 16))+
  theme(panel.background = element_blank(),axis.line = element_line())+   
  stat_compare_means(label="p.format")+              
  stat_boxplot(geom = "errorbar",width=0.5)+                           
  labs(title = "B_C0_IGHA1")+                      
  theme(legend.position = "bottom") +
  ylab('Percentage')+
  scale_fill_manual(values=c('#f16c23','#2b6a99'))



