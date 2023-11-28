
####Acquire the normalized data(celltype percentage)####

library("survival")
library("survminer")
library('ggplot2')
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggprism)
library(ggsci)

cd8= readRDS('...\\cd8_newname.rds')
cd4= readRDS('...\\cd4_newname.rds')
nk = readRDS('...\\nk_newname.rds')
b = readRDS('...\\b_newname.rds')
pc = readRDS('...\\pc_20230312_addname.rds')

cd8 = subset(cd8,seurat_clusters %in% c(0,1,2,3,6,8,9,11,12))
cd4=subset(cd4,seurat_clusters %in% c(0,1,2,3,4,5,7,8))
nk=subset(nk,seurat_clusters %in% c(0,1,2,3,4,6,8,12))
b = subset(b,seurat_clusters %in% c(0,1,3,5,8))
pc = subset(pc,seurat_clusters %in% c(0,1))


df.cd8 = data.frame(patient=cd8$orig.ident,celltype=cd8$celltype)
df.cd8$value = 1
df.cd8=dcast(df.cd8,patient~celltype,fun.aggregate=sum)
rownames(df.cd8)=df.cd8$patient
df.cd8 = df.cd8[,-1]

df.cd4 = data.frame(patient=cd4$orig.ident,celltype=cd4$celltype)
df.cd4$value = 1
df.cd4=dcast(df.cd4,patient~celltype,fun.aggregate=sum)
rownames(df.cd4)=df.cd4$patient
df.cd4 = df.cd4[,-1]

df.nk = data.frame(patient=nk$orig.ident,celltype=nk$celltype)
df.nk$value = 1
df.nk=dcast(df.nk,patient~celltype,fun.aggregate=sum)
rownames(df.nk)=df.nk$patient
df.nk = df.nk[,-1]

df.b = data.frame(patient=b$orig.ident,celltype=b$celltype)
df.b$value = 1
df.b=dcast(df.b,patient~celltype,fun.aggregate=sum)
rownames(df.b)=df.b$patient
df.b = df.b[,-1]

df.pc = data.frame(patient=pc$orig.ident,celltype=pc$celltype)
df.pc$value = 1
df.pc=dcast(df.pc,patient~celltype,fun.aggregate=sum)
rownames(df.pc)=df.pc$patient
df.pc = df.pc[,-1]

df.cd4 = df.cd4[,-9]
df.cd8 = df.cd8[,-10]
df.nk = df.nk[,-9]
df.b = df.b[,-6]

df.cd8$NewSample.ID=rownames(df.cd8)
df.cd4$NewSample.ID=rownames(df.cd4)
df.nk$NewSample.ID=rownames(df.nk)
df.b$NewSample.ID=rownames(df.b)
df.pc$NewSample.ID=rownames(df.pc)

total.count = df.cd4 %>%left_join(df.cd8)%>%left_join(df.nk)%>%
  left_join(df.b)%>%left_join(df.pc)
total.count[is.na(total.count)]=0
rownames(total.count) = total.count$NewSample.ID
total.count = total.count[,-9]
saveRDS(total.count,'...\\count_celltype.rds')

norm=function(x){
  ((x/sum(x))*100)
} 

norm.cd4 = apply(df.cd4, 1 , norm)
norm.cd8 = apply(df.cd8, 1 , norm)
norm.nk = apply(df.nk, 1 , norm)
norm.b = apply(df.b, 1 , norm)

norm.cd4 = as.data.frame(t(norm.cd4))
norm.cd8 = as.data.frame(t(norm.cd8))
norm.nk = as.data.frame(t(norm.nk))
norm.b = as.data.frame(t(norm.b))

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)

norm.cd8$NewSample.ID=rownames(norm.cd8)
norm.cd4$NewSample.ID=rownames(norm.cd4)
norm.nk$NewSample.ID=rownames(norm.nk)
norm.b$NewSample.ID=rownames(norm.b)
#df.pc$NewSample.ID=rownames(df.pc)

norm.total = norm.cd8 %>% left_join(norm.cd4) %>% left_join(norm.nk) %>%
  left_join(norm.b)
rownames(norm.total) = rownames(norm.cd8)
saveRDS(norm.total,'...\\norm_celltype.rds')


####Prognosis analysis####
norm.total =readRDS('...\\02 Lymphoid Celltype Normalized.rds')
info.yuhou=read.table('...\\ls_20220503_yuhou.csv',sep=',',header = TRUE,row.names = 1)

library(dplyr)
yuhou.final.total= norm.total %>% left_join(info.yuhou,by='NewSample.ID')
cellname=colnames(yuhou.final.total)[c(1:64)]
yuhou=c("DFS","OS_month","event")

yuhou.final.total= yuhou.final.total %>% filter(histological.type.short.x=='adeno',Tumors.for.scRNA.seq.short.x=='P')
yuhou.final.total$state = ifelse(yuhou.final.total$event=='dead',"1",'0')

cutoff_yuhou=yuhou.final.total[,c(cellname,yuhou,'state')]
cutoff_yuhou2=cutoff_yuhou[complete.cases(cutoff_yuhou[,'state']),]
cutoff_yuhou2$state=as.numeric(cutoff_yuhou2$state)

##total
plot.list = list()
for (i in c(1:64)){
  cellname2=cellname[i]
  res.cut <- surv_cutpoint(cutoff_yuhou2, time = "OS_month", event = "state",
                           variables =cellname2)
  res.cat <- surv_categorize(res.cut)
  
  yuhou.name=colnames(res.cat)[1:2]
  cell=cellname2
  res.cat2=res.cat[,c(cell,yuhou.name)]
  colnames(res.cat2)[1]='group'
  fit=survfit(Surv(OS_month,state)~group,data=res.cat2)
  a=ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c('#d62728','#1f77b4')
  )+ggtitle(cell)
  plot.list[[i]]=a
}

pdf(file = '...\\total_ct_cutoff.pdf',
    height = 4.5,width = 4)
for (i in c(1:30)){
  print(plot.list[[i]])
}
dev.off()
