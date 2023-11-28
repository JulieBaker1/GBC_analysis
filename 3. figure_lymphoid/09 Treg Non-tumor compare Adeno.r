
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


cd4.sub =subset(cd4,celltype == 'CD4T_C1_FOXP3')
df= data.frame(patient = cd4.sub$orig.ident)
df$group =  case_when(df$patient %in%c(XGC$NewSample.ID,CC$NewSample.ID,HG$NewSample.ID,LG$NewSample.ID)~ 'Normal',
                      df$patient %in%adeno$NewSample.ID~ 'Tumor')

cd4.sub$group = df$group

Idents(cd4.sub) = cd4.sub$group
table(cd4.sub$group)
marker.treg = FindAllMarkers(cd4.sub,only.pos = TRUE)  
marker.treg = filter(marker.treg,p_val_adj<0.05)


library(dplyr)
library(reshape2)
marker.plot = marker.treg[,c('avg_log2FC','cluster','gene')]
marker.plot = marker.plot[order(marker.plot$avg_log2FC,decreasing = TRUE),]
marker.plot = marker.plot[c(1:75),]
marker.plot$gene = factor(marker.plot$gene,levels = marker.plot$gene)
marker.plot$avg_log2FC=ifelse(marker.plot$cluster=='Normal',-marker.plot$avg_log2FC,marker.plot$avg_log2FC)
a=ggplot(marker.plot,aes(gene,avg_log2FC,fill=cluster))+
  geom_bar(stat = 'identity',width = 0.4)+
  theme(axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = "black"))+
  scale_fill_manual(values = c('#6A91B7','#D18796'))



immune.check=c("PDCD1","CD274", "CTLA4","LAG3","PDCD1LG2","BTLA", "HAVCR2" ,"TIGIT","VSIR" ,   
               "C10orf54")
t.activate = c('CD8A','CD8B','CD3E','CD3G','GZMA','GZMB','PRF1','IFNG','TBX21',
               'EOMES','IRF1','IRF8','STAT1','STAT4')
naive = c("CCR7","TCF7","SELL","IL7R" ,"IL2RG",'LEF1')

func.gene = list()
func.gene[['immune.check']]=immune.check
func.gene[['t.activate']]=t.activate
func.gene[['naive']]=naive

sign.name = names(func.gene)
for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(cd4.sub))
  a=list(a)
  cd4.sub<- AddModuleScore(cd4.sub,features = a,
                           ctrl = 100,
                           name = i)
}

df$naive = cd4.sub$naive1
df$immune.check = cd4.sub$immune.check1
df$activated = cd4.sub$t.activate1
df = df[!is.na(df$group),]

group_medians <- df %>%
  group_by(group) %>%
  summarise(median_value = median(naive))

p1=ggplot(df, aes(x = group, y = naive,fill=group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  stat_compare_means(aes(label="p.signif"),
                     method = "wilcox.test",label = "p.signif")+
  ggtitle('Naive Score')+
  scale_fill_manual(values = c( '#6A91B7','#D18796')) +
  theme_pubr()+
  geom_segment(data = group_medians, aes(x = 'Normal', xend = 'Tumor', y = median_value, yend = lead(median_value)),
               color = "black", linetype = "dashed", size = 0.3) +
  geom_point(aes(x=group,y=median_value),data=group_medians)

group_medians <- df %>%
  group_by(group) %>%
  summarise(median_value = median(immune.check))

p2=ggplot(df, aes(x = group, y = immune.check,fill=group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  stat_compare_means(aes(label="p.signif"),
                     method = "wilcox.test",label = "p.signif")+
  ggtitle('Immune Checkpoint Score')+
  scale_fill_manual(values = c( '#6A91B7','#D18796')) +
  theme_pubr()+
  geom_segment(data = group_medians, aes(x = 'Normal', xend = 'Tumor', y = median_value, yend = lead(median_value)),
               color = "black", linetype = "dashed", size =0.5) +
  geom_point(aes(x=group,y=median_value),data=group_medians)

group_medians <- df %>%
  group_by(group) %>%
  summarise(median_value = median(activated))

p3=ggplot(df, aes(x = group, y = activated,fill=group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  labs(title = "Violin Plot", x = "Groups", y = "Value")+
  stat_compare_means(aes(label="p.signif"),
                     method = "wilcox.test",label = "p.signif")+
  ggtitle('Activated Score')+
  scale_fill_manual(values = c( '#6A91B7','#D18796')) +
  theme_pubr()+
  geom_segment(data = group_medians, aes(x = 'Normal', xend = 'Tumor', y = median_value, yend = lead(median_value)),
               color = "black", linetype = "dashed", size = 0.3) +
  geom_point(aes(x=group,y=median_value),data=group_medians)
p=p1+p2+p3




