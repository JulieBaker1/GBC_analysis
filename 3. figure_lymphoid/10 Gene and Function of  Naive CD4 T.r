### Find DEGs of naive CD4 T (P v.s. Pmet)
library(dplyr)
library(reshape2)
library(Seurat)
library(ggplot2)
library(ggpubr)

cd4= readRDS('...\\cd4_newname.rds')
cd4 = subset(cd4,seurat_clusters ==2)

norm.total = readRDS('...\\02 Lymphoid Celltype Normalized.rds')
colnames(norm.total)[42]='NK_C10_HLA_DRA'
P_TOTAL= norm.total %>% filter(histological.type.short=='adeno',metastasis.type=='P' |metastasis.type=='P_LI')
Metastasis= norm.total %>% filter(histological.type.short=='adeno',metastasis.type=='P_LM'|metastasis.type=='P_LN')
sample.p = P_TOTAL$NewSample.ID
Metastasis.p = Metastasis$NewSample.ID

cd4$group = case_when(cd4$orig.ident %in% sample.p~'P',
                      cd4$orig.ident %in% Metastasis.p~'P(Met)')                                                                                                 

cd4 = cd4[,!is.na(cd4$group)]
dim(cd4)
Idents(cd4) = cd4$group

marker = FindAllMarkers(cd4,only.pos = TRUE)
marker$label = ifelse(marker$p_val_adj<0.01,'p_val_adj<0.01','p_val_adj>=0.01')
colnames(marker)
top10sig0 <- filter(marker,cluster=="P(Met)") %>% top_n(10,abs(avg_log2FC))
top10sig1 <- filter(marker,cluster=="P") %>% top_n(10,abs(avg_log2FC))
top10sig = rbind(top10sig0,top10sig1)
marker = marker[!rownames(marker) %in% rownames(top10sig),]

marker$cluster = factor(marker$cluster,levels = c("P",'P(Met)'))
top10sig$cluster = factor(top10sig$cluster,levels = c("P",'P(Met)'))
library(ggrepel)
a=ggplot()+
  geom_jitter(data = marker,
              aes(x = cluster, y = avg_log2FC,color=label),
              size = 0.85,
              width =0.3)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.3)+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=avg_log2FC,label=gene),size=2,
    force = 0.3,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"))+
  ylim(0,1.2)+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black"),
    axis.line.y = element_line(color = "black",
                               size = 0.5),
    axis.line.x =element_line(color = "black",
                              size =0.5),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )+scale_color_manual(values = c("#A55D51","#1F3B54"))


### Compare the function(p vs pmet) of naive cd4 t 
sign1 <- read.delim("...\\immune_sign.csv",sep=',')
colnames(sign1)
sign1 = sign1[,c(4,5,6,7,8,9)]
sign2<- read.delim("...\\nature sig.csv",sep=',')
colnames(sign2)
sign2 = sign2[,c(8,9,10,13,15,16,17)]
sign = cbind(sign1,sign2)

cytotoxic=c("EOMES","TBX21","GZMB","PRF1","FASL","GZMH" ,"GZMA" ,"IFNG"  ,"GZMK",
            "ZAP70","GNLY" ,"FASLG","NKG7",'FCGR3A','KLRC1','NKG2D','KLRC2','KLRD1',
            'STAT1','NCR1','NCR2','NCR3')
immune.check=c("PDCD1","CD274", "CTLA4","LAG3","PDCD1LG2","BTLA", "HAVCR2" ,"TIGIT","VSIR" ,   
               "C10orf54")
naive = c("CCR7","TCF7","SELL","IL7R" ,"IL2RG",'LEF1')

func.gene = list()
for (i in c(1:6)){
  func.gene[[colnames(sign1)[i]]]=sign1[,i]
}

for (i in c(1:7)){
  func.gene[[colnames(sign2)[i]]]=sign2[,i]
}

func.gene[['cytotoxic']]=cytotoxic
func.gene[['immune.check']]=immune.check
func.gene[['naive']]=naive

sign.name = names(func.gene)
for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(cd4))
  a=list(a)
  cd4<- AddModuleScore(cd4,features = a,
                       ctrl = 100,
                       name = i)
}

df =cd4@meta.data
colnames(df)
df = df[,c(82:98)]

group_medians <- df %>%
  group_by(group) %>%
  summarise_each(funs = median)

for (i in c(2:17)){
  func.name = colnames(df)[i]
  group_median = group_medians[,c(func.name,'group')]
  colnames(group_median)[1] = 'median_value'
  a=ggplot(df, aes_string(x = 'group', y = func.name,fill='group')) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
    labs(title = "Violin Plot", x = "Groups", y = "Value")+
    stat_compare_means(aes(label="p.signif"),
                       method = "wilcox.test",label = "p.signif")+
    ggtitle(func.name)+
    scale_fill_manual(values = c( '#6A91B7','#D18796')) +
    theme_pubr()+
    geom_segment(data = group_median, aes(x = 'P', xend = 'P(Met)', y = median_value, yend = lead(median_value)),
                 color = "black", linetype = "dashed", size = 0.3) +
    geom_point(aes(x=group,y=median_value),data=group_median)
  path = paste0('...\\',func.name,'.pdf')
  pdf(file = path,
      height=3.5,width=2.5)
  print(a)
  dev.off()
}




















