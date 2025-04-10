
library(Seurat)
library(dplyr)
library(reshape2)

cd4= readRDS('...\\cd4_newname.rds')
expr = subset(cd4,seurat_clusters==1)

meta.data = data.frame(NewSample.ID = expr$orig.ident)
info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
meta.data = meta.data %>% left_join(info.patient)
sign1 <- read.delim("...\\immune_sign.csv",sep=',')
colnames(sign1)
sign1 = sign1[,c(4,5,6,7,8,9)]

sign2<- read.delim("...\\nature sig.csv",sep=',')
colnames(sign2)
sign2 = sign2[,c(4,6,7,11,13,8,9,10,13,15,16,17)]
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

for (i in c(1:12)){
  func.gene[[colnames(sign2)[i]]]=sign2[,i]
}

func.gene[['cytotoxic']]=cytotoxic
func.gene[['immune.check']]=immune.check
func.gene[['naive']]=naive

sign.name = names(func.gene)
for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(expr))
  a=list(a)
  expr<- AddModuleScore(expr,features = a,
                        ctrl = 100,
                        name = i)
}

expr$clinical.stage = meta.data$Clinical.stage
expr$histological.type.short  = meta.data$histological.type.short
expr$metastasis.type =meta.data$metastasis.type

df =expr@meta.data
colnames(df)
df = df[,c(82:105)]

df = filter(df,histological.type.short == 'adeno')
df = df[!is.na(df$metastasis.type),]

colnames(df)
df_filter = df[,c(1:22)]
group_medians <- df_filter %>%
  group_by(clinical.stage) %>%
  summarise_each(funs = median)


df_filter[which(df_filter$clinical.stage %in% c('IIA','IIB')),'clinical.stage']='II'
colnames(df_filter)
--- # Go to 138 lines

total.sub = df_filter[,c(2,22)]
colnames(total.sub)[1]='cell'
final2= total.sub %>% group_by(clinical.stage) %>% summarise(mean = median(cell),N=length(cell),
                                                             sd = sd(cell),se=sd/sqrt(N))
a=ggplot(final2,aes(clinical.stage,mean,group = 1))+
  geom_ribbon(aes(ymin=mean-se,ymax=mean+se),fill='#E6DCDB')+
  geom_line(color='#C86F76',size=1)+
  geom_point(color='#C86F76',size=2)+
  theme(panel.background = element_blank(),
        panel.grid = element_line(color="white"),
        # axis.title = element_blank(),
        #axis.ticks = element_blank(),
        axis.text.x=element_text(colour='black',size=12,angle = 45,hjust = 1))+
  ggtitle('Anergy1')+theme(axis.line=element_line(colour='black',size=1,lineend = 'square'))+
  theme (axis.text.x = element_text (colour='black', size=12), 
         axis.text.y = element_text (colour='black',size=12))+
  labs(x='Stage',y='Percentage')

for (i in c(1:21)){
  total.sub = df_filter[,c(i,22)]
  colnames(total.sub)[1]='cell'
  final2= total.sub %>% group_by(clinical.stage) %>% summarise(mean = median(cell),N=length(cell),
                                                               sd = sd(cell),se=sd/sqrt(N))
  func.name = colnames(df_filter)[i]
  
  a=ggplot(final2,aes(clinical.stage,mean,group = 1))+
    geom_ribbon(aes(ymin=mean-se,ymax=mean+se),fill='#E6DCDB')+
    geom_line(color='#C86F76',size=1)+
    geom_point(color='#C86F76',size=2)+
    theme(panel.background = element_blank(),
          panel.grid = element_line(color="white"),
          # axis.title = element_blank(),
          #axis.ticks = element_blank(),
          axis.text.x=element_text(colour='black',size=12,angle = 45,hjust = 1))+
    ggtitle(func.name)+theme(axis.line=element_line(colour='black',size=1,lineend = 'square'))+
    theme (axis.text.x = element_text (colour='black', size=12), 
           axis.text.y = element_text (colour='black',size=12))+
    labs(x='Stage',y='Percentage')


###Pheatmap

df_filter = df_filter %>% group_by(clinical.stage) %>% summarise_each(funs = median)
df_filter = as.data.frame(df_filter)
rownames(df_filter) = df_filter$clinical.stage
df_filter = df_filter[,-1]
df_filter = scale(df_filter)
df_filter = df_filter[,-15] 
df_filter = t(df_filter)
library(ComplexHeatmap)
Heatmap(df_filter,cluster_columns = FALSE)

library(dendextend)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#404398", "white", "#F18932"))
row_dend = as.dendrogram(hclust(dist(df_filter)))
row_dend = color_branches(row_dend, k = 5) 
a=Heatmap(df_filter,cluster_columns = FALSE,  cluster_rows = row_dend,
        col = col_fun,name = 'Scale median',
        row_split = 5)




                            
                            
                            
                            
                            
