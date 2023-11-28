
library(Seurat)
library(ggsci)
library(ggplot2)

cd8= readRDS('...\\cd8_newname.rds')
expr = subset(cd8,seurat_clusters==1)

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
df=data.frame(NewSample.ID =  expr$orig.ident)
df$cellid = rownames(df)
df = df %>% left_join(info.patient)
df = filter(df,Tumors.for.scRNA.seq.short=='P',histological.type.short=='adeno')
expr = expr[,df$cellid]

meta.data = data.frame(NewSample.ID = expr$orig.ident)
info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
meta.data = meta.data %>% left_join(info.patient)

expr$clinical.stage = meta.data$Clinical.stage
expr$metastasis.type =meta.data$metastasis.type
expr$Tumors.for.scRNA.seq.short =meta.data$Tumors.for.scRNA.seq.short
expr$histological.type.short =meta.data$histological.type.short
table(expr$clinical.stage)

immune.check=c("PDCD1","CD274", "CTLA4","LAG3","PDCD1LG2","BTLA", "HAVCR2" ,"TIGIT","VSIR" ,   
               "C10orf54")

cytotoxic=c("EOMES","TBX21","GZMB","PRF1","FASL","GZMH" ,"GZMA" ,"IFNG"  ,"GZMK",
            "ZAP70","GNLY" ,"FASLG","NKG7",'FCGR3A','KLRC1','NKG2D','KLRC2','KLRD1',
            'STAT1','NCR1','NCR2','NCR3')

func.gene = list()
func.gene[['immune.check']]=immune.check
func.gene[['cytotoxic']]=cytotoxic

sign.name = names(func.gene)
for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(cd8))
  a=list(a)
  expr<- AddModuleScore(expr,features = a,
                        ctrl = 100,
                        name = i)
}

df = FetchData(expr,vars = c('immune.check1','clinical.stage','metastasis.type','orig.ident'))
df = df[!is.na(df$metastasis.type),]
df = df[,c(1,2)]
colnames(df)=c('score','Clinical.Stage')

df2= FetchData(expr,vars = c('immune.check1','Tumors.for.scRNA.seq.short','histological.type.short'))

total = df
total = subset(total,Clinical.Stage != '0')
total= total %>% group_by(Clinical.Stage) %>% summarise(mean = median(score),N=length(score),
                                                        sd = sd(score),se=sd/sqrt(N))
total$group='Immune.check'
final1 = total
table(final1$Clinical.Stage)

df = FetchData(expr,vars = c('cytotoxic1','clinical.stage','metastasis.type','orig.ident'))
df = df[!is.na(df$metastasis.type),]
df = df[,c(1,2)]
colnames(df)=c('score','Clinical.Stage')

total = df
total = subset(total,Clinical.Stage != '0')
total= total %>% group_by(Clinical.Stage) %>% summarise(mean = median(score),N=length(score),
                                                        sd = sd(score),se=sd/sqrt(N))
total$group='Cytotoxic'
final2 = total

df = rbind(final1,final2)
table(df$Clinical.Stage)
df$Clinical.Stage = factor(df$Clinical.Stage,levels =c("I","IIA","IIB","IIIA","IIIB","IVA","IVB"))

a=ggplot(df,aes(Clinical.Stage,mean,group = group,color=group))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se,color = group), width = 0.1)+
  geom_line(size=0.3)+
  geom_point(size=0.1)+
  scale_color_aaas()+
  theme(panel.background = element_blank(),
        panel.grid = element_line(color="white"),
        axis.text.x=element_text(colour='black',size=12,angle = 45,hjust = 1))+
  ggtitle('')+theme(axis.line=element_line(colour='black',size=1,lineend = 'square'))+
  theme (axis.text.x = element_text (colour='black', size=12, angle=45), 
         axis.text.y = element_text (colour='black',size=12, angle=45))+
  labs(x='Stage',y='Score')+
  ylim(-0.35,0.55)


