
#### Functional Module score####

cd8= readRDS('...\\cd8_newname.rds')
cd4= readRDS('...\\cd4_newname.rds')
nk = readRDS('...\\nk_newname.rds')


cd8 = subset(cd8,seurat_clusters %in% c(0,1,2,3,6,8,9,11,12))
cd4=subset(cd4,seurat_clusters %in% c(0,1,2,3,4,5,7,8))
nk=subset(nk,seurat_clusters %in% c(0,1,2,3,4,6,8,12))

sign <- read.delim("D:\\GBC\\01 common file\\immune_sign.csv",sep=',')
sign = sign %>% select('Proinflammatory','Glycolysis','Glucose_Deprivation','HALLMARK_TGF_BETA_SIGNALING')

cytotoxic=c("EOMES","TBX21","GZMB","PRF1","FASL","GZMH" ,"GZMA" ,"IFNG"  ,"GZMK",
            "ZAP70","GNLY" ,"FASLG","NKG7",'FCGR3A','KLRC1','NKG2D','KLRC2','KLRD1',
            'STAT1','NCR1','NCR2','NCR3')
immune.check=c("PDCD1","CD274", "CTLA4","LAG3","PDCD1LG2","BTLA", "HAVCR2" ,"TIGIT","VSIR" ,   
               "C10orf54")
t.activate = c('CD8A','CD8B','CD3E','CD3G','GZMA','GZMB','PRF1','IFNG','TBX21',
               'EOMES','IRF1','IRF8','STAT1','STAT4')
naive = c("CCR7","TCF7","SELL","IL7R" ,"IL2RG",'LEF1')
Proinflammatory = c("IL1A","IL1B","TNF","IFNG","TBX21","CCL3","CCL4","PRF1","GZMA","GZMB","GZMK","GZMH","CD8A","FASLG", "CCL2" ,
                "CCL20","IL2","IL6","IL12A","IL17A","IL23A","PTGS2","TLR4","TNF")

func.gene = list()
func.gene[['cytotoxic']]=cytotoxic
func.gene[['immune.check']]=immune.check
func.gene[['t.activate']]=t.activate
func.gene[['naive']]=naive
func.gene[[colnames(sign)[1]]]=sign$Proinflammatory


sign.name = names(func.gene)
for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(cd8))
  a=list(a)
  cd8<- AddModuleScore(cd8,features = a,
                       ctrl = 100,
                       name = i)
}

for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(cd4))
  a=list(a)
  cd4<- AddModuleScore(cd4,features = a,
                       ctrl = 100,
                       name = i)
}

for (i in sign.name){
  a=func.gene[[i]]
  a= intersect(a,rownames(nk))
  a=list(a)
  nk<- AddModuleScore(nk,features = a,
                       ctrl = 100,
                       name = i)
}

sign.name = paste0(sign.name,1)
total=cd4@meta.data
total=total[,sign.name]
total$group=cd4$celltype
total2 = total %>% group_by(group) %>% summarise_each(funs = median)
total2=as.data.frame(total2)
rownames(total2)= total2$group
total2=total2[,-1]
total2=scale(total2)
colnames(total2)=sign.name
cd4.final=total2

total=cd8@meta.data
total=total[,sign.name]
total$group=cd8$celltype
total2 = total %>% group_by(group) %>% summarise_each(funs = median)
total2=as.data.frame(total2)
rownames(total2)= total2$group
total2=total2[,-1]
total2=scale(total2)
colnames(total2)=sign.name
cd8.final=total2

total=nk@meta.data
total=total[,sign.name]
total$group=nk$celltype
total2 = total %>% group_by(group) %>% summarise_each(funs = median)
total2=as.data.frame(total2)
rownames(total2)= total2$group
total2=total2[,-1]
total2=scale(total2)
colnames(total2)=sign.name
nk.final=total2

linba.final=rbind(cd8.final,cd4.final,nk.final)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
library(pheatmap)

colnames(linba.final)=c("Cytotoxic" ,"Immune Checkpoint","T Activated",                 
                        "Naive", "Proinflammatory")
a=pheatmap(linba.final,cluster_rows = FALSE,cluster_cols = FALSE,
           gaps_col = c(5),
         gaps_row = c(9,17),
           cellwidth = 15, cellheight = 15,
           color = c(colorRampPalette(colors = c("#3071B1","white"))(length(bk)/2),
                     colorRampPalette(colors = c("white","#BB2C30"))(length(bk)/2)),
           legend_breaks=seq(-2,2,1),
           breaks=bk)

#### Immune checkpoint####

immune.check=c("PDCD1","CD274", "CTLA4","LAG3","PDCD1LG2","BTLA", "HAVCR2" ,"TIGIT","VSIR" ,   
               "C10orf54")
ic.cd8=FetchData(cd8,vars = immune.check,slot = 'data')
ic.nk =FetchData(nk,vars = immune.check,slot = 'data')
ic.cd4 =FetchData(cd4,vars = immune.check,slot = 'data')

ic.cd8$celltype = cd8$celltype
ic.nk$celltype = nk$celltype
ic.cd4$celltype= cd4$celltype

ic.cd8 = ic.cd8 %>% group_by(celltype) %>% summarise_each(funs=mean)
ic.cd4 = ic.cd4 %>% group_by(celltype) %>% summarise_each(funs=mean)
ic.nk = ic.nk %>% group_by(celltype) %>% summarise_each(funs=mean)

ic.t = rbind(ic.cd8,ic.cd4,ic.nk)
ic.t = as.data.frame(ic.t)
rownames(ic.t)=ic.t$celltype
ic.t = ic.t[,-1]
ic.t=scale(ic.t)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
a=pheatmap(ic.t,cluster_rows = FALSE,cluster_cols = FALSE,
         #gaps_col = c(5),
         gaps_row = c(9,17),
         cellwidth = 15, cellheight = 15,
         color = c(colorRampPalette(colors = c("#f0cfac","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#425a04"))(length(bk)/2)),
        #legend_breaks=seq(-1,0,1),
         breaks=bk)

