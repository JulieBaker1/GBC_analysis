
library(reshape2)
total.count.info=readRDS('...\\02 Lymphoid Celltype Normalized.rds')
celltype=colnames(total.count.info)[1:64]

p.sample=filter(total.count.info,Tumors.for.scRNA.seq.short=='P',histological.type.short=='adeno')
p.sample$NewSample.ID

ln.sample=filter(total.count.info,Tumors.for.scRNA.seq.short=='LN')
ln.sample$NewSample.ID

li.sample=filter(total.count.info,Tumors.for.scRNA.seq.short=='LI')
li.sample$NewSample.ID

lm.sample=filter(total.count.info,Tumors.for.scRNA.seq.short=='LM')
lm.sample$NewSample.ID

intersect(substr(p.sample$NewSample.ID,1,7),substr(ln.sample$NewSample.ID,1,7))

intersect(substr(ln.sample$NewSample.ID,1,7),substr(lm.sample$NewSample.ID,1,7))

p.name=paste0(intersect(substr(p.sample$NewSample.ID,1,7),substr(ln.sample$NewSample.ID,1,7)),
              '_P')
ln.name=paste0(intersect(substr(p.sample$NewSample.ID,1,7),substr(ln.sample$NewSample.ID,1,7)),
               '_LN')

rownames(total.count.info)=total.count.info$NewSample.ID
p.sample=total.count.info[p.name,]
ln.sample=total.count.info[ln.name,]

p.sample=p.sample[,celltype]
ln.sample=ln.sample[,celltype]

p.sample$group='p'
p.sample2=melt(p.sample)

ln.sample$group='ln'
ln.sample2=melt(ln.sample)


total.count.info=readRDS('...\\01 Lymphoid Celltype Count.rds')
celltype = colnames(total.count.info)[1:10]
total.count.info = total.count.info[,celltype]

normalize_to1 <- function(x){
  (x/sum(x))*100
} 

total.count.info = apply(total.count.info,1,normalize_to1)
total.count.info = as.data.frame(t(total.count.info))

p.sample=total.count.info[p.name,]
ln.sample=total.count.info[ln.name,]

p.sample=p.sample[,celltype]
ln.sample=ln.sample[,celltype]

p.sample$group='p'
p.sample2=melt(p.sample)

ln.sample$group='ln'
ln.sample2=melt(ln.sample)

plot.list=list()
for (i in c(1:10)){
  test1=subset(p.sample2,variable == celltype[i])
  test2=subset(ln.sample2,variable == celltype[i])
  test=cbind(test1,test2)
  colnames(test)=c( "group1","variable1","value1","group2" ,"variable2","value2")
  
  library(ggpubr)
  a=ggscatter(test, x = 'value1', y = 'value2',
              fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
              add = "reg.line",# Add regressin line
              add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = T, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"))+xlab('P')+ylab('LN')+
    ggtitle(celltype[i])
  plot.list[[i]]=a
}

test1=subset(p.sample2,variable == celltype[3])
test2=subset(ln.sample2,variable == celltype[3])
test=cbind(test1,test2)
colnames(test)=c( "group1","variable1","value1","group2" ,"variable2","value2")

library(ggpubr)
a=ggscatter(test, x = 'value1', y = 'value2',
            fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
            add = "reg.line",# Add regressin line
            add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = T, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"))+xlab('P')+ylab('LN')+
  ggtitle('CD4T_C0_CCR7')


### Compare the P samples with LN samples 

total.count.info=readRDS('...\\01 Lymphoid Celltype Count.rds')
celltype = colnames(total.count.info)[1:10]
total.count.info = total.count.info[,celltype]

normalize_to1 <- function(x){
  (x/sum(x))*100
} 

total.count.info = apply(total.count.info,1,normalize_to1)
total.count.info = as.data.frame(t(total.count.info))

p.sample=total.count.info[p.name,]
ln.sample=total.count.info[ln.name,]

p.sample=p.sample[,celltype]
ln.sample=ln.sample[,celltype]

p.sample$group='p'
p.sample2=melt(p.sample)

ln.sample$group='ln'
ln.sample2=melt(ln.sample)

cd4= readRDS('...\\cd4_newname.rds')
cd4 = subset(cd4,seurat_clusters ==2)
cd4.p = subset(cd4,orig.ident %in% rownames(p.sample))
cd4.p$group = 'P'
cd4.ln = subset(cd4,orig.ident %in% rownames(ln.sample))
cd4.ln$group = 'LN'
cd4 = merge(cd4.p,cd4.ln)
library(Seurat)
Idents(cd4) = cd4$group
marker = FindAllMarkers(cd4,only.pos = TRUE) 

library(ggrepel)
Dat = marker
Dat$Gene = rownames(Dat)
Dat$avg_log2FC = ifelse(Dat$cluster == 'LN',-Dat$avg_log2FC,Dat$avg_log2FC)
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) > 0, ifelse(Dat$avg_log2FC > 0 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))

a=ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  geom_text_repel(
    data = Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>0,],
    aes(label = Gene),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  theme_bw()+
  theme(
    legend.title = element_blank()
  )+
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
  ggtitle("CD4T_C2_CCR7")







