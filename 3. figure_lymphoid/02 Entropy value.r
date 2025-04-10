### This scirpt is used to caculate the entropy value and plot with barplot

library(reshape2)
library(dplyr)
library(seurat)
library(ggplot2)

cd8= readRDS('...\\cd8_newname.rds')
cd4= readRDS('...\\cd4_newname.rds')
nk = readRDS('...\\nk_newname.rds')
b = readRDS('....\\b_newname.rds')
pc = readRDS('...\\pc_20230312_addname.rds')

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


#function about entropy 
fn = function(x){
  return(-x*log(x))
}
softmax=function(x){
  (x/sum(x))
} 
fn2 = function(x){
  return(sum(x)/log(length(x)))
}

a=apply(df.cd8,2,softmax)
a=apply(df.cd4,2,softmax)
a=apply(df.nk,2,softmax)
a=apply(df.b,2,softmax)
a=apply(df.pc,2,softmax)

a=apply(a, 2, fn)
a[is.na(a)]=0
a=apply(a, 2, fn2)
shang.value=as.data.frame(a)
shang.value$Cluster=rownames(shang.value)
shang.value$Cluster=factor(shang.value$Cluster,levels = shang.value$Cluster)


cd8.s=shang.value
cd4.s=shang.value
nk.s=shang.value
b.s=shang.value
pc.s=shang.value

total.entropy = rbind(cd8.s,cd4.s,nk.s,b.s,pc.s)
total.entropy$group=case_when(total.entropy$a >=0.625~'High Shanno Value',
                              total.entropy$a <0.625 ~ 'Low Shanno Value')
total.entropy$color=case_when(total.entropy$a >=0.625~'darkred',
                              total.entropy$a < 0.625 ~ 'gray')

a=ggplot(total.entropy, aes(x=Cluster,y=a,fill=group)) +  
  geom_col()+
  geom_hline(aes(yintercept=0.625), size=1.5,colour="black", linetype="dashed")+
  theme(axis.line=element_line(colour='black',size=1,lineend = 'square'),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,colour = total.entropy$color,angle=60,hjust = 1),
        legend.direction = "vertical",
        legend.position = "left",
        panel.grid = element_blank(),
        #legend.box = "vertical"
  )+scale_fill_manual(values = c('#EA927F','#778899'))+
  xlab('')+ylab('Entropy')

