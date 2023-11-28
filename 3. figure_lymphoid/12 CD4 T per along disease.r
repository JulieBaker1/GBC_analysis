
library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(dplyr)
library(ggplot2)

total.count.info=readRDS('...\\01 Lymphoid Celltype Count.rds')
celltype = colnames(total.count.info)[c(1:10)]
total.count.info = total.count.info[,celltype]

normalize_to1 <- function(x){
  (x/sum(x))*100
} 

total.count.info = apply(total.count.info,1,normalize_to1)
total.count.info = as.data.frame(t(total.count.info))
total.count.info$NewSample.ID = rownames(total.count.info)

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
total.count.info = total.count.info %>% left_join(info.patient)

tumor=total.count.info%>% filter(histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
colnames(tumor)
celltype = c("CD4T_C0_ANXA1","CD4T_C1_FOXP3","CD4T_C2_CCR7","CD4T_C3_CXCL13","CD4T_C4_Unassign","CD4T_C5_MT1X",    
             "CD4T_C7_ISG15","CD4T_C8_MKI67")
tumor = tumor[,c(celltype,'Clinical.stage')]
tumor[tumor$Clinical.stage=='IIA','Clinical.stage']='II'
tumor[tumor$Clinical.stage=='IIB','Clinical.stage']='II'
mean.df2= tumor %>% group_by(Clinical.stage) %>% summarise_each(fun=mean)
mean.df2 = as.data.frame(mean.df2)
rownames(mean.df2) = mean.df2$Clinical.stage
mean.df2 = mean.df2[,-1]

df.cor = cor(mean.df2)
col_fun = colorRamp2(c(-1, 0, 1), c("#404398", "white", "#F18932"))
a=Heatmap(df.cor,col = col_fun)
pdf(file = '...\\Corrlation.pdf',
    height= 4,width=5)
print(a)
dev.off()

total.count.info=readRDS('...\\01 Lymphoid Celltype Count.rds')
celltype = colnames(total.count.info)[c(1:10)]
total.count.info = total.count.info[,celltype]

normalize_to1 <- function(x){
  (x/sum(x))*100
} 

total.count.info = apply(total.count.info,1,normalize_to1)
total.count.info = as.data.frame(t(total.count.info))
total.count.info$NewSample.ID = rownames(total.count.info)

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
total.count.info = total.count.info %>% left_join(info.patient)

tumor=total.count.info%>% filter(histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
tumor = tumor[,c(celltype,'Clinical.stage')]
tumor[tumor$Clinical.stage=='IIA','Clinical.stage']='II'
tumor[tumor$Clinical.stage=='IIB','Clinical.stage']='II'
mean.df2= tumor %>% group_by(Clinical.stage) %>% summarise_each(fun=mean)
mean.df2 = as.data.frame(mean.df2)
rownames(mean.df2) = mean.df2$Clinical.stage
mean.df2 = mean.df2[,-1]

sum(mean.df2[1,])
mean.df2$other = mean.df2$CD4T_C6_AOPE+mean.df2$CD4T_C9_COL1A1
mean.df2 = mean.df2[,-c(7,10)]
mean.df2$stage = rownames(mean.df2)
mean.df3 = melt(mean.df2)

library(ggplot2)
mean.df3$stage = rep(c(1:6),9)
a=ggplot(mean.df3,aes(x=stage,y=value,fill=variable))+
  geom_area()+
  scale_fill_manual(values = c('#99CCB3','#CC99CC','#CCCC99',
                               '#996666','#809966','#CAAFCA',
                               '#FFFF8F','#3399CC','gray'))


total.count.info=readRDS('...\\01 Lymphoid Celltype Count.rds')
celltype = colnames(total.count.info)[c(1:10)]
total.count.info = total.count.info[,celltype]

normalize_to1 <- function(x){
  (x/sum(x))*100
} 

total.count.info = apply(total.count.info,1,normalize_to1)
total.count.info = as.data.frame(t(total.count.info))
total.count.info$NewSample.ID = rownames(total.count.info)

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
total.count.info = total.count.info %>% left_join(info.patient)

tumor=total.count.info%>% filter(histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
colnames(tumor)
celltype = c("CD4T_C0_ANXA1","CD4T_C1_FOXP3","CD4T_C2_CCR7","CD4T_C3_CXCL13","CD4T_C4_Unassign","CD4T_C5_MT1X",    
             "CD4T_C7_ISG15","CD4T_C8_MKI67")
tumor = tumor[,c(celltype,'Clinical.stage')]
tumor[tumor$Clinical.stage=='IIA','Clinical.stage']='II'
tumor[tumor$Clinical.stage=='IIB','Clinical.stage']='II'
mean.df2= tumor %>% group_by(Clinical.stage) %>% summarise_each(fun=mean)
mean.df2 = as.data.frame(mean.df2)
rownames(mean.df2) = mean.df2$Clinical.stage
mean.df2 = mean.df2[,-1]

df.cor = cor(mean.df2)
col_fun = colorRamp2(c(-1, 0, 1), c("#404398", "white", "#F18932"))
a=Heatmap(df.cor,col = col_fun)

total.count.info=readRDS('...\\01 Lymphoid Celltype Count.rds')
celltype = colnames(total.count.info)[c(1:10)]
total.count.info = total.count.info[,celltype]

normalize_to1 <- function(x){
  (x/sum(x))*100
} 

total.count.info = apply(total.count.info,1,normalize_to1)
total.count.info = as.data.frame(t(total.count.info))
total.count.info$NewSample.ID = rownames(total.count.info)

info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
total.count.info = total.count.info %>% left_join(info.patient)

tumor=total.count.info%>% filter(histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
tumor = tumor[,c(celltype,'Clinical.stage')]
tumor[tumor$Clinical.stage=='IIA','Clinical.stage']='II'
tumor[tumor$Clinical.stage=='IIB','Clinical.stage']='II'
mean.df2= tumor %>% group_by(Clinical.stage) %>% summarise_each(fun=mean)
mean.df2 = as.data.frame(mean.df2)
rownames(mean.df2) = mean.df2$Clinical.stage
mean.df2 = mean.df2[,-1]

sum(mean.df2[1,])
mean.df2$other = mean.df2$CD4T_C6_AOPE+mean.df2$CD4T_C9_COL1A1
mean.df2 = mean.df2[,-c(7,10)]
mean.df2$stage = rownames(mean.df2)

mean.df3 = melt(mean.df2)
mean.df3$stage = rep(c(1:6),9)

a=ggplot(mean.df3,aes(x=stage,y=value,fill=variable))+
  geom_area()+
  scale_fill_manual(values = c('#99CCB3','#CC99CC','#CCCC99',
                               '#996666','#809966','#CAAFCA',
                               '#FFFF8F','#3399CC','gray'))


### Liner Regression
test = tumor
library(ggpubr)
colnames(tumor)
tumor$group = case_when(test$Clinical.stage == 'I'~1,
                        test$Clinical.stage == 'II'~2,
                        test$Clinical.stage == 'IIIA'~3,
                        test$Clinical.stage == 'IIIB'~4,
                        test$Clinical.stage == 'IVA'~5,
                        test$Clinical.stage == 'IVB'~6)
model <- lm(tumor$group~ tumor$CD4T_C2_CCR7)
model_summary <- summary(model)
p_value <- coef(model_summary)[2, 4]
k = coef(model_summary)[2,1]
p_value
k

ggplot(data = tumor, aes(x = group, y =CD4T_C2_CCR7)) +
  geom_jitter(size=0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue",size=2) +
  labs(x = "Clinical Stage", y = 'Percentage', 
       title = 'CD4T_C2_CCR7',
       subtitle = paste("p-value =", round(p_value, 3), "  Coefficient =", round(k, 3)))











