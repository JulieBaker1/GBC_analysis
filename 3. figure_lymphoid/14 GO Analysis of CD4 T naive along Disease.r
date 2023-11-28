
###This script is used to analysis the trend of all genes of cd4 t naive along clinical stage

library(Seurat)
library(reshape2)
library(dplyr)
cd4=readRDS('...\\cd4_newname.rds')
cd4 = subset(cd4,seurat_clusters %in% c(2))
cd4 = FindVariableFeatures(cd4)
variable.gene = VariableFeatures(cd4)

meta.data = data.frame(NewSample.ID = cd4$orig.ident)
info.patient=read.csv('...\\xzh_220109_final_wyh_220426.csv',row.names = 1)
meta.data = meta.data %>% left_join(info.patient)
meta.data = meta.data[,c('Clinical.stage','metastasis.type','Tumors.for.scRNA.seq.short','histological.type.short')]

expr.data = cd4@assays$RNA@counts
expr.data = as.data.frame(expr.data)
expr.data = as.data.frame(t(expr.data))
expr.data = expr.data[,variable.gene]
expr.total = cbind(expr.data,meta.data)

expr.filtered = expr.total%>% filter(histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
expr.filtered[expr.filtered$Clinical.stage=='IIA','Clinical.stage']='II'
expr.filtered[expr.filtered$Clinical.stage=='IIB','Clinical.stage']='II'
table(expr.filtered$Clinical.stage)

clinical_stage <- expr.filtered$Clinical.stage
#gene_expression = expr.filtered
gene_expression <- expr.filtered[,c(1:2000)]
threshold = 0.05

expr.filtered$Clinical.stage <- as.numeric(factor(expr.filtered$Clinical.stage, levels = c("I", "II", "IIIA",
                                                                                           "IIIB","IVA",'IVB')))
clinical_stage <- expr.filtered$Clinical.stage
gene_column = colnames(gene_expression)[1]

significant_genes <- c()
p_value_final = c()
k_final = c()
for (gene_column in colnames(gene_expression)) {
  
  model <- lm(gene_expression[, gene_column] ~ clinical_stage)
  model_summary <- summary(model)
  p_value <- coef(model_summary)[2, 4]
  k = coef(model_summary)[2,1]
  if (is.nan(p_value)) {
    p_value <- 1
  }
  
  if (p_value <= 0.05&coef(model_summary)[2,1]>0){
    significant_genes <- c(significant_genes, gene_column)
  }
  p_value_final = c(p_value_final,p_value)
  k_final = c(k_final,k)
}

df.final = data.frame(K=k_final,P.value = p_value_final,gene = colnames(gene_expression))

library(clusterProfiler)
library(org.Hs.eg.db)
gene.id <- mapIds(org.Hs.eg.db, 
                  keys =significant_genes , keytype = "SYMBOL", column="ENTREZID")
GO=enrichGO(gene = gene.id,
            OrgDb = org.Hs.eg.db, 
            pvalueCutoff =0.05,	
            qvalueCutoff = 0.05,	
            ont='BP',	
            readable =T)	
dotplot(GO,showCategory=5)

result=GO@result
interest.pathway=c('regulation of cell-cell adhesion',
                   'positive regulation of response to endoplasmic reticulum stress',
                   'response to reactive oxygen species',
                   'response to temperature stimulus',
                   'positive regulation of cellular catabolic process')

enrich2=filter(result,Description %in% interest.pathway)
dotplot(enrich2)

enrich2.df <- data.frame(enrich2) 
a=enrich2.df$GeneRatio
a <- data.frame(name=a)
b <- apply(a,1,function(x) eval(parse(text=x)))
enrich2.df$GeneRatio=b

enrich2.df$Description = factor(enrich2.df$Description,levels = c('regulation of cell-cell adhesion',
                                                                  'positive regulation of cellular catabolic process',
                                                                  'response to reactive oxygen species',
                                                                  'response to temperature stimulus',
                                                                  'positive regulation of response to endoplasmic reticulum stress'))

library(ggplot2)

a=ggplot(enrich2.df, aes(x = Description, y = GeneRatio,fill=pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(high= "#408591",low="#C6736E")+
  labs(x = "Pathway Description", y = "Gene Count", title = "Enrichment of Interest Pathways") +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


gene_of_interest <- "TPT1"  
gene_index <- match(gene_of_interest, df.final$gene)

p_value <- p_value_final[gene_index]
k_value <- k_final[gene_index]

name.row = rownames(expr.filtered)[sample(1:10107,2000)]
expr.filtered2 = expr.filtered[name.row,]
gene_expression2= gene_expression[name.row,]
a=ggplot(data = expr.filtered2, aes(x = Clinical.stage, y = gene_expression2[, gene_of_interest])) +
  geom_jitter(size=0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue",size=2) +
  labs(x = "Clinical Stage", y = paste(gene_of_interest, "Expression"), 
       title = paste("Regression Trend of", gene_of_interest),
       subtitle = paste("p-value =", round(p_value, 3), "  Coefficient =", round(k_value, 3))) +
  #theme_minimal()+
  ylim(0,10)+
  labs(x = "Clinical Stage", y = paste(gene_of_interest, "Expression"), 
       title = paste("Regression Trend of", gene_of_interest),
       subtitle = paste("p-value =", round(p_value, 3), "  Coefficient =", round(k_value, 3)))


