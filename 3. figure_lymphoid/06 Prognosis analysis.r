
library(dplyr)
library(reshape2)
library("survival")
library("survminer")
library('ggplot2')

norm.total =readRDS('...\\02 Lymphoid Celltype Normalized.rds')
info.yuhou=read.table('...\\GBC Prognosis 0916 simplify.csv',sep=',',header = TRUE,
                     quote = "\"")

yuhou.final.total= norm.total %>% left_join(info.yuhou,by='NewSample.ID')
cellname=colnames(yuhou.final.total)[c(1:64)]
yuhou=c("DFS_month","OS_month","event")
colnames(info.yuhou)
yuhou.final.total= yuhou.final.total %>% filter(histological.type.short=='adeno',Tumors.for.scRNA.seq.short=='P')
yuhou.final.total$state = ifelse(yuhou.final.total$event=='dead',"1",'0')
cutoff_yuhou=yuhou.final.total[,c(cellname,yuhou,'state')]
cutoff_yuhou2=cutoff_yuhou[complete.cases(cutoff_yuhou[,'state']),]
cutoff_yuhou2$state=as.numeric(cutoff_yuhou2$state)

plot_list=list()
j=1

for (cell.name in cellname){
  temp = cutoff_yuhou2[,c(cell.name,yuhou,'state')]
  colnames(temp)[1]='celltype'
  temp$group=ifelse(temp$celltype>=median(temp$celltype),
                    "High",'Low')
  temp$state=as.numeric(temp$state)
  fit=survfit(Surv(OS_month,state)~group,data=temp)
  a=ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               #palette = c("#E7B800", "#2E9FDF")
  )+ggtitle(cell.name)
  plot_list[[j]]=a
  j=j+1
}

colnames(cutoff_yuhou2)
cell.name = 'CD8T_C3_KLRC1'
temp = cutoff_yuhou2[,c(cell.name,yuhou,'state')]
colnames(temp)[1]='celltype'
temp$group=ifelse(temp$celltype>=median(temp$celltype),
                  "High",'Low')
temp$state=as.numeric(temp$state)
fit=survfit(Surv(OS_month,state)~group,data=temp)
a = ggsurvplot(
  fit,                     # survfit object with calculated statistics.
 # data = sub,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  palette = "npg",
  xlab = "Time in months",   # customize X axis label.
  #ggtheme = theme_bw(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
  #tables.y.text = T,
  risk.table.pos = "in",
  risk.table.col = "strata",
  fontsize = 4,
  pval.size = 4,
  #surv.plot.height = 0.8,
  #tables.height = 0.5,
  pval.coord = c(9, 0.9),
  legend = "top",
  title = cell.name
)


a = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               #linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_classic(), # Change ggplot2 theme
               palette = c('#d62728','#1f77b4')
)+ggtitle(cell.name)





