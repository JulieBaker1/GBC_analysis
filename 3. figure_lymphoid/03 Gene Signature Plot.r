#### Gene dotplot of celltypes####
cd8.marker = readRDS('...\\cd8_marker.rds')
cd4.marker = readRDS('...\\cd4_marker.rds')
nk.marker = readRDS('...\\nk_marker.rds')
b.marker = readRDS('...\\b_marker.rds')
pc.marker = readRDS('...\\pc__marker.rds')


top50=cd8.marker%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
cd8.sub=subset(cd8,seurat_clusters %in% c(0,1,2,3,6,8,9,11,12))
top50 = subset(top50,cluster %in%c(0,1,2,3,6,8,9,11,12))
a=DotPlot(cd8.sub, features = unique(top50$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())


top50=cd4.marker%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
cd4.sub=subset(cd4,seurat_clusters %in% c(0,1,2,3,4,5,7,8))
top50 = subset(top50,cluster %in%c(0,1,2,3,4,5,7,8))
a=DotPlot(cd4.sub, features = unique(top50$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())


top50=nk.marker%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
nk.sub=subset(nk,seurat_clusters %in% c(0,1,2,3,4,6,8,12))
top50 = subset(top50,cluster %in%c(0,1,2,3,4,6,8,12))
a=DotPlot(nk.sub, features = unique(top50$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())


top50=b.marker%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
b.sub=subset(b,seurat_clusters %in% c(0,1,3,5,8))
top50 = subset(top50,cluster %in%c(0,1,3,5,8))
a=DotPlot(b.sub, features = unique(top50$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())


top50=pc.marker%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
pc.sub=subset(pc,seurat_clusters %in% c(0))
top50 = subset(top50,cluster %in%c(0))
a=DotPlot(pc.sub, features = unique(top50$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5))+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())


### Figure 4F
cd4= readRDS('D:\\GBC\\03 CD4\\01 矩阵RDS\\cd4_newname.rds')
cd4 = subset(cd4,seurat_clusters == 2)
marker = readRDS('D:\\GBC\\18 Final Final code file\\04 Differential Genes\\02 CD4 markers.rds')
sub.marker = subset(marker,cluster ==2)
features = c(top_n(sub.marker,5,avg_log2FC)$gene,'SELL','TCF7','LEF1','IL7R')
a=DotPlot(cd4,features = features)+scale_size(range = c(3, 10))+
  scale_color_gradient(low = "gray",high = "#BD3438")
