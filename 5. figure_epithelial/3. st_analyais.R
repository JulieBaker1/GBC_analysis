library(Seurat)
setwd("D:/postgraduate/GBC/combined_data/epithelial/spatial/")
sample_name_list = c("GBC_055","GBC_025")
celltype = readRDS("D:/postgraduate/GBC/combined_data/cellchat/result20230312/combined_celltype_include_normal20230312.RDS")
# add the celltype information to single cell data
for(iname in sample_name_list[1:2]){
  sample_RDS = readRDS(paste0("./single_cell/",iname,"_P.RDS"))
  sample_RDS = RenameCells(sample_RDS,new.names = paste0(iname,"_P_",colnames(sample_RDS)))
  sample_RDS$celltype = celltype[colnames(sample_RDS),"celltype"]
  sample_RDS$subtype = celltype[colnames(sample_RDS),"subtype"]
  saveRDS(sample_RDS,paste0("./",iname,"_P.RDS"))
}

# seurat_pipline ----------------------------------------------------------
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
sample_name_list = c("GBC_055","GBC_025")
file_dir = list.files("./spatial/")
for(iindex in 1:2){
  iname = sample_name_list[iindex]
  ifile_dir = file_dir[[iindex]]
  brain = Load10X_Spatial(data.dir = paste0("./spatial/",ifile_dir))
  median(brain$nFeature_Spatial)
  median(brain$nCount_Spatial)
  dim(brain)
  #brain = brain[,brain$nCount_Spatial<25000]
  brain = brain[,brain$nCount_Spatial>1000]
  dim(brain)
  brain = Load10X_Spatial(data.dir = paste0("./spatial/",ifile_dir))
  dim(brain)
  brain = brain[,brain$nFeature_Spatial>200]
  dim(brain)
  plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  #plot1
  plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
  #plot2
  plot3 <- VlnPlot(brain, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  #plot3
  plot4 <- SpatialFeaturePlot(brain, features = "nFeature_Spatial") + theme(legend.position = "right")
  #plot4
  # wrap_plots(plot1, plot2,plot3,plot4,ncol = 2)
  brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
  brain = AddModuleScore(brain,features = as.list(gene_list)[c(1,3,5,6,7,16)],name = paste0("GM",c(1,3,5,6,7,16)),k = FALSE)
  GM = apply(brain@meta.data[,grepl("GM",colnames(brain@meta.data))],1,FUN = function(i){which(i==max(i))})
  brain$epi_state = GM
  library(ggplot2)
  brain$GM1 = brain$GM11
  brain$GM3 = brain$GM32
  brain$GM5 = brain$GM53
  brain$GM6 = brain$GM64
  brain$GM7 = brain$GM75
  brain$GM16 = brain$GM166
  plot <- SpatialFeaturePlot(brain, features = c("GM1","GM3","GM5","GM6","GM7","GM16")) 
  # + 
  #   theme(legend.text = element_text(size = 0),
  #         legend.title = element_text(size = 20), 
  #         legend.key.size = unit(1, "cm"))
  #plot
  png(filename = paste0("./output/images/",iname,"GM_score.png"), height = 455, width = 615)
  print(plot)
  dev.off()
  
  brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
  brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
  brain <- FindClusters(brain, verbose = FALSE,resolution = 0.5)
  brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
  
  p1 <- DimPlot(brain, reduction = "umap", label = TRUE)

  p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
  png(filename = paste0("./output/images/",iname,"_cluster.png"), height = 455, width = 615)
  print(p1 + p2)
  dev.off()
  png(filename = paste0("./output/images/",iname,"_cluster_separate.png"), height = 361, width = 588)
  print(SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(0,1,2,3,4,5)), facet.highlight = TRUE, ncol = 3))
  dev.off()
  #SpatialDimPlot(brain, interactive = TRUE)
  #SpatialFeaturePlot(brain, features = "VEGFA")
  
  de_markers <- FindAllMarkers(brain)
  saveRDS(de_markers,paste0("./spatial/",ifile_dir,"/",iname,"_P_markers.RDS"))
  # SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
  
  de_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  png(filename = paste0("./output/images/",iname,"_DE.png"), height = 712, width = 671)
  print(DoHeatmap(brain, features = top10$gene) + NoLegend())
  dev.off()
  
  brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                         selection.method = "moransi")
  top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
  png(filename = paste0("./output/images/",iname,"_SVG.png"), height = 595, width = 629)
  print(SpatialFeaturePlot(brain, features = top.features, ncol = 3))
  # SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
  dev.off()
  

  allen_reference = readRDS(paste0("./single_cell/",iname,"_P.RDS"))
  allen_reference = allen_reference[,!is.na(allen_reference$celltype)]
  allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
  png(filename = paste0("./output/images/",iname,"celltype.png"), height = 512, width = 630)
  print(DimPlot(allen_reference, group.by = "celltype", label = TRUE))
  dev.off()
  
  
  anchors <- FindTransferAnchors(reference = allen_reference, query = brain, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$celltype, prediction.assay = TRUE,
                                    weight.reduction = brain[["pca"]], dims = 1:30)
  brain[["predictions"]] <- predictions.assay
  DefaultAssay(brain) <- "predictions"
  png(filename = paste0("./output/images/",iname,"celltype_score.png"), height = 595, width = 629)
  print(SpatialFeaturePlot(brain, features = unique(allen_reference$celltype), pt.size.factor = 1.6, ncol = 4, crop = TRUE))
  dev.off()
  saveRDS(brain,paste0("./spatial/",ifile_dir,"/",iname,"_spatial.RDS")) 
}


### plot 
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
setwd("D:/postgraduate/GBC/combined_data/epithelial/spatial/")
GBC_spatial_list = list()
GBC_single_cell_list = list()
celltype_max_list = list()
celltype_max_score_list = list()

sample_name = c("GBC_025","GBC_055")
spatial_dir = c("20230322-V42U13-120-A1-GBC-025-outs",'20230322-V42U13-004-D1-GBC-055-outs')
for(i in 1:2){
  GBC_spatial_list[[paste0(sample_name[i],"_spatial")]] = readRDS(paste0("./spatial/",spatial_dir[i],"/",sample_name[i],"_spatial.RDS"))
  GBC_single_cell_list[[paste0(sample_name[i],"_P")]] = readRDS(paste0("./single_cell/",sample_name[i],"_P.RDS"))
  
  celltype = GBC_spatial_list[[paste0(sample_name[i],"_spatial")]]@assays$predictions@data
  celltype_max_list[[paste0(sample_name[i],"_spatial")]] = rownames(celltype)[apply(celltype[1:(nrow(celltype)-1),],2,FUN = function(i){which(i==max(i))})]
  celltype_max_score_list[[paste0(sample_name[i],"_spatial")]] = apply(celltype[1:(nrow(celltype)-1),],2,FUN = function(i){max(i)})
}

# calculate the correlation ------------------------------------------------------------------

for(i in 1:2){
  index = 1
  ggplot_list = list()
  pdf(paste0("./figures/",sample_name[i],"_GM16_immune_subtype.pdf"),height = 3,width = 15)
  for(itype in c("CD4T_C2_CCR7","DC_C4_cDC2_PPP1R14A","M_C1_S100A8")){
    for(igene in c(200)){
      #pdf(paste0("./immune/figures_GM16/",itype,"_top",igene,".pdf"))
      celltype.markers = read.csv(paste0("./immune/",itype,"_markers.csv"),row.names = 1)
      celltype.markers = celltype.markers[,1][1:igene]
      
      GM6.markers = read.csv(paste0("./immune/","GM16","_markers.csv"),row.names = 1)
      GM6.markers = GM6.markers[,1][GM6.markers$GM16_p<0.05][1:200]
      
      
      GBC_spatial_list[[paste0(sample_name[i],"_spatial")]] = AddModuleScore(GBC_spatial_list[[paste0(sample_name[i],"_spatial")]],
                                                                             features = list(celltype.markers),name = "score")
      GBC_spatial_list[[paste0(sample_name[i],"_spatial")]] = AddModuleScore(GBC_spatial_list[[paste0(sample_name[i],"_spatial")]],
                                                                             features = list(GM6.markers),name = "GM16")
      data_plot = data.frame(score = GBC_spatial_list[[i]]$score1[celltype_max_list[[i]]=="Malignant" & 
                                                                    celltype_max_score_list[[i]]<=upper_bound & celltype_max_score_list[[i]]>0.3],
                             GM6 = GBC_spatial_list[[i]]$GM161[celltype_max_list[[i]]=="Malignant" & 
                                                                 celltype_max_score_list[[i]]<=upper_bound & celltype_max_score_list[[i]]>0.3])
      data_plot$GM6 = (data_plot$GM6 - min(data_plot$GM6))/(max(data_plot$GM6)-min(data_plot$GM6))
      data_plot$score = (data_plot$score - min(data_plot$score))/(max(data_plot$score)-min(data_plot$score))
      print(index)
      ggplot_list[[index]] = ggplot(data_plot, aes(score,GM6))+theme_classic()+
        geom_point(size=1,shape=1)+
        geom_smooth(method = "lm")+stat_cor( method = "pearson")+xlab(paste0(itype,"_score"))+ylab("GM16")+ggtitle(paste0(itype,"top_",igene))
      index = index + 1 
    }
  }
  print(ggplot_list[[1]]+ggplot_list[[2]]+ggplot_list[[3]]+
          ggplot_list[[4]]+ggplot_list[[5]]+
          plot_layout(nrow = 1,guides = 'collect') & theme(legend.position = "bottom"))
  dev.off()
}

# plot the spatial transcriptome ------------------------------------------
library(RColorBrewer)
iindex = 1
ggplot_list_v3 = list()
for(i in c(1,2)){
  for(itype in c("GM16","CD4T_C2_CCR7",'DC_C4_cDC2_PPP1R14A',"M_C1_S100A8")){
    
    plot_spatial = GBC_spatial_list[[i]][,celltype_max_list[[i]]=="Malignant"]
    score = plot_spatial@meta.data[,paste0(itype,"1")]
    score = (score - min(score))/(max(score)-min(score))
    print(max(score))
    plot_spatial@meta.data[,paste0(itype,"1")] = score
    
    
    ggplot_list_v3[[iindex]] = SpatialPlot(plot_spatial,
                                           features = paste0(itype,"1"),image.alpha = c(0.1),stroke = NA,
                                           pt.size.factor = 4,crop = TRUE)+
      theme(legend.title=element_blank())+scale_fill_gradientn(colours = 
                                                                 colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                                                               values = seq(0,1,0.01),limits = c(0,1),breaks = seq(0,1,0.2))
    print(max(ggplot_list_v3[[iindex]]$data[,3]))
    iindex = iindex+1
  }
}




pdf("./figures/combined_spatial_image0.1_size4_V3.pdf",height = 2,width = 4)
print(ggplot_list_v3[[1]]+ggplot_list_v3[[2]]+ggplot_list_v3[[3]]+
        ggplot_list_v3[[4]]+ggplot_list_v3[[5]]+ggplot_list_v3[[6]]+
        ggplot_list_v3[[7]]+ggplot_list_v3[[8]]+
        plot_layout(nrow = 2,guides = 'collect') & theme(legend.position = "bottom"))
dev.off()





## the position of malignant spot
for(i in c(1,2)){
  GBC_spatial_list[[i]]$spot_type = "Non-malignant"
  GBC_spatial_list[[i]]$spot_type[celltype_max_list[[i]]=="Malignant"] = "Malignant"
  Idents(GBC_spatial_list[[i]]) = GBC_spatial_list[[i]]$spot_type
  aa = GBC_spatial_list[[i]][,celltype_max_list[[i]]=="Malignant"]
  pdf(paste0("./figures/",sample_name[i],".pdf"))
  print(SpatialDimPlot(aa, cells.highlight = CellsByIdentities(object = aa, idents = c("Malignant")), facet.highlight = TRUE,
                       stroke = NA))
  dev.off()
}

