## 3-4. Endothelial Cell Subtype ####
Endothelial_RawData = readRDS("./1 - GBC_input data/Endothelial.RDS")
MetaData0609 = read.csv("./1 - GBC_input data/meta_data0609.csv", row.names = 1)

Endothelial_cellid_celltype = rownames(subset(MetaData0609, celltype == "Endothelial"))
Endothelial_cellid_celltype_scanpy = rownames(subset(MetaData0609, celltype_scanpy == "Endothelial"))
Endothelial_cellid_selected = intersect(Endothelial_cellid_celltype, Endothelial_cellid_celltype_scanpy)

ig_genes = c(rownames(Endothelial_RawData)[grep("^IGJ",rownames(Endothelial_RawData))],
             rownames(Endothelial_RawData)[grep("^IGH",rownames(Endothelial_RawData))],
             rownames(Endothelial_RawData)[grep("^IGK",rownames(Endothelial_RawData))],
             rownames(Endothelial_RawData)[grep("^IGL",rownames(Endothelial_RawData))])

Endothelial_filtered = subset(Endothelial_RawData, features = rownames(Endothelial_RawData)[!(rownames(Endothelial_RawData) %in% ig_genes)])
Endothelial_filtered = subset(Endothelial_filtered, cells = Endothelial_cellid_selected)

VlnPlot(Endothelial_filtered, features = c("nFeature_RNA", "nCount_RNA", "mito.percent"), ncol = 3, raster=FALSE) # mito.percent has been subset
plot1 <- FeatureScatter(Endothelial_filtered, feature1 = "nCount_RNA", feature2 = "mito.percent")
plot2 <- FeatureScatter(Endothelial_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Endothelial_filtered <- NormalizeData(Endothelial_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
Endothelial_filtered <- FindVariableFeatures(Endothelial_filtered, selection.method = "vst", nfeatures = 2000)
Endothelial_filtered <- ScaleData(Endothelial_filtered)
Endothelial_filtered <- RunPCA(Endothelial_filtered, features = VariableFeatures(object = Endothelial_filtered))
Endothelial_filtered <- JackStraw(Endothelial_filtered, num.replicate = 100)
Endothelial_filtered <- ScoreJackStraw(Endothelial_filtered, dims = 1:30)
ElbowPlot(Endothelial_filtered)
Endothelial_filtered <- FindNeighbors(Endothelial_filtered, dims = 1:20)
Endothelial_filtered <- FindClusters(Endothelial_filtered, resolution = 0.3)
Endothelial_filtered <- RunUMAP(Endothelial_filtered, dims = 1:20)

DimPlot(Endothelial_filtered, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(Endothelial_filtered, features = "doublet.score")

Endothelial_filtered <- subset(Endothelial_filtered, idents = c(4,5,7,9,12), invert = TRUE)
Endothelial_filtered <- NormalizeData(Endothelial_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
Endothelial_filtered <- FindVariableFeatures(Endothelial_filtered, selection.method = "vst", nfeatures = 2000)
Endothelial_filtered <- ScaleData(Endothelial_filtered)
Endothelial_filtered <- RunPCA(Endothelial_filtered, features = VariableFeatures(object = Endothelial_filtered))
Endothelial_filtered <- JackStraw(Endothelial_filtered, num.replicate = 100)
Endothelial_filtered <- ScoreJackStraw(Endothelial_filtered, dims = 1:30)
ElbowPlot(Endothelial_filtered)
Endothelial_filtered <- FindNeighbors(Endothelial_filtered, dims = 1:20)
Endothelial_filtered <- FindClusters(Endothelial_filtered, resolution = 0.3)
Endothelial_filtered <- RunUMAP(Endothelial_filtered, dims = 1:20)
DimPlot(Endothelial_filtered, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(Endothelial_filtered, features = c('CXCR4', # common tip
                                               'ACKR1', # human vein
                                               'GJA5', # human artery
                                               # 'CA4', # human capillary
                                               'PROX1', # human lymphatic
                                               # 'VWF', # mouse vein
                                               'KDR', # mouse capillary
                                               # 'SOX17', # mouse artery
                                               'MKI67'
))

Endothelial_filtered_DEGs <- FindAllMarkers(Endothelial_filtered, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Endothelial_filtered_DEGs %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DotPlot(Endothelial_filtered, features = unique(top5$gene))  +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5))
new.cluster.ids = c(
  "EC_C0_ACKR1",
  "EC_C1_KDR",
  "EC_C2_CXCR4",
  "EC_C3_GJA5",
  "EC_C4_TMSB4X",
  "EC_C5_PROX1",
  "EC_C6_MKI67",
  "EC_C7_FCN3",
  "EC_C8_ACKR1_MT1X"
)
names(new.cluster.ids) <- levels(Endothelial_filtered)
Endothelial_filtered <- RenameIdents(Endothelial_filtered, new.cluster.ids)

Idents(Endothelial_filtered) = Endothelial_filtered$RNA_snn_res.0.3
saveRDS(Endothelial_filtered, file = './GBC_Endothelium_v3_output/Endothelium.RDS')



# 8. Enrichment of endothelial subtypes: NT vs Adeno ####
# * including SFig7N
Endo_Subtype_counts = Myeloid_Endo_Subtype_counts[grep("EC_C", Myeloid_Endo_Subtype_counts$Subtype),]
cycle = unique(Endo_Subtype_counts$orig.ident)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(Endo_Subtype_counts, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round((table(sub$Subtype)/sum(table(sub$Subtype)))*100,2))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0 # Export to WYH
Freq = Freq[,30:38]
saveRDS(Freq, file = "./3 - GBC_output data_subtype/04_Endothelium/GBC_allsample_EndoPct.RDS")
Freq_melt = Freq
Freq_melt$Patient = rownames(Freq_melt)
Freq_melt$Patient = factor(Freq_melt$Patient, levels = Freq_melt$Patient)
Freq_melt <- reshape2::melt(Freq_melt, id.vars=c("Patient"), variable.name="Subtype", value.name="Percentage")
Freq_melt = Freq_melt[order(Freq_melt$Subtype),]
Freq_melt$Group = "Adeno"
Freq_melt$Group[grep("CC|HG|LG|XGC", Freq_melt$Patient)] = "Non-Tumor"
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_metastasis = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% "adeno",]
PatientInfo_allP = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% "P",]
PatientInfo_nonT = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c('CC','XGC','HG','LG'),]
Freq_melt = Freq_melt[Freq_melt$Patient %in% c(PatientInfo_allP$X,
                                               PatientInfo_nonT$X),]
Freq_melt$Group = factor(Freq_melt$Group, levels = c("Non-Tumor","Adeno"))
Freq_melt$Subtype = as.character(Freq_melt$Subtype)
cycle = unique(Freq_melt$Subtype)
for (i in 1:length(cycle)) {
  temp = Freq_melt[grep(cycle[i],Freq_melt$Subtype),]
  
  p_adj <- get_adj_p(temp,
                     .col = "Percentage", .grp = "Group", p.adjust.method = "BH",
                     comparisons = list(c("Non-Tumor","Adeno")))
  
  p <- ggboxplot(temp, x = "Group", y = "Percentage", color = "Group", palette = "jco",add = "jitter") + 
    stat_pvalue_manual(p_adj, label = "p.adj", hide.ns = T)+
    theme(legend.position = "top") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_blank())  +
    ggtitle(cycle[i])
  pdf(file = paste0("./3 - GBC_output data_subtype/06_EndoSubtype_Enrichment_NTvsAdenoP/", cycle[i], ".pdf"), width =2, height = 3)
  print(p)
  dev.off()
}



# 10. Enrichment of endothelial subtypes: P vs P_Mets ####
# * including SFig7O
Endo_Subtype_counts = Myeloid_Endo_Subtype_counts[grep("EC_C", Myeloid_Endo_Subtype_counts$Subtype),]
cycle = unique(Endo_Subtype_counts$orig.ident)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(Endo_Subtype_counts, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round((table(sub$Subtype)/sum(table(sub$Subtype)))*100,2))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0
Freq = Freq[,30:38]
Freq_melt = Freq
Freq_melt$Patient = rownames(Freq_melt)
Freq_melt$Patient = factor(Freq_melt$Patient, levels = Freq_melt$Patient)
Freq_melt <- reshape2::melt(Freq_melt, id.vars=c("Patient"), variable.name="Subtype", value.name="Percentage")
Freq_melt = Freq_melt[order(Freq_melt$Subtype),]
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_DC = subset(PatientInfo_DC, select = c("NewSample.ID","metastasis.type","histological.type.short","Tumors.for.scRNA.seq.short"))
PatientInfo_DC = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c("adeno"),]
PatientInfo_DC$metastasis.type = ifelse(PatientInfo_DC$metastasis.type == "P_LI", "P", PatientInfo_DC$metastasis.type)
PatientInfo_DC = PatientInfo_DC[PatientInfo_DC$Tumors.for.scRNA.seq.short %in% c("P"),]
PatientInfo_DC$metastasis.type = ifelse(PatientInfo_DC$metastasis.type %in% c("P_LM","P_LN"), "P_Mets", PatientInfo_DC$metastasis.type)
colnames(PatientInfo_DC)[1] = "Patient"
Freq_melt = left_join(Freq_melt,PatientInfo_DC, by="Patient")
Freq_melt = Freq_melt[Freq_melt$metastasis.type %in% c("P","P_Mets"),]
Freq_melt$Subtype = as.character(Freq_melt$Subtype)
colnames(Freq_melt)[4] = "Group1"
Freq_melt$Group1 = factor(Freq_melt$Group1, levels = c("P","P_Mets"))
cycle = unique(Freq_melt$Subtype)
for (i in 1:length(cycle)) {
  temp = Freq_melt[grep(cycle[i],Freq_melt$Subtype),]
  
  p_adj <- get_adj_p(temp,
                     .col = "Percentage", .grp = "Group1", p.adjust.method = "BH",
                     comparisons = list(c("P","P_Mets")))
  
  p <- ggboxplot(temp, x = "Group1", y = "Percentage", color = "Group1", palette = "jco",add = "jitter") + 
    stat_pvalue_manual(p_adj, label = "p.adj", hide.ns = T)+
    theme(legend.position = "top") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_blank())  +
    ggtitle(cycle[i])
  pdf(file = paste0("./3 - GBC_output data_subtype/08_EndoSubtype_Enrichment_AdenoPs/Pmet_", cycle[i], ".pdf"), width =2, height = 3)
  print(p)
  dev.off()
}

# 16. EC ####
## Fig1E ####
DimPlot(Endothelium, reduction = "umap", pt.size = 2, raster = T) + ggplot2::theme(legend.position = "right", panel.background = element_blank())

## FigS2H ####
p = FeaturePlot(Endothelium, features = c("ACKR1","KDR","CXCR4","GJA5","PROX1"),ncol=5, pt.size = 3, raster = T, max.cutoff = c(4,3,3,3,3))
pdf(file = "./3 - GBC_output data_subtype/04_Endothelium/Endothelium_SubtypeDefinition_MarkerExprs.pdf", width =14, height = 2.5)
print(p)
dev.off()

## * Subset common subgroups ####
Endothelium = subset(Endothelium, idents = c("EC_C8_ACKR1_MT1X","EC_C7_FCN3"), invert=T)

## FigS7K ####
Endothelium_DEGs_allpos <- FindAllMarkers(Endothelium, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Endothelium_DEGs_allpos = Endothelium_DEGs_allpos[Endothelium_DEGs_allpos$p_val_adj < 0.05,]
Endothelium_DEGs_allpos %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
p = DotPlot(Endothelium, features = unique(top5$gene))  +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
pdf(file = "./3 - GBC_output data_subtype/04_Endothelium/Endo_SubtypeDEGs_Top5.pdf", width =15, height = 3.5)
print(p)
dev.off()

## FigS7L ####
Endothelium_DEGs_allpos$cluster = as.character(Endothelium_DEGs_allpos$cluster)
Endothelium_DEGs_allpos$cluster = as.factor(Endothelium_DEGs_allpos$cluster)
temp <- bitr(Endothelium_DEGs_allpos$gene,fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb='org.Hs.eg.db')
Endothelium_DEGs_allpos = Endothelium_DEGs_allpos[temp$SYMBOL,]
Endothelium_DEGs_allpos$gene = temp$ENTREZID
formula_res <- compareCluster(gene~cluster, data=Endothelium_DEGs_allpos,
                              fun="enrichGO",
                              OrgDb = org.Hs.eg.db,
                              ont = "BP" ,
                              pAdjustMethod = "BH",
                              readabl = TRUE)
formula_res <- filter(formula_res, cluster %in% unique(Idents(Endothelium)))
dotplot(formula_res, showCategory = 5) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "horizontal")
formula_res <- filter(formula_res, Description %in% c("antigen processing and presentation",
                                                      "receptor-mediated endocytosis", # EC5 top
                                                      "ATP metabolic process", # EC4/6 top
                                                      "oxidative phosphorylation",
                                                      "cell-substrate adhesion",
                                                      "regulation of angiogenesis",
                                                      "extracellular matrix organization"))
p = dotplot(formula_res, showCategory = 5) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 10),
        axis.text.y = element_text(size = 10),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "vertical")
pdf(file = "./3 - GBC_output data_subtype/04_Endothelium/Endo_Sub_Enriched_pathway.pdf", width =6, height = 4.5)
print(p)
dev.off()

## FigS7P ####
Endothelium1 = subset(Endothelium, idents = c("EC_C0_ACKR1","EC_C3_GJA5","EC_C5_PROX1"))
DotPlot(Endothelium1, features = "FLT1")  +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

## FigS7Q ####
Freq <- readRDS("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/04_Endothelium/GBC_allsample_EndoPct.RDS")
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_metastasis = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% "adeno",]
PatientInfo_allP = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% "P",]
Freq = Freq[PatientInfo_allP$NewSample.ID,]
p = ggscatter(Freq, x = 'EC_C2_CXCR4', y = 'EC_C0_ACKR1',
              fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
              add = "reg.line",# Add regressin line
              add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = T, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"))+xlab('EC_C2_CXCR4')+ylab('EC_C0_ACKR1')
pdf(file = "./3 - GBC_output data_subtype/04_Endothelium/EC2_EC0_PctCor_Pearson.pdf", width =4, height = 4)
print(p)
dev.off()
p = ggscatter(Freq, x = 'EC_C2_CXCR4', y = 'EC_C3_GJA5',
              fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
              add = "reg.line",# Add regressin line
              add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = T, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"))+xlab('EC_C2_CXCR4')+ylab('EC_C3_GJA5')
pdf(file = "./3 - GBC_output data_subtype/04_Endothelium/EC2_EC3_PctCor_Pearson.pdf", width =4, height = 4)
print(p)
dev.off()

## FigS7R ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Endothelium@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Endothelium@meta.data$NewSample.ID = Endothelium@meta.data$orig.ident
temp = left_join(Endothelium@meta.data, PatientInfo, by = "NewSample.ID")
identical(Endothelium@meta.data$orig.ident, temp$orig.ident)
Endothelium@meta.data$histological.type.short = temp$histological.type.short
Endothelium@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Endothelium@meta.data$metastasis.type = temp$metastasis.type
Endothelium@meta.data$Clinical.stage = temp$Clinical.stage.x
## EC_C0_ACKR1
Endothelium_0 = subset(Endothelium, idents = "EC_C0_ACKR1")
Endothelium_0 = subset(Endothelium_0, histological.type.short %in% c("adeno","CC","HG","LG","XGC"))
Endothelium_0 = subset(Endothelium_0, Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"))
Endothelium_0$group = ifelse(Endothelium_0$histological.type.short == "adeno", "adeno", "NT")
Idents(Endothelium_0) = Endothelium_0$group
## DEGs: Volcano
Endothelium_0_P_DEGs <- FindMarkers(Endothelium_0, ident.1 = "adeno", only.pos = T,min.pct = 0.25, logfc.threshold = 0.25)
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[Endothelium_0_P_DEGs$p_val_adj < 0.05,]
Endothelium_0_P_DEGs$cluster = "Vein"
ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+
  geom_text_repel(
    data = Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>0,],
    aes(label = Gene),
    size = 3,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+
  theme(
    legend.title = element_blank()
  )+
  ylab('-log10 (p-adj)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
  ggtitle("EC_C0_ACKR1 [Adeno (vs NT)]")

## EC_C3_GJA5
Endothelium_3 = subset(Endothelium, idents = "EC_C3_GJA5")
Endothelium_3 = subset(Endothelium_3, histological.type.short %in% c("adeno","CC","HG","LG","XGC"))
Endothelium_3 = subset(Endothelium_3, Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"))
Endothelium_3$group = ifelse(Endothelium_3$histological.type.short == "adeno", "adeno", "NT")
Idents(Endothelium_3) = Endothelium_3$group
## DEGs: Volcano
Endothelium_3_P_DEGs <- FindMarkers(Endothelium_3, ident.1 = "adeno", only.pos = T,min.pct = 0.25, logfc.threshold = 0.25)
Endothelium_3_P_DEGs = Endothelium_3_P_DEGs[Endothelium_3_P_DEGs$p_val_adj < 0.05,]
Endothelium_3_P_DEGs$cluster = "Artery"

## Vein + Artery GO 
Endothelium_P_DEGs = rbind(Endothelium_0_P_DEGs,Endothelium_3_P_DEGs)
Endothelium_P_DEGs$cluster = factor(Endothelium_P_DEGs$cluster, levels=c("Vein","Artery"))
Endothelium_P_DEGs$gene = rownames(Endothelium_P_DEGs)
temp <- bitr(Endothelium_P_DEGs$gene,fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb='org.Hs.eg.db')
Endothelium_P_DEGs = Endothelium_P_DEGs[temp$SYMBOL,]
Endothelium_P_DEGs$gene = temp$ENTREZID
formula_res <- compareCluster(gene~cluster, data=Endothelium_P_DEGs,
                              fun="enrichGO",
                              OrgDb = org.Hs.eg.db,
                              ont = "BP" ,
                              pAdjustMethod = "BH",
                              readabl = TRUE)
dotplot(formula_res, showCategory = 10) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "horizontal")+
  ggtitle("EC_C0/3_UP in Adeno")
formula_res1 <- filter(formula_res, Description %in% c("endothelial cell migration",
                                                       "regulation of angiogenesis",
                                                       "endothelial cell differentiation",
                                                       "endothelial cell proliferation",
                                                       "vascular endothelial growth factor receptor signaling pathway",
                                                       "ATP metabolic process",
                                                       "glycolytic process",
                                                       "response to hypoxia"))
p = dotplot(formula_res1, showCategory = 10) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "horizontal")
pdf(file = "./3 - GBC_output data_subtype/04_Endothelium/EC0_EC3_GO_Up in Adeno.pdf", width = 4.5, height = 4)
print(p)
dev.off()