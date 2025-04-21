library(Seurat)
library(dplyr)
library(patchwork)
library(plotly)
#library(SingleR)
library(ggplot2)
#library(clusterProfiler) #
#library(org.Hs.eg.db) #
#library(enrichplot) #
library(ggpubr)
library(tidyverse)
#library(monocle) # Need R 4.0
#library(GSVA) #
library(fgsea)
library(msigdbr)
library(survminer)
library(survival)
library(ggrepel)
library(corrplot)
#library(IOBR)
library(circlize)
library(RColorBrewer)
library(GGally)
library(data.table)
#library(ComplexHeatmap) #
library(dendextend)
library(stringr)
library(psych)
library(reshape)
#library(jjPlot)
#library(gg.gap)
library(gmodels)
library(reshape2)
##library(ggradar)
library(pheatmap)

# 1. Function ####
get_adj_p <- function(data, .col, .grp = "Sample", comparisons = NULL,
                      method = "wilcox.test", p.adjust.method = "fdr", p.digits = 3L,symnum.args = NULL, ...) {
  # Compute p-values
  comparison.formula <- paste0(.col, "~", .grp) %>%
    as.formula()
  pvalues <- ggpubr::compare_means(
    formula = comparison.formula, data = data,
    method = method,
    p.adjust.method = p.adjust.method,
    ...
  )
  
  # If a comparison list is provided, extract the comparisons of interest for plotting
  if (!is.null(comparisons)) {
    pvalues <- purrr::map_df(comparisons, ~ pvalues %>% dplyr::filter(group1 == .x[1] & group2 == .x[2]))
  }
  
  # P-value y coordinates
  y.max <- data %>%
    dplyr::pull(.col) %>%
    max(na.rm = TRUE)
  p.value.y.coord <- rep(y.max, nrow(pvalues))
  
  step.increase <- (1:nrow(pvalues)) * (y.max / 10)
  p.value.y.coord <- p.value.y.coord + step.increase
  if (is.null(symnum.args)){
    symnum.args <- list(cutpoints = c(0, 1e-04, 0.001, 0.01, 
                                      0.05, 1), symbols = c("****", "***", "**", "*", 
                                                            "ns"))
  } 
  
  symnum.args$x <- as.numeric(pvalues$p.adj)
  p.adj.signif <- do.call(stats::symnum, symnum.args) %>% 
    as.character()
  pvalues$p.adj.signif = p.adj.signif
  pvalues <- pvalues %>%
    dplyr::mutate(
      y.position = p.value.y.coord,
      p.adj = format.pval(.data$p.adj, digits = p.digits)
      
    )
  
  pvalues
}
# ####

# 2. Color ####
mycolor <- readRDS("./mycolor.RDS")
scales::show_col(mycolor)

# * Summarize patient info. and related prognosis info. ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20240916预后更新.xlsx")
PatientInfo = left_join(PatientInfo, PatientInfo_surv, "NewSample.ID")
write_csv(PatientInfo, file = './1 - GBC_input data/GBC_SampleInfo_230719.csv')

# 3. Subtype Definition ####
## 3-1. Neutrophil Subtype ####
Myeloid_filtered_processed_F = readRDS("./1 - GBC_input data/Myeloid_filtered_processed_F.RDS")
Myeloid_cell_Neu <- subset(Myeloid_filtered_processed_F, idents = "Neu")

KRT_genes = sort(rownames(Myeloid_cell_Neu)[grep("^KRT",rownames(Myeloid_cell_Neu))])[1:53]
Myeloid_cell_Neu = subset(Myeloid_cell_Neu, features = rownames(Myeloid_cell_Neu)[!(rownames(Myeloid_cell_Neu) %in% KRT_genes)])

Myeloid_cell_Neu <- NormalizeData(Myeloid_cell_Neu, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_cell_Neu <- FindVariableFeatures(Myeloid_cell_Neu, selection.method = "vst", nfeatures = 2000)
Myeloid_cell_Neu <- ScaleData(Myeloid_cell_Neu)
Myeloid_cell_Neu <- RunPCA(Myeloid_cell_Neu, features = VariableFeatures(object = Myeloid_cell_Neu))
Myeloid_cell_Neu <- FindNeighbors(Myeloid_cell_Neu, dims = 1:10)
Myeloid_cell_Neu <- FindClusters(Myeloid_cell_Neu, resolution = c(0.5))
Myeloid_cell_Neu <- RunUMAP(Myeloid_cell_Neu, dims = 1:10)
DimPlot(Myeloid_cell_Neu, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

NeuSub.markers <- FindAllMarkers(Myeloid_cell_Neu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NeuSub.markers, file = "./3 - GBC_output data_subtype/01_Neutrophil/DEGs among all neutrophil subtypes_v240422.csv")

NeuSub.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DotPlot(Myeloid_cell_Neu, features = unique(top5$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "vertical")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

new.cluster.ids = c(
  "N_C0_MNDA", # 0
  "N_C1_ISG15", # 1
  "N_C2_CCL3L1", # 2
  "N_C3_CTSD", # 3
  "N_C4_CD44", # 4
  "N_C5_DEFA5", # 5
  "N_C6_ISG15_IFIT1", # 6
  "N_C7_MNDA_S100A8", # 7
  "N_C8_IL1R2", # 8
  "N_C9_S100A10", # 9
  "N_C10_STXBP2", # 10
  "N_C11_MGP" # 11
)
names(new.cluster.ids) <- levels(Myeloid_cell_Neu)
Myeloid_cell_Neu <- RenameIdents(Myeloid_cell_Neu, new.cluster.ids)
DimPlot(Myeloid_cell_Neu, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(Myeloid_cell_Neu, file = "./3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/Myeloid_cell_Neu.RDS")

### NG Appeal - to test per tumor ####
Neu_Shanno = read.csv("./3 - GBC_output data_subtype/01_Neutrophil/Neu_Entropy.csv")
ggplot(Neu_Shanno, aes(x=Cell.subtype, y=Entropy.value, fill = Cell.subtype)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#faae61","#faae61","#c9c6c3","#c9c6c3","#faae61",
                               "#c9c6c3","#c9c6c3","#c9c6c3","#c9c6c3","#faae61",
                               "#c9c6c3","#faae61")) +
  geom_hline(aes(yintercept = 0.625), linetype = "dashed") +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) + NoLegend()

Neu_H = subset(Myeloid_cell_Neu, idents = "N_C2_CCL3L1")
temp = as.data.frame(sort((table(Neu_H$orig.ident)/length(Neu_H$orig.ident)),decreasing = T))
temp$Freq = as.numeric(temp$Freq)
ggplot(temp, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity", fill = "#faae61") +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) + NoLegend()

Neu_L = subset(Myeloid_cell_Neu, idents = "N_C11_MGP")
temp = as.data.frame(sort((table(Neu_L$orig.ident)/length(Neu_L$orig.ident)),decreasing = T))
temp$Freq = as.numeric(temp$Freq)
ggplot(temp, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity", fill = "#c9c6c3") +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) + NoLegend()

## 3-2. Monocyte/Macrophage Subtype ####
Myeloid_filtered_processed_F = readRDS("./1 - GBC_input data/Myeloid_filtered_processed_F.RDS")
Myeloid_cell_MonoMacro <- subset(Myeloid_filtered_processed_F, idents = "Mono_Macro")
Myeloid_cell_MonoMacro <- NormalizeData(Myeloid_cell_MonoMacro, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_cell_MonoMacro <- FindVariableFeatures(Myeloid_cell_MonoMacro, selection.method = "vst", nfeatures = 2000)
Myeloid_cell_MonoMacro <- ScaleData(Myeloid_cell_MonoMacro)
Myeloid_cell_MonoMacro <- RunPCA(Myeloid_cell_MonoMacro, features = VariableFeatures(object = Myeloid_cell_MonoMacro))
Myeloid_cell_MonoMacro <- FindNeighbors(Myeloid_cell_MonoMacro, dims = 1:10)
Myeloid_cell_MonoMacro <- FindClusters(Myeloid_cell_MonoMacro, resolution = 0.3)
Myeloid_cell_MonoMacro <- RunUMAP(Myeloid_cell_MonoMacro, dims = 1:10)
DimPlot(Myeloid_cell_MonoMacro, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

MonoMacroSub.markers <- FindAllMarkers(Myeloid_cell_MonoMacro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(MonoMacroSub.markers, file = "./3 - GBC_output data_subtype/03_MM/DEGs among all mm subtypes_v240422.csv")

MonoMacroSub.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DotPlot(Myeloid_cell_MonoMacro, features = unique(top5$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "vertical")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
new.cluster.ids = c(
  "M_C0_FOLR2",
  "M_C1_S100A8",
  "M_C2_SPP1",
  "M_C3_FCGBP",
  "M_C4_others",
  "M_C5_PCLAF",
  "M_C6_CCL18",
  "M_C7_AGR2",
  "M_C8_MMP9"
)
names(new.cluster.ids) <- levels(Myeloid_cell_MonoMacro)
Myeloid_cell_MonoMacro <- RenameIdents(Myeloid_cell_MonoMacro, new.cluster.ids)
saveRDS(Myeloid_cell_MonoMacro, file = "./3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/Myeloid_cell_MonoMacro.RDS")

## 3-3. DC Subtype ####
Myeloid_filtered_processed_F = readRDS("./1 - GBC_input data/Myeloid_filtered_processed_F.RDS")
Myeloid_cell_DC <- subset(Myeloid_filtered_processed_F, idents = "DC")
Myeloid_cell_DC <- NormalizeData(Myeloid_cell_DC, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_cell_DC <- FindVariableFeatures(Myeloid_cell_DC, selection.method = "vst", nfeatures = 2000)
Myeloid_cell_DC <- ScaleData(Myeloid_cell_DC)
Myeloid_cell_DC <- RunPCA(Myeloid_cell_DC, features = VariableFeatures(object = Myeloid_cell_DC))
Myeloid_cell_DC <- FindNeighbors(Myeloid_cell_DC, dims = 1:10)
Myeloid_cell_DC <- FindClusters(Myeloid_cell_DC, resolution = 0.2)
Myeloid_cell_DC <- RunUMAP(Myeloid_cell_DC, dims = 1:10)
DimPlot(Myeloid_cell_DC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

DCSub.markers <- FindAllMarkers(Myeloid_cell_DC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DCSub.markers, file = "./3 - GBC_output data_subtype/02_DC/DEGs among all DC subtypes_v240422.csv")

DCSub.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DotPlot(Myeloid_cell_DC, features = unique(top5$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "vertical")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
new.cluster.ids = c(
  "DC_C0_IL1B_cDC2", # 0
  "DC_C1_others_cDC2", # 1
  "DC_C2_FSCN1", # 2
  "DC_C3_FCGBP", # 3
  "DC_C4_PPP1R14A_cDC2", # 4
  "DC_C5_CLEC9A_cDC1", # 5
  "DC_C6_GZMB_pDC", # 6
  "DC_C7_ACY3" # 7
)
names(new.cluster.ids) <- levels(Myeloid_cell_DC)
Myeloid_cell_DC <- RenameIdents(Myeloid_cell_DC, new.cluster.ids)
saveRDS(Myeloid_cell_DC, file = "./3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/Myeloid_cell_DC.RDS")

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
#Endothelial_filtered <- JackStraw(Endothelial_filtered, num.replicate = 100)
#Endothelial_filtered <- ScoreJackStraw(Endothelial_filtered, dims = 1:30)
#ElbowPlot(Endothelial_filtered)
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
#Endothelial_filtered <- JackStraw(Endothelial_filtered, num.replicate = 100)
#Endothelial_filtered <- ScoreJackStraw(Endothelial_filtered, dims = 1:30)
#ElbowPlot(Endothelial_filtered)
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
write.csv(Endothelial_filtered_DEGs, file = "./3 - GBC_output data_subtype/04_Endothelium/DEGs among all EC subtypes_v240422.csv")

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

# 4. Entropy input for myeloids and endothelial cells ####
Myeloid_cell_Neu$Subtype = Idents(Myeloid_cell_Neu)
Myeloid_cell_DC$Subtype = Idents(Myeloid_cell_DC)
Myeloid_cell_MonoMacro$Subtype = Idents(Myeloid_cell_MonoMacro)
Endothelium$Subtype = Idents(Endothelium)
Myeloid_cell_Neu_SubMeta = subset(Myeloid_cell_Neu@meta.data, select = c('orig.ident','Subtype'))
colnames(Myeloid_cell_Neu_SubMeta)[2] = "Subtype"
Myeloid_cell_DC_SubMeta = subset(Myeloid_cell_DC@meta.data, select = c('orig.ident','Subtype'))
colnames(Myeloid_cell_DC_SubMeta)[2] = "Subtype"
Myeloid_cell_MonoMacro_SubMeta = subset(Myeloid_cell_MonoMacro@meta.data, select = c('orig.ident','Subtype'))
colnames(Myeloid_cell_MonoMacro_SubMeta)[2] = "Subtype"
Endothelium_SubMeta = subset(Endothelium@meta.data, select = c('orig.ident','Subtype'))
colnames(Endothelium_SubMeta)[2] = "Subtype"
Myeloid_Endo_Subtype_counts = rbind(Myeloid_cell_Neu_SubMeta,Myeloid_cell_DC_SubMeta,Myeloid_cell_MonoMacro_SubMeta,Endothelium_SubMeta)
celltype_percent = Myeloid_Endo_Subtype_counts
cycle = unique(celltype_percent$orig.ident)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent, orig.ident == cycle[i])
  Freq = bind_rows(Freq, table(sub$Subtype))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0
ZTSubtypes_Counts_Shanno = Freq
saveRDS(ZTSubtypes_Counts_Shanno, file = './3 - GBC_output data_subtype/ToLS_Shanno/ZTSubtypes_Counts_Shanno_230605.RDS')

# 5. Roe input of histo-type for subtypes ####
Freq = ZTSubtypes_Counts_Shanno
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_Myeloid = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_adeno = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "adeno",]
PatientInfo_adenoP = PatientInfo_adeno[PatientInfo_adeno$Tumors.for.scRNA.seq.short %in% "P",]
PatientInfo_adenoP$Group = "Adeno"
PatientInfo_AdenoSqua = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "adeno squa",]
PatientInfo_AdenoSquaP = PatientInfo_AdenoSqua[PatientInfo_AdenoSqua$Tumors.for.scRNA.seq.short %in% "P",]
PatientInfo_AdenoSquaP$Group = "Adeno Squa"
PatientInfo_neuro = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "neuro",]
PatientInfo_neuroP = PatientInfo_neuro[PatientInfo_neuro$Tumors.for.scRNA.seq.short %in% "P",]
PatientInfo_neuroP$Group = "Neuro"
PatientInfo_undiff = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "undiff",]
PatientInfo_undiffP = PatientInfo_undiff[PatientInfo_undiff$Tumors.for.scRNA.seq.short %in% "P",]
PatientInfo_undiffP$Group = "Undiff"
PatientInfo_squa = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "squa",]
PatientInfo_squa$Group = "Squa"
PatientInfo_CC = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "CC",]
PatientInfo_CC$Group = "CC"
PatientInfo_HG = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "HG",]
PatientInfo_HG$Group = "HG"
PatientInfo_LG = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "LG",]
PatientInfo_LG$Group = "LG"
PatientInfo_XGC = PatientInfo_Myeloid[PatientInfo_Myeloid$histological.type.short %in% "XGC",]
PatientInfo_XGC$Group = "XGC"
PatientInfo_bind = rbind(PatientInfo_adenoP,PatientInfo_AdenoSquaP,PatientInfo_neuroP,PatientInfo_undiffP,
                         PatientInfo_squa,PatientInfo_CC,PatientInfo_HG,PatientInfo_LG,PatientInfo_XGC)
PatientInfo_bind = subset(PatientInfo_bind, select = c("X","Group"))
Freq = Freq[c(PatientInfo_adenoP$X, PatientInfo_AdenoSquaP$X, PatientInfo_neuroP$X,
              PatientInfo_undiffP$X, PatientInfo_squa$X, PatientInfo_CC$X,
              PatientInfo_HG$X, PatientInfo_LG$X, PatientInfo_XGC$X),]
Freq_melt = Freq
Freq_melt$X = rownames(Freq_melt)
Freq_melt = left_join(Freq_melt,PatientInfo_bind,by="X")
Freq_CC = Freq_melt[Freq_melt$Group %in% c("CC"),]
Freq_CC = subset(Freq_CC, select = colnames(Freq_melt)[1:38])
Freq_XGC = Freq_melt[Freq_melt$Group %in% c("XGC"),]
Freq_XGC = subset(Freq_XGC, select = colnames(Freq_melt)[1:38])
Freq_LG = Freq_melt[Freq_melt$Group %in% c("LG"),]
Freq_LG = subset(Freq_LG, select = colnames(Freq_melt)[1:38])
Freq_HG = Freq_melt[Freq_melt$Group %in% c("HG"),]
Freq_HG = subset(Freq_HG, select = colnames(Freq_melt)[1:38])
Freq_adeno = Freq_melt[Freq_melt$Group %in% c("Adeno"),]
Freq_adeno = subset(Freq_adeno, select = colnames(Freq_melt)[1:38])
Freq_adeno_squa = Freq_melt[Freq_melt$Group %in% c("Adeno Squa"),]
Freq_adeno_squa = subset(Freq_adeno_squa, select = colnames(Freq_melt)[1:38])
Freq_neuro = Freq_melt[Freq_melt$Group %in% c("Neuro"),]
Freq_neuro = subset(Freq_neuro, select = colnames(Freq_melt)[1:38])
Freq_squa = Freq_melt[Freq_melt$Group %in% c("Squa"),]
Freq_squa = subset(Freq_squa, select = colnames(Freq_melt)[1:38])
Freq_undiff = Freq_melt[Freq_melt$Group %in% c("Undiff"),]
Freq_undiff = subset(Freq_undiff, select = colnames(Freq_melt)[1:38])
Freq_NTvsT_detail = bind_rows(colSums(Freq_CC),
                              colSums(Freq_XGC),
                              colSums(Freq_LG),
                              colSums(Freq_HG),
                              colSums(Freq_adeno),
                              colSums(Freq_adeno_squa),
                              colSums(Freq_neuro),
                              colSums(Freq_squa),
                              colSums(Freq_undiff))
rownames(Freq_NTvsT_detail) = c("CC", "XGC", "LG", "HG", "Adeno", "Adeno Squa", "Neuro", "Squa", "Undiff")
ZTSubtypes_Counts_Roe_p = Freq_NTvsT_detail
saveRDS(ZTSubtypes_Counts_Roe_p, file = './3 - GBC_output data_subtype/ToLS_Roe/ZTSubtypes_Counts_Roe_p_230605.RDS')

# 6. Roe input of Adeno-Site for subtypes ####
Freq = ZTSubtypes_Counts_Shanno
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_Myeloid = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_Myeloid = subset(PatientInfo_Myeloid, select = c(X,Tumors.for.scRNA.seq.short,histological.type.short))
Freq_melt = Freq
Freq_melt$X = rownames(Freq_melt)
Freq_melt = left_join(Freq_melt, PatientInfo_Myeloid, by = "X")
Freq_melt = subset(Freq_melt, histological.type.short %in% c('adeno'))
Freq_melt = Freq_melt[,-c(39,41)]
colnames(Freq_melt)[39] = "Group"
Freq_P = Freq_melt[Freq_melt$Group %in% c("P"),]
Freq_P = subset(Freq_P, select = colnames(Freq_melt)[1:38])
Freq_PO = Freq_melt[Freq_melt$Group %in% c("PO"),]
Freq_PO = subset(Freq_PO, select = colnames(Freq_melt)[1:38])
Freq_LN = Freq_melt[Freq_melt$Group %in% c("LN"),]
Freq_LN = subset(Freq_LN, select = colnames(Freq_melt)[1:38])
Freq_LI = Freq_melt[Freq_melt$Group %in% c("LI"),]
Freq_LI = subset(Freq_LI, select = colnames(Freq_melt)[1:38])
Freq_LM = Freq_melt[Freq_melt$Group %in% c("LM"),]
Freq_LM = subset(Freq_LM, select = colnames(Freq_melt)[1:38])
Freq_OM = Freq_melt[Freq_melt$Group %in% c("OM"),]
Freq_OM = subset(Freq_OM, select = colnames(Freq_melt)[1:38])
Freq_NTvsT_detail = bind_rows(colSums(Freq_P),
                              colSums(Freq_PO),
                              colSums(Freq_LN),
                              colSums(Freq_LI),
                              colSums(Freq_LM),
                              colSums(Freq_OM))
rownames(Freq_NTvsT_detail) = c("P", "PO", "LN", "LI", "LM", "OM")
ZTSubtypes_Counts_Roe_site = Freq_NTvsT_detail
saveRDS(ZTSubtypes_Counts_Roe_site, file='./3 - GBC_output data_subtype/ToLS_Roe/ZTSubtypes_Counts_Roe_site_230605.RDS')

# 7. Enrichment of myeloid subtypes: NT vs Adeno ####
# * including Fig3D/SFig4D
Myeloid_Subtype_counts = Myeloid_Endo_Subtype_counts[-grep("EC_C", Myeloid_Endo_Subtype_counts$Subtype),]
myeloid_num = readRDS("./3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/myeloid_num_rmdoublets.RDS")
cycle = unique(Myeloid_Subtype_counts$orig.ident)
myeloid_num = myeloid_num[cycle]
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(Myeloid_Subtype_counts, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round((table(sub$Subtype)/myeloid_num[i])*100,2))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0 # Export to WYH
saveRDS(Freq, file = './3 - GBC_output data_subtype/Myeloid_.RDS')
Freq = Freq[,1:33]
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
    stat_pvalue_manual(p_adj, label = "p.adj", hide.ns = F)+
    theme(legend.position = "top") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    ggtitle(cycle[i])
  pdf(file = paste0("./3 - GBC_output data_subtype/05_MyeloidSubtype_Enrichment_NTvsAdenoP/AddMastVersion_", cycle[i], ".pdf"), width =2, height = 3)
  print(p)
  dev.off()
}

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

# 9. Enrichment of myeloid subtypes: P vs P_Mets ####
# * including Fig2I/Fig3E/SFig4E
Myeloid_Subtype_counts = Myeloid_Endo_Subtype_counts[-grep("EC_C", Myeloid_Endo_Subtype_counts$Subtype),]
myeloid_num = readRDS("./3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/myeloid_num_rmdoublets.RDS")
cycle = unique(Myeloid_Subtype_counts$orig.ident)
myeloid_num = myeloid_num[cycle]
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(Myeloid_Subtype_counts, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round((table(sub$Subtype)/myeloid_num[i])*100,2))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0
Freq = Freq[,1:33]
Freq_melt = Freq
Freq_melt$Patient = rownames(Freq_melt)
Freq_melt$Patient = factor(Freq_melt$Patient, levels = Freq_melt$Patient)
Freq_melt <- reshape2::melt(Freq_melt, id.vars=c("Patient"), variable.name="Subtype", value.name="Percentage")
Freq_melt = Freq_melt[order(Freq_melt$Subtype),]
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_DC = subset(PatientInfo_DC, select = c("NewSample.ID","metastasis.type","histological.type.short","Tumors.for.scRNA.seq.short"))
PatientInfo_DC = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c("adeno"),]
PatientInfo_DC = PatientInfo_DC[PatientInfo_DC$Tumors.for.scRNA.seq.short %in% c("P"),]
PatientInfo_DC$Group = PatientInfo_DC$metastasis.type
PatientInfo_DC$Group = ifelse(PatientInfo_DC$Group == "P_LI", "P", PatientInfo_DC$Group)
PatientInfo_DC$Group = ifelse(PatientInfo_DC$Group %in% c("P_LN","P_LM"), "P_Mets", PatientInfo_DC$Group)
PatientInfo_DC$metastasis.type = PatientInfo_DC$Group
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
  pdf(file = paste0("./3 - GBC_output data_subtype/07_MyeloidSubtype_Enrichment_AdenoPs/AddMastVersion_Pmet_", cycle[i], ".pdf"), width =2, height = 3)
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

# 11. Prognosis of myeloid subtypes in Adeno ####
# * including Fig2G/Fig2M
Myeloid_Subtype_counts = Myeloid_Endo_Subtype_counts[-grep("EC_C", Myeloid_Endo_Subtype_counts$Subtype),]
myeloid_num = readRDS("./3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/myeloid_num_rmdoublets.RDS")
cycle = unique(Myeloid_Subtype_counts$orig.ident)
myeloid_num = myeloid_num[cycle]
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(Myeloid_Subtype_counts, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round((table(sub$Subtype)/myeloid_num[i])*100,2))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0
Freq = Freq[,1:29]
Myeloid_Neu_Freq = Freq
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20240916预后更新.xlsx")
PatientInfo_Neu = PatientInfo[PatientInfo$NewSample.ID %in% unique(rownames(Freq)),]
PatientInfo_Neu = left_join(PatientInfo_Neu,PatientInfo_surv,by="NewSample.ID")
Myeloid_Neu_Freq_surv = Myeloid_Neu_Freq
Myeloid_Neu_Freq_surv$NewSample.ID = rownames(Myeloid_Neu_Freq_surv)
PatientInfo_Neu = left_join(PatientInfo_Neu,Myeloid_Neu_Freq_surv,by="NewSample.ID")
Surv_data_Neu = subset(PatientInfo_Neu, histological.type.short %in% "adeno")
Surv_data_Neu = subset(Surv_data_Neu, Tumors.for.scRNA.seq.short %in% "P")
Surv_data = Surv_data_Neu
Surv_data = Surv_data[!is.na(Surv_data$event),]
Surv_data$event01 = ifelse(Surv_data$event=="dead",1,0)
Surv_data = subset(Surv_data, select = c("OS_month","event01",as.character(sort(unique(Myeloid_Subtype_counts$Subtype)))))
cycle = as.character(sort(unique(Myeloid_Subtype_counts$Subtype)))
for (i in 1:length(cycle)) {
  temp = subset(Surv_data, select = c("OS_month","event01",cycle[i]))
  temp[,3] = as.numeric(temp[,3])
  temp[,3] = ifelse(temp[,3]>median(temp[,3]),"High","Low")
  fit<-survfit(Surv(OS_month, event01)~temp[,3], data=temp)
  p = ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = temp,             # data used to fit survival curves.
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    #conf.int = TRUE,         # show confidence intervals for
    palette = "npg",
    xlab = "Time in months",   # customize X axis label.
    #ggtheme = theme_bw(),
    #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
    #surv.median.line = "hv",  # add the median survival pointer.
    # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
    tables.y.text = T,
    risk.table.pos = "in",
    risk.table.col = "strata",
    fontsize = 4,
    pval.size = 4,
    #surv.plot.height = 0.8,
    #tables.height = 0.2,
    pval.coord = c(9, 0.9),
    legend = "top",
    title = cycle[i]
  )
  pdf(file = paste0("./3 - GBC_output data_subtype/09_MyeloidSubtype_Prognosis_AdenoP/241218_", cycle[i], ".pdf"), width =6, height = 6, onefile = F)
  print(p)
  dev.off()
}

# 11. Prognosis of Endo subtypes in Adeno ####
Myeloid_Subtype_counts = Myeloid_Endo_Subtype_counts[grep("EC_C", Myeloid_Endo_Subtype_counts$Subtype),]
cycle = unique(Myeloid_Subtype_counts$orig.ident)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(Myeloid_Subtype_counts, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round((table(sub$Subtype)/nrow(sub))*100,2))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0
Freq = Freq[,30:38]
Myeloid_Neu_Freq = Freq
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20240916预后更新.xlsx")
PatientInfo_Neu = PatientInfo[PatientInfo$NewSample.ID %in% unique(rownames(Freq)),]
PatientInfo_Neu = left_join(PatientInfo_Neu,PatientInfo_surv,by="NewSample.ID")
Myeloid_Neu_Freq_surv = Myeloid_Neu_Freq
Myeloid_Neu_Freq_surv$NewSample.ID = rownames(Myeloid_Neu_Freq_surv)
PatientInfo_Neu = left_join(PatientInfo_Neu,Myeloid_Neu_Freq_surv,by="NewSample.ID")
Surv_data_Neu = subset(PatientInfo_Neu, histological.type.short %in% "adeno")
Surv_data_Neu = subset(Surv_data_Neu, Tumors.for.scRNA.seq.short %in% "P")
Surv_data = Surv_data_Neu
Surv_data = Surv_data[!is.na(Surv_data$event),]
Surv_data$event01 = ifelse(Surv_data$event=="dead",1,0)
Surv_data = subset(Surv_data, select = c("OS_month","event01",as.character(sort(unique(Myeloid_Subtype_counts$Subtype)))))
cycle = as.character(sort(unique(Myeloid_Subtype_counts$Subtype)))
for (i in 1:length(cycle)) {
  temp = subset(Surv_data, select = c("OS_month","event01",cycle[i]))
  temp[,3] = as.numeric(temp[,3])
  temp[,3] = ifelse(temp[,3]>median(temp[,3]),"High","Low")
  fit<-survfit(Surv(OS_month, event01)~temp[,3], data=temp)
  p = ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = temp,             # data used to fit survival curves.
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    #conf.int = TRUE,         # show confidence intervals for
    palette = "npg",
    xlab = "Time in months",   # customize X axis label.
    #ggtheme = theme_bw(),
    #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
    #surv.median.line = "hv",  # add the median survival pointer.
    # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
    tables.y.text = T,
    risk.table.pos = "in",
    risk.table.col = "strata",
    fontsize = 4,
    pval.size = 4,
    #surv.plot.height = 0.8,
    #tables.height = 0.2,
    pval.coord = c(9, 0.9),
    legend = "top",
    title = cycle[i]
  )
  pdf(file = paste0("./3 - GBC_output data_subtype/09_MyeloidSubtype_Prognosis_AdenoP/241218_", cycle[i], ".pdf"), width = 6, height = 6, onefile = F)
  print(p)
  dev.off()
}

# 12. cellchat input ####
Mast <- readRDS("./3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/Myeloid_cell_Mast.RDS")
Mast$Subtype = factor("Mast")
Myeloid_cell_Neu$Subtype = Idents(Myeloid_cell_Neu)
Myeloid_cell_DC$Subtype = Idents(Myeloid_cell_DC)
Myeloid_cell_MonoMacro$Subtype = Idents(Myeloid_cell_MonoMacro)
Endothelium$Subtype = Idents(Endothelium)
Endothelium$Myeloid_celltype = factor("Endothelium")
cellchat_zt = rbind(subset(Mast@meta.data, select = c("Subtype","Myeloid_celltype")),
                    subset(Myeloid_cell_Neu@meta.data, select = c("Subtype","Myeloid_celltype")),
                    subset(Myeloid_cell_DC@meta.data, select = c("Subtype","Myeloid_celltype")),
                    subset(Myeloid_cell_MonoMacro@meta.data, select = c("Subtype","Myeloid_celltype")),
                    subset(Endothelium@meta.data, select = c("Subtype","Myeloid_celltype")))
cellchat_zt$cellid = rownames(cellchat_zt)
colnames(cellchat_zt)[1:2] = c("subtype","celltype")
cellchat_zt = cellchat_zt[,c("cellid","subtype","celltype")]
cellchat_zt$celltype = as.character(cellchat_zt$celltype)
saveRDS(cellchat_zt, file = "./01_GBC_subtypes_zt/cellchat_zt_230605.RDS")

# 13. Monocyte/Macrophage ####
## Fig1E ####
DimPlot(Myeloid_cell_MonoMacro, reduction = "umap", pt.size = 2, raster = T) + ggplot2::theme(legend.position = "right", panel.background = element_blank())

## FigS4A ####
Myeloid_cell_MonoMacro_DEGs_allpos <- FindAllMarkers(Myeloid_cell_MonoMacro, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Myeloid_cell_MonoMacro_DEGs_allpos = Myeloid_cell_MonoMacro_DEGs_allpos[Myeloid_cell_MonoMacro_DEGs_allpos$p_val_adj < 0.05,]
Myeloid_cell_MonoMacro_DEGs_allpos %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
p = DotPlot(Myeloid_cell_MonoMacro, features = c(unique(top5$gene),"CD14","FCGR3A","CD68","CD163"))  +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
pdf(file = "./3 - GBC_output data_subtype/03_MM/MM_SubtypeDEGs_Top5.pdf", width =15, height = 4)
print(p)
dev.off()

Myeloid_cell_MonoMacro = subset(Myeloid_cell_MonoMacro, idents = c("M_C6_CCL18","M_C8_MMP9"), invert=T)
Myeloid_cell_MonoMacro_DEGs_allpos <- FindAllMarkers(Myeloid_cell_MonoMacro, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Myeloid_cell_MonoMacro_DEGs_allpos = Myeloid_cell_MonoMacro_DEGs_allpos[Myeloid_cell_MonoMacro_DEGs_allpos$p_val_adj < 0.05,]
write_csv(Myeloid_cell_MonoMacro_DEGs_allpos, file = "./3 - GBC_output data_subtype/03_MM/DEGs among common mm subtypes_v240422.csv")

Myeloid_cell_MonoMacro_DEGs_allpos %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DotPlot(Myeloid_cell_MonoMacro, features = c(unique(top5$gene),"CD14","FCGR3A","CD68","CD163"))  +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

## * Subset subgroups ####
Myeloid_cell_MonoMacro = subset(Myeloid_cell_MonoMacro, idents = c("M_C7_AGR2"), invert=T)

## FigS4B ####
IFN_TAM_Sig = list(c("CASP1", "CASP4", "CCL2", "CCL3", "CCL4", "CCL7", "CCL8", "CD274", "CD40", 
                     "CXCL2", "CXCL3", "CXCL9", "CXCL10", "CXCL11", "IDO1", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", 
                     "IRF1", "IRF7", "ISG15", "LAMP3", "PDCD1LG2", "TNFSF10", "C1QA", "C1QC", "CD38", "IL4I1", "IFI44L"))
Inflam_TAM_Sig = list(c("CCL2", "CCL3", "CCL4", "CCL5", "CCL20", "CCL3L1", "CCL3L3", "CCL4L2", "CCL4L4", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", 
                        "G0S2", "IL1B", "IL1RN", "IL6", "INHBA", "KLF2", "KLF6", "NEDD9", "PMAIP1", "S100A8", "S100A9", "SPP1"))
LA_TAM_Sig = list(c("ACP5", "AOPE", "APOC1", "ATF1", "C1QA", "C1QB", "C1QC", "CCL18", "CD163", "CD36", "CD63", "CHI3L1", "CTSB", "CTSD", "CTSL", 
                    "F13A1", "FABP5", "FOLR2", "GPNMB", "IRF3", "LGALS3", "LIPA", "LPL", "MACRO", "MerTK", "MMP7", "MMP9", "MMP12", "MRC1", 
                    "NR1H3", "NRF1", "NUPR1", "PLA2G7", "RNASE1", "SPARC", "SPP1", "TFDP2", "TREM2", "ZEB1"))
Angio_TAM_Sig = list(c("ADAM8", "AREG", "BNIP3", "CCL2", "CCL4", "CCL20", "CD163", "CD300E", "CD44", "CD55", "CEBPB", "CLEC5A", "CTSB", 
                       "EREG", "FCN1", "FLT1", "FN1", "HES1", "IL1B", "IL1RN", "IL8", "MAF", "MIF", "NR1H3", "OLR1", "PPARG", "S100A8", "S100A9", "S100A12", 
                       "SERPINB2", "SLC2A1", "SPIC", "SPP1", "THBS1", "TIMP1", "VCAN", "VEGFA"))
Reg_TAM_Sig = list(c("CCL2", "CD274", "CD40", "CD80", "CD86", "CHIT1", "CX3CR1", "HLA-A", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", 
                     "ICOSLG", "IL-10", "ITGA4", "LGALS9", "MACRO", "MRC1", "TGFB2"))
Prolif_TAM_Sig = list(c("CCNA2", "CDC45", "CDK1", "H2AFC", "HIST1H4C", "HMGB1", "HMGN2", "MKI67", "RRM2", "STMN1", "TOP2A", "TUBA1B", "TUBB", "TYMS"))

IFN_TAM_Sig_genes = as.data.frame(unlist(IFN_TAM_Sig))
write_csv(IFN_TAM_Sig_genes, file = "./3 - GBC_output data_subtype/03_MM/IFN_TAM_Sig_genes_v240422.csv")

Inflam_TAM_Sig_genes = as.data.frame(unlist(Inflam_TAM_Sig))
write_csv(Inflam_TAM_Sig_genes, file = "./3 - GBC_output data_subtype/03_MM/Inflam_TAM_Sig_genes_v240422.csv")

LA_TAM_Sig_genes = as.data.frame(unlist(LA_TAM_Sig))
write_csv(LA_TAM_Sig_genes, file = "./3 - GBC_output data_subtype/03_MM/LA_TAM_Sig_genes_v240422.csv")

Angio_TAM_Sig_genes = as.data.frame(unlist(Angio_TAM_Sig))
write_csv(Angio_TAM_Sig_genes, file = "./3 - GBC_output data_subtype/03_MM/Angio_TAM_Sig_genes_v240422.csv")

Reg_TAM_Sig_genes = as.data.frame(unlist(Reg_TAM_Sig))
write_csv(Reg_TAM_Sig_genes, file = "./3 - GBC_output data_subtype/03_MM/Reg_TAM_Sig_genes_v240422.csv")

Prolif_TAM_Sig_genes = as.data.frame(unlist(Prolif_TAM_Sig))
write_csv(Prolif_TAM_Sig_genes, file = "./3 - GBC_output data_subtype/03_MM/Prolif_TAM_Sig_genes_v240422.csv")

Myeloid_cell_MonoMacro <- AddModuleScore(Myeloid_cell_MonoMacro, features = IFN_TAM_Sig, name = "IFN_TAM_Sig")
Myeloid_cell_MonoMacro <- AddModuleScore(Myeloid_cell_MonoMacro, features = Inflam_TAM_Sig, name = "Inflam_TAM_Sig")
Myeloid_cell_MonoMacro <- AddModuleScore(Myeloid_cell_MonoMacro, features = LA_TAM_Sig, name = "LA_TAM_Sig")
Myeloid_cell_MonoMacro <- AddModuleScore(Myeloid_cell_MonoMacro, features = Angio_TAM_Sig, name = "Angio_TAM_Sig")
Myeloid_cell_MonoMacro <- AddModuleScore(Myeloid_cell_MonoMacro, features = Reg_TAM_Sig, name = "Reg_TAM_Sig")
Myeloid_cell_MonoMacro <- AddModuleScore(Myeloid_cell_MonoMacro, features = Prolif_TAM_Sig, name = "Prolif_TAM_Sig")
p = DotPlot(Myeloid_cell_MonoMacro, features = c("IFN_TAM_Sig1","Inflam_TAM_Sig1","LA_TAM_Sig1","Angio_TAM_Sig1","Reg_TAM_Sig1","Prolif_TAM_Sig1"))  +
  # coord_flip() +
  scale_colour_gradient2(low = "#4287f5", mid = "white", high = "#f5424b")+
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
pdf(file = paste0("./3 - GBC_output data_subtype/03_MM/TAM_Sig_Score", ".pdf"), width =7.5, height = 3.5) #
print(p)
dev.off()

## * Subset subgroups ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_MonoMacro@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_MonoMacro@meta.data$NewSample.ID = Myeloid_cell_MonoMacro@meta.data$orig.ident
temp = left_join(Myeloid_cell_MonoMacro@meta.data, PatientInfo, by = "NewSample.ID")
identical(Myeloid_cell_MonoMacro@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_MonoMacro@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_MonoMacro@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_MonoMacro@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_MonoMacro@meta.data$Clinical.stage = temp$Clinical.stage.x
M0 = subset(Myeloid_cell_MonoMacro, idents = "M_C0_FOLR2")
M1 = subset(Myeloid_cell_MonoMacro, idents = "M_C1_S100A8")
M2 = subset(Myeloid_cell_MonoMacro, idents = "M_C2_SPP1")
M3 = subset(Myeloid_cell_MonoMacro, idents = "M_C3_FCGBP")

## FigS4F ####
Myeloid_cell_MonoMacro <- readRDS("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/03_MM/Myeloid_cell_MonoMacro.RDS")
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_MonoMacro@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_MonoMacro@meta.data$NewSample.ID = Myeloid_cell_MonoMacro@meta.data$orig.ident
temp = left_join(Myeloid_cell_MonoMacro@meta.data, PatientInfo, by = "NewSample.ID")
identical(Myeloid_cell_MonoMacro@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_MonoMacro@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_MonoMacro@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_MonoMacro@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_MonoMacro@meta.data$Clinical.stage = temp$Clinical.stage.x
M1 = subset(Myeloid_cell_MonoMacro, idents = "M_C1_S100A8")
M1 = subset(M1, histological.type.short %in% c("adeno","CC","HG","LG","XGC"))
M1 = subset(M1, Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"))
M1$Group = ifelse(M1$histological.type.short %in% c("adeno"), "Adeno", "Non-tumor")
Idents(M1) = M1$Group
## DEGs: Volcano
Endothelium_0_P_DEGs <- FindAllMarkers(M1, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[Endothelium_0_P_DEGs$p_val_adj < 0.05,]
Endothelium_0_P_DEGs$gene = rownames(Endothelium_0_P_DEGs)
temp <- bitr(Endothelium_0_P_DEGs$gene,fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb='org.Hs.eg.db')
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[temp$SYMBOL,]
Endothelium_0_P_DEGs$gene = temp$ENTREZID
formula_res <- compareCluster(gene~cluster, data=Endothelium_0_P_DEGs,
                              fun="enrichGO",
                              OrgDb = org.Hs.eg.db,
                              ont = "BP" ,
                              pAdjustMethod = "BH",
                              readabl = TRUE)
temp = formula_res@compareClusterResult
write.csv(temp, file = "./3 - GBC_output data_subtype/03_MM/M1 GOBP_NTvsAdeno_v240422.csv")

dotplot(formula_res, showCategory = 10) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "horizontal")

formula_res1 <- filter(formula_res, Description %in% c("positive regulation of cytokine production",
                                                       "antigen processing and presentation",
                                                       "response to hypoxia",
                                                       "regulation of angiogenesis",
                                                       "phagocytosis",
                                                       "regulation of inflammatory response",
                                                       "cell chemotaxis"))
dotplot(formula_res1, showCategory = 10) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "horizontal")

## FigS4H ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_MonoMacro@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Adeno = PatientInfo[PatientInfo$histological.type.short %in% "adeno",]
Adeno = Adeno[Adeno$Tumors.for.scRNA.seq.short %in% c("P","LN"),]
Adeno = Adeno[Adeno$病历号 %in% Adeno$病历号[duplicated(Adeno$病历号)],]
Freq_Cor = Freq[Adeno$NewSample.ID,]
Freq_Cor = subset(Freq_Cor, select = "M_C1_S100A8")
Freq_Cor1 = data.frame(P = Freq_Cor$M_C1_S100A8[grep("_P",rownames(Freq_Cor))],
                       LN = Freq_Cor$M_C1_S100A8[grep("_LN",rownames(Freq_Cor))])
Freq_Cor_all = data.frame(P = c(Freq_Cor1$P.Freq, Freq_Cor2$P.Freq),
                          Mets = c(Freq_Cor1$LN.Freq, Freq_Cor2$LM.Freq))
p = ggscatter(Freq_Cor1, x = 'P.Freq', y = 'LN.Freq',
              fill ='#8491B4FF',color='darkgrey',size = 2,shape = 21, # Points color, shape and size
              add = "reg.line",# Add regressin line
              add.params = list(color = "#7E6148FF", fill = "#B09C85FF"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = T, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson", label.x.npc = "left", label.y.npc = "top",label.sep = "\n"))+xlab('P.Freq')+ylab('LN.Freq')
pdf(file = "./3 - GBC_output data_subtype/03_MM/M1_Freq_Cor_P_LN.pdf", width =3, height = 3)
print(p)
dev.off()

## FigS4I ####
MM = M3
MM$metastasis.type = ifelse(MM$metastasis.type %in% c("P_LM","P_LN"),"P_Met",MM$metastasis.type)
MM$metastasis.type = ifelse(MM$metastasis.type %in% c("P_LI"),"P",MM$metastasis.type)
Idents(MM) = MM$metastasis.type
## DEGs: Volcano
MM_NTvsT_DEGs <- FindMarkers(MM,ident.1 = "P_Met", ident.2 = "P",only.pos = F,min.pct = 0.25, logfc.threshold = 0.25)
MM_NTvsT_DEGs = MM_NTvsT_DEGs[MM_NTvsT_DEGs$p_val_adj < 0.05,]
write.csv(MM_NTvsT_DEGs, file = "./3 - GBC_output data_subtype/03_MM/M3 PvsPMets_v240422.csv")

Dat = MM_NTvsT_DEGs
Dat$Gene = rownames(Dat)
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) > 0, ifelse(Dat$avg_log2FC > 0 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
p=ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
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
  ggtitle("M3")
pdf(file = paste0("./3 - GBC_output data_subtype/03_MM/PvsPmet_DEG_M3", ".pdf"), width = 6.5, height = 5) #
print(p)
dev.off()

## Fig2K ####
Myeloid_cell_MonoMacro <- readRDS("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/03_MM/Myeloid_cell_MonoMacro.RDS")
Myeloid_cell_MonoMacro = subset(Myeloid_cell_MonoMacro, idents = c("M_C7_AGR2"), invert=T)
Myeloid_cell_MonoMacro$Subtype = Idents(Myeloid_cell_MonoMacro)
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_MonoMacro@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_MonoMacro@meta.data$NewSample.ID = Myeloid_cell_MonoMacro@meta.data$orig.ident
temp = left_join(Myeloid_cell_MonoMacro@meta.data, PatientInfo, by = "NewSample.ID")
identical(Myeloid_cell_MonoMacro@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_MonoMacro@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_MonoMacro@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_MonoMacro@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_MonoMacro@meta.data$Clinical.stage = temp$Clinical.stage.x
Idents(Myeloid_cell_MonoMacro) = Myeloid_cell_MonoMacro@meta.data$histological.type.short
Myeloid_cell_MonoMacro = subset(Myeloid_cell_MonoMacro, idents = "adeno")
Idents(Myeloid_cell_MonoMacro) = Myeloid_cell_MonoMacro@meta.data$Tumors.for.scRNA.seq.short
Myeloid_cell_MonoMacro = subset(Myeloid_cell_MonoMacro, idents = "P")
Myeloid_cell_MonoMacro@meta.data$metastasis.type = ifelse(Myeloid_cell_MonoMacro@meta.data$metastasis.type %in% "P_LI", "P", Myeloid_cell_MonoMacro@meta.data$metastasis.type)
Myeloid_cell_MonoMacro@meta.data$metastasis.type = ifelse(Myeloid_cell_MonoMacro@meta.data$metastasis.type %in% c("P_LN","P_LM"), "P_Met", Myeloid_cell_MonoMacro@meta.data$metastasis.type)
Idents(Myeloid_cell_MonoMacro) = Myeloid_cell_MonoMacro$Subtype
Myeloid_cell_MonoMacro = subset(Myeloid_cell_MonoMacro, idents = c("M_C0_FOLR2","M_C1_S100A8","M_C2_SPP1"))
Idents(Myeloid_cell_MonoMacro) = Myeloid_cell_MonoMacro@meta.data$metastasis.type
## M0
vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(subset(Myeloid_cell_MonoMacro, Subtype == "M_C0_FOLR2"), features = signature, pt.size = 0, 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + NoLegend() +
      stat_compare_means(comparisons = test_sign, label = "p") + 
      stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
  }
  map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = ., ncol = 3)
  file_name <- paste0(file_name, ".pdf")
  ggsave(file_name, width = 6, height = 3)
}
gene_sig <- c("CXCL2", "CXCL3", "CXCL8")
comparisons <- list(c("P", "P_Met"))
vp_case1(gene_signature = gene_sig, file_name = "../2 - Analysis_V3_checked/3 - GBC_output data_subtype/03_MM/Cellchat_Ligand_Exprs_M0", test_sign = comparisons, y_max = 7)
## M2
vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(subset(Myeloid_cell_MonoMacro, Subtype == "M_C2_SPP1"), features = signature, pt.size = 0, 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + NoLegend() +
      stat_compare_means(comparisons = test_sign, label = "p") + 
      stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
  }
  map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = ., ncol = 3)
  file_name <- paste0(file_name, ".pdf")
  ggsave(file_name, width = 6, height = 3)
}
gene_sig <- c("CXCL2", "CXCL3", "CXCL8")
comparisons <- list(c("P", "P_Met"))
vp_case1(gene_signature = gene_sig, file_name = "../2 - Analysis_V3_checked/3 - GBC_output data_subtype/03_MM/Cellchat_Ligand_Exprs_M2", test_sign = comparisons, y_max = 8.5)

## ####
## Fig2L ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_MonoMacro@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_MonoMacro@meta.data$NewSample.ID = Myeloid_cell_MonoMacro@meta.data$orig.ident
temp = left_join(Myeloid_cell_MonoMacro@meta.data, PatientInfo, by = "NewSample.ID")
identical(Myeloid_cell_MonoMacro@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_MonoMacro@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_MonoMacro@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_MonoMacro@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_MonoMacro@meta.data$Clinical.stage = temp$Clinical.stage.x
M0 = subset(Myeloid_cell_MonoMacro, idents = "M_C0_FOLR2")
M1 = subset(Myeloid_cell_MonoMacro, idents = "M_C1_S100A8")
M2 = subset(Myeloid_cell_MonoMacro, idents = "M_C2_SPP1")
expr = M0
df = FetchData(expr,vars = c('CXCL2','CXCL3','CXCL8','Clinical.stage','orig.ident','metastasis.type','Tumors.for.scRNA.seq.short','histological.type.short'))
df = df[df$histological.type.short %in% c("adeno","CC","XGC","HG","LG"),]
df = df[df$Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"),]
df = df[!(df$histological.type.short %in% c("CC","XGC","HG","LG")),]
#df$Clinical.stage = ifelse(df$histological.type.short %in% c("CC","XGC","HG","LG"), "NT", df$Clinical.stage)
df$Clinical.stage = ifelse(df$Clinical.stage %in% c('IIA','IIB'),'II',df$Clinical.stage)

final1= df %>% group_by(Clinical.stage) %>% summarise(mean = mean(CXCL2),N=length(CXCL2),
                                                      sd = sd(CXCL2),se=sd/sqrt(N))
final1$group = "CXCL2"
final2= df %>% group_by(Clinical.stage) %>% summarise(mean = mean(CXCL3),N=length(CXCL3),
                                                      sd = sd(CXCL3),se=sd/sqrt(N))
final2$group = "CXCL3"
final3= df %>% group_by(Clinical.stage) %>% summarise(mean = mean(CXCL8),N=length(CXCL8),
                                                      sd = sd(CXCL8),se=sd/sqrt(N))
final3$group = "CXCL8"
final = rbind(final1,final2,final3)
final$Clinical.stage = factor(final$Clinical.stage, levels = c("I","II","IIIA","IIIB","IVA","IVB"))
p=ggplot(final,aes(Clinical.stage,mean,group = group))+
  geom_ribbon(aes(ymin=mean-se,ymax=mean+se),fill='#e3e3e3')+
  #geom_line(color=c('#B7916C'),size=1)+
  #geom_point(color='#B7916C',size=2)+
  geom_line(aes(color=group),size=1)+
  geom_point(color='#B7916C',size=2)+
  facet_grid(.~group)+
  theme(panel.background = element_blank(),
        panel.grid = element_line(color="white"),
        # axis.title = element_blank(),
        #axis.ticks = element_blank(),
        axis.text.x=element_text(colour='black',size=12,angle = 45,hjust = 1))+
  ggtitle('M0')+theme(axis.line=element_line(colour='black',size=1,lineend = 'square'))+
  theme (axis.text.x = element_text (colour='black', size=12, angle=45), 
         axis.text.y = element_text (colour='black',size=12, angle=45))+
  labs(x='Stage',y='Score')
pdf(file = paste0("./3 - GBC_output data_subtype/03_MM/Cellchat_Ligand_Dynamic_M0", ".pdf"), width =7, height = 3.5) #
print(p)
dev.off()
compare_means(CXCL2~Clinical.stage, df, method = "anova") # <2e-16
compare_means(CXCL3~Clinical.stage, df, method = "anova") # <2e-16
compare_means(CXCL8~Clinical.stage, df, method = "anova") # <2e-16

expr = M2
df = FetchData(expr,vars = c('CXCL2','CXCL3','CXCL8','Clinical.stage','orig.ident','metastasis.type','Tumors.for.scRNA.seq.short','histological.type.short'))
df = df[df$histological.type.short %in% c("adeno","CC","XGC","HG","LG"),]
df = df[df$Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"),]
df = df[!(df$histological.type.short %in% c("CC","XGC","HG","LG")),]
#df$Clinical.stage = ifelse(df$histological.type.short %in% c("CC","XGC","HG","LG"), "NT", df$Clinical.stage)
df$Clinical.stage = ifelse(df$Clinical.stage %in% c('IIA','IIB'),'II',df$Clinical.stage)
final1= df %>% group_by(Clinical.stage) %>% summarise(mean = mean(CXCL2),N=length(CXCL2),
                                                      sd = sd(CXCL2),se=sd/sqrt(N))
final1$group = "CXCL2"
final2= df %>% group_by(Clinical.stage) %>% summarise(mean = mean(CXCL3),N=length(CXCL3),
                                                      sd = sd(CXCL3),se=sd/sqrt(N))
final2$group = "CXCL3"
final3= df %>% group_by(Clinical.stage) %>% summarise(mean = mean(CXCL8),N=length(CXCL8),
                                                      sd = sd(CXCL8),se=sd/sqrt(N))
final3$group = "CXCL8"
final = rbind(final1,final2,final3)
final$Clinical.stage = factor(final$Clinical.stage, levels = c("I","II","IIIA","IIIB","IVA","IVB"))
p=ggplot(final,aes(Clinical.stage,mean,group = group))+
  geom_ribbon(aes(ymin=mean-se,ymax=mean+se),fill='#e3e3e3')+
  #geom_line(color=c('#B7916C'),size=1)+
  #geom_point(color='#B7916C',size=2)+
  geom_line(aes(color=group),size=1)+
  geom_point(color='#B7916C',size=2)+
  facet_grid(.~group)+
  theme(panel.background = element_blank(),
        panel.grid = element_line(color="white"),
        # axis.title = element_blank(),
        #axis.ticks = element_blank(),
        axis.text.x=element_text(colour='black',size=12,angle = 45,hjust = 1))+
  ggtitle('M2')+theme(axis.line=element_line(colour='black',size=1,lineend = 'square'))+
  theme (axis.text.x = element_text (colour='black', size=12, angle=45), 
         axis.text.y = element_text (colour='black',size=12, angle=45))+
  labs(x='Stage',y='Score')
pdf(file = paste0("./3 - GBC_output data_subtype/03_MM/Cellchat_Ligand_Dynamic_M2", ".pdf"), width =7, height = 3.5) #
print(p)
dev.off()

compare_means(CXCL2~Clinical.stage, df, method = "anova") # <2e-16
compare_means(CXCL3~Clinical.stage, df, method = "anova") # <2e-16
compare_means(CXCL8~Clinical.stage, df, method = "anova") # <2e-16

# 14. Neutrophil ####
## Fig1E ####
p = DimPlot(Myeloid_cell_Neu, reduction = "umap", pt.size = 2, raster = T, cols = c("N_C0_MNDA" = mycolor[1],
                                                                                         "N_C1_ISG15" = mycolor[2],
                                                                                         "N_C2_CCL3L1" = mycolor[3],
                                                                                         "N_C3_CTSD" = mycolor[4],
                                                                                         "N_C4_CD44" = mycolor[5],
                                                                                         "N_C5_DEFA5" = mycolor[6],
                                                                                         "N_C6_ISG15_IFIT1" = mycolor[7],
                                                                                         "N_C7_MNDA_S100A8" = mycolor[8],
                                                                                         "N_C8_IL1R2" = mycolor[9],
                                                                                         "N_C9_S100A10" = mycolor[10],
                                                                                         "N_C10_STXBP2" = mycolor[11],
                                                                                         "N_C11_MGP" = mycolor[12])) + ggplot2::theme(legend.position = "right", 
                                                                                                                                            panel.background = element_blank())
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/Neutrophil_SubtypeDefinition.pdf", width =5, height = 2.7)
print(p)
dev.off()

# * Subset common subgroups ####
Myeloid_cell_Neu = subset(Myeloid_cell_Neu, idents = c("N_C0_MNDA","N_C1_ISG15","N_C2_CCL3L1","N_C7_MNDA_S100A8","N_C9_S100A10"))
saveRDS(Myeloid_cell_Neu, file = "./3 - GBC_output data_subtype/01_Neutrophil/Myeloid_cell_Neu.RDS")

## Fig2A ####
Myeloid_cell_Neu_DEGs_allpos <- FindAllMarkers(Myeloid_cell_Neu, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Myeloid_cell_Neu_DEGs_allpos = Myeloid_cell_Neu_DEGs_allpos[Myeloid_cell_Neu_DEGs_allpos$p_val_adj < 0.05,]
write.csv(Myeloid_cell_Neu_DEGs_allpos, file = "./3 - GBC_output data_subtype/01_Neutrophil/DEGs among common neutrophil subtypes_v240422.csv")

Myeloid_cell_Neu_DEGs_allpos %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
p = DotPlot(Myeloid_cell_Neu, features = top5$gene)  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/Neutrophil_SubtypeDEGs_Top5.pdf", width =12, height = 2.5)
print(p)
dev.off()

## Fig2B ####
Neu_sig = read.csv("./3 - GBC_output data_subtype/01_Neutrophil/Liver Cancer Atlas_Neu_Signature.csv")
Neu_sig= Neu_sig[Neu_sig$p_val_adj<0.05,]
Neu_sig = Neu_sig[Neu_sig$Gene %in% VariableFeatures(Myeloid_cell_Neu),]
TAN_sig = Neu_sig[Neu_sig$Cell.cluster %in% c("Neu_01_MMP8","Neu_07_APOA2","Neu_08_CD74","Neu_09_IFIT1","Neu_10_SPP1","Neu_11_CCL4"),]$Gene
ALN_sig = Neu_sig[Neu_sig$Cell.cluster %in% c("Neu_05_ELL2","Neu_06_PTGS2"),]$Gene
PBN_sig = Neu_sig[Neu_sig$Cell.cluster %in% c("Neu_02_S100A12","Neu_03_ISG15","Neu_04_TXNIP"),]$Gene
TAN_sig = list(TAN_sig)
ALN_sig = list(ALN_sig)
PBN_sig = list(PBN_sig)

temp = as.data.frame(unlist(TAN_sig))
write.csv(temp, file = "./3 - GBC_output data_subtype/01_Neutrophil/TAN_Sig_v240422.csv")

temp = as.data.frame(unlist(ALN_sig))
write.csv(temp, file = "./3 - GBC_output data_subtype/01_Neutrophil/ALN_sig_v240422.csv")

temp = as.data.frame(unlist(PBN_sig))
write.csv(temp, file = "./3 - GBC_output data_subtype/01_Neutrophil/PBN_sig_v240422.csv")

Myeloid_cell_Neu <- AddModuleScore(Myeloid_cell_Neu, features = PBN_sig, name = "PBN_sig")
Myeloid_cell_Neu <- AddModuleScore(Myeloid_cell_Neu, features = ALN_sig, name = "ALN_sig")
Myeloid_cell_Neu <- AddModuleScore(Myeloid_cell_Neu, features = TAN_sig, name = "TAN_sig")
p = DotPlot(Myeloid_cell_Neu, features = c("PBN_sig1","ALN_sig1","TAN_sig1"))  +
  # coord_flip() +
  scale_colour_gradient2(low = "#4287f5", mid = "white", high = "#f5424b")+
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/Neutrophil_Subtype_addmodule_TAN.pdf", width =7, height = 2.5)
print(p)
dev.off()

## Fig2C & Fig2D ####
Myeloid_cell_Neu$Subtype = Idents(Myeloid_cell_Neu)
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_Neu = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_Neu@meta.data$orig.ident),]
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo_Neu = left_join(PatientInfo_Neu,PatientInfo_surv,by="NewSample.ID")
PatientInfo_Neu$Clinical.stage.y[PatientInfo_Neu$Clinical.stage.y %in% c("IIA","IIB")] = "II"
Myeloid_cell_Neu@meta.data$NewSample.ID = Myeloid_cell_Neu@meta.data$orig.ident
temp = left_join(Myeloid_cell_Neu@meta.data, PatientInfo_Neu, by = "NewSample.ID")
identical(Myeloid_cell_Neu@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_Neu@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_Neu@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_Neu@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_Neu@meta.data$Clinical.stage = temp$Clinical.stage.y
temp = left_join(Myeloid_cell_Neu@meta.data, PatientInfo_surv, by = "NewSample.ID")
identical(Myeloid_cell_Neu@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_Neu@meta.data$OS_month = temp$OS_month
Myeloid_cell_Neu@meta.data$DFS_month = temp$DFS_month
Myeloid_cell_Neu@meta.data$event = temp$event
Idents(Myeloid_cell_Neu) = Myeloid_cell_Neu@meta.data$histological.type.short
Myeloid_cell_Neu = subset(Myeloid_cell_Neu, idents = "adeno")
Idents(Myeloid_cell_Neu) = Myeloid_cell_Neu@meta.data$Tumors.for.scRNA.seq.short
Myeloid_cell_Neu = subset(Myeloid_cell_Neu, idents = "P")
Idents(Myeloid_cell_Neu) = Myeloid_cell_Neu$Subtype
## sample 5000 cells for monocle
set.seed(2023)
pbmc = Myeloid_cell_Neu[,sample(1:ncol(Myeloid_cell_Neu),5000)]
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds = detectGenes(cds, min_expr = 0.1)
expressed_genes = row.names(subset(fData(cds),num_cells_expressed >= 10))
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~Subtype",cores=1) 
saveRDS(diff, file = "./3 - GBC_output data_subtype/01_Neutrophil/Pseudotime_diff.RDS")
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
saveRDS(cds, file = "./3 - GBC_output data_subtype/01_Neutrophil/Pseudotime_obj.RDS")
p = plot_cell_trajectory(cds,color_by="Clinical.stage",cell_size = 0.5,show_branch_points = FALSE) + facet_wrap("~Clinical.stage", nrow = 2) + NoLegend()
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/Neutrophil_Subtype_ClinicalStage.pdf", width =4, height = 3)
print(p)
dev.off()
p = plot_cell_trajectory(cds,color_by="Subtype", show_branch_points = FALSE, cell_size = 0.5)  + 
  ggplot2::theme(legend.position = "right") # + facet_wrap("~Subtype", nrow = 2)
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/Neutrophil_Subtype_OnPseudotime.pdf", width =4, height = 2)
print(p)
dev.off()

## Fig2E ####
Myeloid_cell_Neu_3k = Myeloid_cell_Neu[,sample(1:ncol(Myeloid_cell_Neu),3000)]
exprMat  <-  as.matrix(Myeloid_cell_Neu_3k@assays$RNA@data)
dim(exprMat)
cellInfo <-  subset(Myeloid_cell_Neu_3k@meta.data, select = c("RNA_snn_res.0.2","nFeature_RNA","nCount_RNA"))
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
cellInfo$CellType = paste0("N_",cellInfo$CellType)
head(cellInfo)
table(cellInfo$CellType)
scenicOptions <- initializeScenic(org="hgnc",
                                  dbDir="./cisTarget_databases/", nCores=1)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log )
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
library(Seurat) 
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
load("./3 - GBC_output data_subtype/01_Neutrophil/SCENIC/GBC_subtype_SCENIC_230605.RData")
setwd("./3 - GBC_output data_subtype/01_Neutrophil/SCENIC/")
scenicOptions=readRDS(file="int/scenicOptions.Rds")
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="IRF7"]
viewMotifs(tableSubset) 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="IRF7" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
## cluster-specific regulator 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
colnames(regulonActivity_byCellType_Scaled) = c("N_C0_S100A8/9+","N_C1_MIF+_S100A10-","N_C2_IFIT1/2/3+","N_C5_MIF+_S100A10+")
p = pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                       color=colorRampPalette(c("#3FBF70","white","#785F2A"))(100), breaks=seq(-3, 3, length.out = 100),
                       treeheight_row=10, treeheight_col=10, border_color=NA, cluster_cols = T)
pdf(file = "../Neu_Sub_TF.pdf", width =3.5, height = 4)
print(p)
dev.off()

## Fig2F ####
features = c("STAT2","IRF7","STAT1","ETV7","IRF2","CD274")
p = DotPlot(Myeloid_cell_Neu, features = features)  +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/Neu_Sub_TF_exprs.pdf", width =7, height = 2)
print(p)
dev.off()

## Fig2J ####
library(CellChat)
library(patchwork)
library(Seurat)
library(ggpubr)
library(reshape2)
library(ComplexHeatmap)
library(purrr)

subtype = read.csv("./3 - GBC_output data_subtype/ToLS_Shanno/Subtype_Entropy_230613.csv", row.names = 1)
subtype = subtype[subtype$Cluster > 0.625,]
subtype = subtype$X

netVisual_bubble_replot <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                                    comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                                    min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, color.text = NULL,
                                    title.name = NULL, font.size = 10, font.size.title = 10, show.legend = TRUE,
                                    grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                                    return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
  } else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  
  df = object
  
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
  
  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  #g <- g + scale_radius(range = c(1,3), breaks = values,labels = names(values), name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  
  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5+length(dataset.name[comparison]), length(group.names0)*length(dataset.name[comparison]), by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      } else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      #names(color) <- dataset.name[comparison]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  } else {
    return(g)
  }
  
}

netVisual_chord_gene_zt <- function(object, 
                                    #slot.name = "net", 
                                    color.use = NULL,
                                    signaling = NULL, pairLR.use = NULL, net = NULL,
                                    sources.use = NULL, targets.use = NULL,
                                    lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                    link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                    transparency = 0.4, link.border = NA,
                                    title.name = NULL, legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE,
                                    thresh = 0.05,
                                    ...){
  # if (!is.null(pairLR.use)) {
  #   if (!is.data.frame(pairLR.use)) {
  #     stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
  #   } else if ("pathway_name" %in% colnames(pairLR.use)) {
  #     message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
  #     slot.name = "netP"
  #   }
  # }
  
  # if (!is.null(pairLR.use) & !is.null(signaling)) {
  #   stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  # }
  
  # if (is.null(net)) {
  #   prob <- slot(object, "net")$prob
  #   pval <- slot(object, "net")$pval
  #   prob[pval > thresh] <- 0
  #   net <- reshape2::melt(prob, value.name = "prob")
  #   colnames(net)[1:3] <- c("source","target","interaction_name")
  #   
  #   pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence"))
  #   idx <- match(net$interaction_name, rownames(pairLR))
  #   temp <- pairLR[idx,]
  #   net <- cbind(net, temp)
  # }
  
  # if (!is.null(signaling)) {
  #   pairLR.use <- data.frame()
  #   for (i in 1:length(signaling)) {
  #     pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
  #     pairLR.use <- rbind(pairLR.use, pairLR.use.i)
  #   }
  # }
  
  # if (!is.null(pairLR.use)){
  #   if ("interaction_name" %in% colnames(pairLR.use)) {
  #     net <- subset(net,interaction_name %in% pairLR.use$interaction_name)
  #   } else if ("pathway_name" %in% colnames(pairLR.use)) {
  #     net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
  #   }
  # }
  
  # if (slot.name == "netP") {
  #   net <- dplyr::select(net, c("source","target","pathway_name","prob"))
  #   net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
  #   net <- net %>% dplyr::group_by(source_target, pathway_name) %>% dplyr::summarize(prob = sum(prob))
  #   a <- stringr::str_split(net$source_target, "sourceTotarget", simplify = T)
  #   net$source <- as.character(a[, 1])
  #   net$target <- as.character(a[, 2])
  #   net$ligand <- net$pathway_name
  #   net$receptor <- " "
  # }
  
  # # keep the interactions associated with sources and targets of interest
  # if (!is.null(sources.use)){
  #   if (is.numeric(sources.use)) {
  #     sources.use <- levels(object@idents)[sources.use]
  #   }
  #   net <- subset(net, source %in% sources.use)
  # } else {
  #   sources.use <- levels(object@idents)
  # }
  # if (!is.null(targets.use)){
  #   if (is.numeric(targets.use)) {
  #     targets.use <- levels(object@idents)[targets.use]
  #   }
  #   net <- subset(net, target %in% targets.use)
  # } else {
  #   targets.use <- levels(object@idents)
  # }
  # # remove the interactions with zero values
  # df <- subset(net, prob > 0)
  # 
  # if (nrow(df) == 0) {
  #   stop("No signaling links are inferred! ")
  # }
  # 
  # if (length(unique(net$ligand)) == 1) {
  #   message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  # }
  df = object
  df$id <- 1:nrow(df)
  # deal with duplicated sector names
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }
  
  # cell.order.sources <- levels(object@idents)[levels(object@idents) %in% sources.use]
  # cell.order.targets <- levels(object@idents)[levels(object@idents) %in% targets.use]
  
  cell.order.sources <- sort(unique(df$source))
  cell.order.targets <- sort(unique(df$target))
  object_idents = factor(unique(c(cell.order.sources,cell.order.targets)),levels = unique(c(cell.order.sources,cell.order.targets)))
  
  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  # df.ordered.source <- df[with(df, order(source, target, -prob)), ]
  # df.ordered.target <- df[with(df, order(target, source, -prob)), ]
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]
  
  order.source <- unique(df.ordered.source[ ,c('ligand','source')])
  order.target <- unique(df.ordered.target[ ,c('receptor','target')])
  
  # define sector order
  order.sector <- c(order.source$ligand, order.target$receptor)
  
  # define cell type color
  if (is.null(color.use)){
    color.use = scPalette(nlevels(object_idents))
    names(color.use) <- levels(object_idents)
    color.use <- color.use[levels(object_idents) %in% as.character(union(df$source,df$target))]
  } else if (is.null(names(color.use))) {
    names(color.use) <- levels(object_idents)
    color.use <- color.use[levels(object_idents) %in% as.character(union(df$source,df$target))]
  }
  
  # define edge color
  edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)
  
  # define grid colors
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector
  
  df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]
  
  if (directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }
  circos.clear()
  chordDiagram(df.plot,
               order = order.sector,#
               col = edge.color,
               grid.col = grid.col,#
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type,
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)
  
  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }
  
  circos.clear()
  if(!is.null(title.name)){
    text(-0, 1.02, title.name, cex=1)
  }
  gg <- recordPlot()
  return(gg)
}

cellchat_NT <- readRDS("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/CellChat/cellchat_P.RDS")
cellchat_T <- readRDS("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/CellChat/cellchat_P_Mets.RDS")

group.cellchat_P = unique(union(levels(cellchat_NT@idents),levels(cellchat_T@idents)))
cellchat_T <- liftCellChat(cellchat_T, group.cellchat_P)
cellchat_NT <- liftCellChat(cellchat_NT, group.cellchat_P)
object.list <- list(NT = cellchat_NT, `T` = cellchat_T)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

p = netVisual_bubble(cellchat, 
                     sources.use = group.cellchat_P, 
                     targets.use = group.cellchat_P,  
                     comparison = c(1, 2), 
                     max.dataset = 2, title.name = "Increased signaling in T", angle.x = 45, remove.isolate = T,return.data = T)

object = cellchat
meta <- object@meta
if (is.list(object@idents)) {
  meta$group.cellchat <- object@idents$joint
} else {
  meta$group.cellchat <- object@idents
}
if (!identical(rownames(meta), colnames(object@data.signaling))) {
  cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
  warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
  rownames(meta) <- colnames(object@data.signaling)
}
w10x <- Seurat::CreateSeuratObject(counts = cellchat@data.signaling, meta.data = meta)
Idents(w10x) = w10x$group.cellchat

source = list()
cycle = as.character(unique(p$communication$source))
for (i in 1:length(cycle)) {
  p_sub = p$communication[p$communication$source %in% cycle[i],]
  w10x_sub = subset(w10x, 
                    idents = cycle[i], 
                    features = unique(unlist(stringr::str_split(unlist(stringr::str_split(unique(p_sub$ligand),'_')),':'))))
  Idents(w10x_sub) = w10x_sub$datasets
  
  w10x_toplot = as.data.frame(t(w10x_sub@assays$RNA@data))
  w10x_toplot$bind = rownames(w10x_toplot)
  group_meta = as.data.frame(w10x_sub@active.ident)
  group_meta$bind = rownames(group_meta)
  w10x_toplot = left_join(w10x_toplot,group_meta,by='bind')
  
  if (length(unique(w10x_toplot$`w10x_sub@active.ident`))==1) {
    
    ligand = colnames(w10x_toplot)[1:(length(w10x_toplot)-2)]
    source[[i]] = ligand
    names(source)[i] = as.character(cycle[i])
    
  } else {
    colnames(w10x_toplot) = gsub("-","_",colnames(w10x_toplot))
    ligand = data.frame()
    
    for(j in 1:(ncol(w10x_toplot)-2)){
      variable1 = colnames(w10x_toplot)[j]
      a = compare_means(as.formula(sprintf("%s ~ %s ", variable1,"w10x_sub@active.ident")), data = w10x_toplot)
      ligand = rbind(ligand,a)
    }
    ligand$.y. = gsub("_","-",(ligand$.y.))
    colnames(w10x_toplot) = gsub("_","-",colnames(w10x_toplot))
    
    ligand = ligand[ligand$p.format<0.05,]$.y.
    
    w10x_P = subset(w10x_sub, idents = 'NT')
    w10x_P = as.data.frame(t(w10x_P@assays$RNA@data))
    w10x_P = subset(w10x_P, select = ligand)
    
    w10x_LN = subset(w10x_sub, idents = 'T')
    w10x_LN = as.data.frame(t(w10x_LN@assays$RNA@data))
    w10x_LN = subset(w10x_LN, select = ligand)
    
    identical(colnames(w10x_P),colnames(w10x_LN))
    
    ligand = colnames(w10x_P)[apply(w10x_P,2,mean) < apply(w10x_LN,2,mean)]
    
    source[[i]] = ligand
    names(source)[i] = as.character(cycle[i])
  }
  
  print(paste0(i,' / ',length(cycle)))
  
}

p_id = c()
for (i in 1:length(source)) {
  temp = p$communication[p$communication$source %in% names(source)[i],]
  temp1 = unique(temp$ligand[temp$ligand %in% source[[i]]])
  if(length(grep("_",temp1)) == 0){
    temp2 = temp1
  } else {
    temp2 = temp1[-grep("_",temp1)]
    temp_test = temp1[grep("_",temp1)]
    for (j in 1:length(temp_test)) {
      temp_test_split = unique(unlist(stringr::str_split(unlist(stringr::str_split(temp_test[j],'_')),':')))
      if(length(temp_test_split[temp_test_split %in% source[[i]]]) == length(temp_test_split)){
        temp2 = c(temp2,temp_test[j])
      }
    }
  }
  p_id = c(p_id,rownames(temp[temp$ligand %in% temp2,])) 
  print(paste0(i, " : ", length(temp1), " , ", length(temp2)))
}

pp = p$communication[p_id,]

receptor = unique(unlist(stringr::str_split(unlist(stringr::str_split(unique(pp$receptor),'_')),':')))

object = cellchat
meta <- object@meta
if (is.list(object@idents)) {
  meta$group.cellchat <- object@idents$joint
} else {
  meta$group.cellchat <- object@idents
}
if (!identical(rownames(meta), colnames(object@data.signaling))) {
  cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
  warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
  rownames(meta) <- colnames(object@data.signaling)
}
w10x <- Seurat::CreateSeuratObject(counts = cellchat@data.signaling, meta.data = meta)

receptor = receptor[receptor %in% rownames(w10x)]

Idents(w10x) = w10x$datasets
w10x = subset(w10x, idents = 'T')

Idents(w10x) = w10x$group.cellchat
w10x = subset(w10x, features = receptor)

target = list()
cycle = as.character(unique(pp$target))
for (i in 1:length(cycle)) {
  w10x_target = subset(w10x, idents = cycle[i])
  w10x_target = as.data.frame(t(w10x_target@assays$RNA@data))
  receptor_target = apply(w10x_target,2,function(x) length(x[x>0])/length(x))
  target[[i]] = names(receptor_target[receptor_target > 0.3]) # cutoff set 0.3
  names(target)[i] = as.character(cycle[i])
  print(paste0(i,' / ',length(cycle)))
}

p_id = c()
for (i in 1:length(target)) {
  temp = pp[pp$target %in% names(target)[i],]
  temp1 = unique(temp$receptor[temp$receptor %in% target[[i]]])
  if(length(grep("_",temp1)) == 0){
    temp2 = temp1
  } else {
    temp2 = temp1[-grep("_",temp1)]
    temp_test = temp1[grep("_",temp1)]
    for (j in 1:length(temp_test)) {
      temp_test_split = unique(unlist(stringr::str_split(unlist(stringr::str_split(temp_test[j],'_')),':')))
      if(length(temp_test_split[temp_test_split %in% target[[i]]]) == length(temp_test_split)){
        temp2 = c(temp2,temp_test[j])
      }
    }
  }
  p_id = c(p_id,rownames(temp[temp$receptor %in% temp2,])) 
}

ppp = pp[p_id,]

ligand = unique(unlist(stringr::str_split(unlist(stringr::str_split(unique(ppp$ligand),'_')),':')))

object = cellchat
meta <- object@meta
if (is.list(object@idents)) {
  meta$group.cellchat <- object@idents$joint
} else {
  meta$group.cellchat <- object@idents
}
if (!identical(rownames(meta), colnames(object@data.signaling))) {
  cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
  warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
  rownames(meta) <- colnames(object@data.signaling)
}
w10x <- Seurat::CreateSeuratObject(counts = cellchat@data.signaling, meta.data = meta)

ligand = ligand[ligand %in% rownames(w10x)]

Idents(w10x) = w10x$datasets
w10x = subset(w10x, idents = 'T')

Idents(w10x) = w10x$group.cellchat
w10x = subset(w10x, features = ligand)

target = list()
cycle = as.character(unique(ppp$source))
for (i in 1:length(cycle)) {
  w10x_target = subset(w10x, idents = cycle[i])
  w10x_target = as.data.frame(t(w10x_target@assays$RNA@data))
  receptor_target = apply(w10x_target,2,function(x) length(x[x>0])/length(x))
  target[[i]] = names(receptor_target[receptor_target > 0.3]) # cutoff set 0.3
  names(target)[i] = as.character(cycle[i])
  print(paste0(i,' / ',length(cycle)))
}

p_id = c()
for (i in 1:length(target)) {
  temp = ppp[ppp$source %in% names(target)[i],]
  temp1 = unique(temp$ligand[temp$ligand %in% target[[i]]])
  if(length(grep("_",temp1)) == 0){
    temp2 = temp1
  } else {
    temp2 = temp1[-grep("_",temp1)]
    temp_test = temp1[grep("_",temp1)]
    for (j in 1:length(temp_test)) {
      temp_test_split = unique(unlist(stringr::str_split(unlist(stringr::str_split(temp_test[j],'_')),':')))
      if(length(temp_test_split[temp_test_split %in% target[[i]]]) == length(temp_test_split)){
        temp2 = c(temp2,temp_test[j])
      }
    }
  }
  p_id = c(p_id,rownames(temp[temp$ligand %in% temp2,])) 
}

pppp = ppp[p_id,]

pppp = pppp[pppp$dataset %in% "T",]
pppp = pppp[pppp$source %in% subtype,]
pppp = pppp[pppp$target %in% subtype,]

ppppp = pppp[pppp$target %in% unique(as.character(pppp$target))[grep("N_C0|N_C7",unique(as.character(pppp$target)))],]
ppppp = ppppp[ppppp$pathway_name %in% ppppp$pathway_name[grep("CXCL",ppppp$pathway_name)],]
g=netVisual_chord_gene_zt(ppppp, 
                          title.name = "N0_N7_target_Pmet",
                          link.target.prop = F,
                          annotationTrackHeight = 0.1
)
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N0_N7_target_Pmet_chord.pdf", width =10, height = 10)
print(g)
dev.off()

## ####
## Fig2N ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_Neu@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_Neu@meta.data$NewSample.ID = Myeloid_cell_Neu@meta.data$orig.ident
temp = left_join(Myeloid_cell_Neu@meta.data, PatientInfo, by = "NewSample.ID")
identical(Myeloid_cell_Neu@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_Neu@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_Neu@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_Neu@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_Neu@meta.data$Clinical.stage = temp$Clinical.stage.x
N0 = subset(Myeloid_cell_Neu, idents = "N_C0_MNDA")
N0 = subset(N0, histological.type.short %in% c("adeno","CC","XGC","LG","HG"))
N0 = subset(N0, Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"))
N0$histological.type.short = ifelse(N0$histological.type.short == "adeno", "adeno", "NT")
Idents(N0) = N0$histological.type.short
Endothelium_0_P_DEGs <- FindMarkers(N0, ident.1 = "adeno", min.pct = 0.25, logfc.threshold = 0.25)
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[Endothelium_0_P_DEGs$p_val_adj < 0.05,]
Endothelium_0_P_DEGs$gene = rownames(Endothelium_0_P_DEGs)
temp <- bitr(Endothelium_0_P_DEGs$gene,fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb='org.Hs.eg.db')
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[temp$SYMBOL,]
Endothelium_0_P_DEGs$gene = temp$ENTREZID
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[Endothelium_0_P_DEGs$avg_log2FC > 0,]
Endothelium_0_P_DEGs$cluster = "N0_UP"
Endothelium_0_P_DEGs1 = Endothelium_0_P_DEGs
N7 = subset(Myeloid_cell_Neu, idents = "N_C7_MNDA_S100A8")
N7 = subset(N7, histological.type.short %in% c("adeno","CC","XGC","LG","HG"))
N7 = subset(N7, Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"))
N7$histological.type.short = ifelse(N7$histological.type.short == "adeno", "adeno", "NT")
Idents(N7) = N7$histological.type.short
Endothelium_0_P_DEGs <- FindMarkers(N7, ident.1 = "adeno", min.pct = 0.25, logfc.threshold = 0.25)
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[Endothelium_0_P_DEGs$p_val_adj < 0.05,]
Endothelium_0_P_DEGs$gene = rownames(Endothelium_0_P_DEGs)
temp <- bitr(Endothelium_0_P_DEGs$gene,fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb='org.Hs.eg.db')
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[temp$SYMBOL,]
Endothelium_0_P_DEGs$gene = temp$ENTREZID
Endothelium_0_P_DEGs = Endothelium_0_P_DEGs[Endothelium_0_P_DEGs$avg_log2FC > 0,]
Endothelium_0_P_DEGs$cluster = "N7_UP"
Endothelium_0_P_DEGs = rbind(Endothelium_0_P_DEGs,Endothelium_0_P_DEGs1)
formula_res <- compareCluster(gene~cluster, data=Endothelium_0_P_DEGs,
                              fun="enrichGO",
                              OrgDb = org.Hs.eg.db,
                              ont = "BP" ,
                              pAdjustMethod = "BH",
                              readabl = TRUE)
temp = formula_res@compareClusterResult
write.csv(temp, file = "./3 - GBC_output data_subtype/01_Neutrophil/N0N7up_NTvsAdeno_v240422.csv")

formula_res1 <- filter(formula_res, Description %in% c("actin filament organization",
                                                       "regulation of apoptotic signaling pathway",
                                                       "cell chemotaxis",
                                                       "regulation of inflammatory response",
                                                       "phagocytosis",
                                                       "regulation of cell morphogenesis",
                                                       "regulation of angiogenesis"))
p=dotplot(formula_res1, showCategory = 10) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "horizontal")
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N0_N7_UP-GO in Adeno vs NT.pdf", width =4.5, height = 3)
print(p)
dev.off()

## Fig2O ####
temp = c("CYP4F3","SLC22A4","IL17A","MPO","F3","TNFRSF10C","CREB5","SELP","VNN3","MME","ENTPD4","G0S2","S100A12","BST1","CLEC6A","HPSE","FPR2","CEACAM3","IL8")
write.csv(temp, file = "./3 - GBC_output data_subtype/01_Neutrophil/NET_sig_v240422.csv")

temp = list(temp)
Myeloid_cell_Neu <- AddModuleScore(Myeloid_cell_Neu, features = temp, name = "NET")
p <- VlnPlot(Myeloid_cell_Neu,"NET1",pt.size = 0) +
  stat_compare_means(label = "p", ref.group = ".all.",)+
  geom_hline(yintercept = mean(Myeloid_cell_Neu$NET1),linetype=2)+
  stat_summary(fun = "mean", fill = "white", size = 2, geom = "point", shape = 23) +NoLegend()
pdf(file = paste0("./3 - GBC_output data_subtype/01_Neutrophil/Score_NET1", ".pdf"), width =3, height = 4.5)
print(p)
dev.off()

# 15. DC ####
## Fig1E ####
DimPlot(Myeloid_cell_DC, reduction = "umap", pt.size = 2, raster = T) + ggplot2::theme(legend.position = "right", panel.background = element_blank())

## FigS2G ####
p = FeaturePlot(Myeloid_cell_DC, features = c("CD1E","LAMP3","LILRA4","CLEC9A"),ncol=4, pt.size = 3, raster = T, max.cutoff = c(3,3,3,4))
pdf(file = "./3 - GBC_output data_subtype/02_DC/DC_SubtypeDefinition_MarkerExprs.pdf", width =12, height = 2.5)
print(p)
dev.off()

## * Subset common subgroups ####
Myeloid_cell_DC = subset(Myeloid_cell_DC, idents = c("DC_C3_cDC2_FCGBP"), invert=T)
saveRDS(Myeloid_cell_DC, file = "./01_GBC_subtypes_zt/02_DC/Myeloid_cell_DC.RDS")

## * Subset subgroups ####
Myeloid_cell_DC <- readRDS("./3 - GBC_output data_subtype/02_DC/Myeloid_cell_DC.RDS")
Myeloid_cell_DC2 = subset(Myeloid_cell_DC, idents = c("DC_C0_cDC2_IL1B","DC_C1_cDC2_others","DC_C4_cDC2_PPP1R14A"))

## Fig3A ####
Myeloid_cell_DC_DEGs_allpos <- FindAllMarkers(Myeloid_cell_DC, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
Myeloid_cell_DC_DEGs_allpos = Myeloid_cell_DC_DEGs_allpos[Myeloid_cell_DC_DEGs_allpos$p_val_adj < 0.05,]
write.csv(Myeloid_cell_DC_DEGs_allpos, file = "./3 - GBC_output data_subtype/02_DC/DEGs among common DC subtypes_v240422.csv")

Myeloid_cell_DC_DEGs_allpos %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
p = DotPlot(Myeloid_cell_DC, features = unique(top5$gene))  +
  # coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "vertical")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
pdf(file = "./3 - GBC_output data_subtype/02_DC/DC_SubtypeDEGs_Top5.pdf", width =12, height = 3.5)
print(p)
dev.off()

## Fig3B ####
temp <- bitr(Myeloid_cell_DC_DEGs_allpos$gene,fromType = 'SYMBOL',
             toType = 'ENTREZID',
             OrgDb='org.Hs.eg.db')
Myeloid_cell_DC_DEGs_allpos = Myeloid_cell_DC_DEGs_allpos[temp$SYMBOL,]
Myeloid_cell_DC_DEGs_allpos$gene = temp$ENTREZID
formula_res <- compareCluster(gene~cluster, data=Myeloid_cell_DC_DEGs_allpos,
                              fun="enrichGO",
                              OrgDb = org.Hs.eg.db,
                              ont = "BP" ,
                              pAdjustMethod = "BH",
                              readabl = TRUE)
temp = formula_res@compareClusterResult
write.csv(temp, file = "./3 - GBC_output data_subtype/02_DC/GOBP among common DC subtypes_v240422.csv")

g = dotplot(formula_res, showCategory = 10) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "horizontal")
formula_res <- filter(formula_res, Description %in% c("regulation of cell-cell adhesion", # DC0/2 top
                                                      "generation of precursor metabolites and energy", # DC1 top
                                                      "positive regulation of cellular component biogenesis", # DC4 top
                                                      "aerobic electron transport chain", # DC5 top
                                                      "mRNA processing", # DC6 top
                                                      "intracellular receptor signaling pathway", # DC7 top
                                                      "antigen processing and presentation",
                                                      "ATP metabolic process",
                                                      "oxidative phosphorylation"))
p = dotplot(formula_res, showCategory = 5) +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5,size = 10),
        axis.text.y = element_text(size = 10),
        legend.direction = "vertical",
        legend.position = "right",
        legend.box = "vertical")
pdf(file = "./3 - GBC_output data_subtype/02_DC/DC_SubtypeDEPathways_Selected_230912.pdf", width =6.5, height = 5.5)
print(p)
dev.off()

## Fig3F ####
Myeloid_cell_DC <- readRDS("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/02_DC/Myeloid_cell_DC.RDS")
Myeloid_cell_DC$Subtype = Idents(Myeloid_cell_DC)
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_DC@meta.data$orig.ident),]
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo_DC = left_join(PatientInfo_DC,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_DC@meta.data$NewSample.ID = Myeloid_cell_DC@meta.data$orig.ident
temp = left_join(Myeloid_cell_DC@meta.data, PatientInfo_DC, by = "NewSample.ID")
identical(Myeloid_cell_DC@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_DC@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_DC@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_DC@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_DC@meta.data$Clinical.stage = temp$Clinical.stage.x
temp = left_join(Myeloid_cell_DC@meta.data, PatientInfo_surv, by = "NewSample.ID")
identical(Myeloid_cell_DC@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_DC@meta.data$OS_month = temp$OS_month
Myeloid_cell_DC@meta.data$DFS_month = temp$DFS_month
Myeloid_cell_DC@meta.data$event = temp$event
Idents(Myeloid_cell_DC) = Myeloid_cell_DC@meta.data$histological.type.short
Myeloid_cell_DC = subset(Myeloid_cell_DC, idents = "adeno")
Idents(Myeloid_cell_DC) = Myeloid_cell_DC@meta.data$Tumors.for.scRNA.seq.short
Myeloid_cell_DC = subset(Myeloid_cell_DC, idents = "P")
Idents(Myeloid_cell_DC) = Myeloid_cell_DC$Subtype

Myeloid_cell_DC = subset(Myeloid_cell_DC, idents = c(
  "DC_C0_cDC2_IL1B", 
  "DC_C1_cDC2_others", 
  "DC_C4_cDC2_PPP1R14A"
))
pbmc = Myeloid_cell_DC
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
p_data <- pbmc@meta.data 
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
expressed_genes <- VariableFeatures(pbmc)
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~Subtype",cores=1) 
deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
deg = deg[1:400]
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds,color_by="State", size=1,show_branch_points = FALSE) + facet_wrap("~State", nrow = 1) + NoLegend()
plot_cell_trajectory(cds,color_by="Clinical.stage",size=1,show_branch_points = FALSE) + facet_wrap("~Clinical.stage", nrow = 2) + NoLegend()
cds <- orderCells(cds,root_state = 3)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_branch_points = FALSE)
p = plot_cell_trajectory(cds,color_by="Subtype", show_branch_points = FALSE)# + facet_wrap("~Subtype", nrow = 2) + NoLegend()
pdf(file = "./01_GBC_subtypes_zt/02_DC/DC_cDC2_pseudotime_allDEGs.pdf", width =5, height = 5)
print(p)
dev.off()

## Fig3G & Fig3H ####
Myeloid_cell_DC4 = Myeloid_cell_DC

library(Seurat) 
library(SCENIC)
library(doParallel)
library(SCopeLoomR)

load("./3 - GBC_output data_subtype/02_DC/SCENIC/cDC2/GBC_subtype_SCENIC_230210_cDC2.RData")
setwd("./3 - GBC_output data_subtype/02_DC/SCENIC/cDC2/")
scenicOptions=readRDS(file="int/scenicOptions.Rds")
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="IRF7"]
viewMotifs(tableSubset) 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="IRF7" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
## cluster-specific regulator 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
p = pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                       color=colorRampPalette(c("#3FBF70","white","#785F2A"))(100), breaks=seq(-3, 3, length.out = 100),
                       treeheight_row=10, treeheight_col=10, border_color=NA, cluster_cols = T)
pdf(file = "./3 - GBC_output data_subtype/02_DC/cDC2_SCENIC_230912.pdf", width =3, height = 4.5)
print(p)
dev.off()
p = pheatmap::pheatmap(t(regulonActivity_byCellType_Scaled), #fontsize_row=3,
                       color=colorRampPalette(c("#3FBF70","white","#785F2A"))(100), breaks=seq(-3, 3, length.out = 100),
                       treeheight_row=10, treeheight_col=10, border_color=NA, cluster_cols = T)
pdf(file = "./3 - GBC_output data_subtype/02_DC/cDC2_SCENIC_2309121.pdf", width =6, height = 3)
print(p)
dev.off()
p = DotPlot(Myeloid_cell_cDC2, features = c("STAT1","STAT2","IRF7","IRF9",rownames(Myeloid_cell_cDC2)[grep("HLA-D|CD1A|CD1B|CD1C|CD1D",rownames(Myeloid_cell_cDC2))]),cols = c("blue", "red")) + 
  theme(axis.text.x=element_text(angle = -90,hjust = 0,vjust = 0.5))
pdf(file = "./3 - GBC_output data_subtype/02_DC/cDC2_SCENIC_exprs_2309121.pdf", width =10, height = 3)
print(p)
dev.off()

## Fig3I ####
Myeloid_cell_DC <- readRDS("./3 - GBC_output data_subtype/02_DC/Myeloid_cell_DC.RDS")
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_DC@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_DC@meta.data$NewSample.ID = Myeloid_cell_DC@meta.data$orig.ident
temp = left_join(Myeloid_cell_DC@meta.data, PatientInfo, by = "NewSample.ID")
identical(Myeloid_cell_DC@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_DC@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_DC@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_DC@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_DC@meta.data$Clinical.stage = temp$Clinical.stage.x
Myeloid_cell_DC4 = subset(Myeloid_cell_DC, idents = "DC_C4_cDC2_PPP1R14A")

Idents(Myeloid_cell_DC4) = Myeloid_cell_DC4@meta.data$histological.type.short
Myeloid_cell_DC4 = subset(Myeloid_cell_DC4, idents = c("adeno","CC","LG","XGC")) # No HG
Idents(Myeloid_cell_DC4) = Myeloid_cell_DC4@meta.data$Tumors.for.scRNA.seq.short
Myeloid_cell_DC4 = subset(Myeloid_cell_DC4, idents = c("P","XGC","CC"))
Myeloid_cell_DC4@meta.data$histological.type.short[Myeloid_cell_DC4@meta.data$histological.type.short %in% c("CC","LG","XGC")] = "NT"
Myeloid_cell_DC4@meta.data$histological.type.short[Myeloid_cell_DC4@meta.data$histological.type.short %in% c("adeno")] = "Adeno"
Idents(Myeloid_cell_DC4) = Myeloid_cell_DC4@meta.data$histological.type.short
## DEGs: Volcano
Myeloid_cell_DC_NTvsT_DEGs <- FindMarkers(Myeloid_cell_DC4, ident.1 = "Adeno", min.pct = 0.25, logfc.threshold = 0.25)
Myeloid_cell_DC_NTvsT_DEGs = Myeloid_cell_DC_NTvsT_DEGs[Myeloid_cell_DC_NTvsT_DEGs$p_val_adj < 0.05,]
write.csv(Myeloid_cell_DC_NTvsT_DEGs, file = "./3 - GBC_output data_subtype/02_DC/DEG DC4_NTvsAdeno_v240422.csv")

Dat = Myeloid_cell_DC_NTvsT_DEGs
Dat$Gene = rownames(Dat)
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) > 0, ifelse(Dat$avg_log2FC > 0 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
p = ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
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
  ggtitle("DC_C4_cDC2_PPP1R14A [Adeno (vs NT)]")
pdf(file = "./3 - GBC_output data_subtype/02_DC/DC_DC4_volcano.pdf", width =8, height = 7)
print(p)
dev.off()

## Fig3J ####
Myeloid_cell_DC <- readRDS("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/02_DC/Myeloid_cell_DC.RDS")
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo = PatientInfo[PatientInfo$NewSample.ID %in% unique(Myeloid_cell_DC@meta.data$orig.ident),]
PatientInfo = left_join(PatientInfo,PatientInfo_surv,by="NewSample.ID")
Myeloid_cell_DC@meta.data$NewSample.ID = Myeloid_cell_DC@meta.data$orig.ident
temp = left_join(Myeloid_cell_DC@meta.data, PatientInfo, by = "NewSample.ID")
identical(Myeloid_cell_DC@meta.data$orig.ident, temp$orig.ident)
Myeloid_cell_DC@meta.data$histological.type.short = temp$histological.type.short
Myeloid_cell_DC@meta.data$Tumors.for.scRNA.seq.short = temp$Tumors.for.scRNA.seq.short
Myeloid_cell_DC@meta.data$metastasis.type = temp$metastasis.type
Myeloid_cell_DC@meta.data$Clinical.stage = temp$Clinical.stage.x
DC4 = subset(Myeloid_cell_DC, idents = c("DC_C4_cDC2_PPP1R14A"))
expr = DC4
df = FetchData(expr,vars = c('CD1A','Clinical.stage','orig.ident','metastasis.type','Tumors.for.scRNA.seq.short','histological.type.short'))
df = df[df$histological.type.short %in% c("adeno","CC","XGC","HG","LG"),]
df = df[df$Tumors.for.scRNA.seq.short %in% c("P","CC","XGC"),]
df = df[df$histological.type.short == "adeno",]
#df$Clinical.stage = ifelse(df$histological.type.short %in% c("CC","XGC","HG","LG"), "NT", df$Clinical.stage)
df$Clinical.stage = ifelse(df$Clinical.stage %in% c('IIA','IIB'),'II',df$Clinical.stage)
final2= df %>% group_by(Clinical.stage) %>% summarise(mean = mean(CD1A),N=length(CD1A),
                                                      sd = sd(CD1A),se=sd/sqrt(N))
final2$group = "CD1A"
final = final2
final$Clinical.stage = factor(final$Clinical.stage, levels = c("I","II","IIIA","IIIB","IVA","IVB"))
p = ggplot(final,aes(Clinical.stage,mean,group = 1))+
  geom_ribbon(aes(ymin=mean-se,ymax=mean+se),fill='#F6F1ED')+
  geom_line(color='#B7916C',size=1)+
  geom_point(color='#B7916C',size=2)+
  theme(panel.background = element_blank(),
        panel.grid = element_line(color="white"),
        # axis.title = element_blank(),
        #axis.ticks = element_blank(),
        axis.text.x=element_text(colour='black',size=12,angle = 45,hjust = 1))+
  ggtitle('CD1A')+theme(axis.line=element_line(colour='black',size=1,lineend = 'square'))+
  theme (axis.text.x = element_text (colour='black', size=12, angle=45), 
         axis.text.y = element_text (colour='black',size=12, angle=45))+
  labs(x='Stage',y='Score')
pdf(file = "./3 - GBC_output data_subtype/02_DC/CD1_in DC4_dynamic.pdf", width =3, height = 4)
print(p)
dev.off()

compare_means(CD1A~Clinical.stage, df, method = "anova") # 7.2e-07

## Fig3K ####
ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("_cDC2_",unique(as.character(pppp$source)))],]
ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[grep("ANXA1",ppppp$ligand)],]
g=netVisual_chord_gene_zt(ppppp, 
                          annotationTrackHeight = 0.1,
                          link.target.prop = F,
                          title.name = "DC4_Pmets")
pdf(file = "./3 - GBC_output data_subtype/02_DC/DC4_source_Pmet_chord.pdf", width =10, height = 5)
print(g)
dev.off()

## Fig3L ####
object = cellchat
meta <- object@meta
if (is.list(object@idents)) {
  meta$group.cellchat <- object@idents$joint
} else {
  meta$group.cellchat <- object@idents
}
if (!identical(rownames(meta), colnames(object@data.signaling))) {
  cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
  warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
  rownames(meta) <- colnames(object@data.signaling)
}
w10x <- Seurat::CreateSeuratObject(counts = cellchat@data.signaling, meta.data = meta)
Idents(w10x) = w10x$group.cellchat

cellchat_new = subset(w10x, idents = group.cellchat_P[grep("DC_C4",group.cellchat_P)])
Idents(cellchat_new) = cellchat_new$datasets
vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(cellchat_new, features = signature,pt.size = 0,
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, label = "p.format") + 
      stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
  }
  map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, ".pdf")
  ggsave(file_name, width = 2.5, height = 3)
}
gene_sig <- c("ANXA1")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../2 - Analysis_V3_checked/3 - GBC_output data_subtype/02_DC/PvsPmet_DC4Ligand_ANXA1", test_sign = comparisons, y_max = 5)

## ####
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
write.csv(Endothelium_DEGs_allpos, file = "./3 - GBC_output data_subtype/04_Endothelium/DEG among common EC subtype_v240422.csv")

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
temp = formula_res@compareClusterResult
write.csv(temp, file = "./3 - GBC_output data_subtype/04_Endothelium/GO among common EC subtype_v240422.csv")

formula_res1 <- filter(formula_res, Description %in% c("antigen processing and presentation",
                                                      "receptor-mediated endocytosis", # EC5 top
                                                      "ATP metabolic process", # EC4/6 top
                                                      "oxidative phosphorylation",
                                                      "cell-substrate adhesion",
                                                      "regulation of angiogenesis",
                                                      "extracellular matrix organization"))
p = dotplot(formula_res1, showCategory = 5) +
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
temp = formula_res@compareClusterResult
write.csv(temp, file = "./3 - GBC_output data_subtype/04_Endothelium/EC0EC3_GO_NTvsAdeno_v240422.csv")

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

# 17. Fig7 (Fig5 in the new version) ####
## Fig7A & Fig7F ####
GBC_AllSubtype <- read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
GBC_TIMESubtype = GBC_AllSubtype[!(GBC_AllSubtype$celltype %in% c("tumor","normal")),]
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% unique(GBC_TIMESubtype$orig.ident),]
PatientInfo_metastasis = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c("adeno"),]
PatientInfo_metastasis = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% c("P"),]
PatientInfo_metastasis$Group = PatientInfo_metastasis$metastasis.type
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group == "P_LI", "P", PatientInfo_metastasis$Group)
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group %in% c("P_LN","P_LM"), "P_Mets", PatientInfo_metastasis$Group)
Select_Patient = PatientInfo_metastasis$NewSample.ID
cycle = Select_Patient
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(GBC_TIMESubtype, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$subtype)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0
Freq_adj = scale(Freq)
subtype = c("CD4T_C1_FOXP3",
            "CD4T_C7_ISG15",
            #"CD4T_C8_MKI67",
            "CD8T_C1_CXCL13",
            #"CD8T_C8_MKI67",
            #"CD8T_C11_PCLAF",
            "γdT_C12_TRDC",
            "NKT_C3_CAPG",
            #"NK_C6_MKI67",
            "F_C0_MMP11",
            "Per_C0_RGS5",
            "DC_C0_cDC2_IL1B",
            "DC_C2_cDC3_FSCN1",
            "DC_C5_cDC1",
            "DC_C6_pDC",
            "M_C2_SPP1",
            "M_C7_AGR2",
            "N_C9_Nc",
            "EC_C2_CXCR4",
            "CD4T_C2_CCR7",
            #"CD4T_C4_Unassign",
            "CD8T_C0_CCR7_GZMK",
            "B_C0_IGHA1",
            "F_C1_CFD",
            "Per_C1_MYH11",
            "DC_C4_cDC2_PPP1R14A",
            #"M_C4_others",
            "N_C0_N0",
            "N_C1_N2",
            "N_C2_N1",
            "N_C7_N0",
            "EC_C3_GJA5",
            "EC_C5_PROX1",
            "Per_C3_STEAP4",
            "M_C1_S100A8",
            "CD8T_C9_MT1X_MT1E",
            "EC_C0_ACKR1")
colnames(Freq_adj)[96] = "γdT_C12_TRDC"
Freq_adj = Freq_adj[,subtype]
Freq_adj[Freq_adj > 1] = 1
Freq_adj[Freq_adj < -0.5] = -0.5
Freq_adj = t(Freq_adj)
## annotation
scales::show_col(c("#51574a","#e9d78e",ggsci::pal_d3("category20")(20)[1:20]))
mycol = c("#51574a","#e9d78e",ggsci::pal_d3("category20")(20)[1:20])
annotation = PatientInfo_metastasis
annotation = subset(annotation, select= c("NewSample.ID","Clinical.stage","Group","Sex","Age","Polyps","Gallstones","PBM","Atrophic.cholecystitis",
                                          "Differentiation","liver.invasion","Bile.duct.invasion","Vascular.invasion","Lymph.node.metastasis","Liver.metastasis","Peritoneal..metastasis"))
TMB_TotalMutation = read.csv(file = "../Publication/240923_Nature Genetics_Revise/4 - WYH/MutationMatrix_update_241106/wes_info_TMB-last col.csv", row.names = 1)
colnames(TMB_TotalMutation)[1] = "NewSample.ID"
TMB_TotalMutation = TMB_TotalMutation[,c(1,11)]
TMB_TotalMutation$TMB = as.numeric(TMB_TotalMutation$TMB)
annotation = left_join(annotation, TMB_TotalMutation, by = "NewSample.ID")
Pop_Mutation = read.csv(file = "../Publication/240923_Nature Genetics_Revise/4 - WYH/MutationMatrix_update_241106/onco_matrix_105_patient.csv", row.names = 1)
selected_Mutation = names(sort(apply(Pop_Mutation,1,function(x) length(grep("_",x))),decreasing = T)[1:20])

# selected_Mutation = c("TP53","CDKN2A","ACOT12","E2F8","JUN","HRC","PRKCB","STK11","GBP3","ELF3","WWP1","VSTM2A","MAP3K8","ZNF322",
#                       "AR","METTL7B","CD1A","XCL2","ZNF841","SOX15")

Pop_Mutation = Pop_Mutation[selected_Mutation,]
Pop_Mutation = as.data.frame(t(Pop_Mutation))
Pop_Mutation = Pop_Mutation[rownames(Pop_Mutation) %in% annotation$NewSample.ID,]
Pop_Mutation_MutHeatmap = Pop_Mutation
Pop_Mutation = apply(Pop_Mutation, 2, function(x) ifelse(is.na(x),"No",ifelse(x=="","No","Yes")))
Pop_Mutation = as.data.frame(Pop_Mutation)
Pop_Mutation$NewSample.ID = rownames(Pop_Mutation)
annotation = left_join(annotation, Pop_Mutation, by = "NewSample.ID")
rownames(annotation) = annotation[,1]
annotation = annotation[,-1]
annotation$Clinical.stage = ifelse(annotation$Clinical.stage %in% c("IIA","IIB"), "II", 
                                          annotation$Clinical.stage)
annotation$Age = ifelse(annotation$Age < 50, "<50", ifelse(annotation$Age > 65, ">65", "50~65"))
annotation$Sex = factor(annotation$Sex, levels = c("M","F"))
annotation$Age = factor(annotation$Age, levels = c("<50","50~65",">65"))
annotation_col = list(Clinical.stage = c(I = "#b2f2a5",
                                         II = "#a5eaf2",
                                         IIIA = "#eef2a5",
                                         IIIB = "#edad82",
                                         IVA = "#d754de",
                                         IVB = "#915f19"),
                      Group = c(P = "#09a5ed",
                                P_Mets = "#bbbcbd"),
                      Sex = c(M = mycol[19],
                              F = mycol[12]),
                      
                      #total = colorRampPalette(c("white", "red"))(100),
                      TMB = colorRampPalette(c("white", "purple"))(100),
                      
                      TP53 = c("Yes" = mycol[10],
                                 "No" = mycol[20]),
                      MUC5B = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      ARID1A = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      FLG = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      MUC16 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      TTN = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      CDKN2A = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      MAGI1 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      PLEC = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      MUC12 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      TLE4 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      DISP3 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      CACNA1C = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      FAT2 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      PRRC2B = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      SYNE1 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      ACOT12 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      AHNAK2 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      E2F8 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      FLG2 = c("Yes" = mycol[10],
                               "No" = mycol[20]),
                      
                      Polyps = c("Yes" = mycol[10],
                                   "No" = mycol[20]),
                      Age = c("<50" = mycol[2], 
                                "50~65" = mycol[18], 
                                ">65" = mycol[8]),
                      Gallstones = c("Yes" = mycol[10],
                                 "No" = mycol[20]),
                      PBM = c("Yes" = mycol[10],
                                     "No" = mycol[20]),
                      Atrophic.cholecystitis = c("Yes" = mycol[10],
                                     "No" = mycol[20]),
                      liver.invasion = c("Yes" = mycol[10],
                                     "No" = mycol[20]),
                      Bile.duct.invasion = c("Yes" = mycol[10],
                                     "No" = mycol[20]),
                      Vascular.invasion = c("Yes" = mycol[10],
                                     "No" = mycol[20]),
                      Lymph.node.metastasis = c("Yes" = mycol[10],
                                            "No" = mycol[20]),
                      Liver.metastasis = c("Yes" = mycol[10],
                                            "No" = mycol[20]),
                      Peritoneal..metastasis = c("Yes" = mycol[10],
                                            "No" = mycol[20]),
                      `Differentiation` = c("Good" = mycol[15],
                                            "Moderate" = mycol[22],
                                            "Moderate & poor" = mycol[7],
                                            "Poor" = mycol[16]))
bk <- c(seq(-0.5,0,by=0.01),seq(0.01,1,by=0.01))
p = pheatmap(Freq_adj,
         scale = "none",
         #cluster_rows = F,
         cutree_cols = 5,
         annotation_col = annotation,
         annotation_colors = annotation_col,
         clustering_distance_rows="manhattan",
         clustering_method = "ward.D2",
         clustering_distance_cols="manhattan",
         color = c(colorRampPalette(colors = c("#34ebc3","white","#eb346b"))(150)),
         legend_breaks = seq(-0.5,1,0.5),
         breaks = bk)
# pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/PatientClustering.pdf", width =14, height = 7)
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/PatientClustering_anno_TMB-241123.pdf", width =18, height = 15)
print(p)
dev.off()

library("dendextend")
hclust = hclust(dist(t(Freq_adj), method = "manhattan"),method = "ward.D2")
plot(color_branches(hclust,h=27.3,groupLabels = T))
patient_group = cutree(hclust, h = 27.3, order_clusters_as_data = F)
patient_group = data.frame(NewSample.ID = names(patient_group),
                           group = patient_group)
saveRDS(patient_group, file = "./3 - GBC_output data_subtype/11_PatientClustering/patient_group_5MI.RDS")

Pop_Mutation_MutHeatmap$NewSample.ID = rownames(Pop_Mutation_MutHeatmap)
Pop_Mutation_MutHeatmap = left_join(patient_group, Pop_Mutation_MutHeatmap, by = "NewSample.ID")
Pop_Mutation_MutHeatmap$group = paste0("MI",Pop_Mutation_MutHeatmap$group)
Pop_Mutation_MutHeatmap[is.na(Pop_Mutation_MutHeatmap)] = ""


# patient_order = arrange(Pop_Mutation_MutHeatmap,group,desc(TP53),desc(CDKN2A),desc(ACOT12),desc(E2F8),desc(JUN),desc(HRC),desc(PRKCB)
#                         ,desc(STK11),desc(GBP3),desc(ELF3),desc(WWP1),desc(VSTM2A),desc(MAP3K8),desc(ZNF322),desc(AR),desc(METTL7B),
#                         desc(CD1A),desc(XCL2),desc(ZNF841),desc(SOX15))

patient_order = arrange(Pop_Mutation_MutHeatmap,group,desc(TP53))

rownames(Pop_Mutation_MutHeatmap) = Pop_Mutation_MutHeatmap$NewSample.ID
Pop_Mutation_MutHeatmap = Pop_Mutation_MutHeatmap[patient_order$NewSample.ID,]

library("ComplexHeatmap")
ha = HeatmapAnnotation(group = Pop_Mutation_MutHeatmap$group,
                       col = list(group = c("MI1" = "#915F19",
                                            "MI2" = "#FF7F0E",
                                            "MI3" = "#D754DE",
                                            "MI4" = "#09A5ED",
                                            "MI5" = "#34EBC3")))

color_hp = c("#AA40FC","#D62728","#E377C2","#17BECF","#279E68","#AEC7E8","#1F77B4","#FF7F0E")
names(color_hp) = c("Nonsense_Mutation","Missense_Mutation","In_Frame_Del",
                    "Multi_Hit","Frame_Shift_Ins","Frame_Shift_Del","In_Frame_Ins","Translation_Start_Site")

Pop_Mutation_MutHeatmap_plot = Pop_Mutation_MutHeatmap[,c(3:22)]
Pop_Mutation_MutHeatmap_plot = t(Pop_Mutation_MutHeatmap_plot)
identical(colnames(Pop_Mutation_MutHeatmap_plot), Pop_Mutation_MutHeatmap$NewSample.ID)
p = Heatmap(Pop_Mutation_MutHeatmap_plot,
            col = color_hp,
            top_annotation = ha,
            na_col = "white")
pdf(file = "./3 - GBC_output data_subtype/12_Mutation_Pct_Cor/Mut_Heatmap_New-241212.pdf", width =8, height = 5)
print(p)
dev.off()

write.csv(Pop_Mutation_MutHeatmap_plot, file = "./3 - GBC_output data_subtype/12_Mutation_Pct_Cor/Mut_Heatmap_New-241123.csv")

## * calculate chi square ####
annotation_test = annotation
annotation_test$NewSample.ID = rownames(annotation_test)

temp = subset(patient_group, select = c("NewSample.ID", "group"))
annotation_test = left_join(annotation_test, temp, by = "NewSample.ID")
annotation_test$group1 = annotation_test$group

#### 75 patients
library(gmodels)
CrossTable(annotation_test$group1, annotation_test$Clinical.stage,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Clinical.stage, simulate.p.value=TRUE) 

CrossTable(annotation_test$group1, annotation_test$Group,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Group, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Sex,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Sex, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Age,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Age, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Polyps,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Polyps, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Gallstones,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Gallstones, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$PBM,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$PBM, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Atrophic.cholecystitis,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Atrophic.cholecystitis, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Differentiation,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Differentiation, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$liver.invasion,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$liver.invasion, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Bile.duct.invasion,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Bile.duct.invasion, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Vascular.invasion,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Vascular.invasion, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Lymph.node.metastasis,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Lymph.node.metastasis, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Liver.metastasis,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Liver.metastasis, simulate.p.value=TRUE)

#### 75 patients
CrossTable(annotation_test$group1, annotation_test$TP53,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$TP53, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$MUC5B,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$MUC5B, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$ARID1A,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$ARID1A, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$FLG,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$FLG, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$MUC16,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$MUC16, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$TTN,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$TTN, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$CDKN2A,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$CDKN2A, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$MAGI1,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$MAGI1, simulate.p.value=TRUE)
# annotation_df = subset(annotation_test, select = c("group","MAGI1"))
# annotation_df$MAGI1[is.na(annotation_df$MAGI1)] = "No"
# df_cell = annotation_df %>% group_by(group,MAGI1) %>% summarise(Value = n())
# p = ggplot(data = df_cell, mapping = aes(x = factor(group), y = Value, fill = MAGI1)) +
#   geom_bar(stat = 'identity', position = 'fill') +
#   scale_fill_manual(values=c('grey','#3594f2')) + ggtitle("MAGI1 p-value = 0.04798") + theme_classic()
# pdf(file = "./3 - GBC_output data_subtype/17_Mutation - MI_p/240408_MAGI1.pdf", width = 4, height = 3)
# print(p)
# dev.off()

CrossTable(annotation_test$group1, annotation_test$PLEC,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$PLEC, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$MUC12,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$MUC12, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$TLE4,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$TLE4, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$DISP3,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$DISP3, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$CACNA1C,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$CACNA1C, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$FAT2,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$FAT2, simulate.p.value=TRUE)

# annotation_df = subset(annotation_test, select = c("group","FAT2"))
# annotation_df$FAT2[is.na(annotation_df$FAT2)] = "No"
# df_cell = annotation_df %>% group_by(group,FAT2) %>% summarise(Value = n())
# p = ggplot(data = df_cell, mapping = aes(x = factor(group), y = Value, fill = FAT2)) +
#   geom_bar(stat = 'identity', position = 'fill') +
#   scale_fill_manual(values=c('grey','#3594f2')) + ggtitle("FAT2 p-value = 0.02799") + theme_classic()
# pdf(file = "./3 - GBC_output data_subtype/17_Mutation - MI_p/240408_FAT2.pdf", width = 4, height = 3)
# print(p)
# dev.off()

CrossTable(annotation_test$group1, annotation_test$PRRC2B,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$PRRC2B, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$SYNE1,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$SYNE1, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$ACOT12,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$ACOT12, simulate.p.value=TRUE)
# annotation_df = subset(annotation_test, select = c("group","ACOT12"))
# annotation_df$ACOT12[is.na(annotation_df$ACOT12)] = "No"
# df_cell = annotation_df %>% group_by(group,ACOT12) %>% summarise(Value = n())
# p = ggplot(data = df_cell, mapping = aes(x = factor(group), y = Value, fill = ACOT12)) +
#   geom_bar(stat = 'identity', position = 'fill') +
#   scale_fill_manual(values=c('grey','#3594f2')) + ggtitle("ACOT12 p-value = 0.03748") + theme_classic()
# pdf(file = "./3 - GBC_output data_subtype/17_Mutation - MI_p/240408_ACOT12.pdf", width = 4, height = 3)
# print(p)
# dev.off()

CrossTable(annotation_test$group1, annotation_test$AHNAK2,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$AHNAK2, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$E2F8,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$E2F8, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$FLG2,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$FLG2, simulate.p.value=TRUE)


# ##
# annotation_test$M1 = ifelse(annotation_test$group1 == 1, "M1", "Others")
# annotation_test$M2 = ifelse(annotation_test$group1 == 2, "M2", "Others")
# annotation_test$M3 = ifelse(annotation_test$group1 == 3, "M3", "Others")
# annotation_test$M4 = ifelse(annotation_test$group1 == 4, "M4", "Others")
# annotation_test$M5 = ifelse(annotation_test$group1 == 5, "M5", "Others")
# 
# CrossTable(annotation_test$M1, annotation_test$TP53,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
# set.seed(202308)
# fisher.test(annotation_test$M1, annotation_test$TP53, simulate.p.value=TRUE)
# 
# CrossTable(annotation_test$M2, annotation_test$TP53,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T) # p-value = 0.07644
# set.seed(202308)
# fisher.test(annotation_test$M2, annotation_test$TP53, simulate.p.value=TRUE)
# 
# CrossTable(annotation_test$M3, annotation_test$TP53,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
# set.seed(202308)
# fisher.test(annotation_test$M3, annotation_test$TP53, simulate.p.value=TRUE)
# 
# CrossTable(annotation_test$M4, annotation_test$TP53,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
# set.seed(202308)
# fisher.test(annotation_test$M4, annotation_test$TP53, simulate.p.value=TRUE)
# 
# CrossTable(annotation_test$M5, annotation_test$TP53,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
# set.seed(202308)
# fisher.test(annotation_test$M5, annotation_test$TP53, simulate.p.value=TRUE)

####
compare_means(TMB~group1, annotation_test, method = "anova") # NS
# compare_means(total_perMB~group1, annotation_test, method = "anova")

####
Freq = Freq[,subtype]
Freq$NewSample.ID = rownames(Freq)
Mutation_Pct = left_join(annotation_test, Freq, by = "NewSample.ID")

for (i in 1:20) {
  Mutation_Pct1 = subset(Mutation_Pct, select =c(colnames(Pop_Mutation)[i],subtype))
  Mutation_Pct1 = melt(Mutation_Pct1)
  
  for (j in 1:31) {
    temp = Mutation_Pct1[Mutation_Pct1$variable %in% subtype[j],]
    colnames(temp)[1] = "Group"
    temp = temp[!is.na(temp$Group),]
    aa = compare_means(value~Group,temp)
    if (aa$p.adj < 0.05) {
      p <- ggboxplot(temp, x = "Group", y = "value", color = "Group", palette = "jco",add = "jitter") + 
        stat_compare_means(label = "p.format", hide.ns = F)+
        theme(legend.position = "top") +
        theme(axis.title.x = element_blank()) +
        theme(axis.text.x = element_blank()) +
        ggtitle(paste0(colnames(Pop_Mutation)[i],"_",subtype[j]))
      pdf(file = paste0("./3 - GBC_output data_subtype/12_Mutation_Pct_Cor/", colnames(Pop_Mutation)[i],"_",subtype[j], ".pdf"), width =2, height = 3)
      print(p)
      dev.off()
    }
  }
}

#### summary for Mut-Pct ####
GBC_AllSubtype <- read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
GBC_TIMESubtype = GBC_AllSubtype[!(GBC_AllSubtype$celltype %in% c("tumor","normal")),]
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_DC = PatientInfo[PatientInfo$X %in% unique(GBC_TIMESubtype$orig.ident),]
PatientInfo_metastasis = PatientInfo_DC[PatientInfo_DC$histological.type.short %in% c("adeno"),]
PatientInfo_metastasis = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% c("P"),]
PatientInfo_metastasis$Group = PatientInfo_metastasis$metastasis.type
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group == "P_LI", "P", PatientInfo_metastasis$Group)
PatientInfo_metastasis$Group = ifelse(PatientInfo_metastasis$Group %in% c("P_LN","P_LM"), "P_Mets", PatientInfo_metastasis$Group)
Select_Patient = PatientInfo_metastasis$NewSample.ID
cycle = Select_Patient
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(GBC_TIMESubtype, orig.ident == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$subtype)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

Pop_Mutation = read.csv(file = "../Publication/240923_Nature Genetics_Revise/4 - WYH/MutationMatrix_update_241106/onco_matrix_105_patient.csv", row.names = 1)

selected_Mutation = names(sort(apply(Pop_Mutation,1,function(x) length(grep("_",x))),decreasing = T)[1:20])
# selected_Mutation = c("TP53","CDKN2A","ACOT12","E2F8","JUN","HRC","PRKCB","STK11","GBP3","ELF3","WWP1","VSTM2A","MAP3K8","ZNF322",
#                       "AR","METTL7B","CD1A","XCL2","ZNF841","SOX15")

Pop_Mutation = Pop_Mutation[selected_Mutation,]
Pop_Mutation = as.data.frame(t(Pop_Mutation))
Pop_Mutation = Pop_Mutation[rownames(Pop_Mutation) %in% Select_Patient,]
Pop_Mutation = apply(Pop_Mutation, 2, function(x) ifelse(is.na(x),"No",ifelse(x=="","No","Yes")))
Pop_Mutation = as.data.frame(Pop_Mutation)

yes_count <- apply(Pop_Mutation, 2, function(x) sum(x == "Yes"))
selected_Mutation <- names(yes_count[yes_count >= 3])

subtype = c("CD4T_C1_FOXP3",
            "CD4T_C7_ISG15",
            #"CD4T_C8_MKI67",
            "CD8T_C1_CXCL13",
            #"CD8T_C8_MKI67",
            #"CD8T_C11_PCLAF",
            "γdT_C12_TRDC",
            "NKT_C3_CAPG",
            #"NK_C6_MKI67",
            "F_C0_MMP11",
            "Per_C0_RGS5",
            "DC_C0_cDC2_IL1B",
            "DC_C2_cDC3_FSCN1",
            "DC_C5_cDC1",
            "DC_C6_pDC",
            "M_C2_SPP1",
            "M_C7_AGR2",
            "N_C9_Nc",
            "EC_C2_CXCR4",
            "CD4T_C2_CCR7",
            #"CD4T_C4_Unassign",
            "CD8T_C0_CCR7_GZMK",
            "B_C0_IGHA1",
            "F_C1_CFD",
            "Per_C1_MYH11",
            "DC_C4_cDC2_PPP1R14A",
            #"M_C4_others",
            "N_C0_N0",
            "N_C1_N2",
            "N_C2_N1",
            "N_C7_N0",
            "EC_C3_GJA5",
            "EC_C5_PROX1",
            "Per_C3_STEAP4",
            "M_C1_S100A8",
            "CD8T_C9_MT1X_MT1E",
            "EC_C0_ACKR1")

cell.ct = Freq[,subtype]
cell.ct = cell.ct[Select_Patient,]
cell.ct$PatientID = rownames(cell.ct)
Pop_Mutation$PatientID = rownames(Pop_Mutation)
cell.ct = left_join(cell.ct, Pop_Mutation, by = "PatientID")
cell.ct[is.na(cell.ct)] = "No"

df = data.frame()
final.df = data.frame()
for (i in 1:length(selected_Mutation)){
  cycle.ct = cell.ct[,c(subtype, selected_Mutation[i])]
  colnames(cycle.ct)[ncol(cycle.ct)] = "group"
  
  for (j in 1:31){
    test_result <- compare_means(as.formula(paste(subtype[j], "~ group")), data = cycle.ct, method = 'wilcox.test', ref.group = "No")
    p_value = test_result$p
    
    mean_yes = mean(cycle.ct[cycle.ct$group == "Yes",j])
    mean_no = mean(cycle.ct[cycle.ct$group == "No",j])
    
    df[j,1] = selected_Mutation[i]
    df[j,2] = subtype[j]
    df[j,3] = p_value
    df[j,4] = ifelse(mean_yes > mean_no, "Up", "Down")
  }
  
  final.df = rbind(final.df,df)
  df = data.frame()
  
}

colnames(final.df) = c('Mutation','Subtype','P','Sig')
final.df$Sig = ifelse(final.df$P < 0.05, final.df$Sig, "NS")
final.df$Sig = factor(final.df$Sig, levels = c("Up", "Down", "NS"))
final.df$p_value_log <- -log10(final.df$P)

final.df$Mutation <- factor(final.df$Mutation, levels = names(sort(yes_count)))

p = ggplot(final.df, aes(x = Subtype, y = Mutation, size = p_value_log)) +
  geom_point(aes(color = Sig), alpha = 0.7) +  # 添加气泡
  scale_size_continuous(range = c(1, 10)) +  # 设置气泡大小范围
  labs(title = "MI-specific subtypes") +  # 添加标题和标签
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.3),
        legend.box = "vertical")

pdf(file = "./3 - GBC_output data_subtype/12_Mutation_Pct_Cor/Mut_Pct_Heatmap_New-241212.pdf", width =12, height = 4.5)
print(p)
dev.off()


write.csv(final.df, file = "./3 - GBC_output data_subtype/12_Mutation_Pct_Cor/Mut_Pct_Heatmap_New-241126.csv")


## * add GM annotation for each patient ####
GM_percentage = read.csv("./4 - Cooperation/WYH/GM_patient_ratio.csv")
GM_percentage = melt(GM_percentage)

colnames(GM_percentage)
patient_group$NewSample.ID
GM_percentage$X = factor(GM_percentage$X, levels = patient_group$NewSample.ID)
GM_percentage$variable = factor(GM_percentage$variable, levels = c("GM16","GM1","GM3","GM5","GM6","GM7","GM8"))

p = ggplot(data = GM_percentage, aes(x = X, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme(panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = NULL) +
  #theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1)) + 
  # coord_flip() + 
  scale_fill_manual(values=c("GM16" = "#ff7f0e",
                             "GM1" = "#c7e3a3",
                             "GM3" = "#a3e3ce",
                             "GM5" = "#a3c4e3",
                             "GM6" = "#c3a3e3",
                             "GM7" = "#e3a3c2",
                             "GM8" = "#f7f7b0"))

pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Histo_GM.pdf", width =14, height = 2)
print(p)
dev.off()

## Fig7C ####
library(ggraph)
library(igraph)
library(tidygraph)
library(network)

# jaccard <- function(a, b) {
#   intersection = length(intersect(a, b))
#   union = length(a) + length(b) - intersection
#   return (intersection/union)
# }

##
M1_subtype = c("EC_C3_GJA5","EC_C0_ACKR1","F_C1_CFD","Per_C1_MYH11","Per_C3_STEAP4","DC_C0_cDC2_IL1B")
M2M5_subtype = c("CD4T_C1_FOXP3","CD8T_C1_CXCL13","CD4T_C2_CCR7","B_C0_IGHA1","γdT_C12_TRDC","CD8T_C0_CCR7_GZMK")
M3_subtype = c("N_C1_N2","N_C2_N1","N_C0_N0","N_C7_N0","NKT_C3_CAPG","N_C9_Nc")
M4_subtype = c("M_C7_AGR2","Per_C0_RGS5","EC_C2_CXCR4","M_C2_SPP1","EC_C5_PROX1","F_C0_MMP11")

## MI1
MI_group = annotation_test[annotation_test$group1 == "MI1",]$NewSample.ID

Freq_adj_t = t(Freq_adj)
Freq_adj_t[Freq_adj_t<-0.5] = -0.5
Freq_adj_t[Freq_adj_t>1] = 1

Freq_Jac = Freq_adj_t[MI_group,1:31]
Freq_Jac = round(Freq_Jac,2)

Jaccard_df = data.frame()
for (i in 1:30) {
  for (j in (i+1):31) {
    sub <- Freq_Jac[,c(i,j)]
    sub[sub>=0.5] = 1
    sub[sub<0.5] = 0
    intersection = sum(rowSums(sub)==2)
    union = sum(rowSums(sub)>=1)
    temp = data.frame(Subtype1 = colnames(Freq_Jac)[i],
                      Subtype2 = colnames(Freq_Jac)[j],
                      Weight = round(intersection/union,4))
    Jaccard_df = rbind(Jaccard_df, temp)
  }
}

Jaccard_df = Jaccard_df[Jaccard_df$Subtype1 %in% M1_subtype,]
Jaccard_df = Jaccard_df[Jaccard_df$Subtype2 %in% M1_subtype,]

Jaccard_df$Weight = Jaccard_df$Weight*10
Jaccard_net = as.network(Jaccard_df)
g = ggnet2(Jaccard_net, mode = "circle", label = T, edge.size = "Weight", node.color = "Subtype1", node.size = 30) + NoLegend() + coord_equal()
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Jaccard_MI1.pdf", width =4, height = 4)
print(g)
dev.off()

## MI2
MI_group = annotation_test[annotation_test$group1 == "MI2",]$NewSample.ID

Freq_adj_t = t(Freq_adj)
Freq_adj_t[Freq_adj_t<-0.5] = -0.5
Freq_adj_t[Freq_adj_t>1] = 1

Freq_Jac = Freq_adj_t[MI_group,1:31]
Freq_Jac = round(Freq_Jac,2)

Jaccard_df = data.frame()
for (i in 1:30) {
  for (j in (i+1):31) {
    sub <- Freq_Jac[,c(i,j)]
    sub[sub>=0.5] = 1
    sub[sub<0.5] = 0
    intersection = sum(rowSums(sub)==2)
    union = sum(rowSums(sub)>=1)
    temp = data.frame(Subtype1 = colnames(Freq_Jac)[i],
                      Subtype2 = colnames(Freq_Jac)[j],
                      Weight = round(intersection/union,4))
    Jaccard_df = rbind(Jaccard_df, temp)
  }
}

Jaccard_df = Jaccard_df[Jaccard_df$Subtype1 %in% M2M5_subtype,]
Jaccard_df = Jaccard_df[Jaccard_df$Subtype2 %in% M2M5_subtype,]

Jaccard_df$Weight = Jaccard_df$Weight*10
Jaccard_net = as.network(Jaccard_df)
g = ggnet2(Jaccard_net, mode = "circle", label = T, edge.size = "Weight", node.color = "Subtype1",node.size = 30) + NoLegend() + coord_equal()
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Jaccard_MI2.pdf", width =4, height = 4)
print(g)
dev.off()

## MI3
MI_group = annotation_test[annotation_test$group1 == "MI3",]$NewSample.ID

Freq_adj_t = t(Freq_adj)
Freq_adj_t[Freq_adj_t<-0.5] = -0.5
Freq_adj_t[Freq_adj_t>1] = 1

Freq_Jac = Freq_adj_t[MI_group,1:31]
Freq_Jac = round(Freq_Jac,2)

Jaccard_df = data.frame()
for (i in 1:30) {
  for (j in (i+1):31) {
    sub <- Freq_Jac[,c(i,j)]
    sub[sub>=0.5] = 1
    sub[sub<0.5] = 0
    intersection = sum(rowSums(sub)==2)
    union = sum(rowSums(sub)>=1)
    temp = data.frame(Subtype1 = colnames(Freq_Jac)[i],
                      Subtype2 = colnames(Freq_Jac)[j],
                      Weight = round(intersection/union,4))
    Jaccard_df = rbind(Jaccard_df, temp)
  }
}

Jaccard_df = Jaccard_df[Jaccard_df$Subtype1 %in% M3_subtype,]
Jaccard_df = Jaccard_df[Jaccard_df$Subtype2 %in% M3_subtype,]
Jaccard_df$Weight[Jaccard_df$Weight==0] = 0.001

Jaccard_df$Weight = Jaccard_df$Weight*10
Jaccard_net = as.network(Jaccard_df)
g = ggnet2(Jaccard_net, mode = "circle", label = T, edge.size = "Weight", node.color = "Subtype1",node.size = 30) + NoLegend() + coord_equal()
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Jaccard_MI3.pdf", width =4, height = 4)
print(g)
dev.off()

## MI4
MI_group = annotation_test[annotation_test$group1 == "MI4",]$NewSample.ID

Freq_adj_t = t(Freq_adj)
Freq_adj_t[Freq_adj_t<-0.5] = -0.5
Freq_adj_t[Freq_adj_t>1] = 1

Freq_Jac = Freq_adj_t[MI_group,1:31]
Freq_Jac = round(Freq_Jac,2)

Jaccard_df = data.frame()
for (i in 1:30) {
  for (j in (i+1):31) {
    sub <- Freq_Jac[,c(i,j)]
    sub[sub>=0.5] = 1
    sub[sub<0.5] = 0
    intersection = sum(rowSums(sub)==2)
    union = sum(rowSums(sub)>=1)
    temp = data.frame(Subtype1 = colnames(Freq_Jac)[i],
                      Subtype2 = colnames(Freq_Jac)[j],
                      Weight = round(intersection/union,4))
    Jaccard_df = rbind(Jaccard_df, temp)
  }
}

Jaccard_df = Jaccard_df[Jaccard_df$Subtype1 %in% M4_subtype,]
Jaccard_df = Jaccard_df[Jaccard_df$Subtype2 %in% M4_subtype,]
Jaccard_df$Weight[Jaccard_df$Weight==0] = 0.001

Jaccard_df$Weight = Jaccard_df$Weight*10
Jaccard_net = as.network(Jaccard_df)
g= ggnet2(Jaccard_net, mode = "circle", label = T, edge.size = "Weight", node.color = "Subtype1",node.size = 30) + NoLegend() + coord_equal()
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Jaccard_MI4.pdf", width =4, height = 4)
print(g)
dev.off()

## MI5
MI_group = annotation_test[annotation_test$group1 == "MI5",]$NewSample.ID

Freq_adj_t = t(Freq_adj)
Freq_adj_t[Freq_adj_t<-0.5] = -0.5
Freq_adj_t[Freq_adj_t>1] = 1

Freq_Jac = Freq_adj_t[MI_group,1:31]
Freq_Jac = round(Freq_Jac,2)

Jaccard_df = data.frame()
for (i in 1:30) {
  for (j in (i+1):31) {
    sub <- Freq_Jac[,c(i,j)]
    sub[sub>=0.5] = 1
    sub[sub<0.5] = 0
    intersection = sum(rowSums(sub)==2)
    union = sum(rowSums(sub)>=1)
    temp = data.frame(Subtype1 = colnames(Freq_Jac)[i],
                      Subtype2 = colnames(Freq_Jac)[j],
                      Weight = round(intersection/union,4))
    Jaccard_df = rbind(Jaccard_df, temp)
  }
}

Jaccard_df = Jaccard_df[Jaccard_df$Subtype1 %in% M2M5_subtype,]
Jaccard_df = Jaccard_df[Jaccard_df$Subtype2 %in% M2M5_subtype,]
Jaccard_df$Weight[Jaccard_df$Weight==0] = 0.001

Jaccard_df$Weight = Jaccard_df$Weight*10
Jaccard_net = as.network(Jaccard_df)
g = ggnet2(Jaccard_net, mode = "circle", label = T, edge.size = "Weight", node.color = "Subtype1",node.size = 30) + NoLegend() + coord_equal()
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Jaccard_MI5.pdf", width =4, height = 4)
print(g)
dev.off()

## legend
Jaccard_df = data.frame(Subtype1 = c("A","B","C","D","E"),
           Subtype2 =  c("B","C","D","E","A"),
           Weight = c(0.2,0.4,0.6,0.8,1))
Jaccard_df$Weight = Jaccard_df$Weight*10
Jaccard_net = as.network(Jaccard_df)
g = ggnet2(Jaccard_net, mode = "circle", label = T, edge.size = "Weight", node.color = "Subtype1",node.size = 30) + NoLegend() + coord_equal()
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Jaccard_legend.pdf", width =4, height = 4)
print(g)
dev.off()

## Fig7D ####
library("dendextend")
hclust = hclust(dist(t(Freq_adj), method = "manhattan"),method = "ward.D2")
plot(color_branches(hclust,h=27.3,groupLabels = T))
patient_group = cutree(hclust, h = 27.3, order_clusters_as_data = F)
patient_group = data.frame(NewSample.ID = names(patient_group),
                           group = patient_group)
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20240916预后更新.xlsx")
PatientInfo_surv = subset(PatientInfo_surv, NewSample.ID %in% patient_group$NewSample.ID)
patient_group = left_join(patient_group,PatientInfo_surv,by="NewSample.ID")
patient_group$event[is.na(patient_group$event)] = 0
patient_group$event01 = ifelse(patient_group$event=="dead", 1, 0)
patient_group = left_join(patient_group,PatientInfo_DC,by="NewSample.ID")
PatientInfo_metastasis = patient_group[patient_group$histological.type.short %in% c("adeno"),]
patient_group = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% c("P"),]
patient_group$group1 = patient_group$group
patient_group$group = ifelse(patient_group$group1 ==5,"MI5","other MIs")

c("MI1" = "#915F19",
  "MI2" = "#FF7F0E",
  "MI3" = "#D754DE",
  "MI4" = "#09A5ED",
  "MI5" = "#34EBC3")

fit<-survfit(Surv(OS_month, event01)~group, data=patient_group)
g=ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = patient_group,             # data used to fit survival curves.
  palette = c("#34EBC3","grey"),
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in months",   # customize X axis label.
  risk.table.pos = "in",
  risk.table.col = "strata",
  fontsize = 4,
  pval.size = 4,
  pval.coord = c(9, 0.9),
  legend = "top",
  title = "MI5"
  #ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/New241218_Prognosis_MI5.pdf", width =6, height = 6, onefile = F)
print(g)
dev.off()

## Fig7E ####
base_path <- './4 - Cooperation/WYH/Subtype_signature/subtype_DE/'
dirnames <- stringr::str_sort(list.files(base_path),numeric = TRUE)
dirnames_list <- vector(mode = "list", length = length(dirnames))
i=1
while(i <= length(dirnames)){
  base_path_TMA = file.path(base_path,dirnames[i])
  dirnames_list[i]=base_path_TMA
  i=i+1
}

df_list <- list()
i=1
while(i <= length(dirnames)){
  aa = read.csv(dirnames_list[[i]])
  aa = aa[order(aa$pvals_adj),]
  df_list[[i]] = aa$names[1:30]
  i=i+1
}

dirnames_split = stringr::str_split_fixed(dirnames,"_DE.csv",n=2)
dirnames_split = dirnames_split[,1]
names(df_list) = dirnames_split

downsample_expr <- readRDS("./4 - Cooperation/WYH/Subtype_signature/downsample_expr.RDS")
downsample_expr_sub = downsample_expr[unique(unlist(df_list)),]
downsample_expr_sub = downsample_expr_sub@assays$RNA@data
downsample_expr_sub = as.matrix(downsample_expr_sub)

cellid = read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
cellid = cellid[cellid$cellid %in% colnames(downsample_expr_sub),]
identical(sort(cellid$cellid), sort(colnames(downsample_expr_sub)))

downsample_expr_sub = downsample_expr_sub[,cellid$cellid]
identical(cellid$cellid, colnames(downsample_expr_sub))
colnames(downsample_expr_sub) = cellid$subtype

downsample_expr_sub_avg = data.frame()
for(i in 1:length(unique(cellid$subtype))){
  downsample_expr_sub_avg = rbind(downsample_expr_sub_avg, rowMeans(downsample_expr_sub[,grep(unique(cellid$subtype)[i],colnames(downsample_expr_sub))]))
}
rownames(downsample_expr_sub_avg) = unique(cellid$subtype)
colnames(downsample_expr_sub_avg) = rownames(downsample_expr_sub)
downsample_expr_sub_avg = t(as.matrix(downsample_expr_sub_avg))

write.table(downsample_expr_sub_avg, file = "./downsample_expr_sub.txt" ,sep = '\t')

aa = read.table("./4 - Cooperation/WYH/Subtype_signature/downsample_expr_sub.txt")
aa$GeneSymbol = rownames(aa)
aa = aa[,c(129,1:128)]
aa = as.matrix(aa)
write.table(aa, file = "./4 - Cooperation/WYH/Subtype_signature/downsample_expr_sub.txt" ,sep = '\t')

## CIBERSORTx

RNA_DeconResult = readRDS("./3 - GBC_output data_subtype/13_Deconvolution_RNAseq/CIBERSORTx_Job5_Results.rds")
RNA_DeconResult = RNA_DeconResult[,1:119]
RNA_DeconResult = as.data.frame(t(RNA_DeconResult))
RNA_DeconResult = apply(RNA_DeconResult,2,function(x){x/sum(x)})

ExternalCohort_PctScore = CreateSeuratObject(RNA_DeconResult)

selected_subtypes = list(c("EC-C3-GJA5","EC-C0-ACKR1","F-C1-CFD","Per-C1-MYH11","Per-C3-STEAP4"))
selected_subtypes = list(c("EC-C3-GJA5","EC-C0-ACKR1","F-C1-CFD"))
ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = selected_subtypes, ctrl = 5, name = "MI1_subtype")

selected_subtypes = list(c("CD4T-C2-CCR7","B-C0-IGHA1","γdT-C12-TRDC","CD8T-C0-CCR7-GZMK"))
ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = selected_subtypes, ctrl = 5, name = "MI2_subtype")

selected_subtypes = list(c("CD4T-C1-FOXP3","CD8T-C1-CXCL13"))
ExternalCohort_PctScore = AddModuleScore(ExternalCohort_PctScore, features = selected_subtypes, ctrl = 5, name = "MI5_subtype")

library(survminer)
library("survival")
patient_group = ExternalCohort_PctScore@meta.data
patient_group$NewSample.ID = rownames(patient_group)
PatientInfo_surv = openxlsx::read.xlsx("../8 - GBC_external RNA-seq/NC_gbc_survival_GBCKorea 20210618.xlsx")
PatientInfo_surv = subset(PatientInfo_surv, Patient_Id %in% colnames(ExternalCohort_PctScore))
colnames(PatientInfo_surv)[2] = "NewSample.ID"
patient_group = left_join(PatientInfo_surv,patient_group,by="NewSample.ID")
patient_group$vital_status = ifelse(patient_group$vital_status=="dead", 1, 0)
patient_group$group = ifelse(patient_group$MI1_subtype1 > median(patient_group$MI1_subtype1), "High_MI1", "Low_MI1")
fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=patient_group)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = patient_group,             # data used to fit survival curves.
  palette = c("#34ebc3","#bbbcbd"),
  #risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in months",   # customize X axis label.
  ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)

patient_group$group = ifelse(patient_group$MI2_subtype > median(patient_group$MI2_subtype1), "High_MI2", "Low_MI2")
fit<-survfit(Surv(months_to_last_fu, vital_status)~group, data=patient_group)
g = ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = patient_group,             # data used to fit survival curves.
  palette = c("#ff7f0e","#bbbcbd"),
  #risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in months",   # customize X axis label.
  ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Validate_Surv_MI2.pdf", width =3, height = 3, onefile = F)
print(g)
dev.off()

# Fig 7a ####
Fibroblast_new2 <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/Fibroblast/Fibroblast_new2.RDS")
Fb = subset(Fibroblast_new2, idents = c("F_C1_CFD", "Per_C1_MYH11", "Per_C3_STEAP4"))

Endothelium <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/04_Endothelium/Endothelium.RDS")
Endo = subset(Endothelium, idents = c("EC_C0_ACKR1", "EC_C5_PROX1", "EC_C3_GJA5"))

Myeloid_cell_DC <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/Myeloid_cell_DC.RDS")
DC = subset(Myeloid_cell_DC, idents = c(4))
Idents(DC) = "DC_C4_PPP1R14A_cDC2"

Myeloid_cell_MonoMacro <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/Myeloid_cell_MonoMacro.RDS")
MM = subset(Myeloid_cell_MonoMacro, idents = c(1))
Idents(MM) = "M_C1_S100A8"

Myeloid_cell_Neu <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/Trials/GBC_Myeloid_v3_output/Myeloid_cell_Neu.RDS")
Neu = subset(Myeloid_cell_Neu, idents = c(0,1,2,7))

CD4 <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/18_Fig7a/01 CD4.rds")
CD4 = subset(CD4, idents = c("CD4T_C2_CCR7"))

CD8 <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/18_Fig7a/02 CD8.rds")
CD8 = subset(CD8, idents = c("CD8T_C0_CCR7_GZMK"))

B <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/18_Fig7a/04 B cells.rds")
B = subset(B, idents = c("B_C0_IGHA1"))

# Epithelial <- readRDS("D:/OneDrive/NCLC/1 - Project_220101_GBC/2 - Analysis_V3_checked/1 - GBC_input data/Epithelial.RDS")
# GBC_AllSubtype <- read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
# cell_ids <- GBC_AllSubtype[GBC_AllSubtype$subtype %in% "GM16",]$cellid
# Epithelial <- subset(Epithelial, cells = cell_ids)
# Idents(Epithelial) <- "GM16"

AREG_Subtype = merge(Fb, y = c(Endo, DC, MM, Neu, CD4, CD8, B))
# saveRDS(AREG_Subtype, file = "./3 - GBC_output data_subtype/18_Fig7a/AREG_Subtype.RDS")

# AREG_Subtype <- ScaleData(AREG_Subtype)

p = DotPlot(AREG_Subtype, features = "AREG",  cols = c("lightgrey", "#f0504a"))  +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())
pdf(file = "./3 - GBC_output data_subtype/18_Fig7a/AREG_Subtype.pdf", width =5, height = 5)
print(p)
dev.off()

# # Entropy ####
# GBC_AllSubtype <- read.csv("./4 - Cooperation/WYH/All Subtype Info/celltype_info_all_filter_final.csv")
# GBC_TIMESubtype = GBC_AllSubtype[!(GBC_AllSubtype$celltype %in% c("tumor","normal","Mast")),]
# PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
# PatientInfo_DC = PatientInfo[PatientInfo$X %in% unique(GBC_TIMESubtype$orig.ident),]
# 
# GBC_TIMESubtype = GBC_TIMESubtype[GBC_TIMESubtype$orig.ident %in% PatientInfo_DC$X,]
# Freq = c()
# cycle = unique(GBC_TIMESubtype$subtype)
# for (i in 1:length(cycle)) {
#   sub = subset(GBC_TIMESubtype, subtype == cycle[i])
#   Freq = bind_rows(Freq, round((table(sub$orig.ident)/nrow(sub)),4))
# }
# Freq = as.data.frame(Freq)
# rownames(Freq) = cycle
# Freq[is.na(Freq)]=0
# 
# temp = as.data.frame(sort(apply(Freq, 1, max)))
# temp$X = rownames(temp)
# 
# subtype = read.csv("./3 - GBC_output data_subtype/ToLS_Shanno/Subtype_Entropy_230613.csv", row.names = 1)
# subtype$X = ifelse(subtype$X == "F_C12_ACTG2", "F_C10_ACTG2", 
#                    ifelse(subtype$X == "F_C13_APOD", "F_C11_APOD",
#                           ifelse(subtype$X == "F_C10_MT1X", "F_C8_MT1X",
#                                  ifelse(subtype$X == "F_C11_S100B", "F_C9_S100B",
#                                         ifelse(subtype$X == "F_C14_VCAN", "Per_C2_VCAN",
#                                                ifelse(subtype$X == "F_C15_STEAP4", "Per_C3_STEAP4",
#                                                       ifelse(subtype$X == "F_C2_RGS5", "Per_C0_RGS5",
#                                                              ifelse(subtype$X == "F_C3_COLEC11", "F_C2_COLEC11",
#                                                                     ifelse(subtype$X == "F_C4_MYH11", "Per_C1_MYH11",
#                                                                            ifelse(subtype$X == "F_C5_IGFBP2", "F_C3_IGFBP2",
#                                                                                   ifelse(subtype$X == "F_C6_CCN5", "F_C4_CCN5"）
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
                                                                                         
      