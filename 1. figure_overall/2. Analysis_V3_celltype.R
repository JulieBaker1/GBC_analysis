library(Seurat)
library(dplyr)
library(patchwork)
library(plotly)
# library(SingleR)
library(ggplot2)
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(enrichplot)
library(ggpubr)
library(tidyverse)
library(monocle) # Need R 4.0
library(GSVA)
library(fgsea)
library(msigdbr)
library(survminer)
library(survival)
library(ggrepel)
library(corrplot)
library(IOBR)
library(circlize)
library(RColorBrewer)
library(GGally)
library(data.table)
library(ComplexHeatmap)
library(dendextend)
library(stringr)
library(psych)
library(reshape)
library(jjPlot)
library(gg.gap)
library(gmodels)
library(reshape2)

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
scales::show_col(c("#51574a","#e9d78e",ggsci::pal_d3("category20")(20)[1:20]))
mycol = c("#51574a","#e9d78e",ggsci::pal_d3("category20")(20)[1:20])

# 3. Fig. 1 & Fig. S1 ####
## Fig. 1C & Fig. S1B ####
### Fig. 1C ####
celltype_color = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/palette_celltype_output.csv", row.names = 1)
celltype_percent = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/celltype_combined0823.csv", row.names = 1)
celltype_percent$patient = str_sub(celltype_percent$barcode, end = -18)
celltype_percent$class = ifelse(celltype_percent$class %in% c("CD4+ T cells","CD8+ T cells","NK cells"), "NK & T cells", celltype_percent$class)

cycle = unique(celltype_percent$patient)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent, patient == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$class)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo = PatientInfo[PatientInfo$X %in% rownames(Freq),]

annotation_row = subset(PatientInfo, select = c("NewSample.ID",
                                                "histological.type.short",
                                                "Tumors.for.scRNA.seq.short",
                                                "Sex",
                                                "Age",
                                                "Polyps",
                                                "Gallstones",
                                                "PBM",
                                                "Atrophic.cholecystitis",
                                                "Clinical.stage",
                                                "TNM",
                                                "Differentiation",
                                                "liver.invasion",
                                                "Bile.duct.invasion",
                                                "Vascular.invasion",
                                                "Lymph.node.metastasis",
                                                "Liver.metastasis",
                                                "Peritoneal..metastasis"
))
annotation_row$Tumors.for.scRNA.seq.short = ifelse(annotation_row$Tumors.for.scRNA.seq.short %in% c("CC","P","XGC"), "P", annotation_row$Tumors.for.scRNA.seq.short)
annotation_row$Clinical.stage = ifelse(annotation_row$Clinical.stage %in% c("I","IIIA","IIIB","IVA","IVB"), annotation_row$Clinical.stage,
                                       ifelse(annotation_row$Clinical.stage %in% c("IIA","IIB"), "II", NA))
annotation_row$TNM = ifelse(annotation_row$TNM %in% as.character(unique(annotation_row$TNM)[grep("T1|T2|T3|T4",unique(annotation_row$TNM))]), annotation_row$TNM, NA)
annotation_row$liver.invasion = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$liver.invasion, NA)
annotation_row$Bile.duct.invasion = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Bile.duct.invasion, NA)
annotation_row$Vascular.invasion = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Vascular.invasion, NA)
annotation_row$Lymph.node.metastasis = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Lymph.node.metastasis, NA)
annotation_row$Liver.metastasis = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Liver.metastasis, NA)
annotation_row$Peritoneal..metastasis = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Peritoneal..metastasis, NA)
annotation_row$Differentiation = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Differentiation, NA)
annotation_row$Age = ifelse(annotation_row$Age < 50, "<50", ifelse(annotation_row$Age > 65, ">65", "50~65"))

annotation_row$histological.type.short = factor(annotation_row$histological.type.short, levels = c("CC", "XGC", "LG", "HG", "adeno", "adeno squa", "squa", "neuro", "undiff"))
annotation_row$Tumors.for.scRNA.seq.short = factor(annotation_row$Tumors.for.scRNA.seq.short, levels = c("P", "PO", "LN", "LI", "LM", "OM"))
annotation_row$Clinical.stage = factor(annotation_row$Clinical.stage, levels = c("I","II","IIIA","IIIB","IVA","IVB"))
annotation_row$TNM = factor(annotation_row$TNM, levels = sort(unique(annotation_row$TNM)))
annotation_row$Sex = factor(annotation_row$Sex, levels = c("M","F"))
annotation_row$Age = factor(annotation_row$Age, levels = c("<50","50~65",">65"))

rownames(annotation_row) = annotation_row$NewSample.ID
annotation_row = annotation_row[,-1]

patient_order = rownames(Freq[order(Freq$`Epithelial cells`, decreasing = T),])
annotation_row = annotation_row[patient_order,]
annotation_row$EpiOrder = factor(paste(1:135, rownames(annotation_row), sep = "_"),levels = paste(1:135, rownames(annotation_row), sep = "_"))
patient_order = rownames(arrange(annotation_row,annotation_row[,"histological.type.short"],annotation_row[,"EpiOrder"]))

annotation_row = annotation_row[,1:17]

Freq = Freq[patient_order,]
annotation_row = annotation_row[patient_order,]

annotation_col_ha = ComplexHeatmap::HeatmapAnnotation(`Type` = annotation_row$histological.type.short,
                                                      `Site` = annotation_row$Tumors.for.scRNA.seq.short,
                                                      `Gender` = annotation_row$Sex,
                                                      `Age` = annotation_row$Age,
                                                      `Polyps` = annotation_row$Polyps,
                                                      `Gallstones` = annotation_row$Gallstones,
                                                      `PBM` = annotation_row$PBM,
                                                      `Atrophic cholecystitis` = annotation_row$Atrophic.cholecystitis,
                                                      `Clinical stage` = annotation_row$Clinical.stage,
                                                      `TNM` = annotation_row$TNM,
                                                      `Differentiation` = annotation_row$Differentiation,
                                                      `Liver invasion` = annotation_row$liver.invasion,
                                                      `Bile duct invasion` = annotation_row$Bile.duct.invasion,
                                                      `Vascular invasion` = annotation_row$Vascular.invasion,
                                                      `Lymph node metastasis` = annotation_row$Lymph.node.metastasis,
                                                      `Liver metastasis` = annotation_row$Liver.metastasis,
                                                      `Peritoneal metastasis` = annotation_row$Peritoneal..metastasis,
                                                      
                                                      na_col = "white",
                                                      annotation_name_side = "left",
                                                      
                                                      col = list(`Type` = c("adeno" = mycol[2],
                                                                            "adeno squa" = mycol[3],
                                                                            "CC" = mycol[4],
                                                                            "HG" = mycol[5],
                                                                            "LG" = mycol[6],
                                                                            "neuro" = mycol[7],
                                                                            "squa" = mycol[8],
                                                                            "undiff" = mycol[9],
                                                                            "XGC" = mycol[10]),
                                                                 `Site` = c("LI" = mycol[12],
                                                                            "LM" = mycol[13],
                                                                            "LN" = mycol[11],
                                                                            "OM" = mycol[15],
                                                                            "P" = mycol[16],
                                                                            "PO" = mycol[18]),
                                                                 `Gender` = c("F" = mycol[19],
                                                                              "M" = mycol[12]),
                                                                 `Age` = c("<50" = mycol[2], 
                                                                           "50~65" = mycol[18], 
                                                                           ">65" = mycol[8]),
                                                                 `Polyps` = c("Yes" = mycol[10],
                                                                              "No" = mycol[20]),
                                                                 `Gallstones` = c("Yes" = mycol[10],
                                                                                  "No" = mycol[20]),
                                                                 `PBM` = c("Yes" = mycol[10],
                                                                           "No" = mycol[20]),
                                                                 `Atrophic cholecystitis` = c("Yes" = mycol[10],
                                                                                              "No" = mycol[20]),
                                                                 `Clinical stage` = c("I" = mycol[15],
                                                                                      "II" = mycol[22],
                                                                                      "IIIA" = mycol[3],
                                                                                      "IIIB" = mycol[7],
                                                                                      "IVA" = mycol[9],
                                                                                      "IVB" = mycol[16]),
                                                                 `TNM` = c("T1bN0M0" = mycol[22],
                                                                           "T2aN0M0" = mycol[21],
                                                                           "T2aN1M0" = mycol[20],
                                                                           "T2bN0M0" = mycol[19],
                                                                           "T2bN1M0" = mycol[18],
                                                                           "T2bNxM0" = mycol[17],
                                                                           "T3N0M0" = mycol[16],
                                                                           "T3N0M1" = mycol[15],
                                                                           "T3N1M0" = mycol[14],
                                                                           "T3N1M1" = mycol[13],
                                                                           "T3N2M0" = mycol[12],
                                                                           "T3N2M1" = mycol[11],
                                                                           "T3NxM0" = mycol[10],
                                                                           "T3NxM1" = mycol[9],
                                                                           "T4N0M0" = mycol[8],
                                                                           "T4N1M0" = mycol[7],
                                                                           "T4N1M1" = mycol[6],
                                                                           "T4N2M0" = mycol[5],
                                                                           "T4N2M1" = mycol[4],
                                                                           "T4NxM0" = mycol[3],
                                                                           "T4NxM1" = mycol[2]),
                                                                 `Differentiation` = c("Good" = mycol[15],
                                                                                       "Moderate" = mycol[22],
                                                                                       "Moderate & poor" = mycol[7],
                                                                                       "Poor" = mycol[16]),
                                                                 `Liver invasion` = c("Yes" = mycol[10],
                                                                                      "No" = mycol[20]),
                                                                 `Bile duct invasion` = c("Yes" = mycol[10],
                                                                                          "No" = mycol[20]),
                                                                 `Vascular invasion` = c("Yes" = mycol[10],
                                                                                         "No" = mycol[20]),
                                                                 `Lymph node metastasis` = c("Yes" = mycol[10],
                                                                                             "No" = mycol[20]),
                                                                 `Liver metastasis` = c("Yes" = mycol[10],
                                                                                        "No" = mycol[20]),
                                                                 `Peritoneal metastasis` = c("Yes" = mycol[10],
                                                                                             "No" = mycol[20]))
)

p = ComplexHeatmap::Heatmap(t(as.matrix(Freq)), 
                            bottom_annotation = annotation_col_ha, 
                            cluster_columns = F,
                            cluster_rows = F,
                            show_row_names = F,
                            show_column_names = F,
                            width = ncol(t(as.matrix(Freq)))*unit(2, "mm"))
dev.off()
pdf(file = "./2 - GBC_output data_celltype/01_Output_CellularComponents/TMECellularComponents.pdf", width =20, height = 7)
print(p)
dev.off()

Freq_melt = Freq
Freq_melt$Patient = rownames(Freq_melt)
Freq_melt <- reshape2::melt(Freq_melt, id.vars=c("Patient"), variable.name="Celltype", value.name="Percentage")
Freq_melt = Freq_melt[order(Freq_melt$Patient),]
Freq_melt$Patient = factor(Freq_melt$Patient, levels = patient_order)

p = ggplot(data = Freq_melt, aes(x = Patient, y = Percentage, fill = Celltype)) + 
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
  scale_fill_manual(values=c("B cells" = "#1f77b4",
                             "Dendritic cells" = "#ff7f0e",
                             "Endothelial cells" = "#279e68",
                             "Epithelial cells" = "#d62728",
                             "Fibroblasts" = "#aa40fc",
                             "Mast cells" = "#8c564b",
                             "Monocytes & Macrophages" = "#e377c2",
                             "Neutrophils" = "#aec7e8",
                             "NK & T cells" = "#17becf",
                             "Plasma cells" = "#ffbb78",
                             "NET" = "#b5bd61"))

dev.off()
pdf(file = "./2 - GBC_output data_celltype/01_Output_CellularComponents/TMECellularComponents_barplot.pdf", width =17, height = 3)
print(p)
dev.off()

### Fig. S1B ####
celltype_color = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/palette_celltype_output.csv", row.names = 1)
celltype_percent = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/celltype_combined0823.csv", row.names = 1)
celltype_percent$patient = str_sub(celltype_percent$barcode, end = -18)
celltype_percent$class = ifelse(celltype_percent$class %in% c("CD4+ T cells","CD8+ T cells","NK cells"), "NK & T cells", celltype_percent$class)

celltype_percent = subset(celltype_percent, class != "Epithelial cells")
celltype_percent = subset(celltype_percent, class != "NET")

cycle = unique(celltype_percent$patient)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent, patient == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$class)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

Freq = Freq[patient_order,]

Freq_melt = Freq
Freq_melt$Patient = rownames(Freq_melt)
Freq_melt <- reshape2::melt(Freq_melt, id.vars=c("Patient"), variable.name="Celltype", value.name="Percentage")
Freq_melt = Freq_melt[order(Freq_melt$Patient),]
Freq_melt$Patient = factor(Freq_melt$Patient, levels = patient_order)

p = ggplot(data = Freq_melt, aes(x = Patient, y = Percentage, fill = Celltype)) + 
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
  scale_fill_manual(values=c("B cells" = "#1f77b4",
                             "Dendritic cells" = "#ff7f0e",
                             "Endothelial cells" = "#279e68",
                             "Fibroblasts" = "#aa40fc",
                             "Mast cells" = "#8c564b",
                             "Monocytes & Macrophages" = "#e377c2",
                             "Neutrophils" = "#aec7e8",
                             "NK & T cells" = "#17becf",
                             "Plasma cells" = "#ffbb78"))

dev.off()
pdf(file = "./2 - GBC_output data_celltype/01_Output_CellularComponents/TIMECellularComponents_barplot.pdf", width =17, height = 3)
print(p)
dev.off()

## Fig. S1C-G ####
### Fig. S1C ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo = PatientInfo[PatientInfo$X %in% rownames(Freq),]

annotation_row = subset(PatientInfo, select = c("NewSample.ID",
                                                "metastasis.type",
                                                "histological.type.short",
                                                "Tumors.for.scRNA.seq.short",
                                                "Sex",
                                                "Age",
                                                "Polyps",
                                                "Gallstones",
                                                "PBM",
                                                "Atrophic.cholecystitis",
                                                "Clinical.stage",
                                                "TNM",
                                                "Differentiation",
                                                "liver.invasion",
                                                "Bile.duct.invasion",
                                                "Vascular.invasion",
                                                "Lymph.node.metastasis",
                                                "Liver.metastasis",
                                                "Peritoneal..metastasis"
))

annotation_row$Tumors.for.scRNA.seq.short = ifelse(annotation_row$Tumors.for.scRNA.seq.short %in% c("CC","P","XGC"), "P", annotation_row$Tumors.for.scRNA.seq.short)
annotation_row$Clinical.stage = ifelse(annotation_row$Clinical.stage %in% c("I","IIIA","IIIB","IVA","IVB"), annotation_row$Clinical.stage,
                                       ifelse(annotation_row$Clinical.stage %in% c("IIA","IIB"), "II", NA))
annotation_row$TNM = ifelse(annotation_row$TNM %in% as.character(unique(annotation_row$TNM)[grep("T1|T2|T3|T4",unique(annotation_row$TNM))]), annotation_row$TNM, NA)
annotation_row$liver.invasion = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$liver.invasion, NA)
annotation_row$Bile.duct.invasion = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Bile.duct.invasion, NA)
annotation_row$Vascular.invasion = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Vascular.invasion, NA)
annotation_row$Lymph.node.metastasis = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Lymph.node.metastasis, NA)
annotation_row$Liver.metastasis = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Liver.metastasis, NA)
annotation_row$Peritoneal..metastasis = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Peritoneal..metastasis, NA)
annotation_row$Differentiation = ifelse(annotation_row$histological.type.short %in% c("adeno", "adeno squa", "neuro", "squa", "undiff"), annotation_row$Differentiation, NA)
annotation_row$Age = ifelse(annotation_row$Age < 50, "<50", ifelse(annotation_row$Age > 65, ">65", "50~65"))

annotation_row$histological.type.short = factor(annotation_row$histological.type.short, levels = c("CC", "XGC", "LG", "HG", "adeno", "adeno squa", "squa", "neuro", "undiff"))
annotation_row$Tumors.for.scRNA.seq.short = factor(annotation_row$Tumors.for.scRNA.seq.short, levels = c("P", "PO", "LN", "LI", "LM", "OM"))
annotation_row$Clinical.stage = factor(annotation_row$Clinical.stage, levels = c("I","II","IIIA","IIIB","IVA","IVB"))
annotation_row$TNM = factor(annotation_row$TNM, levels = sort(unique(annotation_row$TNM)))
annotation_row$Sex = factor(annotation_row$Sex, levels = c("M","F"))
annotation_row$Age = factor(annotation_row$Age, levels = c("<50","50~65",">65"))

rownames(annotation_row) = annotation_row$NewSample.ID
annotation_row = annotation_row[,-1]

pie_data = table(annotation_row$histological.type.short)
pdf("./2 - GBC_output data_celltype/02_Output_SamplePercentage/Pie_Histo.pdf", width = 6, height = 4.5)
pie(pie_data, 
    border = "white",
    labels = paste0(names(pie_data), " (", "n=", pie_data, ", ",round((pie_data/sum(pie_data))*100,2), "%)"))
dev.off()

### Fig. S1D ####
pie_data = annotation_row[annotation_row$histological.type.short %in% "adeno",]
pie_data$metastasis.type = ifelse(pie_data$metastasis.type %in% c("P","P_LI"), "P", pie_data$metastasis.type)
pie_data$metastasis.type = ifelse(pie_data$metastasis.type %in% c("P_LN","P_LM"), "P_Mets", pie_data$metastasis.type)
pie_data$Tumors.for.scRNA.seq.short = as.character(pie_data$Tumors.for.scRNA.seq.short)
pie_data$Tumors.for.scRNA.seq.short = ifelse(pie_data$Tumors.for.scRNA.seq.short %in% c("P"), pie_data$metastasis.type,pie_data$Tumors.for.scRNA.seq.short)
pie_data$Tumors.for.scRNA.seq.short = factor(pie_data$Tumors.for.scRNA.seq.short, levels = c("P","P_Mets", "PO", "LN", "LI", "LM", "OM"))
pie_data = table(pie_data$Tumors.for.scRNA.seq.short)

pdf("./2 - GBC_output data_celltype/02_Output_SamplePercentage/Pie_Adeno_Site.pdf", width = 6, height = 4.5)
pie(pie_data, 
    col = c(mycol[16],mycol[19],mycol[18],mycol[11],mycol[12],mycol[13],mycol[15]),
    border = "white",
    main = "adeno",
    labels = paste0(names(pie_data), " (", "n=", pie_data, ", ",round((pie_data/sum(pie_data))*100,2), "%)"))
dev.off()

pie_data = annotation_row[annotation_row$histological.type.short %in% "adeno squa",]
pie_data$metastasis.type = ifelse(pie_data$metastasis.type %in% c("P","P_LI"), "P", pie_data$metastasis.type)
pie_data$metastasis.type = ifelse(pie_data$metastasis.type %in% c("P_LN","P_LM"), "P_Mets", pie_data$metastasis.type)
pie_data$Tumors.for.scRNA.seq.short = as.character(pie_data$Tumors.for.scRNA.seq.short)
pie_data$Tumors.for.scRNA.seq.short = ifelse(pie_data$Tumors.for.scRNA.seq.short %in% c("P"), pie_data$metastasis.type,pie_data$Tumors.for.scRNA.seq.short)
pie_data$Tumors.for.scRNA.seq.short = factor(pie_data$Tumors.for.scRNA.seq.short, levels = c("P","P_Mets", "PO", "LN", "LI", "LM", "OM"))
pie_data = table(pie_data$Tumors.for.scRNA.seq.short)[table(pie_data$Tumors.for.scRNA.seq.short)>0]

pdf("./2 - GBC_output data_celltype/02_Output_SamplePercentage/Pie_AdenoSqua_Site.pdf", width = 6, height = 4.5)
pie(pie_data, 
    col = c(mycol[16],mycol[19],mycol[12],mycol[13]),
    border = "white",
    main = "adeno squa",
    labels = paste0(names(pie_data), " (", "n=", pie_data, ", ",round((pie_data/sum(pie_data))*100,2), "%)"))
dev.off()

pie_data = annotation_row[annotation_row$histological.type.short %in% "neuro",]
pie_data$metastasis.type = ifelse(pie_data$metastasis.type %in% c("P","P_LI"), "P", pie_data$metastasis.type)
pie_data$metastasis.type = ifelse(pie_data$metastasis.type %in% c("P_LN","P_LM"), "P_Mets", pie_data$metastasis.type)
pie_data$Tumors.for.scRNA.seq.short = as.character(pie_data$Tumors.for.scRNA.seq.short)
pie_data$Tumors.for.scRNA.seq.short = ifelse(pie_data$Tumors.for.scRNA.seq.short %in% c("P"), pie_data$metastasis.type,pie_data$Tumors.for.scRNA.seq.short)
pie_data$Tumors.for.scRNA.seq.short = factor(pie_data$Tumors.for.scRNA.seq.short, levels = c("P","P_Mets", "PO", "LN", "LI", "LM", "OM"))
pie_data = table(pie_data$Tumors.for.scRNA.seq.short)[table(pie_data$Tumors.for.scRNA.seq.short)>0]

pdf("./2 - GBC_output data_celltype/02_Output_SamplePercentage/Pie_Neuro_Site.pdf", width = 6, height = 4.5)
pie(pie_data, 
    col = c(mycol[16],mycol[18],mycol[11]),
    border = "white",
    main = "neuro",
    labels = paste0(names(pie_data), " (", "n=", pie_data, ", ",round((pie_data/sum(pie_data))*100,2), "%)"))
dev.off()

### Fig. S1E ####
pie_data = annotation_row[annotation_row$histological.type.short %in% c("adeno","adeno squa","squa","neuro","undiff"),]
pie_data = pie_data[pie_data$Tumors.for.scRNA.seq.short %in% "P",]
pie_data = table(pie_data$Age)

pdf("./2 - GBC_output data_celltype/03_Output_PatientPercentage/Pie_Patient_Age.pdf", width = 6, height = 4.5)
pie(pie_data, 
    col = c(mycol[2], mycol[18], mycol[8]),
    border = "white",
    main = "Age",
    labels = paste0(names(pie_data), " (", "n=", pie_data, ", ",round((pie_data/sum(pie_data))*100,2), "%)"))
dev.off()

### Fig. S1F ####
pie_data = annotation_row[annotation_row$histological.type.short %in% c("adeno","adeno squa","squa","neuro","undiff"),]
pie_data = pie_data[pie_data$Tumors.for.scRNA.seq.short %in% "P",]

pie_data = c(nrow(subset(pie_data, Polyps == "No" & Gallstones == "No")),
             nrow(subset(pie_data, Polyps == "Yes" & Gallstones == "No")),
             nrow(subset(pie_data, Polyps == "No" & Gallstones == "Yes")),
             nrow(subset(pie_data, Polyps == "Yes" & Gallstones == "Yes")))
names(pie_data) = c("No Polyps & Gallstones", "Polyps", "Gallstones", "Polyps & Gallstones")

pdf("./2 - GBC_output data_celltype/03_Output_PatientPercentage/Pie_Patient_PolypsGallstones.pdf", width = 8, height = 4.5)
pie(pie_data, 
    col = c(mycol[15],mycol[14],mycol[18],mycol[16]),
    border = "white",
    main = "Polyps & Gallstones",
    labels = paste0(names(pie_data), " (", "n=", pie_data, ", ",round((pie_data/sum(pie_data))*100,2), "%)"))
dev.off()

### Fig. S1G ####
pie_data = annotation_row[annotation_row$histological.type.short %in% c("adeno","adeno squa","squa","neuro","undiff"),]
pie_data = pie_data[pie_data$Tumors.for.scRNA.seq.short %in% "P",]
pie_data = table(pie_data$Sex)

pdf("./2 - GBC_output data_celltype/03_Output_PatientPercentage/Pie_Patient_Sex.pdf", width = 6, height = 4.5)
pie(pie_data, 
    col = c(mycol[12], mycol[19]),
    border = "white",
    main = "Sex",
    labels = paste0(names(pie_data), " (", "n=", pie_data, ", ",round((pie_data/sum(pie_data))*100,2), "%)"))
dev.off()

# 4. Fig. S2 ####
## Fig. S2A ####
celltype_color = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/palette_celltype_output.csv", row.names = 1)
celltype_percent = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/celltype_combined0823.csv", row.names = 1)
celltype_percent$patient = str_sub(celltype_percent$barcode, end = -18)
celltype_percent$class = ifelse(celltype_percent$class %in% c("CD4+ T cells","CD8+ T cells","NK cells"), "NK & T cells", celltype_percent$class)

celltype_percent = subset(celltype_percent, class != "Epithelial cells")
celltype_percent = subset(celltype_percent, class != "NET")

cycle = unique(celltype_percent$patient)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent, patient == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$class)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

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

Freq = Freq[c(PatientInfo_adenoP$X,
              PatientInfo_AdenoSquaP$X,
              PatientInfo_squa$X,
              PatientInfo_neuroP$X,
              PatientInfo_undiffP$X),]

PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo_surv = subset(PatientInfo_surv, NewSample.ID %in% rownames(Freq))
Freq$NewSample.ID = rownames(Freq)
patient_group = left_join(Freq,PatientInfo_surv,by="NewSample.ID")

patient_group$event01 = ifelse(patient_group$event=="dead", 1, 0)

patient_group = subset(patient_group, select = c(colnames(patient_group)[1:9],"OS_month","event01"))

colnames(patient_group)[1:9] = c("B_cells","Dendritic_cells","Endothelial_cells","Fibroblasts","Mast_cells","Monocytes_Macrophages",
                                 "Neutrophils","NK_T_cells","Plasma_cells")
res.cat = patient_group
res.cat$B_cells = ifelse(res.cat$B_cells > median(res.cat$B_cells), "High", "Low")
res.cat$Dendritic_cells = ifelse(res.cat$Dendritic_cells > median(res.cat$Dendritic_cells), "High", "Low")
res.cat$Endothelial_cells = ifelse(res.cat$Endothelial_cells > median(res.cat$Endothelial_cells), "High", "Low")
res.cat$Fibroblasts = ifelse(res.cat$Fibroblasts > median(res.cat$Fibroblasts), "High", "Low")
res.cat$Mast_cells = ifelse(res.cat$Mast_cells > median(res.cat$Mast_cells), "High", "Low")
res.cat$Monocytes_Macrophages = ifelse(res.cat$Monocytes_Macrophages > median(res.cat$Monocytes_Macrophages), "High", "Low")
res.cat$Neutrophils = ifelse(res.cat$Neutrophils > median(res.cat$Neutrophils), "High", "Low")
res.cat$NK_T_cells = ifelse(res.cat$NK_T_cells > median(res.cat$NK_T_cells), "High", "Low")
res.cat$Plasma_cells = ifelse(res.cat$Plasma_cells > median(res.cat$Plasma_cells), "High", "Low")

for (i in 1:9) {
  sub = subset(res.cat, select = c("OS_month","event01",colnames(patient_group)[i]))
  colnames(sub)[3] = "group"
  fit<-survfit(Surv(OS_month, event01)~group, data=sub)

  p = ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = sub,             # data used to fit survival curves.
    #risk.table = TRUE,       # show risk table.
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
    title = colnames(patient_group)[i]
  )

  pdf(file = paste0("./2 - GBC_output data_celltype/06_Output_CelltypePercentagePrognosis/Celltype_Prognosis_",colnames(patient_group)[i],".pdf"),
      width =4, height = 4.5, onefile = F)
  print(p)
  dev.off()
}

## Fig. S2C ####
celltype_color = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/palette_celltype_output.csv", row.names = 1)
celltype_percent = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/celltype_combined0823.csv", row.names = 1)
celltype_percent$patient = str_sub(celltype_percent$barcode, end = -18)
celltype_percent$class = ifelse(celltype_percent$class %in% c("CD4+ T cells","CD8+ T cells","NK cells"), "NK & T cells", celltype_percent$class)

celltype_percent = subset(celltype_percent, class != "Epithelial cells")
celltype_percent = subset(celltype_percent, class != "NET")

cycle = unique(celltype_percent$patient)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent, patient == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$class)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_Myeloid = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_Myeloid = subset(PatientInfo_Myeloid, select = c(X,Tumors.for.scRNA.seq.short,histological.type.short))

Freq$X = rownames(Freq)
Freq = left_join(Freq, PatientInfo_Myeloid, by = "X")
Freq = subset(Freq, histological.type.short %in% c('adeno'))

Freq = Freq[,-c(10,12)]
colnames(Freq)[10] = "Group"

Freq_melt = Freq
Freq_melt <- reshape2::melt(Freq_melt, id.vars="Group", variable.name="Celltype", value.name="Percentage")

Freq_melt$Percentage = as.numeric(Freq_melt$Percentage)
Freq_melt$Group = factor(Freq_melt$Group, levels = c("P", "PO", "LN", "LI", "LM", "OM"))

for (i in 1:length(unique(Freq_melt$Celltype))) {
  temp = Freq_melt[Freq_melt$Celltype %in% unique(Freq_melt$Celltype)[i],]
  
  p_adj <- get_adj_p(temp,
                     .col = "Percentage", .grp = "Group", p.adjust.method = "BH",
                     comparisons = list(c("P","PO"),c("P","LN"),c("P","LI"),c("P","LM"),c("P","OM"),
                                        c("PO","LN"),c("PO","LI"),c("PO","LM"),c("PO","OM"),
                                        c("LN","LI"),c("LN","LM"),c("LN","OM"),
                                        c("LI","LM"),c("LI","OM"),
                                        c("LM","OM"))
  )
  
  p = ggboxplot(temp, x="Group", y="Percentage",
                color="Group", palette = "jco", add="jitter") +
    stat_pvalue_manual(p_adj, label = "p.adj", hide.ns = T) + 
    NoLegend() +
    rotate_x_text(angle = 45) +
    theme(axis.title.x = element_blank()) +
    ggtitle(paste0("Celltype_",unique(temp$Celltype)))
  
  pdf(file = paste0("./2 - GBC_output data_celltype/07_Output_CellPercentageAdenoPSites/Celltype_BH_",unique(temp$Celltype),".pdf"), width =4, height = 4.5)
  print(p)
  dev.off()
}

## Fig. S2D ####
celltype_color = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/palette_celltype_output.csv", row.names = 1)
celltype_percent = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/celltype_combined0823.csv", row.names = 1)
celltype_percent$patient = str_sub(celltype_percent$barcode, end = -18)
celltype_percent$class = ifelse(celltype_percent$class %in% c("CD4+ T cells","CD8+ T cells","NK cells"), "NK & T cells", celltype_percent$class)

celltype_percent = subset(celltype_percent, class != "Epithelial cells")
celltype_percent = subset(celltype_percent, class != "NET")

cycle = unique(celltype_percent$patient)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent, patient == cycle[i])
  Freq = bind_rows(Freq, round(table(sub$class)/nrow(sub),4))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_Myeloid = PatientInfo[PatientInfo$X %in% rownames(Freq),]
PatientInfo_Myeloid = subset(PatientInfo_Myeloid, select = c(X,Tumors.for.scRNA.seq.short,histological.type.short,Sex))

Freq$X = rownames(Freq)
Freq = left_join(Freq, PatientInfo_Myeloid, by = "X")
Freq = subset(Freq, histological.type.short %in% c('adeno'))
Freq = subset(Freq, Tumors.for.scRNA.seq.short %in% c('P'))

Freq = Freq[,-c(10:12)]
colnames(Freq)[10] = "Group"

Freq_melt = Freq
Freq_melt <- reshape2::melt(Freq_melt, id.vars="Group", variable.name="Celltype", value.name="Percentage")

Freq_melt$Percentage = as.numeric(Freq_melt$Percentage)
Freq_melt$Group = factor(Freq_melt$Group, levels = c("M", "F"))

for (i in 1:length(unique(Freq_melt$Celltype))) {
  temp = Freq_melt[Freq_melt$Celltype %in% unique(Freq_melt$Celltype)[i],]
  
  p_adj <- get_adj_p(temp,
                     .col = "Percentage", .grp = "Group", p.adjust.method = "BH",
                     comparisons = list(c("M","F"))
  )
  
  p = ggboxplot(temp, x="Group", y="Percentage",
                color="Group", palette = "jco", add="jitter") +
    stat_pvalue_manual(p_adj, label = "p.adj", hide.ns = T) + 
    NoLegend() +
    rotate_x_text(angle = 45) +
    theme(axis.title.x = element_blank()) +
    ggtitle(paste0("Celltype_",unique(temp$Celltype)))
  
  pdf(file = paste0("./2 - GBC_output data_celltype/08_Output_CellPercentageAdenoPSex/Celltype_BH_",unique(temp$Celltype),".pdf"), width =4, height = 4.5)
  print(p)
  dev.off()
}

# 5. Roe input ####
## Fig. 1D ####
celltype_color = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/palette_celltype_output.csv", row.names = 1)
celltype_percent = read.csv("./2 - GBC_output data_celltype/Input_FromWYH/celltype_combined0823.csv", row.names = 1)
celltype_percent$patient = str_sub(celltype_percent$barcode, end = -18)
celltype_percent$class = ifelse(celltype_percent$class %in% c("CD4+ T cells","CD8+ T cells","NK cells"), "NK & T cells", celltype_percent$class)

cycle = unique(celltype_percent$patient)
Freq = c()
for (i in 1:length(cycle)) {
  sub = subset(celltype_percent, patient == cycle[i])
  Freq = bind_rows(Freq, table(sub$class))
}
Freq = as.data.frame(Freq)
rownames(Freq) = cycle
Freq[is.na(Freq)]=0

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

Freq_melt = Freq[c(PatientInfo_adenoP$X, PatientInfo_AdenoSquaP$X, PatientInfo_neuroP$X,
              PatientInfo_undiffP$X, PatientInfo_squa$X, PatientInfo_CC$X,
              PatientInfo_HG$X, PatientInfo_LG$X, PatientInfo_XGC$X),]
Freq_melt$X = rownames(Freq_melt)
Freq_melt = left_join(Freq_melt,PatientInfo_bind,by="X")

Freq_CC = Freq_melt[Freq_melt$Group %in% c("CC"),]
Freq_CC = subset(Freq_CC, select = colnames(Freq_melt)[1:11])

Freq_XGC = Freq_melt[Freq_melt$Group %in% c("XGC"),]
Freq_XGC = subset(Freq_XGC, select = colnames(Freq_melt)[1:11])

Freq_LG = Freq_melt[Freq_melt$Group %in% c("LG"),]
Freq_LG = subset(Freq_LG, select = colnames(Freq_melt)[1:11])

Freq_HG = Freq_melt[Freq_melt$Group %in% c("HG"),]
Freq_HG = subset(Freq_HG, select = colnames(Freq_melt)[1:11])

Freq_adeno = Freq_melt[Freq_melt$Group %in% c("Adeno"),]
Freq_adeno = subset(Freq_adeno, select = colnames(Freq_melt)[1:11])

Freq_adeno_squa = Freq_melt[Freq_melt$Group %in% c("Adeno Squa"),]
Freq_adeno_squa = subset(Freq_adeno_squa, select = colnames(Freq_melt)[1:11])

Freq_neuro = Freq_melt[Freq_melt$Group %in% c("Neuro"),]
Freq_neuro = subset(Freq_neuro, select = colnames(Freq_melt)[1:11])

Freq_squa = Freq_melt[Freq_melt$Group %in% c("Squa"),]
Freq_squa = subset(Freq_squa, select = colnames(Freq_melt)[1:11])

Freq_undiff = Freq_melt[Freq_melt$Group %in% c("Undiff"),]
Freq_undiff = subset(Freq_undiff, select = colnames(Freq_melt)[1:11])

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
Celltypes_Counts_Roe_p = Freq_NTvsT_detail
saveRDS(Celltypes_Counts_Roe_p, file = './2 - GBC_output data_celltype/ToLS_Roe/Celltypes_Counts_Roe_p_230117.RDS')

## Fig. S2B ####
PatientInfo = readRDS("./1 - GBC_input data/PatientInfo.Rds")
PatientInfo_Myeloid = PatientInfo[PatientInfo$X %in% rownames(Freq),]

PatientInfo_Myeloid = subset(PatientInfo_Myeloid, select = c(X,Tumors.for.scRNA.seq.short,histological.type.short))

Freq_melt = Freq
Freq_melt$X = rownames(Freq_melt)
Freq_melt = left_join(Freq_melt, PatientInfo_Myeloid, by = "X")
Freq_melt = subset(Freq_melt, histological.type.short %in% c('adeno'))

Freq_melt = Freq_melt[,-c(12,14)]
colnames(Freq_melt)[12] = "Group"

Freq_P = Freq_melt[Freq_melt$Group %in% c("P"),]
Freq_P = subset(Freq_P, select = colnames(Freq_melt)[1:11])

Freq_PO = Freq_melt[Freq_melt$Group %in% c("PO"),]
Freq_PO = subset(Freq_PO, select = colnames(Freq_melt)[1:11])

Freq_LN = Freq_melt[Freq_melt$Group %in% c("LN"),]
Freq_LN = subset(Freq_LN, select = colnames(Freq_melt)[1:11])

Freq_LI = Freq_melt[Freq_melt$Group %in% c("LI"),]
Freq_LI = subset(Freq_LI, select = colnames(Freq_melt)[1:11])

Freq_LM = Freq_melt[Freq_melt$Group %in% c("LM"),]
Freq_LM = subset(Freq_LM, select = colnames(Freq_melt)[1:11])

Freq_OM = Freq_melt[Freq_melt$Group %in% c("OM"),]
Freq_OM = subset(Freq_OM, select = colnames(Freq_melt)[1:11])

Freq_NTvsT_detail = bind_rows(colSums(Freq_P),
                              colSums(Freq_PO),
                              colSums(Freq_LN),
                              colSums(Freq_LI),
                              colSums(Freq_LM),
                              colSums(Freq_OM))
rownames(Freq_NTvsT_detail) = c("P", "PO", "LN", "LI", "LM", "OM")
Celltypes_Counts_Roe_site = Freq_NTvsT_detail
saveRDS(Celltypes_Counts_Roe_site, file='./2 - GBC_output data_celltype/ToLS_Roe/Celltypes_Counts_Roe_site_230117.RDS')

# 6. Myeloid cell type definition ####
Myeloid_RawData = readRDS("./1 - GBC_input data/Myeloid_cell.RDS")
MetaData0609 = read.csv("./1 - GBC_input data/meta_data0609.csv", row.names = 1) ## Finally based on celltype scanpy cells

ig_genes = c(rownames(Myeloid_RawData)[grep("^IGJ",rownames(Myeloid_RawData))],
             rownames(Myeloid_RawData)[grep("^IGH",rownames(Myeloid_RawData))],
             rownames(Myeloid_RawData)[grep("^IGK",rownames(Myeloid_RawData))],
             rownames(Myeloid_RawData)[grep("^IGL",rownames(Myeloid_RawData))])
Myeloid_filtered = subset(Myeloid_RawData, features = rownames(Myeloid_RawData)[!(rownames(Myeloid_RawData) %in% ig_genes)])
VlnPlot(Myeloid_filtered, features = c("nFeature_RNA", "nCount_RNA", "mito.percent"), ncol = 3, raster=FALSE) # mito.percent has been subset
saveRDS(Myeloid_filtered, file = "./GBC_Myeloid_v3_output/Myeloid_filtered.RDS")

Myeloid_filtered <- NormalizeData(Myeloid_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_filtered <- FindVariableFeatures(Myeloid_filtered, selection.method = "vst", nfeatures = 2000)
Myeloid_filtered <- ScaleData(Myeloid_filtered)
Myeloid_filtered <- RunPCA(Myeloid_filtered, features = VariableFeatures(object = Myeloid_filtered))
Myeloid_filtered <- JackStraw(Myeloid_filtered, num.replicate = 100)
Myeloid_filtered <- ScoreJackStraw(Myeloid_filtered, dims = 1:20)
ElbowPlot(Myeloid_filtered)
Myeloid_filtered <- FindNeighbors(Myeloid_filtered, dims = 1:10)
Myeloid_filtered <- FindClusters(Myeloid_filtered, resolution = 0.5)
Myeloid_filtered <- RunUMAP(Myeloid_filtered, dims = 1:10)

Myeloid_filtered_processed = Myeloid_filtered
saveRDS(Myeloid_filtered_processed, file = "./GBC_Myeloid_v3_output/Myeloid_filtered_processed.RDS")

DimPlot(Myeloid_filtered_processed, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(Myeloid_filtered_processed, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "RNA_snn_res.0.5",ncol = 6) + NoLegend()
FeaturePlot(Myeloid_filtered_processed, features = "doublet.score")
table(Idents(Myeloid_filtered_processed))

FeaturePlot(Myeloid_filtered_processed, features = c("TPSAB1","KIT", # Mast
                                                     "CD68","CD163", "CD14","FCGR3A", # Mono/Macro
                                                     "CD1C", # cDC
                                                     "CSF3R","FCGR3B", # Neu
                                                     "CD3D"), ncol = 5)
new.cluster.ids = c(
  "Neu", # 0
  "Mono_Macro", # 1
  "Neu", # 2
  "Mast", # 3
  "Mono_Macro", # 4
  "Mono_Macro", # 5
  "Mono_Macro", # 6
  "Mono_Macro", # 7
  "DC", # 8
  "Doublet", # 9
  "Mono_Macro", # 10
  "Neu", # 11
  "Mono_Macro", # 12
  "Neu", # 13
  "Doublet", # 14
  "Mono_Macro", # 15
  "Doublet", # 16
  "Mono_Macro" # 17
)
names(new.cluster.ids) <- levels(Myeloid_filtered_processed)
Myeloid_filtered_processed <- RenameIdents(Myeloid_filtered_processed, new.cluster.ids)
Myeloid_filtered_processed@meta.data$Myeloid_celltype = Idents(Myeloid_filtered_processed)
DimPlot(Myeloid_filtered_processed, reduction = "umap", label = TRUE, repel = T) + NoLegend()
table(Idents(Myeloid_filtered_processed))

Myeloid_filtered_processed_F <- subset(Myeloid_filtered_processed, idents = "Doublet", invert = TRUE)
Myeloid_filtered_processed_F <- NormalizeData(Myeloid_filtered_processed_F, normalization.method = "LogNormalize", scale.factor = 10000)
Myeloid_filtered_processed_F <- FindVariableFeatures(Myeloid_filtered_processed_F, selection.method = "vst", nfeatures = 2000)
Myeloid_filtered_processed_F <- ScaleData(Myeloid_filtered_processed_F)
Myeloid_filtered_processed_F <- RunPCA(Myeloid_filtered_processed_F, features = VariableFeatures(object = Myeloid_filtered_processed_F))
Myeloid_filtered_processed_F <- FindNeighbors(Myeloid_filtered_processed_F, dims = 1:10)
Myeloid_filtered_processed_F <- FindClusters(Myeloid_filtered_processed_F, resolution = 0.5)
Myeloid_filtered_processed_F <- RunUMAP(Myeloid_filtered_processed_F, dims = 1:10)
DimPlot(Myeloid_filtered_processed_F, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
Idents(Myeloid_filtered_processed_F) = Myeloid_filtered_processed_F@meta.data$Myeloid_celltype
DimPlot(Myeloid_filtered_processed_F, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(Myeloid_filtered_processed_F, file = "./GBC_Myeloid_v3_output/Myeloid_filtered_processed_F.RDS")

cycle = unique(Myeloid_filtered_processed_F$orig.ident)
myeloid_num = c()
for (i in 1:length(cycle)) {
  sub = subset(Myeloid_filtered_processed_F, orig.ident == cycle[i])
  myeloid_num = c(myeloid_num, ncol(sub))
}
names(myeloid_num) = cycle
myeloid_num_rmdoublets = myeloid_num
saveRDS(myeloid_num_rmdoublets, file = "./GBC_Myeloid_v3_output/myeloid_num_rmdoublets.RDS")
