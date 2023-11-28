# 17. Fig7 ####
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
TMB_TotalMutation = read.csv(file = "./4 - Cooperation/WYH/WES/variants_number_patient.csv", row.names = 1)
colnames(TMB_TotalMutation)[1] = "NewSample.ID"
TMB_TotalMutation = TMB_TotalMutation[,c(1:3)]
annotation = left_join(annotation, TMB_TotalMutation, by = "NewSample.ID")
Pop_Mutation = read.csv(file = "./4 - Cooperation/WYH/WES/onco_matrix_105_patient.csv", row.names = 1)
selected_Mutation = names(sort(apply(Pop_Mutation,1,function(x) length(grep("_",x))),decreasing = T)[1:20])
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
                      #total_perMB = colorRampPalette(c("white", "red"))(100),
                      
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
bk <- c(seq(-0.5,0,by=0.01),seq(0,1,by=0.01))
p = pheatmap(Freq_adj,
             scale = "none",
             #cluster_rows = F,
             cutree_cols = 5,
             annotation_col = annotation,
             annotation_colors = annotation_col,
             clustering_distance_rows="manhattan",
             clustering_method = "ward.D2",
             clustering_distance_cols="manhattan",
             color = c(colorRampPalette(colors = c("#34ebc3","white","#eb346b"))(100)),
             legend_breaks = seq(-0.5,1,0.5),
             breaks = bk)
# pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/PatientClustering.pdf", width =14, height = 7)
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/PatientClustering_anno.pdf", width =18, height = 15)
print(p)
dev.off()

library("dendextend")
hclust = hclust(dist(t(Freq_adj), method = "manhattan"),method = "ward.D2")
plot(color_branches(hclust,h=27.3,groupLabels = T))
patient_group = cutree(hclust, h = 27.3, order_clusters_as_data = F)
patient_group = data.frame(NewSample.ID = names(patient_group),
                           group = patient_group)
saveRDS(patient_group, file = "./3 - GBC_output data_subtype/11_PatientClustering/patient_group_5MI.RDS")
patient_group = patient_group[rownames(patient_group) %in% rownames(Pop_Mutation_MutHeatmap),]
patient_group$group = paste0("MI",patient_group$group)
patient_group = patient_group[rownames(Pop_Mutation_MutHeatmap),]
patient_group = patient_group[order(patient_group$group),]
ha = HeatmapAnnotation(group = patient_group$group,
                       col = list(group = c("MI1" = "#915F19",
                                            "MI2" = "#FF7F0E",
                                            "MI3" = "#D754DE",
                                            "MI4" = "#09A5ED",
                                            "MI5" = "#34EBC3")))

color_hp = c("#AA40FC","#D62728","#E377C2","#17BECF","#279E68","#AEC7E8")
names(color_hp) = c("Nonsense_Mutation","Missense_Mutation","In_Frame_Del",
                    "Multi_Hit","Frame_Shift_Ins","Frame_Shift_Del")
Pop_Mutation_MutHeatmap = t(Pop_Mutation_MutHeatmap)
Pop_Mutation_MutHeatmap_t = Pop_Mutation_MutHeatmap[,patient_group$NewSample.ID]
identical(colnames(Pop_Mutation_MutHeatmap_t), patient_group$NewSample.ID)
p = Heatmap(Pop_Mutation_MutHeatmap_t,
            col = color_hp,
            top_annotation = ha,
            na_col = "white")
pdf(file = "./3 - GBC_output data_subtype/12_Mutation_Pct_Cor/Mut_Heatmap.pdf", width =8, height = 5)
print(p)
dev.off()

## * calculate chi square ####
annotation_test = annotation
annotation_test$NewSample.ID = rownames(annotation_test)

temp = subset(patient_group, select = c("NewSample.ID", "group1"))
annotation_test = left_join(annotation_test, temp, by = "NewSample.ID")
annotation_test$group1 = paste0("MI",annotation_test$group1)

CrossTable(annotation_test$group1, annotation_test$Clinical.stage,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Clinical.stage, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$Group,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Sex,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Age,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Polyps,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Gallstones,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$PBM,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Atrophic.cholecystitis,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Differentiation,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
set.seed(202308)
fisher.test(annotation_test$group1, annotation_test$Differentiation, simulate.p.value=TRUE)

CrossTable(annotation_test$group1, annotation_test$liver.invasion,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Bile.duct.invasion,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Vascular.invasion,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Lymph.node.metastasis,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Liver.metastasis,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$Peritoneal..metastasis,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)

####
CrossTable(annotation_test$group1, annotation_test$TP53,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$MUC5B,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$ARID1A,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$FLG,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$MUC16,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)

CrossTable(annotation_test$group1, annotation_test$TTN,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$CDKN2A,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$MAGI1,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T) #
CrossTable(annotation_test$group1, annotation_test$PLEC,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$MUC12,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)

CrossTable(annotation_test$group1, annotation_test$TLE4,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$DISP3,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$CACNA1C,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$FAT2,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T) #
CrossTable(annotation_test$group1, annotation_test$PRRC2B,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)

CrossTable(annotation_test$group1, annotation_test$SYNE1,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$ACOT12,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T) #
CrossTable(annotation_test$group1, annotation_test$AHNAK2,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$E2F8,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)
CrossTable(annotation_test$group1, annotation_test$FLG2,fisher=T,chisq =T,format="SPSS",expected=T,prop.c = T,prop.t = T,prop.chisq = T)

####
compare_means(total~group1, annotation_test, method = "anova")
compare_means(total_perMB~group1, annotation_test, method = "anova")

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
PatientInfo_surv = openxlsx::read.xlsx("./1 - GBC_input data/GBC样本整理信息20220411_预后.xlsx")
PatientInfo_surv = subset(PatientInfo_surv, NewSample.ID %in% patient_group$NewSample.ID)
patient_group = left_join(patient_group,PatientInfo_surv,by="NewSample.ID")
patient_group$event[is.na(patient_group$event)] = 0
patient_group$event01 = ifelse(patient_group$event=="dead", 1, 0)
patient_group = left_join(patient_group,PatientInfo_DC,by="NewSample.ID")
PatientInfo_metastasis = patient_group[patient_group$histological.type.short %in% c("adeno"),]
patient_group = PatientInfo_metastasis[PatientInfo_metastasis$Tumors.for.scRNA.seq.short %in% c("P"),]
patient_group$group1 = patient_group$group
patient_group$group = ifelse(patient_group$group1 ==3,3,"others")
fit<-survfit(Surv(OS_month, event01)~group, data=patient_group)
g=ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = patient_group,             # data used to fit survival curves.
  palette = c("#34ebc3","#bbbcbd"),
  #risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  #conf.int = TRUE,         # show confidence intervals for
  # palette = "npg",
  xlab = "Time in days",   # customize X axis label.
  ggtheme = theme_classic(),
  #conf.int.style = "step",  # customize style of confidence intervals  "ribbon" 'step'
  #surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = c("TP1", "TP2","TP3","TP4","TP5")    # change legend labels.
)
pdf(file = "./3 - GBC_output data_subtype/11_PatientClustering/Prognosis_5_16H59L.pdf", width =3, height = 3, onefile = F)
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