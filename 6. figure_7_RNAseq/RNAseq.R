rm(list = ls())
set.seed(0309)
options(stringsAsFactors = F)
.libPaths("E:/3-Apps/RPackages/")

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(ggrepel)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(stringr)

setwd("E:\\TT\\20231122_RNAseq\\")

cts <- read.csv("cts.csv", row.names = 1, header = T) %>% as.matrix()
cts <- cts[rowMeans(cts)>1,]

condition <- factor(c("AREG_0","AREG_0","AREG_0","AREG_100","AREG_100","AREG_100"))
colData <- data.frame(row.names=colnames(cts), condition)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = colData, design = ~ condition)
dds1 <- DESeq(dds)

resultsNames(dds1)
res <- results(dds1)
summary(res)

res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),] 
res1_total <- rbind(res1_up,res1_down)

genes<- res1
genes <- genes[which(genes$pvalue<0.05 & abs(genes$log2FoldChange)>=1),]

genes$color <- ifelse(genes$pvalue<0.05 & abs(genes$log2FoldChange)>= 1,ifelse(genes$log2FoldChange > 1,'red','blue'),'grey')
color <- c(red = "#DC143C", grey = "grey", blue = "#00008B")

p <- ggplot(
  genes, aes(log2FoldChange, -log10(pvalue), col = color)) +  
  geom_point(size=2.5) +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="Log2 (FoldChange)",y="-Log10 (P)") +
  geom_hline(yintercept = -log10(0.05), lty="dashed", col="black", lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty="dashed", col="black", lwd=0.6) +
  scale_x_continuous(limits=c(-3,3), breaks=seq(-3,3,1.5)) +
  theme(axis.title.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black")) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p

geneList0 <- c('ENSG00000100985',
               'ENSG00000115009',
               'ENSG00000163735',
               'ENSG00000115008',
               'ENSG00000125538',
               'ENSG00000149968',
               'ENSG00000172349',
               'ENSG00000137033'
)

geneList <- genes[geneList0,]

geneList1 <- genes[rownames(genes) %in% geneList0,]
geneList1 <- subset(geneList1, select = -color)
geneList1$label <- c("CXCL5","IL1A","IL1B","CCL20","MMP3","MMP9","IL33","IL16")

p + geom_text_repel(data = geneList1, 
                    aes(x = log2FoldChange, y = -log10(pvalue), label = label),
                    size = 4, color="#DC143C",
                    box.padding = unit(0.75, "lines"),
                    segment.color = "black",
                    segment.size = 0.4
) +
  theme_gray() +
  theme_bw()

pdf("VolcanoPlot_v20231130.pdf", width = 6, height = 5)
p + geom_text_repel(data = geneList1, 
                    aes(x = log2FoldChange, y = -log10(pvalue), label = label),
                    size = 4, color="#DC143C",
                    box.padding = unit(0.75, "lines"),
                    segment.color = "black",
                    segment.size = 0.4
) +
  theme_gray() +
  theme_bw()
dev.off()

go_res <- mutate(go_res, ratio = Count / as.numeric(sub(".*/", "", GeneRatio)))
df <- go_res[go_res$Description %in% c("leukocyte_migration", "regulation_of_inflammatory_response", "cytokine_mediated_signaling_pathway", "ERK1_and_ERK2_cascade", "positive_regulation_of_NF_kappaB_transcription_factor_activity", "type_II_interferon_production", "vascular_endothelial_growth_factor_production"),]
df$Description <- factor(df$Description, levels = df$Description)

p <- ggplot(df)+
  geom_point(aes(ratio, Description,
                 color = p.adjust,
                 size = Count))+
  labs(x="Gene Ratio") + 
  labs(title="")+
  theme_bw()

pdf("GO_plot.pdf", width = 6, height = 5)
print(p)
dev.off()

