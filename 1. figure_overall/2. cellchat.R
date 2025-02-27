suppressWarnings(library(CellChat))
suppressWarnings(library(patchwork))
suppressWarnings(library(Seurat))
library(RhpcBLASctl)
options(future.globals.maxSize= 10485760000)
options(future.fork.multithreading.enable = FALSE)
# future::plan("multicore", workers = 10)
setwd("/public1/home/sch1103/GBC/combined_data")
for(itype in c("Normal","Tumor")){
  expr = readRDS("./RDS_file/combine_final0609.RDS")
  expr$metastasis.type[expr$metastasis.type %in% c("P","P_LI")] = "P"
  celltype_info = readRDS("./meta_data/celltype20230610/combined_celltype_include_normal20230610.RDS")
  expr = expr[,rownames(celltype_info)]
  expr$subtype = celltype_info["subtype"]
  #expr$cell_type = celltype_info["class"]
  
  expr$sample.type = "Tumor"
  expr$sample.type[expr$Tumors.for.scRNA.seq.short %in% c("CC","XGC")] = "Normal"
  expr$sample.type[expr$histological.type.short %in% c("LG","HG")] = "Normal"
  
  expr_P = subset(expr,subset = sample.type %in% c(itype) & Tumors.for.scRNA.seq.short %in% c("CC","XGC","P") & histological.type.short %in% c("CC","XGC","LG","HG","adeno"))
  rm(expr)
  expr_P = NormalizeData(expr_P, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)
  if(itype=="Tumor"){
    expr_P = expr_P[,sample.int(ncol(expr_P),round(ncol(expr_P)/1))]
  }
  meta_data = expr_P@meta.data
  meta_data$subtype = droplevels(meta_data$subtype, exclude = setdiff(levels(meta_data$subtype),unique(meta_data$subtype)))
  cellchat <- createCellChat(object = GetAssayData(expr_P,slot= "data"), meta = meta_data, group.by = "subtype")
  CellChatDB <- CellChatDB.human
  #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,population.size = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  saveRDS(cellchat,paste0("./cellchat/result1201/cellchat_",itype,".RDS"))
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  saveRDS(cellchat,paste0("./cellchat/result20230610/cellchat_",itype,".RDS"))
}