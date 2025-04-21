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

# Modified Functions ####
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


# Load Data ####
cellchat_NT <- readRDS("./cellchat_P.RDS")
cellchat_T <- readRDS("./cellchat_P_Mets.RDS")

group.cellchat_P = unique(union(levels(cellchat_NT@idents),levels(cellchat_T@idents)))
cellchat_T <- liftCellChat(cellchat_T, group.cellchat_P)
cellchat_NT <- liftCellChat(cellchat_NT, group.cellchat_P)

object.list <- list(NT = cellchat_NT, `T` = cellchat_T)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Increased pairs of common subtypes ####
p = netVisual_bubble(cellchat, 
                     sources.use = group.cellchat_P, 
                     targets.use = group.cellchat_P,  
                     comparison = c(1, 2), 
                     max.dataset = 2, title.name = "Increased signaling in T", angle.x = 45, remove.isolate = T,return.data = T)

# Ligands should be up-regulated in source cells ####
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
  
  
  w10x_toplot = as.data.frame(t(w10x_sub@active.assay$)) #@assays$RNA@data
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

# positive cells over 30% ####
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
  w10x_target = as.data.frame(t(w10x_target@assays$RNA@layers$counts)) #@assays$RNA@data
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

# Save Data ####
pppp = pppp[pppp$dataset %in% "T",]
save.image("~/Desktop/PvsPmets.RData")

pppp = pppp[pppp$source %in% subtype,]
pppp = pppp[pppp$target %in% subtype,]

load("~/Desktop/Project/GBC/2 - Analysis_V3_checked/3 - GBC_output data_subtype/CellChat/230612_Pmets.RData")

# Plot ####
# ###### Tumor ####
# ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("M_C",unique(as.character(pppp$source)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[grep("CXCL|CCL|ANXA",ppppp$ligand)],]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = pppp[pppp$target %in% unique(as.character(pppp$target))[grep("EC_C",unique(as.character(pppp$target)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[grep("VEGF",ppppp$ligand)],]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = pppp[pppp$target %in% unique(as.character(pppp$target))[grep("CD4T_C2|CD8T_C0",unique(as.character(pppp$target)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[grep("CXCL|CCL|ANXA1|NECTIN|CD86",ppppp$ligand)],]
# netVisual_bubble_replot(ppppp)

###### !! Neutrophil Plot New ####
## source
ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("N_C2",unique(as.character(pppp$source)))],]
ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[grep("CXCL|CCL|ANXA1",ppppp$ligand)],]
g=netVisual_bubble_replot(ppppp)
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N2_source_Pmet_bubble.pdf", width =7, height = 7)
print(g)
dev.off()

## replot
g=netVisual_chord_gene_zt(ppppp, 
                          annotationTrackHeight = 0.1,
                          link.target.prop = F,
                          title.name = "N2_source_Pmet")
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N2_source_Pmet_chord.pdf", width =7, height = 7)
print(g)
dev.off()

## target
ppppp = pppp[pppp$target %in% unique(as.character(pppp$target))[grep("N_C0|N_C7",unique(as.character(pppp$target)))],]
ppppp = ppppp[ppppp$pathway_name %in% ppppp$pathway_name[grep("CXCL",ppppp$pathway_name)],]
g=netVisual_bubble_replot(ppppp)
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N0_N7_target_Pmet_bubble.pdf", width =7, height = 7)
print(g)
dev.off()


## replot
g=netVisual_chord_gene_zt(ppppp, 
                       title.name = "N0_N7_target_Pmet",
                       link.target.prop = F,
                       annotationTrackHeight = 0.1
                       )
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N0_N7_target_Pmet_chord.pdf", width =10, height = 10)
print(g)
dev.off()

## Ligand expression
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

## N1 ligand
cellchat_new = subset(w10x, idents = group.cellchat_P[grep("N_C1",group.cellchat_P)])
Idents(cellchat_new) = cellchat_new$datasets

vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(cellchat_new, features = signature,pt.size = 0,
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, label = "p.signif") + 
      stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
  }
  map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, ".pdf")
  ggsave(file_name, width = 2.5, height = 3)
}


gene_sig <- c("CCL3L1")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_N1Ligand", test_sign = comparisons, y_max = 8)

## N2 ligand
cellchat_new = subset(w10x, idents = group.cellchat_P[grep("N_C2",group.cellchat_P)])
Idents(cellchat_new) = cellchat_new$datasets

vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(cellchat_new, features = signature,pt.size = 0,
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, label = "p.signif") + 
      stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
  }
  map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, ".pdf")
  ggsave(file_name, width = 2.5, height = 3)
}

gene_sig <- c("CCL3L1")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_N2Ligand_CCL3L1", test_sign = comparisons, y_max = 8)


## Receptor expression for target cells
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

Idents(w10x) = w10x$datasets

cellchat_new = subset(w10x, idents = "T")
Idents(cellchat_new) = cellchat_new$group.cellchat
cellchat_new = subset(cellchat_new, idents = group.cellchat_P[grep("M_C0|M_C2|M_C3|M_C5",group.cellchat_P)])

p = DotPlot(cellchat_new, features = "CCR1")  +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

pdf(file = "../CellChat/FinalVersion/Pmet_target_receptor.pdf", width =6, height = 1.8)
print(p)
dev.off()


###### !! DC Plot ####
ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("_cDC2_",unique(as.character(pppp$source)))],]
ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[grep("ANXA1",ppppp$ligand)],]
g=netVisual_bubble_replot(ppppp)

pdf(file = "../2 - Analysis_V3_checked/3 - GBC_output data_subtype/02_DC/PvsPmet_DC4.pdf", width =6, height = 4)
print(g)
dev.off()

## replot
g=netVisual_chord_gene_zt(ppppp, 
                          annotationTrackHeight = 0.1,
                          link.target.prop = F,
                          title.name = "DC4_Pmets")
pdf(file = "./3 - GBC_output data_subtype/02_DC/DC4_source_Pmet_chord.pdf", width =10, height = 5)
print(g)
dev.off()

## Ligand expression
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

## DC4 ligand
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

# ## Receptor expression for target cells
# object = cellchat
# meta <- object@meta
# if (is.list(object@idents)) {
#   meta$group.cellchat <- object@idents$joint
# } else {
#   meta$group.cellchat <- object@idents
# }
# if (!identical(rownames(meta), colnames(object@data.signaling))) {
#   cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
#   warning("The cell barcodes in 'meta' is different from those in the used data matrix.
#               We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
#   rownames(meta) <- colnames(object@data.signaling)
# }
# w10x <- Seurat::CreateSeuratObject(counts = cellchat@data.signaling, meta.data = meta)
# 
# Idents(w10x) = w10x$datasets
# 
# cellchat_new = subset(w10x, idents = "T")
# Idents(cellchat_new) = cellchat_new$group.cellchat
# 
# celltype_check = sort(as.character(unique(ppppp$target)))
# cellchat_new = subset(cellchat_new, idents = celltype_check[!(celltype_check %in% entropy_removed_subtypes)])
# 
# p = DotPlot(cellchat_new, features = rev(c("INSR","KLRB1","FPR1")))  +
#   coord_flip() +
#   theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_blank())
# 
# pdf(file = "../CellChat/FinalVersion/Pmet_target_receptor.pdf", width =12, height = 3)
# print(p)
# dev.off()

###### !! MM/Mast Plot ####
## source
ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("N_C2",unique(as.character(pppp$source)))],]
ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[grep("CXCL|CCL|ANXA1",ppppp$ligand)],]
g=netVisual_bubble_replot(ppppp)
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N2_source_Pmet_bubble.pdf", width =7, height = 7)
print(g)
dev.off()

## replot
g=netVisual_chord_gene_zt(ppppp, 
                          annotationTrackHeight = 0.1,
                          link.target.prop = F,
                          title.name = "N2_source_Pmet")
pdf(file = "./3 - GBC_output data_subtype/01_Neutrophil/N2_source_Pmet_chord.pdf", width =7, height = 7)
print(g)
dev.off()

## S100 Macro
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

cellchat_new = subset(w10x, idents = group.cellchat_P[grep("M_C1",group.cellchat_P)])
Idents(cellchat_new) = cellchat_new$datasets

vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(cellchat_new, features = signature,pt.size = 0,
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, label = "p.signif") + 
      stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
  }
  map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, ".pdf")
  ggsave(file_name, width = 2.5, height = 3)
}

gene_sig <- c("CXCL8")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "./3 - GBC_output data_subtype/CellChat/PvsPmet_M_S100_CXCL3", test_sign = comparisons, y_max = 7)

gene_sig <- c("CD86")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_M_S100_CD86", test_sign = comparisons, y_max = 4)

## Receptor expression for target cells
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

Idents(w10x) = w10x$datasets

cellchat_new = subset(w10x, idents = "T")
Idents(cellchat_new) = cellchat_new$group.cellchat

celltype_check = sort(as.character(unique(ppppp$target)))
cellchat_new = subset(cellchat_new, idents = celltype_check[!(celltype_check %in% entropy_removed_subtypes)])

p = DotPlot(cellchat_new, features = rev(c("CXCR2","ACKR1","CTLA4")))  +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

pdf(file = "../CellChat/FinalVersion/Pmet_M_S100_target_receptor.pdf", width =7.5, height = 3)
print(p)
dev.off()


ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("M_C2",unique(as.character(pppp$source)))],]
ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("TNFSF|RETN|APP|ADM|ANXA1",ppppp$ligand)],]
ppppp = ppppp[!(ppppp$target %in% entropy_removed_subtypes),]
g=netVisual_bubble_replot(ppppp)

pdf(file = "../CellChat/FinalVersion/PvsPmet_M_SPP1_comm.pdf", width =3, height = 5)
print(g)
dev.off()

## SPP1 Macro
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

cellchat_new = subset(w10x, idents = group.cellchat_P[grep("M_C2",group.cellchat_P)])
Idents(cellchat_new) = cellchat_new$datasets

vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(cellchat_new, features = signature,pt.size = 0,
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, label = "p.signif") + 
      stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
  }
  map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, ".pdf")
  ggsave(file_name, width = 2.5, height = 3)
}

gene_sig <- c("CXCL8")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_M_SPP1_CXCL8", test_sign = comparisons, y_max = 9)

gene_sig <- c("CXCL3")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_M_SPP1_CXCL3", test_sign = comparisons, y_max = 7)

gene_sig <- c("CXCL2")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_M_SPP1_CXCL2", test_sign = comparisons, y_max = 7)

gene_sig <- c("CCL20")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_M_SPP1_CCL20", test_sign = comparisons, y_max = 8)

gene_sig <- c("CCL2")
comparisons <- list(c("NT", "T"))
vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_M_SPP1_CCL2", test_sign = comparisons, y_max = 7)

## Receptor expression for target cells
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

Idents(w10x) = w10x$datasets

cellchat_new = subset(w10x, idents = "T")
Idents(cellchat_new) = cellchat_new$group.cellchat

celltype_check = sort(as.character(unique(ppppp$target)))
cellchat_new = subset(cellchat_new, idents = celltype_check[!(celltype_check %in% entropy_removed_subtypes)])

p = DotPlot(cellchat_new, features = rev(c("CXCR2","ACKR1","CCR6")))  +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank())

pdf(file = "../CellChat/FinalVersion/Pmet_M_SPP1_target_receptor.pdf", width =6, height = 3.5)
print(p)
dev.off()

ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("CD4T_C2_",unique(as.character(pppp$source)))],]
netVisual_bubble_replot(ppppp)


###### !! CD8T ####
## target
ppppp = pppp[pppp$target %in% unique(as.character(pppp$target))[grep("CD8T_C1_",unique(as.character(pppp$target)))],]
ppppp = ppppp[ppppp$receptor %in% ppppp$receptor[grep("CTLA4",ppppp$receptor)],]
g=netVisual_bubble_replot(ppppp)
pdf(file = "./3 - GBC_output data_subtype/CellChat/CD8T_source_Pmet_bubble.pdf", width =3, height = 3)
print(g)
dev.off()














ppppp = pppp_PLN_spesific[pppp_PLN_spesific$source %in% 
                            unique(as.character(pppp_PLN_spesific$source))[grep("T_CD8",unique(as.character(pppp_PLN_spesific$source)))],]
ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
ppppp_PLN = ppppp

ppppp = pppp_PLM_spesific[pppp_PLM_spesific$source %in% 
                            unique(as.character(pppp_PLM_spesific$source))[grep("T_CD8",unique(as.character(pppp_PLM_spesific$source)))],]
ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
ppppp_PLM = ppppp

ppppp_PLN$filter = paste0(ppppp_PLN$source,"_",ppppp_PLN$target,"_",ppppp_PLN$interaction_name)
ppppp_PLM$filter = paste0(ppppp_PLM$source,"_",ppppp_PLM$target,"_",ppppp_PLM$interaction_name)

shared = intersect(x = ppppp_PLN$filter, y = ppppp_PLM$filter)

## PLN
ppppp = ppppp_PLN[ppppp_PLN$filter %in% shared,]
netVisual_bubble_replot(ppppp)

ppppp = ppppp_PLN[!(ppppp_PLN$filter %in% shared),]
netVisual_bubble_replot(ppppp)

## PLM
ppppp = ppppp_PLM[ppppp_PLM$filter %in% shared,]
netVisual_bubble_replot(ppppp)

ppppp = ppppp_PLM[!(ppppp_PLM$filter %in% shared),]
netVisual_bubble_replot(ppppp)

# ###### CD4T ####
# ppppp = pppp_PLN_spesific[pppp_PLN_spesific$source %in% 
#                             unique(as.character(pppp_PLN_spesific$source))[grep("T_CD4",unique(as.character(pppp_PLN_spesific$source)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
# ppppp_PLN = ppppp
# 
# ppppp = pppp_PLM_spesific[pppp_PLM_spesific$source %in% 
#                             unique(as.character(pppp_PLM_spesific$source))[grep("T_CD4",unique(as.character(pppp_PLM_spesific$source)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
# ppppp_PLM = ppppp
# 
# ppppp_PLN$filter = paste0(ppppp_PLN$source,"_",ppppp_PLN$target,"_",ppppp_PLN$interaction_name)
# ppppp_PLM$filter = paste0(ppppp_PLM$source,"_",ppppp_PLM$target,"_",ppppp_PLM$interaction_name)
# 
# shared = intersect(x = ppppp_PLN$filter, y = ppppp_PLM$filter)
# 
# ## PLN
# ppppp = ppppp_PLN[ppppp_PLN$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLN[!(ppppp_PLN$filter %in% shared),]
# netVisual_bubble_replot(ppppp)
# 
# ## PLM
# ppppp = ppppp_PLM[ppppp_PLM$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLM[!(ppppp_PLM$filter %in% shared),]
# netVisual_bubble_replot(ppppp)
# 
# ###### NK ####
# ppppp = pppp_PLN_spesific[pppp_PLN_spesific$source %in% 
#                             unique(as.character(pppp_PLN_spesific$source))[grep("NK_|NKT_",unique(as.character(pppp_PLN_spesific$source)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
# ppppp_PLN = ppppp
# 
# ppppp = pppp_PLM_spesific[pppp_PLM_spesific$source %in% 
#                             unique(as.character(pppp_PLM_spesific$source))[grep("NK_|NKT_",unique(as.character(pppp_PLM_spesific$source)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
# ppppp_PLM = ppppp
# 
# ppppp_PLN$filter = paste0(ppppp_PLN$source,"_",ppppp_PLN$target,"_",ppppp_PLN$interaction_name)
# ppppp_PLM$filter = paste0(ppppp_PLM$source,"_",ppppp_PLM$target,"_",ppppp_PLM$interaction_name)
# 
# shared = intersect(x = ppppp_PLN$filter, y = ppppp_PLM$filter)
# 
# ## PLN
# ppppp = ppppp_PLN[ppppp_PLN$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLN[!(ppppp_PLN$filter %in% shared),]
# netVisual_bubble_replot(ppppp)
# 
# ## PLM
# ppppp = ppppp_PLM[ppppp_PLM$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLM[!(ppppp_PLM$filter %in% shared),]
# netVisual_bubble_replot(ppppp)
# 
# ###### B/Plasma ####
# ppppp = pppp_PLN_spesific[pppp_PLN_spesific$source %in% 
#                             unique(as.character(pppp_PLN_spesific$source))[grep("B_|Plasma_",unique(as.character(pppp_PLN_spesific$source)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
# ppppp_PLN = ppppp
# 
# ppppp = pppp_PLM_spesific[pppp_PLM_spesific$source %in% 
#                             unique(as.character(pppp_PLM_spesific$source))[grep("B_|Plasma_",unique(as.character(pppp_PLM_spesific$source)))],]
# ppppp = ppppp[ppppp$ligand %in% ppppp$ligand[-grep("HLA|COL|FN1|LAM",ppppp$ligand)],]
# ppppp_PLM = ppppp
# 
# ppppp_PLN$filter = paste0(ppppp_PLN$source,"_",ppppp_PLN$target,"_",ppppp_PLN$interaction_name)
# ppppp_PLM$filter = paste0(ppppp_PLM$source,"_",ppppp_PLM$target,"_",ppppp_PLM$interaction_name)
# 
# shared = intersect(x = ppppp_PLN$filter, y = ppppp_PLM$filter)
# 
# ## PLN
# ppppp = ppppp_PLN[ppppp_PLN$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLN[!(ppppp_PLN$filter %in% shared),]
# netVisual_bubble_replot(ppppp)
# 
# ## PLM
# ppppp = ppppp_PLM[ppppp_PLM$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLM[!(ppppp_PLM$filter %in% shared),]
# netVisual_bubble_replot(ppppp)
# 
# ###### EC ####
# ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("EC_C",unique(as.character(pppp$source)))],]
# entropy_removed_subtypes = readRDS("/Users/taozhou/Desktop/Project/GBC/01_Bioinfo_Analysis_V3/CheckedCodes/01_GBC_subtypes_zt/entropy_removed_subtypes.RDS")
# ppppp = ppppp[!(ppppp$target %in% entropy_removed_subtypes),]
# ppppp = ppppp[ppppp$ligand %in% "CXCL12",]
# ppppp = ppppp[ppppp$receptor %in% "CXCR4",]
# g=netVisual_bubble_replot(ppppp)
# 
# pdf(file = "../CellChat/FinalVersion/PvsPmet_EC_comm.pdf", width =15, height = 3)
# print(g)
# dev.off()
# 
# ppppp = pppp[pppp$target %in% unique(as.character(pppp$target))[grep("EC_C",unique(as.character(pppp$target)))],]
# netVisual_bubble_replot(ppppp)
# 
# ## Ligand expression
# object = cellchat
# meta <- object@meta
# if (is.list(object@idents)) {
#   meta$group.cellchat <- object@idents$joint
# } else {
#   meta$group.cellchat <- object@idents
# }
# if (!identical(rownames(meta), colnames(object@data.signaling))) {
#   cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
#   warning("The cell barcodes in 'meta' is different from those in the used data matrix.
#               We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
#   rownames(meta) <- colnames(object@data.signaling)
# }
# w10x <- Seurat::CreateSeuratObject(counts = cellchat@data.signaling, meta.data = meta)
# Idents(w10x) = w10x$group.cellchat
# 
# ## ligand
# cellchat_new = subset(w10x, idents = group.cellchat_P[grep("EC_C3",group.cellchat_P)])
# Idents(cellchat_new) = cellchat_new$datasets
# 
# vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
#   plot_case1 <- function(signature){
#     VlnPlot(cellchat_new, features = signature,pt.size = 0,
#             y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
#     ) + 
#       stat_compare_means(comparisons = test_sign, label = "p.signif") + 
#       stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white")
#   }
#   map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
#   file_name <- paste0(file_name, ".pdf")
#   ggsave(file_name, width = 2.5, height = 3)
# }
# 
# 
# gene_sig <- c("CXCL12")
# comparisons <- list(c("NT", "T"))
# vp_case1(gene_signature = gene_sig, file_name = "../CellChat/FinalVersion/PvsPmet_EC3Ligand", test_sign = comparisons, y_max = 6)
# 
# ## Receptor expression for target cells
# object = cellchat
# meta <- object@meta
# if (is.list(object@idents)) {
#   meta$group.cellchat <- object@idents$joint
# } else {
#   meta$group.cellchat <- object@idents
# }
# if (!identical(rownames(meta), colnames(object@data.signaling))) {
#   cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
#   warning("The cell barcodes in 'meta' is different from those in the used data matrix.
#               We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
#   rownames(meta) <- colnames(object@data.signaling)
# }
# w10x <- Seurat::CreateSeuratObject(counts = cellchat@data.signaling, meta.data = meta)
# 
# Idents(w10x) = w10x$datasets
# 
# cellchat_new = subset(w10x, idents = "T")
# Idents(cellchat_new) = cellchat_new$group.cellchat
# cellchat_new = subset(cellchat_new, idents = as.character(unique(ppppp$target)))
# 
# p = DotPlot(cellchat_new, features = "CXCR4")  +
#   coord_flip() +
#   theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),legend.box = "horizontal")+
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_blank())
# 
# pdf(file = "../CellChat/FinalVersion/Pmet_CXCR4_target_receptor.pdf", width =18, height = 3)
# print(p)
# dev.off()
# 
# ###### F ####
# ppppp = pppp[pppp$source %in% unique(as.character(pppp$source))[grep("F_|Per_",unique(as.character(pppp$source)))],]
# g=netVisual_bubble_replot(ppppp)
# 
# ppppp = pppp[pppp$target %in% unique(as.character(pppp$target))[grep("F_C1",unique(as.character(pppp$target)))],]
# g=netVisual_bubble_replot(ppppp)
# 
# pdf(file = "~/Desktop/F_PvsPmets.pdf", width =100, height = 10)
# print(g)
# dev.off()
# 
# ## PLN
# ppppp = ppppp_PLN[ppppp_PLN$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLN[!(ppppp_PLN$filter %in% shared),]
# netVisual_bubble_replot(ppppp)
# 
# ## PLM
# ppppp = ppppp_PLM[ppppp_PLM$filter %in% shared,]
# netVisual_bubble_replot(ppppp)
# 
# ppppp = ppppp_PLM[!(ppppp_PLM$filter %in% shared),]
# netVisual_bubble_replot(ppppp)




