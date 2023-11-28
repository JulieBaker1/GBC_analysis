

cd4.fuji=readRDS('...\\Fig6.Tcell.rds')
cd4.fuji$celltype=cd4.fuji@active.ident
cd4.naive= subset(cd4.fuji,celltype=='CD4 TN')

expr=cd4.naive
expr <- FindVariableFeatures(expr, selection.method = "vst", nfeatures = 2000)
expr <- ScaleData(object = expr, verbose = T)
PCS = 30
k.param = 20
expr = RunPCA(expr,features = VariableFeatures(expr),verbose = FALSE)
expr=FindNeighbors(expr,dims = 1:PCS,verbose = FALSE,k.param = k.param)
expr=RunUMAP(expr,dims = 1:PCS,verbose = FALSE) 
expr=FindClusters(expr,resolution = 0.5
                  ,verbose = FALSE)
DimPlot(expr)
a=FeaturePlot(expr,features = 'AREG',pt.size = 0.3,max.cutoff =2)




