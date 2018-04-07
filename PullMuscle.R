library(Seurat)
library(ggplot2)

Muscle1.data<-Read10X(data.dir = "Muscle-10X_P7_14/")
Muscle1<-CreateSeuratObject(raw.data = Muscle1.data, project = "Muscle1")
Muscle2.data<-Read10X(data.dir = "Muscle-10X_P7_15/")
Muscle2<-CreateSeuratObject(raw.data = Muscle2.data, project = "Muscle2")

Muscle.combined<-MergeSeurat(object1 = Muscle1,object2 = Muscle2, add.cell.id1 = "Muscle1", add.cell.id2 = "Muscle2", project = "MuscleAll")
Muscle.combined
tail(x = Muscle.combined@meta.data)


ERCC.genes <- grep(pattern = "^ERCC-", x = rownames(x = Muscle.combined@data), value = TRUE)
mito.genes <- grep(pattern = "^mt", x = rownames(x = Muscle.combined@data), value = TRUE)
percent.ERCC <- Matrix::colSums(Muscle.combined@raw.data[ERCC.genes, ])/Matrix::colSums(Muscle.combined@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
Muscle.combined <- AddMetaData(object = Muscle.combined, metadata = percent.ERCC, col.name = "percent.ERCC")
VlnPlot(object = Muscle.combined, features.plot = c("nGene", "nUMI"), nCol = 2)

# Filter cells.
#23433 genes across 4537 samples.
Muscle.combined.Trim <- FilterCells(object = Muscle.combined, subset.names = "nUMI",low.thresholds = 5000, high.thresholds = 20000)
Muscle.combined.Trim
#Now left with  23433 genes across 1888 samples.

#Normalize 
Muscle.combined.Trim <- NormalizeData(object = Muscle.combined.Trim, normalization.method = "LogNormalize", scale.factor = 10000)
Muscle.combined.Trim <- FindVariableGenes(object = Muscle.combined.Trim, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3.5, y.cutoff = 0.5)
length(x = Muscle.combined.Trim@var.genes)

#
Muscle.combined.Trim <- ScaleData(object = Muscle.combined.Trim, vars.to.regress = c("nUMI"))
#PCA
Muscle.combined.Trim <- RunPCA(object = Muscle.combined.Trim, pc.genes = Muscle.combined.Trim@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5)
VizPCA(object = Muscle.combined.Trim, pcs.use = 1:2)

#Clustering cells
Muscle.combined.Trim <- FindClusters(object = Muscle.combined.Trim, reduction.type = "pca", dims.use = 1:5, resolution = 0.2, print.output = 0, save.SNN = TRUE)

Muscle.combined.Trim <- RunTSNE(object = Muscle.combined.Trim, dims.use = 1:5, do.fast = TRUE)

TSNEPlot(object = Muscle.combined.Trim)
#
VlnPlot(object = Muscle.combined, features.plot = c("Myod1","Cd3g","Sdc4","Cdh5","Ptprc","Dll1","Dll4", "Jag1","Jag2","Pdgfra","Nos2"))
Dotplot<-DotPlot(object = Muscle.combined, genes.plot = c("Dll4","Myod1","Cd3g","Cdh5","Pdgfra","Cd19","Kera","Cd68","Ckm","Mylk"),
                 plot.legend = TRUE,
                 x.lab.rot = TRUE,
                 do.return = TRUE,
                 cols.use = c("green","red"))
#DoHeatmap(object = SubsetData(object = Muscle.combined, max.cells.per.ident = 100), genes.use = c("Myod1","Cd3g","Cdh5","Pdgfra","Cd19","Kera","Cd68", "Dll4","Ck"), slim.col.label = TRUE, group.label.rot = TRUE)

#Find markers of clusters
cluster7.markers <- FindMarkers(object = Muscle.combined, ident.1 = "Myogenic cell", min.pct = 0.25)
head(cluster7.markers)
#Coexpression"
FeaturePlot(Muscle.combined, features.plot = c("Dll4"),
            # cols.use = c("grey","red","green","blue"),
            overlay = FALSE, no.legend = FALSE)

#Name the clusters
current.cluster.ids<-c("Endothelial cells","FAP", "B&T cells", "Satellite cell","Fibroblast","Monocyte","Myogenic cell","Myonuclei") 
new.cluster.ids<-c("Endothelial cells","FAP", "B&T cells", "Satellite cell","Fibroblast","Monocyte","Smooth Muscle","Myonuclei") 
Muscle.combined@ident <- plyr::mapvalues(x = Muscle.combined@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = Muscle.combined, do.label = TRUE, pt.size = 0.5)