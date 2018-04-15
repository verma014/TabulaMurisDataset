

library(dplyr)
library(data.table)
library(Seurat)
library(ggplot2)
library(cowplot)
#Load the annotation
setwd("S:/Mayank/singlecellseq/Tabulamuris/SmartSeq")
#Annotation<-read.csv("annotations_FACS.csv")
#metadata<-read.csv("metadata_FACS.csv")
#list all the files in the folder
setwd("S:/Mayank/singlecellseq/Tabulamuris/SmartSeq/FACS")
xlist<-list.files(pattern = "*.csv")
#Open all the CSVs with data.table package and fredd
xlist<-xlist[5:7]

for(i in xlist)  { 
x <- data.frame(fread((i)),row.names = 1)
x <- CreateSeuratObject(raw.data = x,project = gsub("-counts.csv","",i), normalization.method = "LogNormalize",do.scale = TRUE, do.center = TRUE,scale.factor = 100000,min.cells = 0,min.genes = 500)
x@meta.data$tissue<- gsub("-counts.csv","",i)
x <- FindVariableGenes(x,do.plot = FALSE)
#x<-RunPCA(x,do.print = FALSE)
#x<-FindClusters(x, reduction.type = "pca", dims.use = 1:15,resolution = 0.9,print.output = FALSE)
#x<-RunTSNE(x, dims.use = 1:15)
#plotgenes=c("Myf5","Foxg1","Meis2")
#plotgenes=c("Myf5","Cdh5","Ptprc","Cd3g","Cd68","Pdgfra","Foxg1","Meis2","Col19a1","Col8a1")
#temp1<-TSNEPlot(x,do.label = TRUE, do.return =TRUE)
#temp2<-DotPlot(x,genes.plot = rev(plotgenes),plot.legend = TRUE,x.lab.rot = TRUE, do.return = TRUE)
#temp2<-temp2+labs(title=i)
#print(plot_grid(temp1,temp2))
assign(i,x)

}
#

plotgenes<-c("Cdh5","Myf5","Foxg1","Meis2","Pdgfra")
DotPlot(`Muscle-counts.csv`,genes.plot = rev(plotgenes),plot.legend = TRUE,x.lab.rot = TRUE, do.return = TRUE)

 
tabulamuris <- MergeSeurat(object1 = `Bladder-counts.csv`,
                           object2 = `Kidney-counts.csv`,
                           add.cell.id1 = "1",
                           add.cell.id2 = "2",
                           project = "TabulaMuris")
`Kidney-counts.csv`<-rm()
`Bladder-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Liver-counts.csv`,
                          add.cell.id1 = "3",
                          add.cell.id2 = "4",
                          project = "TabulaMuris")
`Liver-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Lung-counts.csv`,
                          add.cell.id1 = "5",
                          add.cell.id2 = "6",
                          project = "TabulaMuris")
`Lung-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Muscle-counts.csv`,
                          add.cell.id1 = "tabulamuris",
                          add.cell.id2 = "7",
                          project = "8")
`Muscle-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Spleen-counts.csv`,
                          add.cell.id1 = "9",
                          add.cell.id2 = "10",
                          project = "TabulaMuris")
`Spleen-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Pancreas-counts.csv`,
                          add.cell.id1 = "11",
                          add.cell.id2 = "12",
                          project = "TabulaMuris")
`Pancreas-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Skin-counts.csv`,
                          add.cell.id1 = "13",
                          add.cell.id2 = "14",
                          project = "TabulaMuris")
`Skin-counts.csv`<-rm()

tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Mammary-counts.csv`,
                          add.cell.id1 = "15",
                          add.cell.id2 = "16",
                          project = "TabulaMuris")
`Mammary-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Brain_Microglia-counts.csv`,
                          add.cell.id1 = "17",
                          add.cell.id2 = "18",
                          project = "TabulaMuris")
`Brain_Microglia-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Brain_Neurons-counts.csv`,
                          add.cell.id1 = "19",
                          add.cell.id2 = "20",
                          project = "TabulaMuris")
`Brain_Neurons-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Colon-counts.csv`,
                          add.cell.id1 = "21",
                          add.cell.id2 = "22",
                          project = "TabulaMuris")
`Colon-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Fat-counts.csv`,
                          add.cell.id1 = "23",
                          add.cell.id2 = "24",
                          project = "TabulaMuris")
`Fat-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Marrow-counts.csv`,
                          add.cell.id1 = "25",
                          add.cell.id2 = "26",
                          project = "TabulaMuris")
`Marrow-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Thymus-counts.csv`,
                          add.cell.id1 = "27",
                          add.cell.id2 = "28",
                          project = "TabulaMuris")
`Thymus-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Tongue-counts.csv`,
                          add.cell.id1 = "29",
                          add.cell.id2 = "30",
                          project = "TabulaMuris")
`Tongue-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Trachea-counts.csv`,
                          add.cell.id1 = "31",
                          add.cell.id2 = "32",
                          project = "TabulaMuris")
`Trachea-counts.csv`<-rm()
tabulamuris<- MergeSeurat(object1 = tabulamuris,
                          object2 = `Heart-counts.csv`,
                          add.cell.id1 = "33",
                          add.cell.id2 = "34",
                          project = "TabulaMuris")
`Heart-counts.csv`<-rm()

`Kidney-counts.csv`-rm()
`Skin-counts.csv`<-rm()
save()


table(tabulamuris@meta.data$tissue)
tabulamuris

tabulamuris<-FilterCells(tabulamuris,subset.names = "nUMI",low.thresholds = 5.0e+04, high.thresholds = 4.0e+06)
VlnPlot(tabulamuris,features.plot = "nUMI",group.by = "tissue")
VlnPlot(tabulamuris,features.plot = "Foxg1",group.by = "tissue")
hist(log10(tabulamuris@meta.data$nUMI))

#Rescale and do TSNE
tabulamuris<-NormalizeData(tabulamuris,display.progress = TRUE)
tabulamuris<-FindVariableGenes(tabulamuris,do.plot = FALSE)
ERCC.genes <- grep(pattern = "^ERCC-", x = rownames(x = tabulamuris@data), value = TRUE)
percent.ERCC <- Matrix::colSums(tabulamuris@raw.data[ERCC.genes, ])/Matrix::colSums(tabulamuris@raw.data)
tabulamuris <- AddMetaData(object = tabulamuris, metadata = percent.ERCC, col.name = "percent.ERCC")
VlnPlot(tabulamuris,features.plot = "percent.ERCC",group.by = "tissue")
hv.genes <- head(rownames(tabulamuris@hvg.info), 1000)
tabulamuris<-ScaleData(tabulamuris,genes.use = hv.genes, num.cores = 4,do.par = TRUE)
tabulamuris<-RunPCA(tabulamuris,pc.genes = hv.genes,pcs.compute = 100,do.print = FALSE)
PCElbowPlot(tabulamuris,num.pc = 100)
tabulamuris <- FindClusters(object = tabulamuris, reduction.type = "pca", dims.use = 1:75, resolution = 3, 
                    save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)
tabulamuris<-RunTSNE(tabulamuris,dims.use = 1:75)
t1<-TSNEPlot(tabulamuris,do.label = TRUE, do.return =TRUE,do.legend = FALSE)
t1
TSNEPlot(tabulamuris)
FeaturePlot(tabulamuris,features.plot = "tissue")
FeaturePlot(tabulamuris,features.plot = c("Myf5","Foxg1"),cols.use = c("grey","red"))
DotPlot(tabulamuris,genes.plot = c("Pax3","Foxg1"))
t1<-DimPlot(tabulamuris,reduction.use = "tsne",group.by = "res.3", vector.friendly = TRUE, do.return = TRUE, no.legend = T,do.label = T)
t2<-DimPlot(tabulamuris,reduction.use = "tsne",group.by = "tissue", vector.friendly = TRUE, do.return = TRUE, no.legend = T, do.label = FALSE)
t1+theme(text = element_text(size = 12))
t2+theme(text = element_text(size = 12))
plot_grid(t1,t2)
#Use the ordering of the cells to bring the + cells to the front 
x<-FeatureHeatmap(tabulamuris,features.plot = c("Foxg1"),
               group.by = "tissue",pt.size = .1,sep.scale = F,
               do.return = TRUE,cols.use = c("grey86","red"),pch.use = 16,rotate.key = FALSE,key.position = "bottom"
               )
x$data<-x$data[order(x$data$expression),]
x+theme(text = element_text(size = 12))
DotPlot(tabulamuris,genes.plot = c("Foxg1","Myf5"),group.by = "tissue")
DotPlot(tabulamuris,genes.plot = c("Foxg1","Myf5"),group.by = "tissue",plot.legend = TRUE, do.return = TRUE)

#Violiin plot with different tissue 
Foxg1<-as.matrix(Foxg1)
Foxg1<-(`Muscle-counts.csv`@data[grep(pattern = "Foxg1", x = (x = `Muscle-counts.csv`@data), value = TRUE),1:2102])
Foxg1[1:5,1:5]
Foxg1$Bladder<-(`Bladder-counts.csv`@data[grep(pattern = "Foxg1", x = rownames(x = `Bladder-counts.csv`@data), value = TRUE),1:1638])

violin<-as.data.frame(1)
violin<-melt(Foxg1)
violin$tissue<-"Muslce"


ggplot(violin,aes(tissue,(log10(value))))+geom_violin()
ggplot(iris,aes(Species,Sepal.Length))+
  geom_violin()
#Probe the marrow data further 
geneplots<-c("Myf5","Foxg1","Meis2")
DotPlot(`Marrow-counts.csv`,geneplots)
VlnPlot(`Thymus-counts.csv`,geneplots)
