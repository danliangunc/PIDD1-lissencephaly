library(Seurat);
library(ggplot2);
library(dplyr);
library(RColorBrewer);
library(pheatmap);
library(SingleCellExperiment);
library(DoubletFinder);
options(stringsAsfactors=F);
set.seed(20220504);
## meta data 
metadata = read.csv("/home/dl2248/project/CoCeZhang/Merged_SeuratAnalysis/RDSfiles.csv");
## read rds files
P1srat = readRDS(metadata$RDSFile[1]);
P2srat = readRDS(metadata$RDSFile[2]);
P3srat = readRDS(metadata$RDSFile[3]);
R1srat = readRDS(metadata$RDSFile[4]);
R2srat = readRDS(metadata$RDSFile[5]);
P4srat = readRDS(metadata$RDSFile[6]);
LFsrat = readRDS(metadata$RDSFile[7]);
C1srat = readRDS(metadata$RDSFile[8]);
C2srat = readRDS(metadata$RDSFile[9]);
C3srat = readRDS(metadata$RDSFile[10]);
C4srat = readRDS(metadata$RDSFile[11]);
CKI1srat = readRDS(metadata$RDSFile[12]);
CKI2srat = readRDS(metadata$RDSFile[13]);
## data list
data.list = list(P1srat,P2srat,P3srat,P4srat,C1srat,C2srat,C3srat,C4srat,CKI1srat,CKI2srat,LFsrat,R1srat,R2srat);
## Normalize and find Variable Features
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures=5000, verbose = FALSE)
})
## features
features <- SelectIntegrationFeatures(object.list = data.list);
## scale and PCA data
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
## integrated using anchors
anchors <- FindIntegrationAnchors(object.list = data.list, reduction = "cca",dims = 1:50);
data.integrated <- IntegrateData(anchorset = anchors, dims = 1:50);
## save data
saveRDS(data.integrated,file="data.integrated.rds");
## integrated data
DefaultAssay(data.integrated) <- "integrated";
## genes need expressed in at least 100 cells
selected_gene <- rownames(data.integrated[Matrix::rowSums(data.integrated) > 100]);
# Run the standard workflow for visualization and clustering
data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, npcs = 55, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:55)
data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:55)
data.integrated <- FindClusters(data.integrated, resolution = 0.8);
## save data
saveRDS(data.integrated,file="data.integrated.rds");
#plot
pdf("MergedData.pdf");
DefaultAssay(data.integrated)="RNA";
DimPlot(data.integrated, group.by = "orig.ident");
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE);

## Radial Glia
FeaturePlot(data.integrated,features ="SOX2") + ggtitle("SOX2: Radial Glia");
FeaturePlot(data.integrated,features ="VIM") + ggtitle("VIM: Radial Glia");
## oRG
FeaturePlot(data.integrated,features ="HOPX") + ggtitle("HOPX: oRG");
FeaturePlot(data.integrated,features ="FAM107A") + ggtitle("FAM107A: oRG");
FeaturePlot(data.integrated,features ="LIFR") + ggtitle("LIFR: oRG");
## Neuron
FeaturePlot(data.integrated,features ="DCX") + ggtitle("DCX: Neuron");
FeaturePlot(data.integrated,features ="NEUROD2") + ggtitle("NEUROD2: Neuron");
FeaturePlot(data.integrated,features ="RBFOX3") + ggtitle("RBFOX3: Neuron");
FeaturePlot(data.integrated,features ="SLC17A7") + ggtitle("SLC17A7: Neuron");
## Deep layer excitatory neuron
FeaturePlot(data.integrated,features ="BCL11B") + ggtitle("BCL11B: Deep layer excitatory neuron");
FeaturePlot(data.integrated,features ="TBR1") + ggtitle("TBR1: Deep layer excitatory neuron");
FeaturePlot(data.integrated,features ="NEUROD2") + ggtitle("NEUROD2: Deep layer excitatory neuron");
## Upper layer excitatory neuron
FeaturePlot(data.integrated,features ="SATB2") + ggtitle("SATB2: Upper layer excitatory neuron");
FeaturePlot(data.integrated,features ="FRMD4B") + ggtitle("FRMD4B: Upper layer excitatory neuron");
FeaturePlot(data.integrated,features ="NEUROD2") + ggtitle("NEUROD2: Upper layer excitatory neuron");
## Immature neuron
FeaturePlot(data.integrated,features ="DCX") + ggtitle("DCX: Immature neuron");
FeaturePlot(data.integrated,features ="STMN2") + ggtitle("STMN2: Immature neuron");
FeaturePlot(data.integrated,features ="GAP43") + ggtitle("GAP43: Immature neuron");
FeaturePlot(data.integrated,features ="NEUROD2") + ggtitle("NEUROD2: Immature neuron");
## Astrocyte
FeaturePlot(data.integrated,features ="S100B") + ggtitle("S100B: Astrocyte");
FeaturePlot(data.integrated,features ="GFAP") + ggtitle("GFAP: Astrocyte");
## Intermediate progenitor
FeaturePlot(data.integrated,features ="EOMES") + ggtitle("EOMES: Intermediate progenitor");
FeaturePlot(data.integrated,features ="TMEM158") + ggtitle("TMEM158: Intermediate progenitor");
## Proliferating IPC
FeaturePlot(data.integrated,features ="EOMES") + ggtitle("EOMES: Proliferating IPC");
FeaturePlot(data.integrated,features ="TMEM158") + ggtitle("TMEM158: Proliferating IPC");
FeaturePlot(data.integrated,features ="MKI67") + ggtitle("MKI67: Proliferating IPC");
FeaturePlot(data.integrated,features ="TOP2A") + ggtitle("TOP2A: Proliferating IPC");

markergenes = c("SOX2","HOPX","FAM107A","LIFR","DCX","STMN2","GAP43","BCL11B","TBR1","NEUROD2","SATB2","FRMD4B","RBFOX3","SLC17A7","S100B","GFAP","EOMES","TMEM158","MKI67","TOP2A");

VlnPlot(data.integrated, markergenes, stack = TRUE, sort = TRUE, flip = TRUE) + theme(legend.position = "none") + ggtitle("identities on x-axis");

dev.off();

pdf("CellType_anno.pdf",height=8,width=12);

data.integrated$CellType = "NPCtoN cells";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(9,10,19,21))] = "Radial Glia";
#data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(10,19))] = "Outer Radial Glia";
#data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(21))] = "Dividing/Outer Radial Glia";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(0,1,2,3,4,5,6,8,11,12,13,15,16,25))] = "Neuron";
#data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(0,5,8))] = "Deep layer excitatory neuron";
#data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(2,4,6,16))] = "Upper layer excitatory neuron";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(24))] = "Astrocyte";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(22,23))] = "Intermediate progenitor";
#data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(22))] = "Proliferating IPC";

DimPlot(data.integrated, group.by = "CellType");

DimPlot(data.integrated, group.by = "CellType", label=T);

DimPlot(data.integrated, group.by = "seurat_clusters", label=T);

DotPlot(data.integrated, features=markergenes, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("identities on x-axis");

DotPlot(data.integrated, features=markergenes, group.by = "CellType") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("identities on x-axis");

RG.data.integrated = subset(data.integrated, subset = CellType=="Radial Glia");
RG.data.integrated$CellType[which(RG.data.integrated$seurat_clusters %in% c(10,19))] = "Outer Radial Glia";
RG.data.integrated$CellType[which(RG.data.integrated$seurat_clusters %in% c(21))] = "Dividing/Outer Radial Glia"; 
DotPlot(RG.data.integrated, features=markergenes, group.by = "CellType") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("RG clusters");
DotPlot(RG.data.integrated, features=markergenes, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("RG clusters");

Neuron.data.integrated = subset(data.integrated, subset = CellType=="Neuron");
Neuron.data.integrated$CellType[which(Neuron.data.integrated$seurat_clusters %in% c(0,5,8))] = "Deep layer excitatory neuron";
Neuron.data.integrated$CellType[which(Neuron.data.integrated$seurat_clusters %in% c(2,4,6,16))] = "Upper layer excitatory neuron";
DotPlot(Neuron.data.integrated, features=markergenes, group.by = "CellType") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Neuron clusters");
DotPlot(Neuron.data.integrated, features=markergenes, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Neuron clusters");

#ImNeuron.data.integrated = subset(data.integrated, subset = CellType=="Immature neuron");
#DotPlot(ImNeuron.data.integrated, features=markergenes, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Immature Neuron clusters");

IPC.data.integrated = subset(data.integrated, subset = CellType=="Intermediate progenitor");
DotPlot(IPC.data.integrated, features=markergenes, group.by = "seurat_clusters") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Intermediate progenitor clusters");

data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(10,19))] = "Outer Radial Glia";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(21))] = "Dividing/Outer Radial Glia";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(0,5,8))] = "Deep layer excitatory neuron";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(2,4,6,16))] = "Upper layer excitatory neuron";
data.integrated$CellType[which(data.integrated$seurat_clusters %in% c(22))] = "Proliferating IPC";

DotPlot(data.integrated, features=markergenes, group.by = "CellType") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("identities on x-axis");

DimPlot(data.integrated, group.by = "CellType", label=T);

dev.off();

#cellmarkers = FindAllMarkers(data.integrated,logfc.threshold = 0.25);
#write.csv(cellmarkers,file="Merged_Data_cellmarkers.csv",quote=F,row.names=F);

saveRDS(data.integrated, file="data.integrated.labeled.rds");

## cell type comp
P_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3","YB712522-CZ24");
CKI_samples = c("YB7HOMO12522-CZ22", "YB7HomoOrg2-CZ33");
PR_samples = c("C9HomoOrg1-CZ32","C9HomoOrg2");
LF_samples = c("Y5LOF12522-CZ23");

data.integrated$Condition = "LF";
data.integrated$Condition[which(data.integrated$orig.ident %in% P_samples)] = "P";
data.integrated$Condition[which(data.integrated$orig.ident %in% C_samples)] = "C";
data.integrated$Condition[which(data.integrated$orig.ident %in% CKI_samples)] = "CKI";
data.integrated$Condition[which(data.integrated$orig.ident %in% PR_samples)] = "PR";

library(RColorBrewer);
cellcomp = table(data.integrated$Condition,data.integrated$CellType);
cellcompper = t(apply(cellcomp,1,function(x) x/sum(x)));
pdf("Cell_comp.pdf");
mycols = brewer.pal(5, "Set3");
barplot(cellcompper,las=2,col=mycols,legend=T, beside=T);
dev.off();








