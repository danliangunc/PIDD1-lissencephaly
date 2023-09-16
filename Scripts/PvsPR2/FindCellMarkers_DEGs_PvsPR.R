library(Seurat);
library(ggplot2);
library(dplyr);
library(RColorBrewer);
library(pheatmap);
library(SingleCellExperiment);
library(DoubletFinder);
options(stringsAsfactors=F);
set.seed(20220504);
## read data
data.integrated = readRDS("../Merged_SeuratAnalysis/data.integrated.labeled.rds");
table(data.integrated$orig.ident);
## Condition 2- Patient vs. Control
## Patient (4): 3751D70Org1_CZ15, 3751D70Org2_CZ16, 3751D70Org3_CZ17, and C9Org12522−CZ21
## Control (4): Y6D70Org1_CZ1, Y6D70Org2_CZ2, Y6D70Org3_CZ3, and YB712522−CZ24

table(data.integrated$orig.ident,data.integrated$CellType);
table(data.integrated$orig.ident,data.integrated$seurat_clusters);
Idents(data.integrated) = "CellType";
## only save P and C
## add conditions
table(data.integrated$orig.ident,data.integrated$CellType);
table(data.integrated$orig.ident,data.integrated$seurat_clusters);
P_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
PR_samples = c("C9HomoOrg2");
data.integrated = subset(data.integrated, subset = orig.ident %in% c(P_samples,PR_samples));

## cell type markers
cellmarkers = FindAllMarkers(data.integrated,logfc.threshold = 0.25);
write.csv(cellmarkers,file="Merged_Data_cellmarkers_PvsPR.csv",quote=F,row.names=F);

## add conditions
data.integrated$PvsPR = "NA";
data.integrated$PvsPR[which(data.integrated$orig.ident %in% P_samples)] = "P";
data.integrated$PvsPR[which(data.integrated$orig.ident %in% PR_samples)] = "PR";

## cell types
## Astrocyte
Astro.integrated = subset(data.integrated, subset = CellType=="Astrocyte");
table(Astro.integrated$PvsPR);
Idents(Astro.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(Astro.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Astrocyte_DEG_PvsPR.csv");

## neuron
Neuron.integrated = subset(data.integrated, subset = CellType=="Neuron");
table(Neuron.integrated$PvsPR);
Idents(Neuron.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(Neuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Neuron_DEG_PvsPR.csv");

## Intermediate progenitor
IPC.integrated = subset(data.integrated, subset = CellType=="Intermediate progenitor");
table(IPC.integrated$PvsPR);
Idents(IPC.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(IPC.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "IntermediateProgenitor_DEG_PvsPR.csv");

## Deep layer excitatory Neuron
DLNeuron.integrated = subset(data.integrated, subset = CellType=="Deep layer excitatory neuron");
table(DLNeuron.integrated$PvsPR);
Idents(DLNeuron.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(DLNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DLNeuron_DEG_PvsPR.csv");

## NPCtoN transition cells
NPCtoN.integrated = subset(data.integrated, subset = CellType=="NPCtoN cells");
table(NPCtoN.integrated$PvsPR);
Idents(NPCtoN.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(NPCtoN.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "NPCtoNcells_DEG_PvsPR.csv");

## Upper layer excitatory Neuron
ULNeuron.integrated = subset(data.integrated, subset = CellType=="Upper layer excitatory neuron");
table(ULNeuron.integrated$PvsPR);
Idents(ULNeuron.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(ULNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "ULNeuron_DEG_PvsPR.csv");

## Radial Glia
RG.integrated = subset(data.integrated, subset = CellType=="Radial Glia");
table(RG.integrated$PvsPR);
Idents(RG.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(RG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "RadialGlia_DEG_PvsPR.csv");

## Outer Radial Glia
ORG.integrated = subset(data.integrated, subset = CellType=="Outer Radial Glia");
table(ORG.integrated$PvsPR);
Idents(ORG.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(ORG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "OuterRadialGlia_DEG_PvsPR.csv");

## Dividing Radial Glia/Outer Radial Glia
DRG.integrated = subset(data.integrated, subset = CellType=="Dividing/Outer Radial Glia");
table(DRG.integrated$PvsPR);
Idents(DRG.integrated) = "PvsPR";
cellmarkers = FindAllMarkers(DRG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DividingRadialGlia_DEG_PvsPR.csv");




