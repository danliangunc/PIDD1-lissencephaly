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
data.integrated = readRDS("data.integrated.labeled.rds");
table(data.integrated$orig.ident);
## Condition 2- Patient vs. Control
## Patient (4): 3751D70Org1_CZ15, 3751D70Org2_CZ16, 3751D70Org3_CZ17, and C9Org12522−CZ21
## Control (4): Y6D70Org1_CZ1, Y6D70Org2_CZ2, Y6D70Org3_CZ3, and YB712522−CZ24

table(data.integrated$orig.ident,data.integrated$CellType);

table(data.integrated$orig.ident,data.integrated$seurat_clusters);

Idents(data.integrated) = "CellType";

## cell type markers
cellmarkers = FindAllMarkers(data.integrated,logfc.threshold = 0.25);
write.csv(cellmarkers,file="Merged_Data_cellmarkers.csv",quote=F,row.names=F);

## only save P and C
## add conditions
P_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3","YB712522-CZ24");
data.integrated = subset(data.integrated, subset = orig.ident %in% c(P_samples,C_samples));

## cell type markers
cellmarkers = FindAllMarkers(data.integrated,logfc.threshold = 0.25);
write.csv(cellmarkers,file="Merged_Data_cellmarkers_Condition1.csv",quote=F,row.names=F);

## add conditions
data.integrated$Condition1 = "NA";
data.integrated$Condition1[which(data.integrated$orig.ident %in% P_samples)] = "P";
data.integrated$Condition1[which(data.integrated$orig.ident %in% C_samples)] = "C";

## cell types
## Astrocyte
Astro.integrated = subset(data.integrated, subset = CellType=="Astrocyte");
table(Astro.integrated$Condition1);
Idents(Astro.integrated) = "Condition1";
cellmarkers = FindAllMarkers(Astro.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Astrocyte_DEG_Condition1.csv");

## neuron
Neuron.integrated = subset(data.integrated, subset = CellType=="Neuron");
table(Neuron.integrated$Condition1);
Idents(Neuron.integrated) = "Condition1";
cellmarkers = FindAllMarkers(Neuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Neuron_DEG_Condition1.csv");

## Intermediate progenitor
IPC.integrated = subset(data.integrated, subset = CellType=="Intermediate progenitor");
table(IPC.integrated$Condition1);
Idents(IPC.integrated) = "Condition1";
cellmarkers = FindAllMarkers(IPC.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "IntermediateProgenitor_DEG_Condition1.csv");

## Deep layer excitatory Neuron
DLNeuron.integrated = subset(data.integrated, subset = CellType=="Deep layer excitatory neuron");
table(DLNeuron.integrated$Condition1);
Idents(DLNeuron.integrated) = "Condition1";
cellmarkers = FindAllMarkers(DLNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DLNeuron_DEG_Condition1.csv");

## NPCtoN transition cells
NPCtoN.integrated = subset(data.integrated, subset = CellType=="NPCtoN cells");
table(NPCtoN.integrated$Condition1);
Idents(NPCtoN.integrated) = "Condition1";
cellmarkers = FindAllMarkers(NPCtoN.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "NPCtoNcells_DEG_Condition1.csv");

## Upper layer excitatory Neuron
ULNeuron.integrated = subset(data.integrated, subset = CellType=="Upper layer excitatory neuron");
table(ULNeuron.integrated$Condition1);
Idents(ULNeuron.integrated) = "Condition1";
cellmarkers = FindAllMarkers(ULNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "ULNeuron_DEG_Condition1.csv");

## Radial Glia
RG.integrated = subset(data.integrated, subset = CellType=="Radial Glia");
table(RG.integrated$Condition1);
Idents(RG.integrated) = "Condition1";
cellmarkers = FindAllMarkers(RG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "RadialGlia_DEG_Condition1.csv");

## Outer Radial Glia
ORG.integrated = subset(data.integrated, subset = CellType=="Outer Radial Glia");
table(ORG.integrated$Condition1);
Idents(ORG.integrated) = "Condition1";
cellmarkers = FindAllMarkers(ORG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "OuterRadialGlia_DEG_Condition1.csv");

## Dividing Radial Glia/Outer Radial Glia
DRG.integrated = subset(data.integrated, subset = CellType=="Dividing/Outer Radial Glia");
table(DRG.integrated$Condition1);
Idents(DRG.integrated) = "Condition1";
cellmarkers = FindAllMarkers(DRG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DividingRadialGlia_DEG_Condition1.csv");




