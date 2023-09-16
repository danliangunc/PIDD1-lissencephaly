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
## Condition 2- Patient + Control KI vs. Control + Patient Rescue
## Patient (4): 3751D70Org1_CZ15, 3751D70Org2_CZ16, 3751D70Org3_CZ17, and C9Org12522−CZ21
## Control (4): Y6D70Org1_CZ1, Y6D70Org2_CZ2, Y6D70Org3_CZ3, and YB712522−CZ24
## Control KI (2): YB7HOMO12522−CZ22, YB7HomoOrg2−CZ33 
## Patient Rescue (1): C9HomoOrg1−CZ32, C9HomoOrg2

table(data.integrated$orig.ident,data.integrated$CellType);

table(data.integrated$orig.ident,data.integrated$seurat_clusters);

Idents(data.integrated) = "CellType";

## add conditions
P_CKI_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21","YB7HOMO12522-CZ22", "YB7HomoOrg2-CZ33");
C_PR_samples=c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3","YB712522-CZ24","C9HomoOrg1-CZ32","C9HomoOrg2");

data.integrated = subset(data.integrated, subset = orig.ident %in% c(P_CKI_samples,C_PR_samples));
data.integrated$Condition2 = "NA";
data.integrated$Condition2[which(data.integrated$orig.ident %in% P_CKI_samples)] = "P_CKI";
data.integrated$Condition2[which(data.integrated$orig.ident %in% C_PR_samples)] = "C_PR";

## cell type markers
cellmarkers = FindAllMarkers(data.integrated,logfc.threshold = 0.25);
write.csv(cellmarkers,file="Merged_Data_cellmarkers_Condition2.csv",quote=F,row.names=F);

## cell types
## Astrocyte
Astro.integrated = subset(data.integrated, subset = CellType=="Astrocyte");
table(Astro.integrated$Condition2);
Idents(Astro.integrated) = "Condition2";
cellmarkers = FindAllMarkers(Astro.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Astrocyte_DEG_Condition2.csv");

## neuron
Neuron.integrated = subset(data.integrated, subset = CellType=="Neuron");
table(Neuron.integrated$Condition2);
Idents(Neuron.integrated) = "Condition2";
cellmarkers = FindAllMarkers(Neuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Neuron_DEG_Condition2.csv");

## Intermediate progenitor
IPC.integrated = subset(data.integrated, subset = CellType=="Intermediate progenitor");
table(IPC.integrated$Condition2);
Idents(IPC.integrated) = "Condition2";
cellmarkers = FindAllMarkers(IPC.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "IntermediateProgenitor_DEG_Condition2.csv");

## Deep layer excitatory Neuron
DLNeuron.integrated = subset(data.integrated, subset = CellType=="Deep layer excitatory neuron");
table(DLNeuron.integrated$Condition2);
Idents(DLNeuron.integrated) = "Condition2";
cellmarkers = FindAllMarkers(DLNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DLNeuron_DEG_Condition2.csv");

## NPCtoN transition cells
NPCtoN.integrated = subset(data.integrated, subset = CellType=="NPCtoN cells");
table(NPCtoN.integrated$Condition2);
Idents(NPCtoN.integrated) = "Condition2";
cellmarkers = FindAllMarkers(NPCtoN.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "NPCtoNcells_DEG_Condition2.csv");

## Upper layer excitatory Neuron
ULNeuron.integrated = subset(data.integrated, subset = CellType=="Upper layer excitatory neuron");
table(ULNeuron.integrated$Condition2);
Idents(ULNeuron.integrated) = "Condition2";
cellmarkers = FindAllMarkers(ULNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "ULNeuron_DEG_Condition2.csv");

## Radial Glia
RG.integrated = subset(data.integrated, subset = CellType=="Radial Glia");
table(RG.integrated$Condition2);
Idents(RG.integrated) = "Condition2";
cellmarkers = FindAllMarkers(RG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "RadialGlia_DEG_Condition2.csv");

## Outer Radial Glia
ORG.integrated = subset(data.integrated, subset = CellType=="Outer Radial Glia");
table(ORG.integrated$Condition2);
Idents(ORG.integrated) = "Condition2";
cellmarkers = FindAllMarkers(ORG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "OuterRadialGlia_DEG_Condition2.csv");

## Dividing Radial Glia/Outer Radial Glia
DRG.integrated = subset(data.integrated, subset = CellType=="Dividing/Outer Radial Glia");
table(DRG.integrated$Condition2);
Idents(DRG.integrated) = "Condition2";
cellmarkers = FindAllMarkers(DRG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DividingRadialGlia_DEG_Condition2.csv");







