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
## Condition 3.1 - Control vs CKI
## Control (4): Y6D70Org1_CZ1, Y6D70Org2_CZ2, Y6D70Org3_CZ3, and YB712522−CZ24
##  CKI_samples = c("YB7HOMO12522-CZ22", "YB7HomoOrg2-CZ33");
## Condition 3.2 - Patient vs Rescue
## Patient (4): 3751D70Org1_CZ15, 3751D70Org2_CZ16, 3751D70Org3_CZ17, and C9Org12522−CZ21
##  PR_samples = c("C9HomoOrg1-CZ32","C9HomoOrg2");

## only save P, C, PR, CKI
## add conditions
P_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3","YB712522-CZ24");
CKI_samples = c("YB7HOMO12522-CZ22", "YB7HomoOrg2-CZ33");
PR_samples = c("C9HomoOrg1-CZ32","C9HomoOrg2");
data.integrated = subset(data.integrated, subset = orig.ident %in% c(P_samples,C_samples,CKI_samples,PR_samples));

## add conditions
data.integrated$Condition3 = "NA";
data.integrated$Condition3[which(data.integrated$orig.ident %in% P_samples)] = "P";
data.integrated$Condition3[which(data.integrated$orig.ident %in% C_samples)] = "C";
data.integrated$Condition3[which(data.integrated$orig.ident %in% CKI_samples)] = "CKI";
data.integrated$Condition3[which(data.integrated$orig.ident %in% PR_samples)] = "PR";


## cell types
## Astrocyte
Astro.integrated = subset(data.integrated, subset = CellType=="Astrocyte" & Condition3 %in% c("C","CKI"));
table(Astro.integrated$Condition3);
Idents(Astro.integrated) = "Condition3";
cellmarkers = FindAllMarkers(Astro.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Astrocyte_DEG_Condition3_C_CKI.csv");

Astro.integrated = subset(data.integrated, subset = CellType=="Astrocyte" & Condition3 %in% c("P","PR"));
table(Astro.integrated$Condition3);
Idents(Astro.integrated) = "Condition3";
cellmarkers = FindAllMarkers(Astro.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Astrocyte_DEG_Condition3_P_PR.csv");

## neuron
Neuron.integrated = subset(data.integrated, subset = CellType=="Neuron" & Condition3 %in% c("C","CKI"));
table(Neuron.integrated$Condition3);
Idents(Neuron.integrated) = "Condition3";
cellmarkers = FindAllMarkers(Neuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Neuron_DEG_Condition3_C_CKI.csv");

Neuron.integrated = subset(data.integrated, subset = CellType=="Neuron" & Condition3 %in% c("P","PR"));
table(Neuron.integrated$Condition3);
Idents(Neuron.integrated) = "Condition3";
cellmarkers = FindAllMarkers(Neuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "Neuron_DEG_Condition3_P_PR.csv");

## Intermediate progenitor
IPC.integrated = subset(data.integrated, subset = CellType=="Intermediate progenitor" & Condition3 %in% c("C","CKI"));
table(IPC.integrated$Condition3);
Idents(IPC.integrated) = "Condition3";
cellmarkers = FindAllMarkers(IPC.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "IntermediateProgenitor_DEG_Condition3_C_CKI.csv");

IPC.integrated = subset(data.integrated, subset = CellType=="Intermediate progenitor" & Condition3 %in% c("P","PR"));
table(IPC.integrated$Condition3);
Idents(IPC.integrated) = "Condition3";
cellmarkers = FindAllMarkers(IPC.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "IntermediateProgenitor_DEG_Condition3_P_PR.csv");

## Deep layer excitatory Neuron
DLNeuron.integrated = subset(data.integrated, subset = CellType=="Deep layer excitatory neuron" & Condition3 %in% c("C","CKI"));
table(DLNeuron.integrated$Condition3);
Idents(DLNeuron.integrated) = "Condition3";
cellmarkers = FindAllMarkers(DLNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DLNeuron_DEG_Condition3_C_CKI.csv");

DLNeuron.integrated = subset(data.integrated, subset = CellType=="Deep layer excitatory neuron" & Condition3 %in% c("P","PR"));
table(DLNeuron.integrated$Condition3);
Idents(DLNeuron.integrated) = "Condition3";
cellmarkers = FindAllMarkers(DLNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DLNeuron_DEG_Condition3_P_PR.csv",row.names=F);

## NPCtoN transition cells
NPCtoN.integrated = subset(data.integrated, subset = CellType=="NPCtoN cells" & Condition3 %in% c("C","CKI"));
table(NPCtoN.integrated$Condition3);
Idents(NPCtoN.integrated) = "Condition3";
cellmarkers = FindAllMarkers(NPCtoN.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "NPCtoNcells_DEG_Condition3_C_CKI.csv");

NPCtoN.integrated = subset(data.integrated, subset = CellType=="NPCtoN cells" & Condition3 %in% c("P","PR"));
table(NPCtoN.integrated$Condition3);
Idents(NPCtoN.integrated) = "Condition3";
cellmarkers = FindAllMarkers(NPCtoN.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "NPCtoNcells_DEG_Condition3_P_PR.csv");

## Upper layer excitatory Neuron
ULNeuron.integrated = subset(data.integrated, subset = CellType=="Upper layer excitatory neuron" & Condition3 %in% c("C","CKI"));
table(ULNeuron.integrated$Condition3);
Idents(ULNeuron.integrated) = "Condition3";
cellmarkers = FindAllMarkers(ULNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "ULNeuron_DEG_Condition3_C_CKI.csv");

ULNeuron.integrated = subset(data.integrated, subset = CellType=="Upper layer excitatory neuron" & Condition3 %in% c("P","PR"));
table(ULNeuron.integrated$Condition3);
Idents(ULNeuron.integrated) = "Condition3";
cellmarkers = FindAllMarkers(ULNeuron.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "ULNeuron_DEG_Condition3_P_PR.csv");

## Radial Glia
RG.integrated = subset(data.integrated, subset = CellType=="Radial Glia" & Condition3 %in% c("C","CKI"));
table(RG.integrated$Condition3);
Idents(RG.integrated) = "Condition3";
cellmarkers = FindAllMarkers(RG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "RadialGlia_DEG_Condition3_C_CKI.csv",);

RG.integrated = subset(data.integrated, subset = CellType=="Radial Glia" & Condition3 %in% c("P","PR"));
table(RG.integrated$Condition3);
Idents(RG.integrated) = "Condition3";
cellmarkers = FindAllMarkers(RG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "RadialGlia_DEG_Condition3_P_PR.csv");

## Outer Radial Glia
ORG.integrated = subset(data.integrated, subset = CellType=="Outer Radial Glia" & Condition3 %in% c("C","CKI"));
table(ORG.integrated$Condition3);
Idents(ORG.integrated) = "Condition3";
cellmarkers = FindAllMarkers(ORG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "OuterRadialGlia_DEG_Condition3_C_CKI.csv");

ORG.integrated = subset(data.integrated, subset = CellType=="Outer Radial Glia" & Condition3 %in% c("P","PR"));
table(ORG.integrated$Condition3);
Idents(ORG.integrated) = "Condition3";
cellmarkers = FindAllMarkers(ORG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "OuterRadialGlia_DEG_Condition3_P_PR.csv");

## Dividing Radial Glia/Outer Radial Glia
DRG.integrated = subset(data.integrated, subset = CellType=="Dividing/Outer Radial Glia" & Condition3 %in% c("C","CKI"));
table(DRG.integrated$Condition3);
Idents(DRG.integrated) = "Condition3";
cellmarkers = FindAllMarkers(DRG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DividingRadialGlia_DEG_Condition3_C_CKI.csv");

DRG.integrated = subset(data.integrated, subset = CellType=="Dividing/Outer Radial Glia" & Condition3 %in% c("P","PR"));
table(DRG.integrated$Condition3);
Idents(DRG.integrated) = "Condition3";
cellmarkers = FindAllMarkers(DRG.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DividingRadialGlia_DEG_Condition3_P_PR.csv",row.names=F);





