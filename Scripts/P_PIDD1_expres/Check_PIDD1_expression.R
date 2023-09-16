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
data.integrated = readRDS("/home/dl2248/project/CoCeZhang/Merged_SeuratAnalysis/data.integrated.labeled.rds");
table(data.integrated$orig.ident);
## Patient (4): 3751D70Org1_CZ15, 3751D70Org2_CZ16, 3751D70Org3_CZ17, and C9Org12522−CZ21

## only save P 
P_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
data.integrated = subset(data.integrated, subset = orig.ident %in% c(P_samples));

## check PIDD1 expression level
pdf("PIDD1_patient.pdf")
DefaultAssay(data.integrated)="RNA";
Idents(data.integrated) = "CellType";
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE);
FeaturePlot(data.integrated,features ="PIDD1") + ggtitle("PIDD1");
dev.off();

## keep the top PIDD1 and tail PIDD1
PIDD1_expression = data.integrated@assays[["RNA"]]@counts["PIDD1",];
index = which(PIDD1_expression>0);

PIDD1_cells = names(PIDD1_expression)[index];
nonPIDD1_cells = sample(names(PIDD1_expression)[-index], length(PIDD1_cells));

## new data
pidd1_data = subset(data.integrated, cells = c(PIDD1_cells, nonPIDD1_cells));
pidd1_data@meta.data$condition = "PIDD1_cells";
pidd1_data@meta.data$condition[which(rownames(pidd1_data@meta.data) %in% nonPIDD1_cells)] = "nonPIDD1_cells";
## cell type markers
DefaultAssay(pidd1_data)="RNA";
Idents(pidd1_data) = "condition";
cellmarkers = FindAllMarkers(pidd1_data, logfc.threshold = 0);
write.csv(cellmarkers,file="Patient_pidd1_data_cellmarkers.csv",quote=F,row.names=F);
## find cell markers
PIDD1_cellmarkers = cellmarkers[which(cellmarkers$cluster=="PIDD1_cells" & cellmarkers$avg_log2FC>0 & cellmarkers$p_val_adj < 0.05),];
nonPIDD1_cellmarkers = cellmarkers[which(cellmarkers$cluster=="nonPIDD1_cells" & cellmarkers$avg_log2FC>0 & cellmarkers$p_val_adj < 0.05),];
write.csv(PIDD1_cellmarkers, file="Patient_PIDD1_cellmarkers.csv",quote=F,row.names=F);
write.csv(nonPIDD1_cellmarkers, file="Patient_nonPIDD1_cellmarkers.csv",quote=F,row.names=F);

## read data
data.integrated = readRDS("/home/dl2248/project/CoCeZhang/Merged_SeuratAnalysis/data.integrated.labeled.rds");
## only save C
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3","YB712522-CZ24");
data.integrated = subset(data.integrated, subset = orig.ident %in% c(C_samples));

## check PIDD1 expression level
pdf("PIDD1_control.pdf")
DefaultAssay(data.integrated)="RNA";
Idents(data.integrated) = "CellType";
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE);
FeaturePlot(data.integrated,features ="PIDD1") + ggtitle("PIDD1");
dev.off();

## keep the top PIDD1 and tail PIDD1
PIDD1_expression = data.integrated@assays[["RNA"]]@counts["PIDD1",];
index = which(PIDD1_expression>0);

PIDD1_cells = names(PIDD1_expression)[index];
nonPIDD1_cells = sample(names(PIDD1_expression)[-index], length(PIDD1_cells));

## new data
pidd1_data = subset(data.integrated, cells = c(PIDD1_cells, nonPIDD1_cells));
pidd1_data@meta.data$condition = "PIDD1_cells";
pidd1_data@meta.data$condition[which(rownames(pidd1_data@meta.data) %in% nonPIDD1_cells)] = "nonPIDD1_cells";
## cell type markers
DefaultAssay(pidd1_data)="RNA";
Idents(pidd1_data) = "condition";
cellmarkers = FindAllMarkers(pidd1_data, logfc.threshold = 0);
write.csv(cellmarkers,file="Control_pidd1_data_cellmarkers.csv",quote=F,row.names=F);
## find cell markers
PIDD1_cellmarkers = cellmarkers[which(cellmarkers$cluster=="PIDD1_cells" & cellmarkers$avg_log2FC>0 & cellmarkers$p_val_adj < 0.05),];
nonPIDD1_cellmarkers = cellmarkers[which(cellmarkers$cluster=="nonPIDD1_cells" & cellmarkers$avg_log2FC>0 & cellmarkers$p_val_adj < 0.05),];
write.csv(PIDD1_cellmarkers, file="Control_PIDD1_cellmarkers.csv",quote=F,row.names=F);
write.csv(nonPIDD1_cellmarkers, file="Control_nonPIDD1_cellmarkers.csv",quote=F,row.names=F);



## read data
data.integrated = readRDS("/home/dl2248/project/CoCeZhang/Merged_SeuratAnalysis/data.integrated.labeled.rds");
## only save C
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3","YB712522-CZ24");
control.integrated = subset(data.integrated, subset = orig.ident %in% c(C_samples));

## keep the top PIDD1
PIDD1_expression = control.integrated@assays[["RNA"]]@counts["PIDD1",];
index = which(PIDD1_expression>0);
controlPIDD1_cells = names(PIDD1_expression)[index];

## only save P
P_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
pat.integrated = subset(data.integrated, subset = orig.ident %in% c(P_samples));

## keep the top PIDD1 
PIDD1_expression = pat.integrated@assays[["RNA"]]@counts["PIDD1",];
index = which(PIDD1_expression>0);
patPIDD1_cells = names(PIDD1_expression)[index];

## new data
pidd1_data = subset(data.integrated, cells = c(controlPIDD1_cells, patPIDD1_cells));
pidd1_data@meta.data$condition = "controlPIDD1_cells";
pidd1_data@meta.data$condition[which(rownames(pidd1_data@meta.data) %in% patPIDD1_cells)] = "patientPIDD1_cells";
## cell type markers
DefaultAssay(pidd1_data)="RNA";
Idents(pidd1_data) = "condition";
cellmarkers = FindAllMarkers(pidd1_data, logfc.threshold = 0);
write.csv(cellmarkers,file="ControlPatient_pidd1_data_cellmarkers.csv",quote=F,row.names=F);
## find cell markers
control_cellmarkers = cellmarkers[which(cellmarkers$cluster=="controlPIDD1_cells" & cellmarkers$avg_log2FC>0 & cellmarkers$p_val_adj < 0.05),];
patient_cellmarkers = cellmarkers[which(cellmarkers$cluster=="patientPIDD1_cells" & cellmarkers$avg_log2FC>0 & cellmarkers$p_val_adj < 0.05),];
write.csv(control_cellmarkers, file="PIDD1_Control_cellmarkers.csv",quote=F,row.names=F);
write.csv(patient_cellmarkers, file="PIDD1_Patient_cellmarkers.csv",quote=F,row.names=F);














