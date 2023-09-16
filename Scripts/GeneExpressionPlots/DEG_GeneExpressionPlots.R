library(Seurat);
library(ggplot2);
library(dplyr);
library(RColorBrewer);
library(pheatmap);
library(ggvenn);
options(stringsAsfactors=F);
set.seed(20220504);
## read data
data.integrated = readRDS("../Merged_SeuratAnalysis/data.integrated.labeled.rds");
DefaultAssay(data.integrated) = "RNA";
## conditions
P_samples = c("3751D70Org1_CZ15", "3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3", "YB712522-CZ24");
CKI_samples = c("YB7HOMO12522-CZ22", "YB7HomoOrg2-CZ33");
PR_samples = c("C9HomoOrg1-CZ32");
## assign conditions
data.integrated$Condition = NA;
data.integrated$Condition[which(data.integrated$orig.ident %in% P_samples)]="Patient";
data.integrated$Condition[which(data.integrated$orig.ident %in% C_samples)]="Control";
data.integrated$Condition[which(data.integrated$orig.ident %in% CKI_samples)]="Control_KI";
data.integrated$Condition[which(data.integrated$orig.ident %in% PR_samples)]="Patient_Rescue";
data.integrated = subset(data.integrated, subset = orig.ident %in% c(P_samples,C_samples,CKI_samples,PR_samples));
## add conditions to cell types
data.integrated$Cell_Condition = paste(data.integrated$FinerCellType,data.integrated$Condition);
Idents(data.integrated)="Cell_Condition";

## cell type 
celltypes = unique(data.integrated$FinerCellType);
files = c("ULNeuron","NPCtoNcells","OuterRadialGlia","Neuron","DLNeuron","Astrocyte","RadialGlia","DividingRadialGlia","IntermediateProgenitor");
pdf("DEG_overlap_celltype.pdf");
for (i in 1:length(celltypes)) {
    DGEs_1 = read.csv(paste0("../Merged_SeuratAnalysis/",files[i],"_DEG_Condition1.csv"));
    DGEs_2 = read.csv(paste0("../Merged_SeuratAnalysis/",files[i],"_DEG_Condition2.csv"));
    DGEs_3_1 = read.csv(paste0("../Merged_SeuratAnalysis/",files[i],"_DEG_Condition3_C_CKI.csv"));
    DGEs_3_2 = read.csv(paste0("../Merged_SeuratAnalysis/",files[i],"_DEG_Condition3_P_PR.csv"));

    DGEs_1_C = DGEs_1[which(DGEs_1$cluster=="C"),];
    DGEs_1_C = DGEs_1_C[which(DGEs_1_C$p_val_adj<0.05 & DGEs_1_C$avg_log2FC>0.5 & DGEs_1_C$pct.1 >0.25),];
    DGEs_1_P = DGEs_1[which(DGEs_1$cluster=="P"),];
    DGEs_1_P = DGEs_1_P[which(DGEs_1_P$p_val_adj<0.05 & DGEs_1_P$avg_log2FC>0.5 & DGEs_1_P$pct.1 >0.25),];

    DGEs_2_C_PR = DGEs_2[which(DGEs_2$cluster=="C_PR"),];
    DGEs_2_C_PR = DGEs_2_C_PR[which(DGEs_2_C_PR$p_val_adj<0.05 & DGEs_2_C_PR$avg_log2FC>0.5 & DGEs_2_C_PR$pct.1 >0.25),];
    DGEs_2_P_CKI = DGEs_2[which(DGEs_2$cluster=="P_CKI"),];
    DGEs_2_P_CKI = DGEs_2_P_CKI[which(DGEs_2_P_CKI$p_val_adj<0.05 & DGEs_2_P_CKI$avg_log2FC>0.5 & DGEs_2_P_CKI$pct.1 >0.25),];

    DGEs_3_1_C = DGEs_3_1[which(DGEs_3_1$cluster=="C"),];
    DGEs_3_1_C = DGEs_3_1_C[which(DGEs_3_1_C$p_val_adj<0.05 & DGEs_3_1_C$avg_log2FC>0.5 & DGEs_3_1_C$pct.1 >0.25),];
    DGEs_3_1_CKI = DGEs_3_1[which(DGEs_3_1$cluster=="CKI"),];
    DGEs_3_1_CKI = DGEs_3_1_CKI[which(DGEs_3_1_CKI$p_val_adj<0.05 & DGEs_3_1_CKI$avg_log2FC>0.5 & DGEs_3_1_CKI$pct.1 >0.25),];

    DGEs_3_2_PR = DGEs_3_2[which(DGEs_3_2$cluster=="PR"),];
    DGEs_3_2_PR = DGEs_3_2_PR[which(DGEs_3_2_PR$p_val_adj<0.05 & DGEs_3_2_PR$avg_log2FC>0.5 & DGEs_3_2_PR$pct.1 >0.25),];
    DGEs_3_2_P = DGEs_3_2[which(DGEs_3_2$cluster=="P"),];
    DGEs_3_2_P = DGEs_3_2_P[which(DGEs_3_2_P$p_val_adj<0.05 & DGEs_3_2_P$avg_log2FC>0.5 & DGEs_3_2_P$pct.1 >0.25),];
    
    plotlist_c = list('Upregulated_C'=DGEs_1_C$gene,'Upregulated_CvsCKI'=DGEs_3_1_C$gene,'Upregulated_PR'=DGEs_3_2_PR$gene);
    p = ggvenn(plotlist_c,show_percentage=F) + ggtitle(paste(celltypes[i],"Control Like"));
    print(p);

    plotlist_p = list('Upregulated_P'=DGEs_1_P$gene,'Upregulated_CKI'=DGEs_3_1_CKI$gene,'Upregulated_PvsPR'=DGEs_3_2_P$gene);
    p = ggvenn(plotlist_p,show_percentage=F) + ggtitle(paste(celltypes[i],"Patient Like"));
    print(p);
}

dev.off();






