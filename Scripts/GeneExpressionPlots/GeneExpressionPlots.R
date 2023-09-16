library(Seurat);
library(ggplot2);
library(dplyr);
library(RColorBrewer);
library(pheatmap);
library(SingleCellExperiment);
library(DoubletFinder);
library(xlsx)
options(stringsAsfactors=F);
set.seed(20220504);
## read data
data.integrated = readRDS("../Merged_SeuratAnalysis/data.integrated.labeled.rds");
## conditions
P_samples = c("3751D70Org1_CZ15", "3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3", "YB712522-CZ24");
CKI_samples = c("YB7HOMO12522-CZ22", "YB7HomoOrg2-CZ33");
PR_samples = c("C9HomoOrg1-CZ32","C9HomoOrg2");
## assign conditions
data.integrated$Condition = NA;
data.integrated$Condition[which(data.integrated$orig.ident %in% P_samples)]="Patient";
data.integrated$Condition[which(data.integrated$orig.ident %in% C_samples)]="Control";
data.integrated$Condition[which(data.integrated$orig.ident %in% CKI_samples)]="Control_KI";
data.integrated$Condition[which(data.integrated$orig.ident %in% PR_samples)]="Patient_Rescue";
data.integrated = subset(data.integrated, subset = orig.ident %in% c(P_samples,C_samples,CKI_samples,PR_samples));
## add conditions to cell types
data.integrated$Cell_Condition = paste(data.integrated$CellType,data.integrated$Condition);
Idents(data.integrated)="Cell_Condition";
## re-order the cell types
myLevels = unique(data.integrated$Cell_Condition)[c(8,18,29,39,4,15,21,32,1,13,25,33,2,12,22,31,5,11,23,36,6,16,26,34,3,14,27,35,7,19,24,38,10,17,30,37,9,20,28,40)];
Idents(data.integrated) =  factor(Idents(data.integrated), levels= myLevels);
## genes
AKT_MTOR_genes = read.xlsx("/home/dl2248/project/CoCeZhang/GeneExpressionPlots/Listofgenes.xlsx",1,header=F);
LISSENCEPHALY_genes = read.xlsx("/home/dl2248/project/CoCeZhang/GeneExpressionPlots/Listofgenes.xlsx",2,header=F);
CELLDEATH_genes = read.xlsx("/home/dl2248/project/CoCeZhang/GeneExpressionPlots/Listofgenes.xlsx",3,header=F);
P53_genes = read.xlsx("/home/dl2248/project/CoCeZhang/GeneExpressionPlots/Listofgenes.xlsx",4,header=F);
## plots
pdf("Selected_genes_plot.pdf",height=12,width=18);
DotPlot(data.integrated, features= unique(AKT_MTOR_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("AKT-MTOR");

DotPlot(data.integrated, features= c(LISSENCEPHALY_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("LISSENCEPHALY");

DotPlot(data.integrated, features= unique(CELLDEATH_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("CELL DEATH");

DotPlot(data.integrated, features= unique(P53_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("P53");

celltypes = unique(data.integrated$CellType);
for (i in 1:length(celltypes)) {
    p = DotPlot(subset(data.integrated,subset = CellType==celltypes[i]), group.by="Cell_Condition", features= unique(AKT_MTOR_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("AKT-MTOR");
    print(p);

    p = DotPlot(subset(data.integrated,subset = CellType==celltypes[i]), group.by="Cell_Condition", features= unique(LISSENCEPHALY_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("LISSENCEPHALY");
    print(p);

    p = DotPlot(subset(data.integrated,subset = CellType==celltypes[i]), group.by="Cell_Condition", features= unique(CELLDEATH_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("CELL DEATH");
    print(p);

    p = DotPlot(subset(data.integrated,subset = CellType==celltypes[i]), group.by="Cell_Condition", features= unique(P53_genes$X1)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("P53");
    print(p);
}
dev.off();

pdf("AKT_MTOR_genes_plot.pdf",height=12,width=18);
AKT_MTOR_genes = read.table("/home/dl2248/project/CoCeZhang/GeneExpressionPlots/AKT_MTOR_genes.txt",header=1);
p = DotPlot(data.integrated, group.by="Condition", features= unique(AKT_MTOR_genes$Gene)) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("AKT MTOR genes");
print(p);
dev.off();

pdf("Heatmap_Oxidative_Phosphorylation_Translation_ribosome.pdf",height=20,width=10);
Translation_ribosome = read.xlsx("/home/dl2248/project/CoCeZhang/GeneExpressionPlots/Heatmapgenes.xls",1);
Oxidative_Phosphorylation = read.xlsx("/home/dl2248/project/CoCeZhang/GeneExpressionPlots/Heatmapgenes.xls",2,header=F);
## average expresion
plotmat=AverageExpression(data.integrated,features =  Oxidative_Phosphorylation$X1, group.by = "Condition")$RNA;
pheatmap(plotmat, cluster_cols=F,cluster_rows=T,fontsize = 5,scale="row",main = "Oxidative Phosphorylation");

plotmat=AverageExpression(data.integrated,features = Translation_ribosome$Symbol, group.by = "Condition")$RNA;
pheatmap(plotmat, cluster_cols=F,cluster_rows=T,fontsize = 5,scale="row",main = "Translation ribosome");

celltypes = unique(data.integrated$CellType);
for (i in 1:length(celltypes)) {
    plotmat=AverageExpression(subset(data.integrated,subset = CellType==celltypes[i]),features =  Oxidative_Phosphorylation$X1, group.by = "Cell_Condition")$RNA;
    p = pheatmap(plotmat, cluster_cols=F,cluster_rows=T,fontsize = 5,scale="row",main = "Oxidative Phosphorylation");
    print(p);
    
    plotmat=AverageExpression(subset(data.integrated,subset = CellType==celltypes[i]),features =  Translation_ribosome$Symbol, group.by = "Cell_Condition")$RNA;
    p = pheatmap(plotmat, cluster_cols=F,cluster_rows=T,fontsize = 5,scale="row",main = "Translation ribosome");
    print(p);
}

dev.off();


