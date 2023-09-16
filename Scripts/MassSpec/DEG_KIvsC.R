library(Seurat);
#library(ggplot2);
#library(dplyr);
#library(RColorBrewer);
#library(pheatmap);
#library(SingleCellExperiment);
#library(DoubletFinder);
options(stringsAsfactors=F);
set.seed(20220504);
## read data
data.integrated=readRDS("../Merged_SeuratAnalysis/data.integrated.labeled.rds");
## cell type comp
P_samples = c("3751D70Org1_CZ15","3751D70Org2_CZ16","3751D70Org3_CZ17","C9Org12522-CZ21");
C_samples = c("Y6D70Org1_CZ1", "Y6D70Org2_CZ2", "Y6D70Org3_CZ3","YB712522-CZ24");
CKI_samples = c("YB7HOMO12522-CZ22", "YB7HomoOrg2-CZ33");
PR_samples = c("C9HomoOrg1-CZ32","C9HomoOrg2");
LF_samples = c("Y5LOF12522-CZ23");
# set conditions
data.integrated$Condition = "LF";
data.integrated$Condition[which(data.integrated$orig.ident %in% P_samples)] = "P";
data.integrated$Condition[which(data.integrated$orig.ident %in% C_samples)] = "C";
data.integrated$Condition[which(data.integrated$orig.ident %in% CKI_samples)] = "CKI";
data.integrated$Condition[which(data.integrated$orig.ident %in% PR_samples)] = "PR";
## subset
data.integrated = subset(data.integrated, subset = orig.ident %in% c(CKI_samples,C_samples));
## cell types
Idents(data.integrated) = "Condition";
cellmarkers = FindAllMarkers(data.integrated,logfc.threshold = 0);
write.csv(cellmarkers, file = "DEG_CKIvsC_s.csv");


