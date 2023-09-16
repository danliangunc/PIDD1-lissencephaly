library(Seurat);
library(ggplot2);
library(dplyr);
library(RColorBrewer);
library(pheatmap);
library(SingleCellExperiment);
library(DoubletFinder);
options(stringsAsfactors=F);
args = commandArgs(trailingOnly=TRUE);
sampleID = args[1];
## PCA dim numbers
PCA_ndims = 1:20;
STC_ndims=1:30;
## selected number of features
n_sfeatures = 2000;
## read cellranger filtered count matrix
filedir = paste0("/home/dl2248/project/CoCeZhang/CellRanger/",sampleID,"_HHT_cellranger/filtered_feature_bc_matrix");
if (!file.exists(filedir)) {
   filedir = paste0("/home/dl2248/project/CoCeZhang/CellRanger/",sampleID,"_HHT_cellranger/outs/filtered_feature_bc_matrix");
}

adj.matrix <- Read10X(filedir);
srat <- CreateSeuratObject(adj.matrix,project = sampleID); 
## label Michochondrial genes
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-");
## label ribosomal proteins (their names begin with RPS or RPL)
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]");
## output pdf
pdf(paste0(sampleID,"_Seurat.pdf"));
## violin plots of the selected metadata features
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
## metadata features correlation
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt") & labs(subtitle="Raw data");
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") & labs(subtitle="Raw data");
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb") & labs(subtitle="Raw data");
FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt") & labs(subtitle="Raw data");
## Normalization and PCA/UMAP
srat <- NormalizeData(srat);
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = n_sfeatures);
all.genes <- rownames(srat);
srat <- ScaleData(srat, features = all.genes);
srat <- RunPCA(srat, features = VariableFeatures(object = srat),npcs = 100);
ElbowPlot(srat,ndims = 100) & labs(subtitle="Raw data");
srat <- FindNeighbors(srat, dims = PCA_ndims);
srat <- FindClusters(srat);
srat <- RunUMAP(srat, dims = PCA_ndims, verbose = F);
## cell cycles
s.genes <- cc.genes.updated.2019$s.genes;
g2m.genes <- cc.genes.updated.2019$g2m.genes;
srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes);
## UMAP
DimPlot(srat,label.size = 4,repel = T,label = T) & labs(subtitle="Raw data");
FeaturePlot(srat, features = "percent.mt") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
FeaturePlot(srat, features = "percent.rb") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
FeaturePlot(srat, features = "nFeature_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
FeaturePlot(srat, features = "nCount_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
FeaturePlot(srat,features =c("G2M.Score"),label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
FeaturePlot(srat,features =c("S.Score"),label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
VlnPlot(srat,features = "percent.mt") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
VlnPlot(srat,features = "percent.rb") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
VlnPlot(srat,features = "nFeature_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
VlnPlot(srat,features = "nCount_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
VlnPlot(srat,features = "G2M.Score") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");
VlnPlot(srat,features = "S.Score") & theme(plot.title = element_text(size=10)) & labs(subtitle="Raw data");

## only keep the Pass cells
## set QC threshold
srat <- subset(x = srat, subset = nFeature_RNA > 500 & percent.mt < 10 );
selected_gene <- rownames(srat[Matrix::rowSums(srat) > 3]);
srat = subset(x = srat, features = selected_gene);
srat;

srat <- NormalizeData(object = srat, normalization.method = "LogNormalize", scale.factor = median(srat$nCount_RNA));
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = n_sfeatures);
all.genes <- rownames(srat);
srat <- ScaleData(srat, features = all.genes);
srat <- RunPCA(srat);
ElbowPlot(srat, n=30);
srat <- RunUMAP(srat, dims = PCA_ndims);

# define the expected number of doublet cellscells.
nExp <- round(ncol(srat) * 0.075)  # expect 7.5% doublets
srat <- doubletFinder_v3(srat, pN = 0.25, pK = 0.09, nExp = nExp, PCs = PCA_ndims);
# name of the DF prediction can change, so extract the correct column name.
colnames(srat@meta.data)[grepl("DF.classification", colnames(srat@meta.data))] = "Doublet";
colnames(srat@meta.data)[grepl("pANN", colnames(srat@meta.data))] = "DoubletScore";
DimPlot(srat, group.by = "Doublet");
srat = subset(x = srat, subset = Doublet == "Singlet");
srat;
saveRDS(srat, file = paste0(sampleID,"_filtered_srat.rds"));
## keep cells passed QC
FeaturePlot(srat, features = "DoubletScore") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
FeaturePlot(srat, features = "percent.mt") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
FeaturePlot(srat, features = "percent.rb") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
FeaturePlot(srat, features = "nFeature_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
FeaturePlot(srat, features = "nCount_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
FeaturePlot(srat, features = "G2M.Score") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
FeaturePlot(srat, features = "S.Score") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
VlnPlot(srat,features = "DoubletScore") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
VlnPlot(srat,features = "percent.mt") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
VlnPlot(srat,features = "percent.rb") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
VlnPlot(srat,features = "nFeature_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
VlnPlot(srat,features = "nCount_RNA") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
VlnPlot(srat,features = "G2M.Score") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
VlnPlot(srat,features = "S.Score") & theme(plot.title = element_text(size=10)) & labs(subtitle="QCed data");
## Normalization and PCA/UMAP
srat <- NormalizeData(srat);
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = n_sfeatures);
all.genes <- rownames(srat);
srat <- ScaleData(srat, features = all.genes);
srat <- RunPCA(srat, features = VariableFeatures(object = srat));
ElbowPlot(srat,ndims = 30);
srat <- FindNeighbors(srat, dims = PCA_ndims);
srat <- FindClusters(srat);
srat <- RunUMAP(srat, dims = PCA_ndims, verbose = F);
## UMAP
DimPlot(srat,label.size = 4,repel = T,label = T);
table(srat[[]]$seurat_clusters);
dev.off();
## options for clusters
pdf(paste0(sampleID,"_Seurat_celltype_manual.pdf"));
## Radial Glia
FeaturePlot(srat,features ="SOX2") + ggtitle("SOX2: Radial Glia");
FeaturePlot(srat,features ="HES1") + ggtitle("HES1: Radial Glia");
FeaturePlot(srat,features ="VIM") + ggtitle("VIM: Radial Glia");
## oRG
FeaturePlot(srat,features ="HOPX") + ggtitle("HOPX: oRG");
FeaturePlot(srat,features ="PTPRZ1") + ggtitle("PTPRZ1: oRG");
FeaturePlot(srat,features ="TNC") + ggtitle("TNC: oRG");
FeaturePlot(srat,features ="LGALS3") + ggtitle("LGALS3: oRG");
FeaturePlot(srat,features ="FAM107A") + ggtitle("FAM107A: oRG");
## Neuron
FeaturePlot(srat,features ="DCX") + ggtitle("DCX: Neuron");
FeaturePlot(srat,features ="STMN2") + ggtitle("STMN2: Neuron");
FeaturePlot(srat,features ="GAP43") + ggtitle("GAP43: Neuron");
## Deep layer excitatory neuron
FeaturePlot(srat,features ="BCL11B") + ggtitle("BCL11B: Deep layer excitatory neuron");
FeaturePlot(srat,features ="TBR1") + ggtitle("TBR1: Deep layer excitatory neuron");
FeaturePlot(srat,features ="NEUROD2") + ggtitle("NEUROD2: Deep layer excitatory neuron");
## Upper layer excitatory neuron
FeaturePlot(srat,features ="SATB2") + ggtitle("SATB2: Upper layer excitatory neuron");
FeaturePlot(srat,features ="FRMD4B") + ggtitle("FRMD4B: Upper layer excitatory neuron");
FeaturePlot(srat,features ="BTG1") + ggtitle("BTG1: Upper layer excitatory neuron");
FeaturePlot(srat,features ="NEUROD2") + ggtitle("NEUROD2: Upper layer excitatory neuron");
## Layer 1 neuron
FeaturePlot(srat,features ="RELN") + ggtitle("RELN: Upper layer excitatory neuron");
## Immature neuron
FeaturePlot(srat,features ="BCL11B") + ggtitle("BCL11B: Immature neuron");
FeaturePlot(srat,features ="TBR1") + ggtitle("TBR1: Immature neuron");
FeaturePlot(srat,features ="SATB2") + ggtitle("SATB2: Immature neuron");
## Astrocyte
FeaturePlot(srat,features ="S100B") + ggtitle("S100B: Astrocyte");
FeaturePlot(srat,features ="GFAP") + ggtitle("GFAP: Astrocyte");
## Intermediate progenitor
FeaturePlot(srat,features ="EOMES") + ggtitle("EOMES: Intermediate progenitor");
FeaturePlot(srat,features ="PPP1R17") + ggtitle("PPP1R17: Intermediate progenitor");
FeaturePlot(srat,features ="TMEM158") + ggtitle("TMEM158: Intermediate progenitor");
## Proliferation cell
FeaturePlot(srat,features ="MKI67") + ggtitle("MKI67: Proliferation cell");
FeaturePlot(srat,features ="TOP2A") + ggtitle("TOP2A: Proliferation cell");

markergenes = c("SOX2","HES1","VIM","HOPX","PTPRZ1","TNC","LGALS3","FAM107A","DCX","STMN2","GAP43","BCL11B","TBR1","NEUROD2","SATB2","FRMD4B","BTG1","RELN","S100B","GFAP","EOMES","PPP1R17","TMEM158","MKI67","TOP2A");

VlnPlot(srat, markergenes, stack = TRUE, sort = TRUE, flip = TRUE) + theme(legend.position = "none") + ggtitle("identities on x-axis");

DimPlot(srat, label = T , repel = T, label.size = 3);
DimPlot(srat, label = T , repel = T, label.size = 3) + NoLegend();
dev.off();




