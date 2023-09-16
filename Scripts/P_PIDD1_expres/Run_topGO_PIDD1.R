##Set up the input files for topGO
library(GenomicRanges);
library(GenomicFeatures);
library(GenomicAlignments);
library(biomaRt);
library(AnnotationDbi);
library(org.Hs.eg.db);
library(topGO);
options(stringsAsFactors=FALSE);
## DEG genes
DEG = read.csv("Control_pidd1_data_cellmarkers.csv");
## C gene list
DEG.C = DEG[which(DEG$cluster=="nonPIDD1_cells"),];
DEG.C = DEG.C[which(DEG.C$p_val_adj < 0.05 & DEG.C$avg_log2FC > 0),];
## P gene list
DEG.P = DEG[which(DEG$cluster=="PIDD1_cells"),];
DEG.P = DEG.P[which(DEG.P$p_val_adj < 0.05 & DEG.P$avg_log2FC > 0),];
## find protein coding promoters and non-protein coding promoters
load("/home/dl2248/project/CoWithDuyPhan/GeneOntology/hgncdata+hg38.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
##All genes
geneNames= unique(codgene$entrezgene_id);
##C gene gene
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.C$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
C_GOdata <- new("topGOdata",ontology = "BP",description = "Control nonPIDD1_cells Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(C_GOdata, algorithm = "classic", statistic = "fisher");
CallRes <- GenTable(C_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## P genes
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.P$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
P_GOdata <- new("topGOdata",ontology = "BP",description = "Control PIDD1_cells Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(P_GOdata, algorithm = "classic", statistic = "fisher");
PallRes <- GenTable(P_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## save resutls
CallRes$Condition = "C";
PallRes$Condition = "P";
output = rbind(CallRes,PallRes);
write.csv(output, file = "Top_GO_Control_pidd1_results.csv",quote=F,row.names=F);
## plot
pdf(file="TopGO_Control_pidd1_results.pdf",height=8,width=6);

par(oma=c(2,10,1,1));
    bp = barplot(-log10(as.numeric(CallRes$classicFisher)),main="Control nonPIDD1_cells Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=CallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

    bp = barplot(-log10(as.numeric(PallRes$classicFisher)),main="Control PIDD1_cells Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=PallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

dev.off();



## DEG genes
DEG = read.csv("Patient_pidd1_data_cellmarkers.csv");
## C gene list
DEG.C = DEG[which(DEG$cluster=="nonPIDD1_cells"),];
DEG.C = DEG.C[which(DEG.C$p_val_adj < 0.05 & DEG.C$avg_log2FC > 0),];
## P gene list
DEG.P = DEG[which(DEG$cluster=="PIDD1_cells"),];
DEG.P = DEG.P[which(DEG.P$p_val_adj < 0.05 & DEG.P$avg_log2FC > 0),];
## find protein coding promoters and non-protein coding promoters
load("/home/dl2248/project/CoWithDuyPhan/GeneOntology/hgncdata+hg38.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
##All genes
geneNames= unique(codgene$entrezgene_id);
##C gene gene
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.C$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
C_GOdata <- new("topGOdata",ontology = "BP",description = "Patient nonPIDD1_cells Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(C_GOdata, algorithm = "classic", statistic = "fisher");
CallRes <- GenTable(C_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## P genes
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.P$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
P_GOdata <- new("topGOdata",ontology = "BP",description = "Patient PIDD1_cells Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(P_GOdata, algorithm = "classic", statistic = "fisher");
PallRes <- GenTable(P_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## save resutls
CallRes$Condition = "C";
PallRes$Condition = "P";
output = rbind(CallRes,PallRes);
write.csv(output, file = "Top_GO_Patient_pidd1_results.csv",quote=F,row.names=F);
## plot
pdf(file="TopGO_Patient_pidd1_results.pdf",height=8,width=6);

par(oma=c(2,10,1,1));
    bp = barplot(-log10(as.numeric(CallRes$classicFisher)),main="Patient nonPIDD1_cells Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=CallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

    bp = barplot(-log10(as.numeric(PallRes$classicFisher)),main="Patient PIDD1_cells Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=PallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

dev.off();





## DEG genes
DEG = read.csv("ControlPatient_pidd1_data_cellmarkers.csv");
## C gene list
DEG.C = DEG[which(DEG$cluster=="controlPIDD1_cells"),];
DEG.C = DEG.C[which(DEG.C$p_val_adj < 0.05 & DEG.C$avg_log2FC > 0),];
## P gene list
DEG.P = DEG[which(DEG$cluster=="patientPIDD1_cells"),];
DEG.P = DEG.P[which(DEG.P$p_val_adj < 0.05 & DEG.P$avg_log2FC > 0),];
## find protein coding promoters and non-protein coding promoters
load("/home/dl2248/project/CoWithDuyPhan/GeneOntology/hgncdata+hg38.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
##All genes
geneNames= unique(codgene$entrezgene_id);
##C gene gene
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.C$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
C_GOdata <- new("topGOdata",ontology = "BP",description = "Control PIDD1_cells Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(C_GOdata, algorithm = "classic", statistic = "fisher");
CallRes <- GenTable(C_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## P genes
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.P$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
P_GOdata <- new("topGOdata",ontology = "BP",description = "Patient PIDD1_cells Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(P_GOdata, algorithm = "classic", statistic = "fisher");
PallRes <- GenTable(P_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## save resutls
CallRes$Condition = "C";
PallRes$Condition = "P";
output = rbind(CallRes,PallRes);
write.csv(output, file = "Top_GO_ControlPatient_pidd1_results.csv",quote=F,row.names=F);
## plot
pdf(file="TopGO_ControlPatient_pidd1_results.pdf",height=8,width=6);

par(oma=c(2,10,1,1));
    bp = barplot(-log10(as.numeric(CallRes$classicFisher)),main="Control PIDD1_cells Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=CallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

    bp = barplot(-log10(as.numeric(PallRes$classicFisher)),main="Patient PIDD1_cells Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=PallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

dev.off();








