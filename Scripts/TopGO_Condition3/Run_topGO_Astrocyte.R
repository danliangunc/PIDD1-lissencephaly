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
AstrDEG = read.csv("../Merged_SeuratAnalysis/Astrocyte_DEG_Condition3_C_CKI.csv",row.names=1);
## C gene list
AstrDEG.C = AstrDEG[which(AstrDEG$cluster=="C"),];
AstrDEG.C = AstrDEG.C[which(AstrDEG.C$p_val_adj < 0.05 & AstrDEG.C$avg_log2FC > 0.5),];
## CKI gene list
AstrDEG.CKI = AstrDEG[which(AstrDEG$cluster=="CKI"),];
AstrDEG.CKI = AstrDEG.CKI[which(AstrDEG.CKI$p_val_adj < 0.05 & AstrDEG.CKI$avg_log2FC > 0.5),];
## find protein coding promoters and non-protein coding promoters
load("/home/dl2248/project/CoWithDuyPhan/GeneOntology/hgncdata+hg38.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
##All genes
geneNames= unique(codgene$entrezgene_id);
##C gene gene
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% AstrDEG.C$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
C_GOdata <- new("topGOdata",ontology = "BP",description = "Astrocytes Control Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(C_GOdata, algorithm = "classic", statistic = "fisher");
CallRes <- GenTable(C_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## CKI genes
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% AstrDEG.CKI$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
CKI_GOdata <- new("topGOdata",ontology = "BP",description = "Astrocytes Control KI Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(CKI_GOdata, algorithm = "classic", statistic = "fisher");
CKIallRes <- GenTable(CKI_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## save resutls
CallRes$Condition = "C";
CKIallRes$Condition = "CKI";
output = rbind(CallRes,CKIallRes);
write.csv(output, file = "Top_GO_Astrocytes_results_Condition3_C_CKI.csv",quote=F,row.names=F);
## plot
pdf(file="TopGO_Astrocytes_DEGene_Condition3_C_CKI.pdf",height=8,width=6);

par(oma=c(2,10,1,1));
    bp = barplot(-log10(as.numeric(CallRes$classicFisher)),main="Astrocytes Control Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=CallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

    bp = barplot(-log10(as.numeric(CKIallRes$classicFisher)),main="Astrocytes Control KI Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=CKIallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

dev.off();

## DEG genes
AstrDEG = read.csv("../Merged_SeuratAnalysis/Astrocyte_DEG_Condition3_P_PR.csv",row.names=1);
## P gene list
AstrDEG.P = AstrDEG[which(AstrDEG$cluster=="P"),];
AstrDEG.P = AstrDEG.P[which(AstrDEG.P$p_val_adj < 0.05 & AstrDEG.P$avg_log2FC > 0.5),];
## PR gene list
AstrDEG.PR = AstrDEG[which(AstrDEG$cluster=="PR"),];
AstrDEG.PR = AstrDEG.PR[which(AstrDEG.PR$p_val_adj < 0.05 & AstrDEG.PR$avg_log2FC > 0.5),];
## find protein coding promoters and non-protein coding promoters
load("/home/dl2248/project/CoWithDuyPhan/GeneOntology/hgncdata+hg38.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
##All genes
geneNames= unique(codgene$entrezgene_id);
##P gene gene
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% AstrDEG.P$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
P_GOdata <- new("topGOdata",ontology = "BP",description = "Astrocytes Patient Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(P_GOdata, algorithm = "classic", statistic = "fisher");
PallRes <- GenTable(P_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## PR genes
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% AstrDEG.PR$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
PR_GOdata <- new("topGOdata",ontology = "BP",description = "Astrocytes Patient Rescue Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(PR_GOdata, algorithm = "classic", statistic = "fisher");
PRallRes <- GenTable(PR_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## save resutls
PallRes$Condition3 = "P";
PRallRes$Condition3 = "PR";
output = rbind(PallRes,PRallRes);
write.csv(output, file = "Top_GO_Astrocytes_results_Condition3_P_PR.csv",quote=F,row.names=F);
## plot
pdf(file="TopGO_Astrocytes_DEGene_Condition3_P_PR.pdf",height=8,width=6);
par(oma=c(2,10,1,1));
    bp = barplot(-log10(as.numeric(PallRes$classicFisher)),main="Astrocytes Patient Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=PallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

    bp = barplot(-log10(as.numeric(PRallRes$classicFisher)),main="Astrocytes Patient Rescue Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=PRallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

dev.off();



