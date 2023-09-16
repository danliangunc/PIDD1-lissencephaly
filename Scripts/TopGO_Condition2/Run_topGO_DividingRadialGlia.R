##Set up the input files for topGO for Hydrocephalus Risk gene and Microcephaly Risk gene
library(GenomicRanges);
library(GenomicFeatures);
library(GenomicAlignments);
library(biomaRt);
library(AnnotationDbi);
library(org.Hs.eg.db);
library(topGO);
options(stringsAsFactors=FALSE);
## DEG genes
DEG = read.csv("../Merged_SeuratAnalysis/DividingRadialGlia_DEG_Condition2.csv",row.names=1);
## C gene list
DEG.C_PR = DEG[which(DEG$cluster=="C_PR"),];
DEG.C_PR = DEG.C_PR[which(DEG.C_PR$p_val_adj < 0.05 & DEG.C_PR$avg_log2FC > 0.3),];
## P gene list
DEG.P_CKI = DEG[which(DEG$cluster=="P_CKI"),];
DEG.P_CKI = DEG.P_CKI[which(DEG.P_CKI$p_val_adj < 0.05 & DEG.P_CKI$avg_log2FC > 0.3),];
## find protein coding promoters and non-protein coding promoters
load("/home/dl2248/project/CoWithDuyPhan/GeneOntology/hgncdata+hg38.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
##All genes
geneNames= unique(codgene$entrezgene_id);
##C gene gene
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.C_PR$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
C_GOdata <- new("topGOdata",ontology = "BP",description = "Dividing Radial Glia Control Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(C_GOdata, algorithm = "classic", statistic = "fisher");
CallRes <- GenTable(C_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## P genes
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.P_CKI$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
P_GOdata <- new("topGOdata",ontology = "BP",description = "Dividing Radial Glia Patient Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(P_GOdata, algorithm = "classic", statistic = "fisher");
PallRes <- GenTable(P_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## save resutls
CallRes$Condition2 = "C_PR";
PallRes$Condition2 = "P_CKI";
output = rbind(CallRes,PallRes);
write.csv(output, file = "Top_GO_DividingRadialGlia_results_Condition2.csv",quote=F,row.names=F);
## plot
pdf(file="TopGO_DividingRadialGlia_DEGene_Condition2.pdf",height=8,width=6);
par(oma=c(2,10,1,1));
    bp = barplot(-log10(as.numeric(CallRes$classicFisher)),main="Dividing Radial Glia Control Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=CallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

    bp = barplot(-log10(as.numeric(PallRes$classicFisher)),main="Dividing Radial Glia Patient Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=PallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

dev.off();

