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
DEG = read.csv("../Merged_SeuratAnalysis/Neuron_DEG_Condition2.csv",row.names=1);
## C_PR gene list
DEG.C_PR = DEG[which(DEG$cluster=="C_PR"),];
DEG.C_PR = DEG.C_PR[which(DEG.C_PR$p_val_adj < 0.05 & DEG.C_PR$avg_log2FC > 0.3),];
## P_CKI gene list
DEG.P_CKI = DEG[which(DEG$cluster=="P_CKI"),];
DEG.P_CKI = DEG.P_CKI[which(DEG.P_CKI$p_val_adj < 0.05 & DEG.P_CKI$avg_log2FC > 0.3),];
## find protein coding promoters and non-protein coding promoters
load("/home/dl2248/project/CoWithDuyPhan/GeneOntology/hgncdata+hg38.Rdata");
codgene=hgncdata[which(hgncdata$gene_biotype=="protein_coding"),];
##All genes
geneNames= unique(codgene$entrezgene_id);
##C_PR gene gene
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.C_PR$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
C_PR_GOdata <- new("topGOdata",ontology = "BP",description = "Neuron Control + Patient Rescue Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(C_PR_GOdata, algorithm = "classic", statistic = "fisher");
C_PRallRes <- GenTable(C_PR_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## P_CKI genes
output=unique(codgene$entrezgene_id[which(codgene$hgnc_symbol %in% DEG.P_CKI$gene)]);
myInterestingGenes = output;
geneList <- factor(as.integer(geneNames %in% myInterestingGenes));
names(geneList) <- geneNames;
P_CKI_GOdata <- new("topGOdata",ontology = "BP",description = "Neuron Patient + Control KI Rescue Up Regulated Genes",allGenes = geneList,nodeSize = 10,annotationFun = annFUN.org,mapping = "org.Hs.eg.db");
resultFisher <- runTest(P_CKI_GOdata, algorithm = "classic", statistic = "fisher");
P_CKIallRes <- GenTable(P_CKI_GOdata, classicFisher = resultFisher,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 25);
## save resutls
C_PRallRes$Condition2 = "C_PR";
P_CKIallRes$Condition2 = "P_CKI";
output = rbind(C_PRallRes,P_CKIallRes);
write.csv(output, file = "Top_GO_Neuron_results_Condition2.csv",quote=F,row.names=F);
## plot
pdf(file="TopGO_Neuron_DEGene_Condition2.pdf",height=8,width=6);
par(oma=c(2,10,1,1));
    bp = barplot(-log10(as.numeric(C_PRallRes$classicFisher)),main="Neuron Control + Patient Rescue Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=C_PRallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

    bp = barplot(-log10(as.numeric(P_CKIallRes$classicFisher)),main="Neuron Patient + Control KI Rescue Up Regulated Genes",horiz=TRUE,yaxt='n',yaxt='n',col="blue",xlab="-log10(P-value)",cex.main=0.5);
    axis(2,at=bp,labels=P_CKIallRes$Term,tick=FALSE,las=2,cex.axis=0.5);
    abline(v=-log(0.01),col="red",lwd=2,lty=1);

dev.off();

