##load package
library(VennDiagram);
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2");
library(pheatmap);
options(stringsAsFactors=F);
## cell types
celltypes = c("Astrocyte", "DividingRadialGlia", "DLNeuron", "IntermediateProgenitor", "Neuron", "NPCtoNcells", "OuterRadialGlia", "RadialGlia", "ULNeuron");
## plot
pdf("GraphDEGoverlap.pdf");
for (i in 1:length(celltypes)) {
    ## DEG and conditions
    DEG1 = read.csv(paste0(celltypes[i],"_DEG_Condition1.csv"));
    DEG2 = read.csv(paste0(celltypes[i],"_DEG_Condition2.csv"));
    DEG3_1 = read.csv(paste0(celltypes[i],"_DEG_Condition3_C_CKI.csv"));
    DEG3_2 = read.csv(paste0(celltypes[i],"_DEG_Condition3_P_PR.csv"));
    ## C gene list
    DEG.C = DEG1[which(DEG1$cluster=="C"),];
    DEG.C = DEG.C[which(DEG.C$p_val_adj < 0.05 & DEG.C$avg_log2FC > 0.3),];
    ## P gene list
    DEG.P = DEG1[which(DEG1$cluster=="P"),];
    DEG.P = DEG.P[which(DEG.P$p_val_adj < 0.05 & DEG.P$avg_log2FC > 0.3),]; 
    ## C_PR gene list
    DEG.C_PR = DEG2[which(DEG2$cluster=="C_PR"),];
    DEG.C_PR = DEG.C_PR[which(DEG.C_PR$p_val_adj < 0.05 & DEG.C_PR$avg_log2FC > 0.3),];
    ## P_CKI gene list
    DEG.P_CKI = DEG2[which(DEG2$cluster=="P_CKI"),];
    DEG.P_CKI = DEG.P_CKI[which(DEG.P_CKI$p_val_adj < 0.05 & DEG.P_CKI$avg_log2FC > 0.3),];
    ## C > CKI gene list
    DEG.CKI_C.C = DEG3_1[which(DEG3_1$cluster=="C"),];
    DEG.CKI_C.C = DEG.CKI_C.C[which(DEG.CKI_C.C$p_val_adj < 0.05 & DEG.CKI_C.C$avg_log2FC > 0.3),];
    ## CKI >C gene list
    DEG.CKI_C.CKI = DEG3_1[which(DEG3_1$cluster=="CKI"),];
    DEG.CKI_C.CKI = DEG.CKI_C.CKI[which(DEG.CKI_C.CKI$p_val_adj < 0.05 & DEG.CKI_C.CKI$avg_log2FC > 0.3),];   
    ## P > PR gene list
    DEG.P_PR.P = DEG3_2[which(DEG3_2$cluster=="P"),];
    DEG.P_PR.P = DEG.P_PR.P[which(DEG.P_PR.P$p_val_adj < 0.05 & DEG.P_PR.P$avg_log2FC > 0.3),];
    ## PR > P gene list
    DEG.P_PR.PR = DEG3_2[which(DEG3_2$cluster=="PR"),];
    DEG.P_PR.PR = DEG.P_PR.PR[which(DEG.P_PR.PR$p_val_adj < 0.05 & DEG.P_PR.PR$avg_log2FC > 0.3),];
    ## plot
    DEG.list = list(DEG.C$gene, DEG.C_PR$gene, DEG.CKI_C.C$gene, DEG.P_PR.PR$gene);
    DEG.names = c("C > P" , "C_PR > P_CKI" , "C > CKI", "PR > P");
    venn.diagram(
        x = DEG.list,
        category.names = DEG.names,
        filename = paste0(celltypes[i],'_control_like_DEG_venn.png'),
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 960 , 
        width = 960 , 
        resolution = 300,

        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol);
    ## heatmap
    output = matrix(NA,4,4);
    colnames(output) = DEG.names;
    rownames(output) = DEG.names;
    for (a in 1:4) {
        for (b in 1:4) {
        output[a,b] = length(intersect(DEG.list[[a]],DEG.list[[b]]))/length(unique(c(DEG.list[[a]],DEG.list[[b]])));
        }
    }
    pheatmap(output,display_numbers=T,main=paste0(celltypes[i]," control-like Jaccard index"));

    ## plot
    DEG.list = list(DEG.P$gene, DEG.P_CKI$gene, DEG.CKI_C.CKI$gene, DEG.P_PR.P$gene);
    DEG.names = c("P > C" , "P_CKI > C_PR" , "CKI > C", "P > PR");
    venn.diagram(
        x = DEG.list,
        category.names = DEG.names,
        filename = paste0(celltypes[i],'_patient_like_DEG_venn.png'),
        output=TRUE,

        # Output features
        imagetype="png" ,
        height = 960 ,
        width = 960 ,
        resolution = 300,

        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol);
    ## heatmap
    output = matrix(NA,4,4);
    colnames(output) = DEG.names;
    rownames(output) = DEG.names;
    for (a in 1:4) {
        for (b in 1:4) {
        output[a,b] = length(intersect(DEG.list[[a]],DEG.list[[b]]))/length(unique(c(DEG.list[[a]],DEG.list[[b]])));
        }
    }
    pheatmap(output,display_numbers=T,main=paste0(celltypes[i]," patient-like Jaccard index"));

}

dev.off();

