pdf('Shared_GO_plots.pdf',width = 12, height = 6)
for (i in 1:length(files)){
  if (grepl("_up", files[i])) {
    par(mar=c(3,25,3,2))
    thiscol = 'red'
    godata = read.csv(paste0('~/Desktop/Zhang_MS2/GO_shared/',files[i]))
    vmax = max(-log10(godata$Adjusted.P.value[20:1]),-log10(0.05))
    barplot(-log10(godata$Adjusted.P.value[20:1]), names.arg = godata$Term[20:1], horiz = T, las=2, col = thiscol, xlim = c(0, vmax*1.5), main=files[i], cex.names = 0.7)
    abline(v=-log10(0.05))
    } else {
      par(mar=c(3,25,3,2))
      thiscol = 'blue'
      godata = read.csv(paste0('~/Desktop/Zhang_MS2/GO_shared/',files[i]))
      vmax = max(-log10(godata$Adjusted.P.value[20:1]),-log10(0.05))
      barplot(-log10(godata$Adjusted.P.value[20:1]), names.arg = godata$Term[20:1], horiz = T, las=2, col = thiscol, xlim = c(0, vmax*1.5), main=files[i], cex.names = 0.7)
      abline(v=-log10(0.05))
      }
}
dev.off()


ms2_PvsC = read.csv("/Users/liangdan/Library/CloudStorage/Box-Box/PIDD1\ Data/CeZhang_MS1_2/Otherfiles/Patient_vs_Control.csv")
ms2_PvsC$log2MeanRatio_C_YB = -ms2_PvsC$log2MeanRatio_C_P
mtorc1_genes = read.csv("/Users/liangdan/Library/CloudStorage/Box-Box/PIDD1\ Data/CeZhang_MS1_2/Otherfiles/mTORC1_Signaling_from_MSigDB_Hallmark_2020.csv")
# P vs C
ms2_PvsC$minuslog_Anova_p = -log10(ms2_PvsC$anova_p)
ms2_PvsC$GO = NA
ms2_PvsC$GO[which(ms2_PvsC$Gene.Names %in% mtorc1_genes$mTORC1_Signaling & ms2_PvsC$fdr < 0.05 )] = 'mTORC1 Signaling'
ms2_PvsC$GO[which(ms2_PvsC$Gene.Names %in% c('CYB5A','NDUFB8','NDUFB5','MRPS22','MFN2','ATP6V1H','NDUFS2','ACADSB','ATP6V1C1'))] = 'Oxidative Phosphorylation'
ms2_PvsC$GO[which(ms2_PvsC$Gene.Names %in% c('PTCD3','GPAT4','MCCC1','ATL2','PREB','GHITM','SLC1A5','HIBCH'))] = 'Adipogenesis'
ms2_PvsC$GO[which(ms2_PvsC$Gene.Names %in% c('LMAN1','DST','AP1G1','COPB1','GBF1','PPT1','AP3B1','AP2B1','ARFGAP3','CTSC','COPE','SEC31A'))] = 'Protein Secretion'
ms2_PvsC$GO[which(ms2_PvsC$Gene.Names %in% c('PDHX','GPI','NQO2','ECHS1','GPX4','ABCB7','IDH3G','TIMM9','IMMT','TIMM10','TIMM50','RHOT2','SUCLA2','ACADM','ACO2'))] = 'Oxidative Phosphorylation'
ms2_PvsC$GO[which(ms2_PvsC$Gene.Names %in% c('SMARCC1','NCBP1','EIF1AX','CAD','MRPS18B','PA2G4','AIMP2','YWHAQ','HPRT1','PABPC1','KPNA2','ACP1','EIF3B'))] = 'Myc Targets V1'


library(RColorBrewer)
library(ggplot2)
myColors <- brewer.pal(7,"Set1")
colScale <- scale_colour_manual(name = "condition",values = myColors)
myalpha = rep(1,dim(ms2_PvsC)[1])
myalpha[which(is.na(ms2_PvsC$GO))] = 0.2
ggplot(ms2_PvsC, aes(x=log2MeanRatio_C_YB, y=minuslog_Anova_p, color = GO)) + geom_point(alpha=myalpha) + colScale + theme_bw() + xlim(-12, 6) + ylim(0, 8)


ms2_KIvsC = read.csv("/Users/liangdan/Library/CloudStorage/Box-Box/PIDD1\ Data/CeZhang_MS1_2/Otherfiles/Ki_vs_Control.csv")
ms2_KIvsC$log2MeanRatio_KI_YB = -ms2_KIvsC$log2MeanRatio_C_Ki
# KI vs C
ms2_KIvsC$minuslog_Anova_p = -log10(ms2_KIvsC$anova_p)
ms2_KIvsC$GO = NA
ms2_KIvsC$GO[which(ms2_KIvsC$Gene.Names %in% mtorc1_genes$mTORC1_Signaling & ms2_KIvsC$fdr < 0.05 )] = 'mTORC1 Signaling'
ms2_KIvsC$GO[which(ms2_KIvsC$Gene.Names %in% c('DDX18','TCOF1','PES1','DCTPP1','SRM'))] = 'Myc Targets V2'
ms2_KIvsC$GO[which(ms2_KIvsC$Gene.Names %in% c('PDHX','AFG3L2','RHOT2','NQO2','GPX4','NDUFB3','NDUFA2','OGDH','NDUFS3','ACO2','NDUFV1','OPA1','ABCB7','MGST3','NDUFS2','MFN2','IDH3A'))] = 'Oxidative Phosphorylation'
ms2_KIvsC$GO[which(ms2_KIvsC$Gene.Names %in% c('NOP56','SMARCC1','MCM7','EIF1AX','EIF2S1','TUFM','AIMP2','POLD2','PSMD1','SNRPA1','HPRT1','ACP1','MCM2','EIF3B'))] = 'Myc Targets V1'
ms2_KIvsC$GO[which(ms2_KIvsC$Gene.Names %in% c('GPX4','POLR2A','RFC2','MPC2','POLR2E','AK3','HPRT1','GTF2F1','SUPT5H','ITPA'))] = 'DNA Repair'
ms2_KIvsC$GO[which(ms2_KIvsC$Gene.Names %in% c('TOP2A','SMARCC1','ACSL1','TSPO','ACSL4','DHCR24','PEX14'))] = 'Pperoxisome'

library(RColorBrewer)
myColors <- brewer.pal(7,"Set1")
colScale <- scale_colour_manual(name = "condition",values = myColors)
myalpha = rep(1,dim(ms2_KIvsC)[1])
myalpha[which(is.na(ms2_KIvsC$GO))] = 0.2
ggplot(ms2_KIvsC, aes(x=log2MeanRatio_KI_YB, y=minuslog_Anova_p, color = GO)) + geom_point(alpha=myalpha) + colScale + theme_bw() + xlim(-12, 6) + ylim(0, 8)


ms2_MDSvsC = read.csv("/Users/liangdan/Library/CloudStorage/Box-Box/PIDD1\ Data/CeZhang_MS1_2/Otherfiles/MDS_vs_Control.csv")
ms2_MDSvsC$log2MeanRatio_MDS_YB = - ms2_MDSvsC$log2MeanRatio_C_MDS
# MDS vs C
ms2_MDSvsC$minuslog_Anova_p = -log10(ms2_MDSvsC$anova_p)
ms2_MDSvsC$GO = NA
ms2_MDSvsC$GO[which(ms2_MDSvsC$Gene.Names %in% mtorc1_genes$mTORC1_Signaling & ms2_MDSvsC$fdr < 0.05 )] = 'mTORC1 Signaling'
ms2_MDSvsC$GO[which(ms2_MDSvsC$Gene.Names %in% c('SQLE','GSTM2','MVK','ACSS2','HMGCS1','SCD','CYP51A1'))] = 'Cholesterol Homeostasis'
ms2_MDSvsC$GO[which(ms2_MDSvsC$Gene.Names %in% c('DDX18','FBL','GLO1','PGK1','EIF4H','PRPF31','SNRPB2','HPRT1','EIF2S2','YWHAE','HDAC2','MCM7','PPM1G','SYNCRIP','SNRPD2','YWHAQ','POLD2','RPS3','RPL14','CYC1',
                                                 'KPNA2','EIF4E','NOP56','SMARCC1','RRM1','NCBP2','RPS5','GOT2','CAD','TUFM','AIMP2','EIF3J','UBA2','MCM4','SNRPA1','MCM5','PABPC1','MCM6','KPNB1','MCM2','EIF3B'))] = 'Myc Targets V1'
ms2_MDSvsC$GO[which(ms2_MDSvsC$Gene.Names %in% c('NDUFB6','COX15','ECI1','TIMM9','NDUFB4','NDUFB1','ACAT1','RHOT2','ACADM','CYC1','NDUFV2','NDUFV1','ATP6V1F','ATP5ME','NDUFA9',
                                                 'IDH2','GOT2','NDUFA2','IMMT','PDP1','UQCRQ','NDUFS6','UQCRC1','NDUFS1','UQCRC2','SLC25A12'))] = 'Oxidative Phosphorylation'
ms2_MDSvsC$GO[which(ms2_MDSvsC$Gene.Names %in% c('TOP2A','NOP56','MCM7','PRKDC','RFC2','RPA1','CCP110','NUP153','CTCF','SSRP1','TUBG1','SYNCRIP','NASP','POLD2',
                                                 'MCM3','MCM4','MCM5','MCM6','KPNA2','LYAR','LBR','TMPO','SNRPB','MCM2'))] = 'E2F Targets'
ms2_MDSvsC$GO[which(ms2_MDSvsC$Gene.Names %in% c('TOP2A','UPF1','TLE3','SMARCC1','DR1','NUMA1','CTCF','PML','SYNCRIP','MEIS1','NASP','MCM3','MCM5','NUP98',
                                                 'MCM6','KPNA2','LBR','KPNB1','TMPO','PAFAH1B1','MCM2'))] = 'G2-M Checkpoint'

ms2_MDSvsC$GO = factor(ms2_MDSvsC$GO, levels = c("Cholesterol Homeostasis","mTORC1 Signaling","E2F Targets","G2-M Checkpoint","Myc Targets V1","Oxidative Phosphorylation"))

library(RColorBrewer)
myColors <- brewer.pal(7,"Set1")
colScale <- scale_colour_manual(name = "condition",values = myColors)
myalpha = rep(1,dim(ms2_MDSvsC)[1])
myalpha[which(is.na(ms2_MDSvsC$GO))] = 0.2
ggplot(ms2_MDSvsC, aes(x=log2MeanRatio_MDS_YB, y=minuslog_Anova_p, color = GO)) + geom_point(alpha=myalpha) + colScale + theme_bw() + xlim(-12, 6) + ylim(0, 8)



draw.triple.venn(252, 139, 258, 58, 49,  61, 26, category = c("P_up","KI_up","MDS_up"))
draw.triple.venn(401, 338, 553, 142, 136, 131, 70, category = c("P_down","KI_down","MDS_down"))





