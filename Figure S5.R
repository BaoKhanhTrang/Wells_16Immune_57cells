library("PRROC")
library("precrec")
library(ggrepel)

allv2g=readRDS("allv2g.rds")

burren=read.table("Burren.v2g_v2g.txt",header=T)
javrie=read.table("jav.v2g_v2g.txt",header=T)
yang=NULL
for(t in c("T0" ,"T20", "T1H", "T4H" ,"T24H")){
  y=read.table(paste0("Yang2020_",t,".v2g_v2g.txt"),header=T)
  y$study = paste0("Yang2020_",t)
  y=y[-grep("close",y$type),]
  yang = rbind(yang,y)
}

# ======= Figure S5 A ======= 
VennDiagram::venn.diagram(x = list(" "=unique(allv2g$gene_name),
                                   " " = unique(javrie$gene_name)),
                          filename="comparewithJavierre2016_pchic_genes.jpg",
                          height = 800,width = 800)

VennDiagram::venn.diagram(x = list(" "=unique(allv2g$gene_name[allv2g$system=="immune"]),
                                   " "=unique(burren$gene_name)),
                          filename="comparewithBurren2017_phic_genes.jpg",
                          height = 800,width = 800)

VennDiagram::venn.diagram(x = list(" "=unique(allv2g$gene_name[allv2g$system=="immune"]),
                                   " "= unique(yang$gene_name)),
                          filename="comparewithYang2020_genes.jpg",
                          height = 800,width = 800)
        
# ======= Figure S5 B =======
hiei = c("IL2RG","JAK3","IL7R","PTPRC","CD3D","CD3E","CD3Z","CORO1A","LAT","LCP2","RAG1","RAG2","DCLRE1C","PRKDC","NHEJ1","LIG4","ADA","AK2","RAC2","CD40","CD40LG","ICOS","ICOSL","CD3G","CD8A","ZAP70","TAP1","TAP2","TAPBP","B2M","CIITA","RFXANK","RFX5","RFXAP","DOCK8","DOCK2","POLD1","POLD2","RHOH","STK4","TRAC","LCK","ITK","MALT1","CARD11","BCL10","IL21","IL21R","TNFRSF4","IKBKB","MAP3K14","RELB","RELA","MSN","TFRC","REL","FCHO1","PAX1","ITPKB","SASH3","MAN2B2","COPG1","IKZF2","CHUK","WAS","WIPF1","ARPC1B","ATM","NBS1","BLM","DNMT3B","ZBTB24","CDCA7","HELLS","PMS2","RNF168","MCM4","POLA1","POLE1","POLE2","LIG1","NSMCE3","ERCC6L2","GINS1","MCM10","TBX1","CHD7","SEMA3E","FOXN1","RMRP","SMARCAL1","MYSM1","RNU4ATAC","EXTL3","STAT3","IL6R","IL6ST","ZNF341","ERBB2IP","TGFBR1","TGFBR2","SPINK5","PGM3","TCN2","SLC46A1","MTHFD1","IKBKG","NFKBIA","ORAI","STIM","CRACR2A","PNP","TTC7A","TTC37","SKIV2L","SP110","BCL11B","EPG5","RBCK1","RNF31","CCBE1","FAT4","NFE2L2","STAT5B","KMT2D","KDM6A","KMT2A","DIAPH1","IKZF3","CD28","BTK","IGHM","IHLL1","CK79A","CD79B","BLNK","PIK3CD","PIK3R1","TCF3","SLC39A7","TOP2B","FNIP1","SPI1","PTEN","CD19","CD81","CD20","CD21","TNFRSF13B","TNFRSF13C","TNFSF12","TRNT1","NFKB1","NFKB2","IKZF1","IRF2BP2","ATP6AP1","ARHGEF1","SH3KBP1","SEC61A1","RAC2","MOGS","PIK3CG","POU2AF1","AICDA","UNG","INO80","MSH6","CTNNBL1","TNFSF13","IGKC","PRF1","UNC13D","STX11","STXBP2","FAAP24","SLC7A7","RHOG","LYST","RAB27A","AP3B1","AP3D1","CEBPE","FOXP3","IL2RA","IL2RB","CTLA4","LRBA","DEF6","STAT3","BACH2","FERMT1","AIRE","ITCH","TPP2","JAK1","PEPD","SOCS1","PDCD1","IL10","IL10RA","IL10RB","NFAT5","TFGB1","RIPK1","ELF4","TNFRSF6","TFNSF6","CASP10","CASP8","FADD","SH2D1A","XIAP","CD27","CD70","CTPS1","TNFRSF9","RASGRP1","CARMIL2","MAGT1","PRKCD","TET2","ELANE","GFI1","HAX1","G6PC3","VPS45","G6PT1","WAS","LAMTOR2","TAZ","VPS13B","USB1","JAGN1","CLPB","CSF3R","SMARCD2","CEBPE","SBDS","DNAJC21","EFL1","HYOU1","SRP54","CXCR2","ITGB2","SLC35C1","FERMT3","RAC2","ACTB","FPR1","CTSC","WDR1","CFTR","MKL1","CYBB","CYBA","CYBC1","NCF1","NCF2","NCF4","GATA2","CSF2RA","CSFR2B","IL12RB1","IL12B","IL12RB2","IL23R","IFNGR1","IFNGR2","STAT1","IRF8","SPPL2A","TYK2","ISG15","RORC","TBX21","IFNG","TMC6","TMC8","CIB1","CXCR4","STAT1","STAT2","IRF9","IRF7","IFNAR1","IFNAR2","FCGR3A","IFIH1","NOS2","ZNFX1","POLR3A","POLR3C","POLR3F","TLR3","UNC93B1","TRAF3","TIRAP","TLR3","IRF3","DBR1","SNORA31","ATC4","MAP1LC3B2","CARD9","IL17RA","IL17RC","IL17F","STAT1","TRAF3IP2","MAPK8","IRAK4","MYD88","IRAK1","TLR7","TLR8","RPSA","HMOX","APOL1","NBAS","RANBP2","CLCN7","SNX10","OSTM1","PLEKHM1","TCIRG1","TNFRSF11A","TNFSF11","NCSTN","PSEN","PSENEN","IRF4","IL18BP","TMEM173","ADA2","TREX1","RNASEH2B","RNASEH2C","RNASEH2A","SAMHD1","ADAR1","IFIH1","DNASE2","LSM11","RNU7-1","DNASE1L3","ACP5","POLA1","USP18","OAS1","CDC42","STAT2","ATAD3A","MEFV","MVK","NLRP3","NLRP12","NLRC4","PLCG2","NLRP1","RIPK1","TNFRSF1A","PSTPIP1","NOD2","ADAM17","LPIN2","IL1RN","IL36RN","SLC29A3","CARD14","SHEBP2","PSMB8","PSMG2","COPA","OTULIN","TNFAIP3","AP1S3","ALPI","TRIM22","HAVCR2","C2ORF69","NCKAP1L","SYK","HCK","PSMB9","IKBKG","TBK1","C1QA","C1QB","C1QC","C1R","C1S","C4A","C4B","C2","C3","C5","C6","C7","C8A","C8G","C8B","C9","MASP2","FCN3","SERPING1","CFB","CFD","CFP","CFI","CFH","CFHR1","CFHR2","CFHR3","CFHR4","CFHR5","THBD","CD46","CD59","CD55","FANCA","FANCB","FANCC","BRCA2","FANCD2","FANCE","FANCF","XRCC9","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","ERCC4","RAD51","BRCA1","UBE2T","XRCC2","MAD2L2","RFWD3","SAMD9","SAMD9L","DKC1","TERC","TERT","TINF2","RTEL1","ACD","NOLA3","NOLA2","WRAP53","PARN","SRP72","ERCC6L2","TP53","STN1","CTC1","MECOM","TNFRSF6","KRAS","NRAS","UBA1")
pp=data.frame(Name=unique(allv2g$Cellline[which(allv2g$system=="immune")]),Precision=0,Recall=0)
for(i in 1:nrow(pp)){
  tmp = unique(allv2g$gene_name[allv2g$Cellline == pp$Name[i]])
  pp$Precision[i] = length(intersect(tmp,hiei))/length(tmp)
  pp$Recall[i] = length(intersect(tmp,hiei))/length(hiei)
  pp$numgene=length(tmp)
}
pp$mystu=TRUE
ppp=data.frame(Name=unique(allv2g$system),Precision=0,Recall=0)
for(i in 1:nrow(ppp)){
  tmp = unique(allv2g$gene_name[allv2g$system == ppp$Name[i]])
  ppp$Precision[i] = length(intersect(tmp,hiei))/length(tmp)
  ppp$Recall[i] = length(intersect(tmp,hiei))/length(hiei)
  ppp$numgene=length(tmp)
}
ppp$mystu=TRUE
pp = rbind(pp,ppp)
rm(ppp,tmp,i)

tmp = pp[1,]
tmp$Name="Javierre2016_pchic"
tmp$Precision = length(intersect(unique(javrie$gene_name),hiei))/length(unique(javrie$gene_name))
tmp$Recall = length(intersect(unique(javrie$gene_name),hiei))/length(hiei)
tmp$mystu=F
tmp$numgene = length(unique(javrie$gene_name))
pp = rbind(pp,tmp)

tmp = pp[1,]
tmp$Name="Burren2017_phic"
tmp$Precision = length(intersect(unique(burren$gene_name),hiei))/length(unique(burren$gene_name))
tmp$Recall = length(intersect(unique(burren$gene_name),hiei))/length(hiei)
tmp$mystu=F
tmp$numgene = length(unique(burren$gene_name))
pp = rbind(pp,tmp)

tmp = pp[1,]
tmp$Name="Yang2020_phic"
tmp$Precision = length(intersect(unique(yang$gene_name),hiei))/length(unique(yang$gene_name))
tmp$Recall = length(intersect(unique(yang$gene_name),hiei))/length(hiei)
tmp$mystu=F
tmp$numgene = length(unique(yang$gene_name))
pp = rbind(pp,tmp)


pp$Name[which(pp$Name=="immune")] = "Total immune cells"
pp$color = pp$Name
pp$color[which(pp$Name %in% allv2g$Cellline)] = "Total immune cells"
pp$color[grep("Yang",pp$Name)] = "Yang2020"

pdf('Precision_recell.pdf',height = 5,width = 7)
ggplot(pp,aes(Recall,Precision,label=Name,color=color))+
  geom_point(aes(shape=mystu,size=numgene))+
  geom_label_repel()+
  scale_size_binned(name="Number of \nimplicated genes",breaks = seq(500,2500,by=500) )+
  scale_color_discrete(guide = "none")+
  scale_shape_manual(values = c(15,16) , guide = "none")+
  guides(size = guide_bins(show.limits = TRUE))+
  theme_bw()
dev.off()