library(ggplot2)
library(ggsci)
library(egg)
library(scales)
library(reshape2)
library(patchwork)
library(paletteer)
library(colorblindr)
options(scipen = 100)

trait_ord = readRDS("trait_ord.rds")
allstat=readRDS("allstat.rds")
g1=ggplot(allstat)+geom_bar(aes(atac_peaks_N,Name,fill=material),stat = "identity",width = .8)+
  geom_bar(aes(ocr_N,Name,fill=Platform),stat = "identity",width = .7)+
  facet_grid(rows=vars(system),scale = "free",space = "free",drop = TRUE)+
  scale_fill_manual(name="",values =c("Bulk ATAC-seq"="#0072B2FF","Single-cell ATAC-seq"="#C0392BFF","Hi-C"="#009E73FF","Capture-C"="#E69F00FF"))+
  theme_bw()+
  scale_y_discrete(name="Cell types",breaks = levels(allstat$Name))+
  scale_x_continuous(expand = c(0,0),name = "#OCRs x1K",
                     breaks = seq(0,300000,by=100000),
                     labels = label_number(scale = 1e-3),
                     limits = c(0,300000))+
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_blank(),
        legend.position = c(0.75,0.9),legend.background = element_blank())+
  guides(fill=guide_legend(nrow=2,byrow=F))

trait=read.delim("trait_gwas.txt",header=T)
nar=NULL
scales <- list()
for(i in 1:nrow(trait)){
  tmp = read.delim(paste0("narrow/summary.lia_",trait$Name[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(trait$Name[i])
  scales[[i]] = scale_x_continuous(breaks = c(1,max(ceiling(tmp$Enrichment))),limits = c(min(floor(tmp$Enrichment))-1,max(ceiling(tmp$Enrichment))+1))

  tmp2 = read.delim(paste0("narrow/h2.lia_",trait$Name[i],".txt"),header=F,sep=":",fill = TRUE,strip.white = T)
  colnames(tmp2)=c("Celltype","h2")
  tmp=merge(tmp,tmp2,by="Celltype")
  nar=rbind(nar,tmp)
}
rm(i,tmp,tmp2)
nar = merge(nar,allstat[,c('Celltype','system','Platform','Name')],by='Celltype')
nar$material="Bulk ATAC-seq"
nar$material[grep("klaus_pancreas_cells",nar$Celltype)]="Single-cell ATAC-seq"
nar$Name = factor(nar$Name,levels = levels(allstat$Name))

sd=read.table("narrow_sd.txt",header=F,sep=":",strip.white = T)
colnames(sd)=c("Celltype","sdc")
nar=merge(nar,sd,by="Celltype")
nar=merge(nar,trait[,c("Name","RemainSNPs")],by.x="trait",by.y="Name")
nar$trait = factor(nar$trait,levels = trait_ord)

# ======= Figure 1 A ======= 
nar$tau = nar$RemainSNPs*nar$sdc*nar$Coefficient/nar$h2
nar$tau_se = nar$RemainSNPs*nar$sdc*nar$Coefficient_std_error/nar$h2
nar$tau_p <- 2 * (1 - pnorm(abs(nar$Coefficient_z.score)))
nar$sig = nar$tau_p < 0.05 
NAR=NULL
for(i in 1:nrow(trait)){
  tmp=nar[nar$trait==trait$Name[i],]
  tmp$FDR = p.adjust(tmp$tau_p,"fdr")
  NAR=rbind(NAR,tmp)
}
NAR$FDRsig = NAR$FDR < 0.05 
NAR$sig = ifelse(NAR$tau_p < 0.05,"p-value < 0.05",NA)
NAR$fdr = p.adjust(NAR$tau_p,"fdr")
NAR$fdrsig = NAR$fdr < 0.05 

minbin = round(min(nar$Prop._SNPs),2)
maxbin = round(max(nar$Prop._SNPs),2)
maxenr = max(ceiling(nar$tau+nar$tau_se)) 
minenr = min(floor(min(nar$tau-nar$tau_se))) 
minp = floor(min(-log10(nar$tau_p)))
maxp = ceiling(max(-log10(nar$tau_p)))

g3=ggplot(NAR)+
  geom_point(aes(trait,Name,fill=tau,size=-log10(tau_p),colour = factor(sig)),alpha=0.7,stat="identity",shape=21)+
  theme_bw()+
  scale_y_discrete(breaks = levels(allstat$Name))+
  geom_point(data=NAR[which(NAR$fdrsig),],aes(trait,Name),color="#019875FF",size=1,shape=8)+
  scale_fill_gradient2(name=paste0("\U1D749","*"),high = "#96281BFF",mid = "white",low = "#0072B2FF",limits = c(minenr,maxenr),guide=guide_colorbar(show.limits = TRUE))+
  scale_color_manual(na.value = "white",values = c("p-value < 0.05"="#CC79A7FF"),na.translate = F,name="",guide = 'none')+
  xlab(paste0("Conditional annotation effect size (", "\U1D749","*) in putative cREs of each cell type"))+
  facet_grid(rows=vars(system),space = "free_y",scales = "free")+
  scale_size_binned(name = "-log10 P-value",breaks = seq(0,maxp,by=2),limits = c(0,maxp))+
  guides(size = guide_bins(show.limits = TRUE))+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),plot.caption = element_text(hjust = 0.5,vjust = 0),
        axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        legend.title.align=0.5)
cairo_pdf("LDSC_taustar_filtered.pdf",height = 8,width = 8)
g1 + g3 + plot_layout(widths = c(1,6),guides = "collect") & theme(legend.position = 'bottom')
dev.off()

# ======= Figure 1 B ======= 
a = all[which(all$Celltype=="CD1cPluscDC"),]
span=max(a$Enrichment+a$Enrichment_std_error) - min(a$Enrichment-a$Enrichment_std_error)
g5=ggplot(a,aes(x=trait,y=Enrichment,color=cate,group=cate))+
    geom_errorbar(aes(ymax=Enrichment+Enrichment_std_error,ymin=Enrichment-Enrichment_std_error),position = position_dodge(width = .75))+
    geom_point(aes(size=`Prop._SNPs`),position = position_dodge(width=.75))+
    geom_point(aes(shape=sig),color="white",position = position_dodge(width=.75))+
    geom_point(aes(shape=tsig,y=Enrichment+Enrichment_std_error+span/20),color="red",position = position_dodge(width=.75))+
    scale_size_binned(name = "%SNP",breaks = pretty_breaks(n=4),
                    labels = function(x) paste0(as.character(round(x,3)*100)),limits = c(0,max(all$Prop._SNPs)))+
    scale_shape_manual(values = c(NA,"p-value < 0.05"=8))+
    theme_bw()+labs(x="")+
    guides(size = guide_bins(show.limits = TRUE) )

  
span=max(a$tau+a$tau_se) - min(a$tau-a$tau_se)
g6=ggplot(a,aes(x=trait,y=tau,color=cate,group=cate))+
    geom_errorbar(aes(ymax=tau+tau_se,ymin=tau-tau_se),position = position_dodge(width = .75))+
    geom_point(aes(size=`Prop._SNPs`),position = position_dodge(width=.75))+
    geom_point(aes(shape=FDRsig),color="white",position = position_dodge(width=.75))+
    geom_point(aes(shape=tsig_tau,y=tau+tau_se+span/20),color="red",position = position_dodge(width=.75))+
    scale_size_binned(name = "%SNP",breaks = pretty_breaks(n=4),
                    labels = function(x) paste0(as.character(round(x,3)*100)),limits = c(0,max(all$Prop._SNPs)))+
    scale_shape_manual(values = c(NA,"p-value < 0.05"=8,"perTrait FDR < 0.05"=8))+
    theme_bw()+labs(x="")+
    guides(size = guide_bins(show.limits = TRUE) )
pdf("CD1cPluscDC.LDSC panels 5 categories",width = 11,height = 4.5)
print(g5+g6+plot_layout(nrow = 2,guides = "collect") )
dev.off()

# ======= Figure S1 A ======= 
nar$sig = nar$Enrichment_p < 0.05 & nar$Enrichment>0
minbin = round(min(nar$Prop._SNPs),2)
maxbin = round(max(nar$Prop._SNPs),2)
maxenr = max(nar$Enrichment[which(nar$system=="immune")])
minenr = min(nar$Enrichment[which(nar$system=="immune")])
minp = floor(min(-log10(nar$Enrichment_p)))
maxp = ceiling(max(-log10(nar$Enrichment_p)))


g2=ggplot(nar)+geom_vline(xintercept = 1, linetype='dashed', linewidth=.5)+
  geom_point(aes(Enrichment,Name,size=`Prop._SNPs`,color=-log10(Enrichment_p)),alpha=0.7)+
  geom_errorbarh(aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#E69F00FF",height=.5,linewidth=.4,position = "dodge")+
  geom_errorbarh(data=nar[nar$Platform=="Hi-C",],aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#009E73FF",height=.5,linewidth=.4,position = "dodge")+
  theme_bw()+
  scale_y_discrete(breaks = levels(allstat$Name))+
  geom_point(data=nar[which(nar$sig),],aes(Enrichment,Name),color="white",size=1,shape=8)+
  scale_size_binned(name = "Proportion\nof SNP(%)",breaks = seq(0,maxbin,by=0.005),
                    labels = function(x) paste0(as.character(round(x,3)*100)),limits = c(0,maxbin))+
  coord_cartesian(xlim=c(minenr,maxenr))+
  scale_color_gradient(name="-log10 p-value",high = "#96281BFF",low = "#0072B2FF",limits = c(minp,maxp))+
  facet_grid(rows=vars(system),cols = vars(trait),space = "free_y",scales = "free")+
  xlab("Enrichment of heritability in putative cREs")+  #,n.breaks = 3, breaks = waiver())+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),plot.caption = element_text(hjust = 0.5,vjust = 0),
        axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust=0.5))+
  guides(size = guide_bins(show.limits = TRUE,title.position="top"),color=guide_colourbar(title.position="top") )

pdf("all traits panel - all cells _ lia_freexaxis.pdf",height = 8.5,width = 12)
g1 + g2 + plot_layout(widths = c(1,12),guides = "collect")  & theme(legend.position = 'bottom')
dev.off()

# ======= Figure S1 B ======= 
g4=ggplot(NAR)+
  geom_bar(aes(tau,Name,fill=-log10(tau_p),colour = factor(sig)),alpha=0.7,stat="identity")+
  geom_errorbarh(aes(y=Name,xmax=tau+tau_se,xmin=tau-tau_se),color="#E69F00FF",height=.5,linewidth=.4,position = "dodge")+
  geom_errorbarh(data=NAR[NAR$Platform=="Hi-C",],aes(y=Name,xmax=tau+tau_se,xmin=tau-tau_se),color="#009E73FF",height=.5,linewidth=.4,position = "dodge")+
  theme_bw()+
  scale_y_discrete(breaks = levels(allstat$Name))+
  geom_point(data=NAR[which(NAR$fdrsig),],aes(tau+tau_se*1.5*tau/abs(tau),Name),color="red",size=1,shape=8)+
  scale_fill_paletteer_c(name="-log10 p-value","ggthemes::Classic Gray",limits = c(minp,maxp))+
  scale_color_manual(na.value = "white",values = c("p-value < 0.05"="#5B3794"),na.translate = F,name="")+
  xlab(paste0("Conditional annotation effect size (", "\U1D749","*) in putative cREs of each cell type"))+
  facet_grid(cols=vars(trait),rows=vars(system),space = "free_y",scales = "free")+
  guides(size = guide_bins(show.limits = TRUE) )+
  ggh4x::facetted_pos_scales(x = scales)+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),plot.caption = element_text(hjust = 0.5,vjust = 0),
        axis.title.x.top = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        legend.title.align=0.5)+
  guides(size = guide_bins(show.limits = TRUE,title.position="top",),fill=guide_colourbar(title.position="top") )

cairo_pdf("LDSC_taustar_bar.pdf",height = 8,width = 7.5)
g1 + g4 + plot_layout(widths = c(1,12),guides = "collect") & theme(legend.position = 'bottom')
dev.off()
