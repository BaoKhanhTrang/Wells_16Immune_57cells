# ======= Figure 4 A created using pyGenomeTracks with correspoding ATAC-seq Bigwig files 
# and chromatin contacts BEDPE files of indivuidual cell types at FDFT1 locus =======

library(data.table)
library(ggtern) # v3.4.0 with ggplot2 v3.4.0 devtools::install_version("ggplot2","3.4.0")
library(RColorBrewer)
library(ggstatsplot)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(EnsDb.Hsapiens.v86)
library(eulerr)
library(dplyr)
library(nVennR)
library(broom)
library(tidyr) # for `unnest`
library(purrr)

match=read.table("matching_celltypes.txt",sep="\t",header=T)
trait_ord = readRDS("trait_ord.rds")
allv2g=readRDS("allv2g.rds")
allv2g = allv2g[which(allv2g$system=="immune"),]

# ======= Figure 4 B ======= 
fdft1=read.table("chr8.9921161.11972641_FDFT1/ENSG00000079459.coloc_results_all_summaryspread_condPP4.txt",header=T)
fdft1 = unique(merge(fdft1,match,by=c("source","filename")))

fdft1 = merge(fdft1,unique(allv2g[,c('Cellline','Celltype')]),by="Celltype")
fdft1$Cell = factor(fdft1$Cell,levels = c("NaiveB","B mem/GCB","plasma/PGF","NaiveT",
                                        "Stimulated T-cells","T follicular helper cells","Th1" ,"Th2","Th17","Treg" ,"Treg memory/activated" , 
                                        "Dendritic cells" ,"Monocytes","Monocytes-nonclassical" ,"Natural killer" ))
fdft1$Trait = factor(fdft1$Trait,levels = trait_ord)


pdf("chr8.9921161.11972641_FDFT1.coloc_perTrait.noImmuNexUT.SLE_RA_MS_T1D.pdf",width = 16,height = 8)
ggtern(fdft1[which(fdft1$Trait %in% c("SLE","MS","RA","T1D")),],aes(x=PP.H3.abf,y=PP.H4.abf,z=PP.H0.abf+PP.H1.abf+PP.H2.abf))+ 
  geom_point(aes(color=Trait,shape=Cellline),size=2) +
  facet_grid(cols=vars(Trait))+
  scale_color_manual(values=trait_colors[which(names(trait_colors) %in% c("SLE","MS","RA","T1D"))],guide="none")+
  theme_bw()+theme_hidetitles()+theme_showarrows()+theme_showgrid_major()+
  guides(shape=guide_legend(ncol=5,position = "bottom"))
dev.off()

# ======= Figure 4 C ======= 
trait="SLE"
trait_reg=read.table(paste0("chr8.9921161.11972641_FDFT1/",trait,".txt"),header=T)
trait_reg=trait_reg[,c('chromosome','position','BETA','P','SNP')]
colnames(trait_reg)=c("CHR","BP","BETA","P","SNP")
LD = read.table("chr8.9921161.11972641_FDFT1/chr8.9921161.11972641.ld2",header = T)
colnames(LD)=c("SNP_A","SNP_B","R2")
ld=merge(LD,trait_reg[,c("CHR","BP","SNP")],by.x="SNP_A",by="SNP")
colnames(ld)[4:5]=c("CHR_A","BP_A")
ld=merge(ld,trait_reg[,c("CHR","BP","SNP")],by.x="SNP_B",by="SNP")
colnames(ld)[6:7]=c("CHR_B","BP_B")

traitleads=intersect(unique(allv2g$index_rsid[which(allv2g$trait==trait)]),leads)
traitproxy = unique(allv2g$proxy_rsid[which(allv2g$index_rsid %in% traitleads & allv2g$gene_name=="FDFT1")])

filename="NaiveT"
b=read.table(paste0("chr8.9921161.11972641_FDFT1/",filename,".ENSG00000079459.txt"),sep="\t",header=F) # b-cell eQTL_cata
header=readLines(paste0(unique(match$source[which(match$filename == filename)]),".header"))
colnames(b)=header
b=b[,c('rsid','Pvalue','Beta')]
colnames(b)=c('SNP.Id',"P.Value",'NES')
b = b[order(b$SNP.Id,b$P.Value,decreasing = F),]
if(any(duplicated(b$SNP.Id))){b = b[-which(duplicated(b$SNP.Id)),]}

sle=trait_reg[which(trait_reg$SNP %in% intersect(trait_reg$SNP,b$SNP.Id)),] 
b=b[which(b$SNP.Id %in% intersect(trait_reg$SNP,b$SNP.Id)),] 
trait2 = paste0("FDFT1 in ",unique(match$Cell[which(match$filename==filename)])," (",unique(match$source[which(match$filename==filename)]),")")
trait2 = sub("/"," or ",trait2)

source("scripts/pairGWASplot.R")
p=pairGWASplot(GWAS.df = sle, eQTL.df = b, R2min=0.8,gbuild = "hg38",
               sigpvalue_GWAS =5e-8,sigpvalue_eQTL = 0.05,leadSNP=unique(c(leads,traitleads,traitproxy)),
               chromosome=8, startpos=9921161, stoppos=11972641,
               LD.df	=ld,  trait1 = trait, trait2 =  trait2, congruence=TRUE,getplot = TRUE)

pdf("8.9921161.11972641.SLE.FDFT1 in NaiveT.pdf",width = 11,height = 11)
p$p1+ p$genetracks +p$ld + 
  plot_layout(nrow = 3, byrow = FALSE, heights = c(8,1,10),guides = "collect") &
  theme(legend.direction = "horizontal",legend.justification = "top", 
        legend.key = element_rect(fill = NA, colour = NA, size = 0.25) ) 
dev.off()


# ======= Figure 4 D & Figure S14 ======= 
trait="SLE"
trait_reg=read.table(paste0("chr8.9921161.11972641_FDFT1/",trait,".txt"),header=T)
trait_reg=trait_reg[,c('chromosome','position','BETA','P','SNP')]
colnames(trait_reg)=c("CHR","BP","BETA","P","SNP")
cell=unique(fdft1$filename[which(round(fdft1$PP.H3.abf,digits = 1)>=0.8 & fdft1$Trait==trait)])
traitleads=intersect(unique(allv2g$index_rsid[which(allv2g$trait==trait)]),leads)
traitproxy = unique(allv2g$proxy_rsid[which(allv2g$index_rsid %in% traitleads & allv2g$gene_name=="FDFT1")])
beta=NULL
for(filename in cell){
  b=read.table(paste0("chr8.9921161.11972641_FDFT1/",filename,".ENSG00000079459.txt"),sep="\t",header=F) 
  header=readLines(paste0(unique(match$source[which(match$filename == filename)]),".header"))
  colnames(b)=header
  b=b[,c('rsid','Pvalue','Beta')]
  colnames(b)=c('SNP.Id',"P.Value",'NES')
  b = b[order(b$SNP.Id,b$P.Value,decreasing = F),]
  if(any(duplicated(b$SNP.Id))){b = b[-which(duplicated(b$SNP.Id)),]}
  b$filename = filename
  b$set = ifelse(b$SNP %in% traitleads,"leads",ifelse(b$SNP %in% traitproxy,"proxy","none"))
  beta = rbind(beta,b)
}
beta = merge(beta,unique(fdft1[which(round(fdft1$PP.H3.abf,digits = 1)>=0.8 & fdft1$Trait==trait),c('filename','source','Cell')]),by="filename")
beta$set=factor(beta$set,levels = c("leads",'proxy','none'))
beta$trait=trait

trait="RA"
trait_reg=read.table(paste0("chr8.9921161.11972641_FDFT1/",trait,".txt"),header=T)
trait_reg=trait_reg[,c('chromosome','base_pair_location','beta','p_value','hm_rsid')]
colnames(trait_reg)=c("CHR","BP","BETA","P","SNP")
cell=unique(fdft1$filename[which(round(fdft1$PP.H3.abf,digits = 1)>=0.8 & fdft1$Trait==trait)])
traitleads=intersect(unique(allv2g$index_rsid[which(allv2g$trait==trait)]),leads)
traitproxy = unique(allv2g$proxy_rsid[which(allv2g$index_rsid %in% traitleads & allv2g$gene_name=="FDFT1")])
for(filename in cell){
  b=read.table(paste0("chr8.9921161.11972641_FDFT1/",filename,".ENSG00000079459.txt"),sep="\t",header=F) 
  header=readLines(paste0(unique(match$source[which(match$filename == filename)]),".header"))
  colnames(b)=header
  b=b[,c('rsid','Pvalue','Beta')]
  colnames(b)=c('SNP.Id',"P.Value",'NES')
  b = b[order(b$SNP.Id,b$P.Value,decreasing = F),]
  if(any(duplicated(b$SNP.Id))){b = b[-which(duplicated(b$SNP.Id)),]}
  b$filename = filename
  b$set = ifelse(b$SNP %in% traitleads,"leads",ifelse(b$SNP %in% traitproxy,"proxy","none"))
  b$trait=trait
  b = merge(b,unique(fdft1[which(round(fdft1$PP.H3.abf,digits = 1)>=0.8 & fdft1$Trait==trait),c('filename','source','Cell')]),by="filename")
  beta = rbind(beta,b)
}
beta$set=factor(beta$set,levels = c("leads",'proxy','none'))
beta$OR=exp(beta$NES)

results_set_comp <- beta %>%
  filter(set %in% c("none", "leads", "proxy")) %>%
  mutate(set_group = ifelse(set == "none", "none", "lead_or_proxy")) %>%
  group_by(Cell, trait, source) %>%
  summarise(
    t_test = list(t.test(NES ~ set_group)),
    .groups = "drop"
  ) %>%
  mutate(tidy_result = map(t_test, tidy)) %>%
  unnest(tidy_result) %>%
  select(Cell, trait, source, estimate, estimate1, estimate2, statistic, p.value, conf.low, conf.high)


plot_data <- beta %>%
  filter(set %in% c("none", "leads", "proxy")) %>%
  mutate(set_group = ifelse(set == "none", "none", "lead_or_proxy")) %>%
  left_join(results_set_comp %>% select(Cell, trait, source, p.value), 
            by = c("Cell", "trait", "source")) %>%
  mutate(p_label = ifelse(p.value < 0.001, "***",
                          ifelse(p.value < 0.01, "**",
                                 ifelse(p.value < 0.05, "*", paste0(signif(p.value, 2))) )))
plot_data <- plot_data %>%
  mutate(set = factor(set, levels = c("none", "proxy", "leads")))


pdf("chr8.9921161.11972641.RA_SLE_eQTL.pdf",height =5,width = 7)
ggplot() +
  geom_jitter(data = plot_data ,
              aes(x = Cell, y = NES),color="grey90",
              position = position_jitter(width = 0.1, height = 0),size=0.3,
              alpha = 0.1) +
  
  geom_boxplot(data = plot_data, aes(x = Cell, y = NES, color=set_group),
               width = 0.5, outlier.shape = NA, lwd = 0.5,alpha=0.5) +
  
  geom_text(data = summary_data %>% distinct(Cell, trait, source, p_label,mean_NES , sd_NES),
            aes(x = Cell, y = max(mean_NES + sd_NES, na.rm = TRUE)*2, label = p_label),
            inherit.aes = FALSE, size = 3,color="red") +
  
  scale_color_manual(name="",values = c("none" = "grey40", "lead_or_proxy" = "red")) +
  facet_grid(rows = vars(trait),cols=vars(source),drop=TRUE,scales="free_x",space="free")+
  theme_minimal() +
  labs(title = "",x = "", y = "effect size (beta)") +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,size=8),
    strip.text = element_text(size = 10),legend.position = "bottom"
)
dev.off()


# ======= Figure 4 E ======= 
a=read_h5ad("GSE138266/MS_CSF.h5ad")
meta=a$obs

pbmc = CreateSeuratObject(count=t(a$X),meta.data = meta)
LayerData(pbmc,"scale.data") <- LayerData(pbmc,"counts")
Idents(pbmc)="labels"
pbmc <- subset(pbmc, subset = labels %in% c("B1","B2","CD4","CD8a","CD8n",
                                            "Mono","ncMono","NK1","NK2","Tdg","Tregs",
                                            "mDC1","mDC2","pDC","plasma") & CSF=="False")
pdf("MS.FDFT1.pdf",height = 4,width = 10)
VlnPlot(pbmc, features = "FDFT1", split.by = "MS")+
  stat_compare_means(method = "wilcox.test",hide.ns = TRUE,label="p.signif",tip.length = 0.05)+
  scale_fill_manual(labels = c("Healthy","MS"),values = c("False" = "#00BFC4", "True" = "#F8766D"), name = "Condition") 
dev.off() 

# ======= Figure 4 F ======= 
library(data.table)
meta=read.delim("E-GEAD-397.processed/E-GEAD-397.sdrf.txt",sep="\t")
rownames(meta)= meta$sample = gsub("-Blood","",meta$Source.Name)
colnames(meta)[24]="disease"
meta=meta[which(meta$disease %in% c("healthy control","Systemic Lupus Erythematosus","Rheumatoid Arthritis")),]
cells=list.files("count/")
cells = gsub("_count.txt","",cells)

# extract genes id at FDFT1 locus#########
tmp=unique(allv2g$proxy_rsid[which(allv2g$gene_name=="FDFT1")])
genes=unique(allv2g[which(allv2g$proxy_rsid %in% tmp ),c("ensembl_gene_id","gene_name")])

fdft1=NULL
gatherfdft1=NULL
for(i in c(1:12,14:length(cells))){
  exp=read.delim(paste0("E-GEAD-397.processed/count/",cells[i],"_count.txt"))
  tmp=unlist(strsplit(paste(exp$Gene_id),split="\\."))
  rownames(exp) = paste(tmp[seq(1,length(tmp),by=2)])
  sam=intersect(rownames(meta),colnames(exp))
  dds <- DESeqDataSetFromMatrix(countData = exp[,sam],
                                colData = meta[sam,],
                                design = ~ disease)
  dds$disease <- relevel(dds$disease, ref = "healthy control")
  dds <- DESeq(dds)
  sam=intersect(rownames(meta),colnames(dds))
  fdftgenes <- as.data.frame(dds["ENSG00000079459"]@assays@data$counts )
  gatherFDFT1genes = as.data.frame(fdftgenes) %>%
    tibble::rownames_to_column(var = "ensembl_gene_id") %>%
    tidyr::pivot_longer(cols = colnames(fdftgenes),names_to = "sample",values_to="count")
  gatherFDFT1genes <- dplyr::inner_join(meta[sam,],gatherFDFT1genes)
  gatherFDFT1genes$disease=factor(gatherFDFT1genes$disease,levels = c("healthy control","Rheumatoid Arthritis","Systemic Lupus Erythematosus"))
  gatherFDFT1genes$celltype = cells[i]
  gatherfdft1 = rbind(gatherfdft1,gatherFDFT1genes)
  
  RA <- results(dds,name = "disease_Rheumatoid.Arthritis_vs_healthy.control")
  SLE <- results(dds,name ="disease_Systemic.Lupus.Erythematosus_vs_healthy.control")

    tmp <- rbind (RA["ENSG00000079459",],
                  SLE["ENSG00000079459",])
    tmp$group1="healthy control"
    tmp$group2 = c("Rheumatoid Arthritis","Systemic Lupus Erythematosus")
    tmp$p.format = formatC(tmp$padj,digits = 1) 
    mx=max(gatherFDFT1genes$count)
    mn=min(gatherFDFT1genes$count)
    tmp$y.position = c(mx + (mx-mn)/10,mx + (mx-mn)/5)
    tmp=as_tibble(tmp)
    tmp$gene_name = "FDFT1"
    tmp$celltype = cells[i]
    fdft1 = rbind(fdft1,tmp)
}   
fdft1=fdft1[which(fdft1$padj <= 0.05),]
fdft1$logFC = formatC(fdft1$log2FoldChange,digits = 2)
fdft1$p.signif = symnum(fdft1$padj,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
pdf("SLE_RA.FDFT1.pdf",width = 20,height = 20)
ggboxplot(gatherfdft1,x="disease",y="count",color="disease",add="jitter",palette="lancet")+
    facet_wrap(~celltype,scale="free")+
    stat_pvalue_manual(fdft1,label = "p={p.format}",  tip.length = 0.01,vjust=1.5)+
    theme(axis.text.x = element_blank(),legend.position = "right")
dev.off()

# ======= Figure S15 ======= 
tmp=unique(allv2g$proxy_rsid[which(allv2g$gene_name=="BLK")])
genes=unique(allv2g[which(allv2g$proxy_rsid %in% tmp ),c("ensembl_gene_id","gene_name")])
blkgenes=NULL
blk=NULL
for(i in c(1:13,14:length(cells))){
  dds <- readRDS(paste0("DE.",cells[i],".rds"))
  sam=intersect(rownames(meta),colnames(dds))
  fdftgenes <- as.data.frame(dds["ENSG00000136573"]@assays@data$counts )
  gatherFDFT1genes = as.data.frame(fdftgenes) %>%
    tibble::rownames_to_column(var = "ensembl_gene_id") %>%
    tidyr::pivot_longer(cols = colnames(fdftgenes),names_to = "sample",values_to="count")
  gatherFDFT1genes <- dplyr::inner_join(meta[sam,],gatherFDFT1genes)
  gatherFDFT1genes$disease=factor(gatherFDFT1genes$disease,levels = c("healthy control","Rheumatoid Arthritis","Systemic Lupus Erythematosus"))
  gatherFDFT1genes$celltype = cells[i]
  blkgenes = rbind(blkgenes,gatherFDFT1genes)
  
  RA <- results(dds,name = "disease_Rheumatoid.Arthritis_vs_healthy.control")
  SLE <- results(dds,name ="disease_Systemic.Lupus.Erythematosus_vs_healthy.control")
  
  tmp <- rbind (RA["ENSG00000136573",],
                SLE["ENSG00000136573",])
  tmp$group1="healthy control"
  tmp$group2 = c("Rheumatoid Arthritis","Systemic Lupus Erythematosus")
  tmp$p.format = formatC(tmp$padj,digits = 1) 
  mx=max(gatherFDFT1genes$count)
  mn=min(gatherFDFT1genes$count)
  tmp$y.position = c(mx + (mx-mn)/10,mx + (mx-mn)/5)
  tmp=as_tibble(tmp)
  tmp$gene_name = "BLK"
  tmp$celltype = cells[i]
  blk = rbind(blk,tmp)
}   
blk=blk[which(blk$padj <= 0.05),]
blk$logFC = formatC(blk$log2FoldChange,digits = 2)
blk$p.signif = symnum(blk$padj,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
pdf("SLE_RA.BLK.pdf",width = 20,height = 20)
ggboxplot(blkgenes,x="disease",y="count",color="disease",add="jitter",palette="lancet")+
  facet_wrap(~celltype,scale="free_y",drop = F)+
  stat_pvalue_manual(blk,label = "p={p.format}",  tip.length = 0.01,vjust=1.5)+
  theme(axis.text.x = element_blank(),legend.position = "right")
dev.off()