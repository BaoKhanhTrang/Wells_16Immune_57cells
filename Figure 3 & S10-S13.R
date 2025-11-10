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


match=read.table("matching_celltypes.txt",sep="\t",header=T)
done=readRDS("258loci.rds")
trait_ord = readRDS("trait_ord.rds")
allv2g=readRDS("allv2g.rds")
allv2g = allv2g[which(allv2g$system=="immune"),]
gtex=readRDS("GTEx.eGenes_ld.rds")
gtex$keep <- ifelse(gtex$shared == "shared", "overlapped",
                   ifelse(gtex$shared == "V2G", "V2G only",
                          ifelse(gtex$shared == "outsideLD", "outsideLD", "Diverged")))

dice=readRDS("DICE.eGenes_ld.rds")
dice$keep <- case_when(
  dice$shared == "shared"      ~ "overlapped",
  dice$shared == "V2G"         ~ "V2G only",
  dice$shared == "outsideLD"   ~ "outsideLD",
  TRUE                        ~ "Diverged"
)

onek=readRDS("OneK1K.eGenes_ld.rds")
onek$keep <- case_when(
  onek$shared == "shared"      ~ "overlapped",
  onek$shared == "V2G"         ~ "V2G only",
  onek$shared == "outsideLD"   ~ "outsideLD",
  TRUE                        ~ "Diverged"
)

eqtlc=readRDS("eQTL_catalogue.eGenes_ld.rds")
eqtlc$keep <- case_when(
  eqtlc$shared == "shared"      ~ "overlapped",
  eqtlc$shared == "V2G"         ~ "V2G only",
  eqtlc$shared == "outsideLD"   ~ "outsideLD",
  TRUE                        ~ "Diverged"
)


# ======= Figure S10 ======= 
genes = list(
  V2G = unique(allv2g$ensembl_gene_id),
  GTEx_WholeBlood = unique(gtex$ensembl_gene_id[which(gtex$keep!="V2G only")]),
  DICE = unique(dice$ensembl_gene_id[which(dice$keep!="V2G only")]),
  OneK1K = unique(onek$ensembl_gene_id[which(onek$keep!="V2G only")]),
  eQTL_catalogue = unique(eqtlc$ensembl_gene_id[which(eqtlc$keep!="V2G only")])#,
)
venn <- Venn(genes)
data <- process_data(venn)

library(colorjam)

base_colors=c("#e0d3ec","#4d99ff","#FF0000","#FFC000","#70AD47")
names(base_colors)=seq(1,5)
# Generate all possible color intersections (1 to 5 colors)
all_combinations <- unlist(lapply(1:length(base_colors), function(n) {
  combn(base_colors, n, simplify = FALSE)
}), recursive = FALSE)
# Create a palette based on these intersections
palette_intersections <- sapply(all_combinations, blend_colors)
names(palette_intersections) = names(all_combinations) <- sapply(all_combinations, function(x) paste(names(x), collapse = "/"))

data$regionEdge$col <- plyr::mapvalues(
  data$regionEdge$id,
  unique(data$regionEdge$id),
  palette_intersections
)
zer=data$regionLabel$id[which(data$regionData$count==0)]
if(length(zer)>0){data$regionEdge$col[which(data$regionEdge$id %in% zer)] = "#ffffff"}

pdf("shared_genes.pdf")
ggplot() +
  geom_polygon(aes(X, Y, fill = I(col), group = id), data = venn_regionedge(data), show.legend = F) +
  geom_path(aes(X, Y, color = id, group = id), data = venn_setedge(data), show.legend = F) +
  geom_text(aes(X, Y, label = name), fontface = "bold",data = venn_setlabel(data)) +
  geom_label(aes(X, Y, label = count), data = venn_regionlabel(data),alpha = 0,label.size = 0) +
  coord_equal() +
  theme_void() + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2))
dev.off()


# ======= Figure 3 A ======= 
# filter diff-crhom and outsideLD #######
gtex=gtex[-which(gtex$cis=="diff_chrom" | gtex$keep=="outsideLD" ),]
dice=dice[-which(dice$cis=="diff_chrom" | dice$keep=="outsideLD" ),]
onek=onek[-which(onek$cis=="diff_chrom" | onek$keep=="outsideLD" ),]
eqtlc=eqtlc[-which(eqtlc$cis=="diff_chrom" | eqtlc$keep=="outsideLD" ),]
imu=imu[-which(imu$cis=="diff_chrom" | imu$keep=="outsideLD" ),]

genes = list(
  V2G = unique(allv2g$ensembl_gene_id),
  GTEx_WholeBlood = unique(gtex$ensembl_gene_id[which(gtex$keep!="V2G only")]),
  DICE = unique(dice$ensembl_gene_id[which(dice$keep!="V2G only")]),
  OneK1K = unique(onek$ensembl_gene_id[which(onek$keep!="V2G only")]),
  eQTL_catalogue = unique(eqtlc$ensembl_gene_id[which(eqtlc$keep!="V2G only")])
)

venn <- Venn(genes)
data <- process_data(venn)

library(colorjam)

base_colors=c("#e0d3ec","#4d99ff","#FF0000","#FFC000","#70AD47")
names(base_colors)=seq(1,5)
# Generate all possible color intersections (1 to 5 colors)
all_combinations <- unlist(lapply(1:length(base_colors), function(n) {
  combn(base_colors, n, simplify = FALSE)
}), recursive = FALSE)
# Create a palette based on these intersections
palette_intersections <- sapply(all_combinations, blend_colors)
names(palette_intersections) = names(all_combinations) <- sapply(all_combinations, function(x) paste(names(x), collapse = "/"))

data$regionEdge$col <- plyr::mapvalues(
  data$regionEdge$id,
  unique(data$regionEdge$id),
  palette_intersections
)
zer=data$regionLabel$id[which(data$regionData$count==0)]
if(length(zer)>0){data$regionEdge$col[which(data$regionEdge$id %in% zer)] = "#ffffff"}

pdf("shared_genes_exclude_chromdiff_andOutsideLD.noImmuNexUT.pdf")
ggplot() +
  geom_polygon(aes(X, Y, fill = I(col), group = id), data = venn_regionedge(data), show.legend = F) +
  geom_path(aes(X, Y, color = id, group = id), data = venn_setedge(data), show.legend = F) +
  geom_text(aes(X, Y, label = name), fontface = "bold",data = venn_setlabel(data)) +
  geom_label(aes(X, Y, label = count), data = venn_regionlabel(data),alpha = 0,label.size = 0) +
  coord_equal() +
  theme_void() + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2))
dev.off()

# ======= Figure 3 B ======= 
Cellline=data.frame(Cellline=unique(allv2g$Cellline[which(allv2g$system=="immune")]))
for(i in 1:nrow(Cellline)){
  v2g = unique(allv2g$ensembl_gene_id[which( allv2g$Cellline ==Cellline$Cellline[i] )])
  index = unique(allv2g$index_rsid[which( allv2g$Cellline ==Cellline$Cellline[i])])
  locus=NULL
  for(e in 1:length(index)){
    locus = union(locus,grep(index[e],done$rsid))
  }
  cells= unique(match$filename[which(match$Celltype %in% unique(allv2g$Celltype[which(allv2g$Cellline ==Cellline$Cellline[i])]) )])
  eq = which(eGenes$filename %in% cells &
               eGenes$locus %in% done$id[locus] & eGenes$keep!="V2G only")
  
  Cellline$`V2G only`[i] = length(setdiff(v2g,eGenes$ensembl_gene_id[eq]))
  Cellline$overlapped[i] = length(unique(eGenes$ensembl_gene_id[intersect(eq,which(eGenes$keep=="overlapped"))]))
  Cellline$`eQTL only`[i] = length(unique(eGenes$ensembl_gene_id[intersect(eq,which(eGenes$keep=="Diverged"))]))
  Cellline$total[i] = length(unique(eGenes$ensembl_gene_id[eq]))
}
tmp2=reshape2::melt(Cellline[,1:4])
pdf("NumerGene_perCellLine_v4.noImmuNexUT.pdf",height = 5,width = 6)
ggplot(tmp2)+
  geom_bar(aes(1,y=value,fill=variable),position = position_fill(),stat = "identity")+
  geom_text(
    data = Cellline,                   # Use the totals data frame
    aes(x = 0, y = 0, label = total),   # Position text at x=1, y=0 (center) and use 'total' column for label
    inherit.aes = FALSE,                # VERY Important: Prevent inheriting aesthetics (like fill, y) from ggplot()
    size = 3.5,                         # Adjust text size as needed
    color = "black"                     # Adjust text color as needed (e.g., "white")
  ) +
  coord_polar(theta = "y")+
  theme_void()+
  facet_grid(cols=vars(Cellline),scales = "free")+
  scale_fill_manual(values = c("eQTL only"="#826aed","overlapped"="#f6cd30","V2G only"="grey80"))+
  guides(fill=guide_legend(nrow=1,title = ""))+
  theme(aspect.ratio = 1,legend.position = "bottom",strip.text.y = element_text(hjust = 0))
dev.off()


# ======= Figure S11 ======= 
pdf("shared_genes_perLoci_exclude_chromdiff_andOutsideLD.noImmuNexUT.pdf")
for(i in 1:nrow(done)){
  genes = list(
    V2G = unique(allv2g$ensembl_gene_id[which(allv2g$index_rsid %in% unlist(done$rsid[i]))]),
    GTEx_WholeBlood = unique(gtex$ensembl_gene_id[which(gtex$locus == paste(done$id[i]) & gtex$keep!="V2G only")]),
    DICE = unique(dice$ensembl_gene_id[which(dice$locus == paste(done$id[i]) & dice$keep!="V2G only")]),
    OneK1K = unique(onek$ensembl_gene_id[which(onek$locus == paste(done$id[i]) & onek$keep!="V2G only")]),
    eQTL_catalogue = unique(eqtlc$ensembl_gene_id[which(eqtlc$locus == paste(done$id[i]) & eqtlc$keep!="V2G only")])
  )
  venn <- Venn(genes)
  data <- process_data(venn)
  data$regionEdge$col <- plyr::mapvalues(
    data$regionEdge$id,
    unique(data$regionEdge$id),
    palette_intersections
  )
  zer=data$regionLabel$id[which(data$regionData$count==0)]
  if(length(zer)>0){data$regionEdge$col[which(data$regionEdge$id %in% zer)] = "#ffffff"}
  g=ggplot() +
    geom_polygon(aes(X, Y, fill = I(col), group = id), data = venn_regionedge(data), show.legend = F) +
    geom_path(aes(X, Y, color = id, group = id), data = venn_setedge(data), show.legend = F) +
    geom_text(aes(X, Y, label = name), fontface = "bold",data = venn_setlabel(data)) +
    geom_label(aes(X, Y, label = count), data = venn_regionlabel(data),alpha = 0,label.size = 0) +
    coord_equal() +
    ggtitle(paste("locus:",done$id[i]),)+theme_void()+
    theme(plot.title = element_text(hjust = 0.5))+
     ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2))
  print(g)
}
dev.off()


# ======= Figure S12 ======= 
dice$source="DICE"
eqtlc$source="eQTL_catalogue"
onek$source="OneK1K"
eGenes = rbind(dice,eqtlc,onek)

pair=expand.grid(unique(match$Cell),trait_ord)
colnames(pair)=c("Cell","trait")
for(i in 1:nrow(pair)){
  v2g = unique(allv2g$ensembl_gene_id[which(allv2g$trait == pair$trait[i] & allv2g$Celltype %in% match$Celltype[which(match$Cell == pair$Cell[i])])])
  index = unique(allv2g$index_rsid[which(allv2g$trait == pair$trait[i] & allv2g$Celltype %in% match$Celltype[which(match$Cell == pair$Cell[i])])])
  if(length(index)==0){
    index = unique(allv2g$index_rsid[which(allv2g$Celltype %in% match$Celltype[which(match$Cell == pair$Cell[i])])])
  }
  locus=NULL
  for(e in 1:length(index)){
    locus = union(locus,grep(index[e],done$rsid))
  }
  eq = which(eGenes$filename %in% match$filename[which(match$Cell == pair$Cell[i])] &
                 eGenes$locus %in% done$id[locus] & eGenes$keep!="V2G only")
  
  pair$`V2G only`[i] = length(setdiff(v2g,eGenes$ensembl_gene_id[eq]))
  pair$overlapped[i] = length(unique(eGenes$ensembl_gene_id[intersect(eq,which(eGenes$keep=="overlapped"))]))
  pair$`eQTL only`[i] = length(unique(eGenes$ensembl_gene_id[intersect(eq,which(eGenes$keep=="Diverged"))]))
}
tmp2=reshape2::melt(pair)
tmp2$Cell = factor(tmp2$Cell,levels = c("NaiveB","B mem/GCB","plasma/PGF","NaiveT",
                                      "Stimulated T-cells","T follicular helper cells","Th1" ,"Th2","Th17","Treg" ,"Treg memory/activated" , 
                                      "Dendritic cells" ,"Monocytes","Monocytes-nonclassical" ,"Natural killer" ))
tmp2$trait=factor(tmp2$trait,levels = trait_ord)

g1=ggplot(tmp2)+
  geom_bar(aes(1,value,fill=variable) ,stat = 'identity',position = position_fill())+
  coord_polar(theta = "y")+
  facet_grid(cols = vars(trait),rows = vars(Cell))+
  theme_void()+
  scale_fill_manual(values = c("eQTL only"="#826aed","overlapped"="#f6cd30","V2G only"="grey80"))+
  guides(fill=guide_legend(ncol=1,title = ""))+
  theme(aspect.ratio = 1,legend.position = "bottom",strip.text.y = element_text(hjust = 0))

percell = data.frame(Cell=unique(match$Cell))
percell$Cell = factor(percell$Cell,levels = c("NaiveB","B mem/GCB","plasma/PGF","NaiveT",
                                              "Stimulated T-cells","T follicular helper cells","Th1" ,"Th2","Th17","Treg" ,"Treg memory/activated" , 
                                              "Dendritic cells" ,"Monocytes","Monocytes-nonclassical" ,"Natural killer" ))
for(i in 1:nrow(percell)){
  v2g = unique(allv2g$ensembl_gene_id[which( allv2g$Celltype %in% match$Celltype[which(match$Cell == percell$Cell[i])])])
  index = unique(allv2g$index_rsid[which( allv2g$Celltype %in% match$Celltype[which(match$Cell == percell$Cell[i])])])
  locus=NULL
  for(e in 1:length(index)){
    locus = union(locus,grep(index[e],done$rsid))
  }
  eq = which(tmp$filename %in% match$filename[which(match$Cell == percell$Cell[i])] &
               tmp$locus %in% done$id[locus] & tmp$keep!="V2G only")
  
  percell$`V2G only`[i] = length(setdiff(v2g,tmp$ensembl_gene_id[eq]))
  percell$overlapped[i] = length(unique(tmp$ensembl_gene_id[intersect(eq,which(tmp$keep=="overlapped"))]))
  percell$`eQTL only`[i] = length(unique(tmp$ensembl_gene_id[intersect(eq,which(tmp$keep=="Diverged"))]))
  
}
tmp2=reshape2::melt(percell)
tmp2$Cell = factor(tmp2$Cell,levels = c("NaiveB","B mem/GCB","plasma/PGF","NaiveT",
                                        "Stimulated T-cells","T follicular helper cells","Th1" ,"Th2","Th17","Treg" ,"Treg memory/activated" , 
                                        "Dendritic cells" ,"Monocytes","Monocytes-nonclassical" ,"Natural killer" ))

g2=ggplot(tmp2)+
  geom_bar(aes(y=variable,x=-value,fill=variable),position = "stack",stat = "identity")+
  facet_grid(rows=vars(Cell),scales = "free")+
  scale_fill_manual(values = c("eQTL only"="#826aed","overlapped"="#f6cd30","V2G only"="grey80"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1),
        strip.text = element_blank(),
        legend.title = element_blank(),legend.position = "none",axis.title = element_blank())


pertrait = data.frame(trait=trait_ord)
for(i in 1:nrow(pertrait)){
  v2g = unique(allv2g$ensembl_gene_id[which(allv2g$trait == pertrait$trait[i] )])
  cells=unique(allv2g$Celltype[which(allv2g$trait == pertrait$trait[i] )])
  index = unique(allv2g$index_rsid[which(allv2g$trait == pertrait$trait[i] )])
  locus=NULL
  for(e in 1:length(index)){
    locus = union(locus,grep(index[e],done$rsid))
  }
  eq = which(tmp$filename %in% match$filename[which(match$Celltype %in% cells)] &
               tmp$locus %in% done$id[locus] & tmp$keep!="V2G only")
  
  pertrait$`V2G only`[i] = length(setdiff(v2g,tmp$ensembl_gene_id[eq]))
  pertrait$overlapped[i] = length(unique(tmp$ensembl_gene_id[intersect(eq,which(tmp$keep=="overlapped"))]))
  pertrait$`eQTL only`[i] = length(unique(tmp$ensembl_gene_id[intersect(eq,which(tmp$keep=="Diverged"))]))
}
tmp2=reshape2::melt(pertrait)
tmp2$trait = factor(tmp2$trait,levels = trait_ord)
g3=ggplot(tmp2)+
  geom_bar(aes(variable,y=-value,fill=variable),position = "stack",stat = "identity")+
  facet_grid(cols=vars(trait),scales = "free")+
  scale_fill_manual(values = c("eQTL only"="#826aed","overlapped"="#f6cd30","V2G only"="grey80"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1),
        strip.text = element_blank(),
        legend.title = element_blank(),legend.position = "none",axis.title = element_blank())


pdf("NumerGene_perCell_perTrait_v4.noImmuNexUT.pdf",height = 10,width = 12)
g2+g1+plot_spacer()+g3+plot_layout(widths = c(1,8),nrow = 2,heights = c(8,1),guides = "collect") & theme(legend.position = 'bottom')
dev.off()

# ======= Figure 3 C ======= 
colo=NULL
for(t in trait_ord){
  for(db in c("DICE","OneK1K","eQTL_catalogue")){
    tmp=read.table(paste0(t,".",db,".coloc.txt"),header=T)
    if(nrow(tmp)>0){
      tmp$trait=t
      tmp$source=db
      colo = rbind(colo,tmp)
    }
  }
}
colnames(colo)[which(colnames(colo)=="Cell")]="filename"
colo=merge(colo,match[,c('source','filename','Cell')],by=c('source','filename'))
colo = unique(colo)
colo$Cell = factor(colo$Cell,levels = c("NaiveB","B mem/GCB","plasma/PGF","NaiveT",
                                                         "Stimulated T-cells","T follicular helper cells","Th1" ,"Th2","Th17","Treg" ,"Treg memory/activated" , 
                                                         "Dendritic cells" ,"Monocytes","Monocytes-nonclassical" ,"Natural killer" ))
colo$trait = factor(colo$trait,levels = trait_ord)

eGenes = unique(eGenes[-which(is.na(eGenes$trait)),c("locus","filename","ensembl_gene_id","Celltype","trait","keep","source")])
eGenes$eg="eg"
eGenes = merge(eGenes,coloc[,c(1:16)],by= intersect(colnames(eGenes),colnames(coloc)[1:16]),all=TRUE)

gen <- getBM(
  attributes = c("ensembl_gene_id", "uniprot_gn_symbol","hgnc_symbol"), 
  values = eGenes$ensembl_gene_id,filters = "ensembl_gene_id",
  mart = ensembl
)

for(i in 1:nrow(eGenes)){
  v2g=unique(allv2g$proxy_rsid[which(allv2g$trait==eGenes$trait[i] & allv2g$ensembl_gene_id==eGenes$ensembl_gene_id[i] & allv2g$Celltype == eGenes$Celltype[i])] )
  eGenes$numSnps[i] = length(v2g)
  if(is.na(eGenes$nsnps[i])){
    if(is.na(eGenes$gene_name[i])){
      eGenes$gene_name[i] = unique(allv2g$gene_name[which(allv2g$ensembl_gene_id==eGenes$ensembl_gene_id[i])])
    }
    
    name = unique(c(coloc$gene_name[which(coloc$ensembl_gene_id == eGenes$gene_name[i])],eGenes$gene_name[i],
                    unique(allv2g$gene_name[which(allv2g$ensembl_gene_id==eGenes$ensembl_gene_id[i])])))
    file=paste(eGenes$locus[i],eGenes$trait[i],eGenes$source[i],eGenes$filename[i],name,"",sep=".")
    l=list.files(path="RDS_result_files/",pattern = file)
    if(length(l)>0){
      file=paste0("RDS_result_files/" ,l)
      if(length(file)>1){
        truefile = paste(eGenes$locus[i],eGenes$trait[i],eGenes$source[i],eGenes$filename[i],eGenes$gene_name[i],"coloc.rds",sep=".")
        truefile=paste0("RDS_result_files/" ,truefile)
        tmp=readRDS(truefile)
      }else{
        tmp=readRDS(file)
      }
      
      if(purrr::is_empty(tmp)){
        eGenes[i,c("nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","cond.PP.H4.abf")]=0
      }else{
        eGenes[i,names(tmp$summary)] = tmp$summary
        eGenes$cond.PP.H4.abf[i] = tmp$summary["PP.H4.abf"] / sum(tmp$summary[c("PP.H3.abf","PP.H4.abf")])
        eGenes$Cell[i] = paste(unique(match$Cell[which(match$filename==eGenes$filename[i] & match$source==eGenes$source[i])]))
      }
    }
  }
}

eGenes$trait = factor(eGenes$trait,levels = trait_list$TRAIT)
eGenes = merge(eGenes,unique(allv2g[,c('Cellline','Celltype')]),by="Celltype")
pdf("coloc_hilightP4.pdf",width = 10,height = 6)
ggtern(eGenes,aes(x=PP.H3.abf,y=PP.H4.abf,z=PP.H0.abf+PP.H1.abf+PP.H2.abf))+ 
  geom_point(alpha = 0.02,size = 1)+
  geom_point(data=eGenes[which(eGenes$PP.H4.abf>=0.8 & eGenes$keep=="overlapped" ),],aes(color=trait,shape=Cellline),size=2)+
  geom_point(data=eGenes[which(eGenes$PP.H3.abf>=0.8 & eGenes$keep=="overlapped" ),],aes(color=trait,shape=Cellline),size=2)+
  theme_rgbw()+
  scale_color_manual(values=trait_colors)+
  theme_hidetitles()+
  theme_hidegrid_minor()+
  guides(color=guide_legend(ncol=2))
dev.off()


# ======= Figure 3 D ======= 
p4=unique(eGenes[which(eGenes$PP.H4.abf>=0.8 & eGenes$keep=="overlapped" ),])
p4$v2g=""
p4$coloc=""
p4$matchSNP=""
inp=NULL
library(biomaRt)
ensembl <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
for(i in 1:nrow(p4)){
  file=paste0("coloc_results/",p4$trait[i],"/",p4$source[i],"/",p4$locus[i],".",p4$ensembl_gene_id[i],".",p4$gene_name[i],".",p4$filename[i],".coloc_results_full.txt")
  if(file.exists(file)){
    tmp=read.table(file,header=T)
  }else{
    file=paste(p4$locus[i],p4$trait[i],p4$source[i],p4$filename[i],p4$gene_name[i],"coloc.rds",sep=".")
    file=paste0( "RDS_result_files/",file)
    tmp=readRDS(file)
    tmp=tmp$results
  }
  
  
  fm=unique(tmp$snp[tmp$SNP.PP.H4 >= 0.8])
  p4$coloc[i]=ifelse(length(fm) > 0, paste(fm, collapse = ","),"")
  
  vg = unique(allv2g$proxy_rsid[which(allv2g$trait==p4$trait[i] & allv2g$ensembl_gene_id==p4$ensembl_gene_id[i] & allv2g$Celltype == p4$Celltype[i])])
  inp = rbind(inp,
              data.frame(LOCUS_ID=p4$locus[i],CATEGORY=trait_list$CATEGORY[which(trait_list$TRAIT==p4$trait[i])],
                         TRAIT = p4$trait[i],
                         CHR=sub("chr","",unique(allv2g$chr[which(allv2g$ensembl_gene_id==p4$ensembl_gene_id[i])])),
                         BP=unique(allv2g$variant_pos[which(allv2g$proxy_rsid %in% vg)]),
                         MARKER=vg,GENE=p4$gene_name[i],Type="V2G"))
  
  if(length(fm)>0){
    snp_info <- getBM(
      attributes = c("refsnp_id", "chr_name", "chrom_start"),  # rsID, chromosome, position
      filters = "snp_filter",
      values = fm,
      mart = ensembl
    )
    inp = rbind(inp,
                data.frame(LOCUS_ID=p4$locus[i],CATEGORY=trait_list$CATEGORY[which(trait_list$TRAIT==p4$trait[i])],
                           TRAIT = p4$trait[i],
                           CHR=sub("chr","",unique(allv2g$chr[which(allv2g$ensembl_gene_id==p4$ensembl_gene_id[i])])),
                           BP=snp_info$chrom_start,
                           MARKER=snp_info$refsnp_id,
                           GENE=p4$gene_name[i],Type="eQTL")
    )
  }
  tmp=tmp[which(tmp$snp %in% vg),]
  if(nrow(tmp)>0){
    tmp$vg = paste0(tmp$snp," (",round(tmp$SNP.PP.H4,digits = 2),")")
    p4$v2g[i]=paste(tmp$vg,collapse = ",")
  }
  if(length(fm)>0 & length(vg)>0){
    tmp = intersect(fm,vg)
    if(length(tmp)>0){
      p4$matchSNP[i] = paste(tmp)
    }
  }
  rm(vg,fm,tmp,file)
}
inp=unique(inp)
write.table(inp,file="Immune.p4.input.txt",col.names=T,row.names = F,sep="\t",quote=F)

p4=merge(p4,unique(allv2g[,c("Celltype","Cellline")]),by="Celltype")


stat=reshape2::dcast(trait+Cellline~"num",data=unique(p4[,c('trait','Cellline','gene_name')]),value.var = "gene_name")
stat$trait = factor(stat$trait,levels = rev(trait_list$TRAIT))

# Create the topleft dotplot of Fig3-D
pdf("dotplot.pdf",height = 3,width = 4.5)
ggplot(stat)+geom_point(aes(trait,num,shape=Cellline,color=trait),size=2)+
  scale_color_manual(values = trait_colors,guide="none")+
  coord_flip()+
  scale_y_reverse(name="#genes",expand = c(0,1),breaks=seq(1,19,by=3)) + 
  scale_x_discrete(name="") +
  theme_bw()+
  theme(legend.position = "left", 
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill='transparent'))
dev.off()

# The main circos plot of Fig3-D is created using the input file "Immune.p4.input.txt" and
# Fiju plot https://github.com/mkanai/fujiplot, a circos representation of multiple GWAS results based on Circos software (http://circos.ca/).
Rscript fujiplot/fujiplot.R Immune.p4.input.txt Immune.traitlist.txt Immune.p4
circos -conf fujiplot/config/circos.conf -debug_group summary


# ======= Figure 3 E ======= 
pdf("p4_highH4_8SNP.pdf",height = 3.5,width = 6)
ggplot(p4)+
  geom_point(aes(x=BETA,y=Zeqtl,color=trait,shape=Cellline))+
  geom_text_repel(aes(x=BETA,y=Zeqtl,label = gene_name,color=trait))+
  scale_shape_manual(values = Cellline_shapes)+
  theme_bw()+
  scale_color_manual(values=trait_colors)+
  labs(x="Effect size on immune trait",y="Effect size on gene expression")
dev.off()  

# ======= Figure S13 ======= 
library(motifbreakR)

mb=readRDS("allSNPs_motifBreakR_jaspar_v2g.rds") # from TF motifbreakR analysis
mb4=mb[mb$SNP_id %in% unique(p4$matchSNP) & mb$effect=="strong"]
pdf("rs1131017.pdf",height = 5,width = 4)
plotMB(mb,rsid = "rs1131017",altAllele = "G")
dev.off()
pdf("rs113075206.pdf",height = 5,width = 4)
plotMB(mb4,rsid = "rs113075206",altAllele = "G")
dev.off()
pdf("rs17323934.pdf",height = 12,width = 4)
plotMB(mb,rsid = "rs17323934",altAllele = "G")
dev.off()
pdf("rs7528684.pdf",height = 5,width = 4)
plotMB(mb4,rsid = "rs7528684",altAllele = "T",reverseMotif = T)
dev.off()
pdf("rs479777.pdf",height = 12,width = 5)
plotMB(mb,rsid = "rs479777")
dev.off()
pdf("rs7731626.pdf",height = 11,width = 6)
plotMB(mb4,rsid = "rs7731626")
dev.off()
pdf("rs6908626.pdf",height = 5,width = 6)
plotMB(mb4,rsid = "rs6908626",altAllele = "T",reverseMotif = T)
dev.off()
pdf("rs72928038.pdf",height = 7,width = 5)
plotMB(mb4,rsid = "rs72928038")
dev.off()