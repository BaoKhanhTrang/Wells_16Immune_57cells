library(tibble)
library(reshape2)
library(ggplot2)
library(patchwork)
library(nVennR)
library(UpSetR)
library(eulerr)
library(mixOmics)
library(pathfindR)
library(clusterProfiler)
library(DESeq2)
library(ggrepel)
library(ggpubr)
library(ggnewscale)
library(clustifyr)

trait_ord = readRDS("trait_ord.rds")
allstat=readRDS("allstat.rds")
allv2g = readRDS("allv2g.rds")

# ======= Figure 2 B =======
genes = dcast(gene_id~Name,data=unique(allv2g[,c('gene_id','Name')]))
rownames(genes)=genes$gene_id
genes$numCell = apply(genes[,-1],1,function(x) length(which(x!=0)))
genes$Cell = apply(genes[,2:58],1,function(x) paste(sort(na.exclude(unlist(x))),collapse=","))
genes$Cell[which(genes$numCell != 1)] = "shared"

snps = dcast(proxy_rsid~Name,data=unique(allv2g[,c('proxy_rsid','Name')]))
rownames(snps)=snps$proxy_rsid
snps$numCell = apply(snps[,-1],1,function(x) length(which(x!=0)))
snps$Cell = apply(snps[,2:58],1,function(x) paste(sort(na.exclude(unlist(x))),collapse=","))
snps$Cell[which(snps$numCell != 1)] = "shared"

snpgene2 = dcast(proxy_rsid+gene_id~"a",
                data=unique(allv2g[,c('proxy_rsid','gene_id')]))
snpgene2 = merge(snpgene2[,1:2],snps[,c('proxy_rsid','Cell')],by="proxy_rsid")
snpgene2 = merge(snpgene2,genes[,c('gene_id','Cell')],by="gene_id")
colnames(snpgene2)[3:4]=c("variants_inCell","target_genes_inCell")
snpgene = merge(snpgene,snpgene2,by=c('proxy_rsid','gene_id'))
snpgene2=snpgene
snpgene2$variants[which(snpgene2$variants!="shared")]="trait-unique"
snpgene2$variants_inCell[which(snpgene2$variants_inCell!="shared")]="Celltype-specific"
snpgene2$target_genes[which(snpgene2$target_genes != "shared")] = "trait-unique"
snpgene2$target_genes_inCell[which(snpgene2$target_genes_inCell != "shared")] = "Celltype-specific"
snpgene2$variants[which(snpgene2$variants=="shared")]="across-trait"
snpgene2$variants_inCell[which(snpgene2$variants_inCell=="shared")]="across-Celltype"
snpgene2$target_genes[which(snpgene2$target_genes == "shared")] = "across-trait"
snpgene2$target_genes_inCell[which(snpgene2$target_genes_inCell == "shared")] = "across-Celltype"

stat = dcast(variants+variants_inCell+target_genes+target_genes_inCell~"num",data=snpgene2)
stat = stat[order(stat$num,decreasing = T),]
stat$cate=seq(1,nrow(stat))
tmp=melt(stat,id.vars = "cate",measure.vars = c("variants", "target_genes",   "variants_inCell" ,"target_genes_inCell"))

g1=ggplot(tmp)+geom_tile(aes(x=factor(variable),y=factor(cate),group=factor(variable),fill=value),stat = 'identity',position = "dodge")+theme_minimal()+
  theme(legend.position = "bottom",legend.title = element_blank(),axis.text = element_blank())+
  labs(x="",y="")
g2=ggplot(stat)+geom_bar(aes(y=factor(cate),x=num),stat = "identity")+labs(x="x1K SNP-gene pairs",y="")+
  scale_x_continuous(labels = function(x) as.character(x/1000))+
  theme_minimal()+theme(axis.text.y = element_blank())
pdf("how_SNPgene_shared_across_traits_and_cells.pdf",height = 5,width = 5)
g1+g2+plot_layout(ncol=2,widths = c(1,1),guides = "collect")&theme(legend.position = "bottom")
dev.off()

# ======= Figure 2 C created using pyGenomeTracks with correspoding ATAC-seq Bigwig files 
# and chromatin contacts BEDPE files of indivuidual cell types at XPO1 and MTF1 locus =======


# ======= Figure 2 D,E created using Excel with the numbers =======

# ======= Figure 2 F =======
genes = dcast(gene_name+gene_id~trait,data=unique(allv2g[,c('gene_name','gene_id','trait')]))
rownames(genes)=genes$gene_id
genes$numTrait = apply(genes[,-c(1,2)],1,function(x) length(which(x!=0)))
genes$traits = apply(genes[,3:18],1,function(x) paste(sort(na.exclude(unlist(x))),collapse=","))
tmp=data.frame(sort(table(genes$traits[which(genes$numTrait==1)])))
colnames(tmp)=c("trait","unique_genes")
tmp2=dcast(trait~'number of target genes',data=unique(allv2g[,c('trait','gene_id')]))
tmp=merge(tmp,tmp2,by="trait")
tmp$ratio = paste0(round(tmp$unique_genes/tmp$`number of target genes`,digits = 2)*100,"%")
tmp$shared = tmp$`number of target genes` - tmp$unique_genes
tmp = tmp[order(tmp$`number of target genes`,decreasing = T),]
trait_ord=paste(tmp$trait)
tmp$trait = factor(tmp$trait,levels = paste(tmp$trait))
tmp2=melt(tmp[,c('trait','unique_genes','shared')],id.vars = "trait",value.name = "#genes",variable.name = "category")
tmp2$trait=factor(tmp2$trait,levels = levels(tmp$trait))
pdf("numGene_perTrait_sharedvsunique.pdf",height = 5,width = 5)
ggplot(tmp2)+
  geom_bar(aes(`#genes`,trait,fill=category),stat = "identity")+
  geom_text(data=tmp, aes(`number of target genes`+100,trait,label=ratio))+
  scale_fill_grey(name="")+
  theme_bw()+theme(legend.position = c(0.8,0.8))
dev.off()


# ======= Figure 2 G =======
snp=reshape2::dcast(proxy_rsid~trait,data=unique(allv2g[,c('trait','proxy_rsid')]),value.var = "trait")
snp$trait=apply(snp[,-1],1,function(x) paste(na.exclude(x),collapse = ","))
snp$trait[grep(",",snp$trait)] = "shared"

gene=reshape2::dcast(ensembl_gene_id~trait,data=unique(allv2g[,c('trait','ensembl_gene_id')]),value.var = "trait")
gene$trait=apply(gene[,-1],1,function(x) paste(na.exclude(x),collapse = ","))
gene$trait[grep(",",gene$trait)] = "shared"

pair=unique(allv2g[,c('proxy_rsid','ensembl_gene_id')])
pair=merge(pair,snp[,c('proxy_rsid','trait')],by='proxy_rsid')
pair=merge(pair,gene[,c('ensembl_gene_id','trait')],by='ensembl_gene_id')
colnames(pair)[3:4]=c("snp","gene")
rm(snp,gene)
pair$snp=factor(pair$snp,levels = c(trait_ord,"shared"))
pair$gene=factor(pair$gene,levels = c(trait_ord,"shared"))
m=reshape2::dcast(snp~gene,data=pair)
rownames(m)=m$snp
m=m[,-1]
colnames(m)=paste0("v.",colnames(m))
rownames(m)=paste0("g.",rownames(m))
grid.col = structure(rep("",length(unlist(dimnames(m)))),names = unlist(dimnames(m)))
for(i in 1:length(trait_colors)){
  grid.col[grep(names(trait_colors)[i],names(grid.col))] = trait_colors[i]
}
grid.col[grep("shared",names(grid.col))]="grey75"

library(circlize)
pdf("SNP2gene_shared_and_unique_acrossTrait.pdf",height=8,width=8)
chordDiagram(as.matrix(m),annotationTrack = "grid",grid.col = grid.col,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(m))))/1.5))
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1]+.1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  #circos.axis(h = "bottom", labels.cex = 0.5, major.tick.percentage = 0.05, sector.index = sector.name, track.index = 1)
}, bg.border = NA) 
dev.off()

# ======= Figure 2 H =======
gene=dcast(trait~'number of target genes',data=unique(allv2g[,c('trait','ensembl_gene_id')]))
gene=gene[order(gene$`number of target genes`,decreasing = T),]
gene$trait = factor(gene$trait,levels = paste(gene$trait))
gene=merge(gene,snp,by='trait')
library(ggrepel)
pdf("numgene_numsnp_pertrait.pdf",height=5,width = 3)
ggplot(gene,aes(num,`number of target genes`))+
  geom_point()+
  geom_label_repel(aes(label=trait),box.padding   = 0.5, 
                   point.padding = 0.1,
                   segment.color = 'grey50')+
  theme_minimal()+
  labs(x='number of variants\nfound contacting gene')
dev.off()

# ======= Figure 2 I =======
gene=dcast(Name+system~'number of target genes',data=unique(allv2g[,c('Name','gene_id','system')]))
gene=merge(gene,snp,by=c('Name','system'))
gene=gene[order(gene$`number of target genes`,decreasing = T),]
gene$Name = factor(gene$Name,levels = paste(gene$Name))
gene$sec=1
gene$sec[which(gene$`number of target genes`<1000)]=2
gene$sec[which(gene$num>300)]=3

pdf("numgene_numsnp_perCells.pdf",height=5,width = 9)
ggplot(gene,aes(num,`number of target genes`))+
  geom_point()+
  geom_label_repel(aes(label=Name,color=system), 
                   max.overlaps = 55,
                   box.padding   = 0.01,
                   point.padding = 0.1,
                   segment.color = 'grey50',min.segment.length=0.1)+
  theme_minimal()+theme(legend.position = "none")+
  labs(x='number of variants\nfound contacting gene')
dev.off()

# ======= Figure S4 A =======
gene = dcast(gene_id~trait,data=unique(allv2g[,c('gene_id','trait')]),fun.aggregate = length)
rownames(gene)=gene$gene_id
gene=gene[,-1]
pdf("gene_per_trait_upset.pdf",height = 5,width = 22)
upset(gene,nsets = 16,nintersects = 334,mb.ratio = c(0.5,0.5),order.by = 'degree')
dev.off()

# ======= Figure S4 B =======
gene=dcast(ensembl_gene_id~Cellline,data=unique(allv2g[which(allv2g$system=="immune"),c('ensembl_gene_id','Cellline')]),fun.aggregate = length)
rownames(gene)=paste(gene$ensembl_gene_id)
gene=gene[,-1]
pdf("ImmuneGenes_across_cellLines.pdf",height = 4,width = 6)
upset(gene,nintersects = NA,mb.ratio = c(0.7,0.3))
dev.off()

# ======= Figure S4 C =======
gene = dcast(gene_id+gene_name~'num',data=unique(allv2g[,c('gene_id','gene_name','Name')]),fun.aggregate = length)
pdf("gene_sharedCell_barplot.pdf",width = 5,height = 8)
ggplot(as.data.frame(table(gene$num)))+ 
  geom_bar(aes(Freq,Var1),stat = 'identity')+theme_minimal()+
  geom_text(aes(Freq+150,Var1,label=Freq))+
  labs(x="#genes",y="#cell types shared")
dev.off()

# ======= Figure S4 D =======
proxy = dcast(proxy_rsid~Name,data=unique(allv2g[,c('proxy_rsid','Name')]),fun.aggregate = length)
rownames(proxy)=proxy$proxy_rsid
proxy=proxy[,-1]
proxy_cell = apply(proxy,1,function(x) paste(colnames(proxy)[which(x!=0)],collapse = ","))
proxy_cellnum = apply(proxy,1,function(x) length(which(x!=0)))
pdf("SNP_sharedCell_barplot.pdf",width = 5,height = 8)
ggplot(as.data.frame(table(proxy_cellnum)))+ 
  geom_bar(aes(Freq,proxy_cellnum),stat = 'identity')+
  geom_text(aes(Freq+75,proxy_cellnum,label=Freq))+
  theme_minimal()+
  labs(x="#variants",y="#cell types shared")
dev.off()

# ======= Figure S8 A-B-C =======
tmp=unique(allv2g[which(allv2g$Name=="Enteroid" & allv2g$trait %in% c("CRO","IBD","UC")),c('gene_name','trait')])
library(clusterProfiler)
re=compareCluster(gene_name~trait,data=tmp,keyType="SYMBOL" ,fun = enrichGO,OrgDb = "org.Hs.eg.db",pAdjustMethod = "BH",pvalueCutoff  = 1,qvalueCutoff  = 1)

l=list()
for(i in 1:3){
  l[[i]] = unique(tmp$gene_name[which(tmp$trait == unique(tmp$trait)[i])])
}
names(l) = unique(tmp$trait)
pdf("UC_IBD_CRO_in_enteroid.pdf",height = 3,width = 3)
upset(fromList(l),nintersects = NA,nsets = 7)
dev.off()
pdf("UC_IBD_CRO_in_enteroid_path.pdf",height = 7,width = 5)
dotplot(re)
dev.off()
pdf("UC_IBD_CRO_in_enteroid_network.pdf",height = 7,width = 7)
cnetplot(re)
dev.off()