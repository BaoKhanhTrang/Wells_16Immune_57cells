library(reshape2)
library(ggplot2)
library(ggsci)
library(paletteer)


allv2g=readRDS("allv2g.rds")
allstat=readRDS("allstat.rds")
metaData=readRDS("metaData.rds")
allExp=readRDS(file="allExp.rds")

pir = resmelt(dcast(data=unique(allv2g[which(allv2g$OCRtype=="PIR"),c('proxy_rsid','gene_name','Celltype')]),gene_name~Celltype,fun.aggregate = length))
colnames(pir)[3] = "PIR"
pir = pir[-which(pir$PIR==0),]
prom = melt(dcast(data=unique(allv2g[which(allv2g$OCRtype=="Prom"),c('proxy_rsid','gene_name','Celltype')]),gene_name~Celltype,fun.aggregate = length))
colnames(prom)[3] = "Prom"
prom = prom[-which(prom$Prom==0),]
colnames(prom)[2]=colnames(pir)[2]="Celltype"
prom$pair=paste(prom$gene_name,prom$Celltype,sep = "_")
pir$pair=paste(pir$gene_name,pir$Celltype,sep = "_")
prom$PIR=0
pir$Prom=0
exp_num = merge(pir,prom,by="pair",all=TRUE)
exp_num$gene_name = exp_num$gene_name.x
exp_num$gene_name[which(is.na(exp_num$gene_name.x))] = exp_num$gene_name.y[which(is.na(exp_num$gene_name.x))]
exp_num$Celltype = apply(exp_num[,c('Celltype.x','Celltype.y')],1,function(x) paste(unique(na.exclude(x))))
exp_num$PIR=rowSums(exp_num[,c('PIR.x','PIR.y')],na.rm = TRUE)
exp_num$Prom=rowSums(exp_num[,c('Prom.x','Prom.y')],na.rm = TRUE)
exp_num = exp_num[,10:13]
rm(pir,prom)


scaled = allExp
scaled[is.na(scaled)]=0
for(i in 2:ncol(scaled)){
  scaled[,i] = scales::rescale(x = log2(unlist(allExp[,i])),to = c(0, 100))
}
scaled = scaled[-which(duplicated(scaled$gene_name)),]
rownames(scaled)=scaled$gene_name
scaled=scaled[which(scaled$gene_name %in% unique(allv2g$gene_name)),]
tmp = melt(scaled,id.vars = 'gene_name')
colnames(tmp)[2]="id"
tmp = merge(tmp,metaData[,c('id','Celltype')],by="id")
tmp = tmp[which(tmp$Celltype %in% allstat$Celltype),-1]
scaled = dcast(gene_name~Celltype,data=tmp,value.var = 'value',fun.aggregate = function(x) mean(x))
tmp = melt(scaled,id.vars = 'gene_name',variable.name = "Celltype")

exp_num = merge(exp_num,tmp,by=c("gene_name","Celltype"),all.x=TRUE)
colnames(exp_num)[ncol(exp_num)]="Scaled\nexp.level"
exp_num$gene_name=factor(exp_num$gene_name,levels = levels(allv2g$gene_name))
exp_num$`Scaled\nexp.level`=ceiling(exp_num$`Scaled\nexp.level`)
exp_num = merge(exp_num,allv2g[,c('gene_name','chr')],by="gene_name")
exp_num=unique(exp_num)
exp_num = merge(exp_num,unique(allstat[,c('Celltype','Name','system')]),by="Celltype")
exp_num=unique(exp_num)
exp_num$chr = factor(exp_num$chr,levels = unique(allv2g$chr))
exp_num$`Scaled\nexp.level`[which(is.na(exp_num$`Scaled\nexp.level`))]=0
exp_num$`Scaled\nexp.level`[which(is.infinite(exp_num$`Scaled\nexp.level`))]=0
exp_num$color=paste(exp_num$PIR,exp_num$Prom,sep=",")
exp_num = exp_num[order(exp_num$PIR,exp_num$Prom,decreasing = F),]
exp_num$shape=2
exp_num$shape[which(exp_num$PIR==0)]=15
exp_num$shape[which(exp_num$Prom==0)]=19
exp_num$shape[which(exp_num$PIR!=0 & exp_num$Prom !=0)]=18
exp_num$shape[which(exp_num$`Scaled\nexp.level`==0)] = 2

mycol=sample(paletteer_d("ggsci::default_igv"),size = length(unique(exp_num$color)),replace = F)
names(mycol)=paste(unique(exp_num$color))
for(i in 1:nrow(exp_num)){
  exp_num$mycolor[i]=mycol[[exp_num$color[i]]] 
}
exp_num=merge(exp_num,unique(allv2g[,c('gene_name','Name','locus')]),by=c("gene_name","Name"))
exp_num$locus = factor(exp_num$locus, levels = levels(proxy$locus))

saveRDS(exp_num,file="exp_num.rds")