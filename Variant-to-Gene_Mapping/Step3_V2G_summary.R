trait_ord = readRDS("trait_ord.rds")
allstat=readRDS("allstat.rds")

allv2g = NULL
for(t in 1:length(trait_ord)){
  for(i in 1:nrow(allstat)){
    tmp =read.delim(paste0("v2g/",trait_ord[t],".topLD_v2g/",trait_ord[t],".topLD_",allstat$Celltype[i],".cis_v2g.txt"),header=T)
    if(nrow(tmp)!=0){
      tmp$trait=trait_ord[t]
    }
    colnames(tmp)[11]="Celltype"
    allv2g = rbind(allv2g,tmp)
  }
}

rm(tmp,i,t)
prom=read.table("gencode.v40.promoter.bed",header=F)
colnames(prom)=c('pro_chr','pro_start','pro_end','pro_anno','gene_id','strand')
prom=prom[order(as.numeric(paste(gsub("chr","",prom$pro_chr))),as.numeric(prom$pro_start),decreasing = F),]
allv2g$gene_id=factor(allv2g$gene_id,level=unique(prom$gene_id)) 
allv2g = allv2g[order(allv2g$gene_id,decreasing = FALSE),]


tmp=unlist(strsplit(allv2g$variant_id,split = "\\|"))
allv2g$index_rsid = sub("index_rsid:","",paste(tmp[grep("index",tmp)]))
allv2g$proxy_rsid = sub("proxy_rsid:","",paste(tmp[grep("proxy",tmp)]))

hla=read.table("hla.txt",header=F)
hla=paste(hla$V1)
files=list.files(pattern = "_filter.txt")
filter=NULL # filter out suggestive variants
for(i in 1:length(files)){
  tmp=read.delim(file=files[i],header=F)
  tmp=paste(tmp$V1)
  filter = c(filter,tmp)
}
filter = unique(filter)

if(any(c(filter,hla)%in% allv2g$proxy_rsid)){allv2g=allv2g[-which(allv2g$proxy_rsid %in% c(filter,hla)),]}
rm(files,filter,hla,i,tmp,prom)

allv2g = merge(allv2g,allstat[,c('Celltype','Name','system',"Cellline")],by="Celltype")
tmp=unlist(strsplit(paste(allv2g$gene_id),split="\\."))
allv2g$ensembl_gene_id = paste(tmp[seq(1,length(tmp),by=2)])
tmp=unlist(strsplit(paste(allv2g$tx_id),split="\\."))
allv2g$ensembl_tx_id = paste(tmp[seq(1,length(tmp),by=2)])


genes=dcast(gene_name~trait,data=unique(allv2g[,c('gene_name','Name','trait')]),
            value.var = "Name" ,fun.aggregate = function(x) paste(unique(x),collapse = ",") )
genes$sum=apply(genes[,-1],1,function(x) length(which(x!="")))
write.table(genes,file="total_genes.txt",col.names = T,row.names = T,sep="\t",quote=F)
for(i in 1:length(trait_ord)){
  id = genes$gene_name[which(genes$sum==1 & genes[,which(colnames(genes)==trait_ord[i])]!="")]
  gene=dcast(gene_name~Name,data=unique(allv2g[which(allv2g$trait==trait_ord[i]),c('proxy_rsid','gene_name','Name')]),
              value.var = "proxy_rsid",fun.aggregate = function(x) paste(unique(x),collapse = ", "))
  gene$Num.cell=apply(gene[,-1],1,function(x) length(which(x!="")))
  gene = gene[which(gene$gene_name %in% id),]
  write.table(gene,file=paste("unique_genes",trait_ord[i],".txt"),col.names = T,row.names = F,sep="\t",quote=F)
  
}

genes=dcast(gene_name~Name,data=unique(allv2g[,c('gene_name','Name','trait')]),
            value.var = "trait" ,fun.aggregate = function(x) paste(unique(x),collapse = ",") )
genes$sum=apply(genes[,-1],1,function(x) length(which(x!="")))
write.table(genes,file="total_genes_across_cell.txt",col.names = T,row.names = F,sep="\t",quote=F)

immune=allstat$Name[which(allstat$system=="immune")]
for(i in 1:length(immune)){
  id = genes$gene_name[which(genes$sum==1 & genes[,which(colnames(genes)==immune[i])]!="")]
  gene=dcast(gene_name~trait,data=unique(allv2g[which(allv2g$Name==immune[i]),c('proxy_rsid','gene_name','trait')]),
             value.var = "proxy_rsid",fun.aggregate = function(x) paste(unique(x),collapse = ", "))
  gene$Num.trait=apply(gene[,-1],1,function(x) length(which(x!="")))
  gene = gene[which(gene$gene_name %in% id),]
  write.table(gene,file=paste("unique_genes",immune[i],".txt"),col.names = T,row.names = F,sep="\t",quote=F)
  
}

saveRDS(allv2g,file="allv2g.rds")
