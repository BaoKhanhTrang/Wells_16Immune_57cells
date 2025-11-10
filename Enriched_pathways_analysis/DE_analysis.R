args = commandArgs(trailingOnly=TRUE)
cell=args[1]

library(DESeq2)
library(scales)

dir="pathway/"
metaData=readRDS("metaData.rds")
allExp=readRDS("allExp.rds")

scaled = allExp
# rescale concatinated values from single cell RNA-seq data to match bulk data
for(e in grep("klaus",colnames(scaled))){
  scaled[,e] = scales::rescale(x = allExp[,e],to = c(0, mean(apply(allExp[,-c(1,2)],2,function(x) max(x,na.rm=TRUE)))))
}
scaled=scaled[,-2]
scaled = scaled[-which(duplicated(scaled$gene_id)),]
rownames(scaled)=paste(scaled$gene_id)
scaled[is.na(scaled)]=0
scaled = scaled[-which(rowSums(scaled[,-1])==0),]
scaled[,-1] = round(scaled[,-1]+1)
scaled = scaled[,c('gene_id',metaData$id)]
metaData$Celltype=paste(metaData$Celltype)


md = metaData
md$Celltype[which(md$Celltype != cell)]="other"
scl = scaled
ddsMat <- DESeqDataSetFromMatrix(scaled,colData=md, design=~Celltype+system,tidy = TRUE)
dds <- estimateSizeFactors(ddsMat)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- DESeq(dds)
saveRDS(dds,file=paste0(dir,"DEanalysis/",cell,"_overAll.rds"))

# shrink logFC
a=lfcShrink(dds,coef=2,apeAdapt=T,apeMethod="nbinomC")
saveRDS(a,file=paste0(dir,"DEanalysis/",cell,"_overAll.apeglm.rds"))

