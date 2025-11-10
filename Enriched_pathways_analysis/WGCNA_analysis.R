library(DESeq2)
library(scales)
library(WGCNA)
library(tidyverse)
library(magrittr) 

dir="pathway/"
metaData=readRDS("metaData.rds")
allExp=readRDS("allExp.rds")
allv2g=readRDS("allv2g.rds")

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

ddsMat <- DESeqDataSetFromMatrix(countData = scaled,colData=metaData, design=~Celltype+system,ignoreRank = T)
dds <- estimateSizeFactors(ddsMat)
dds <- estimateDispersionsGeneEst(dds,)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- DESeq(dds)

wpn_vsd <- getVarianceStabilizedData(dds)
expr_normalized <- wpn_vsd[ which(rownames(wpn_vsd) %in% unique(allv2g$gene_name)),] # only play with V2G genes
expr_normalized = expr_normalized[,paste(metaData$id)]

input_mat = t(expr_normalized)
allowWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 10 # or 10 or 18
temp_cor <- cor       
cor <- WGCNA::cor 

netwk <- blockwiseModules(input_mat, 
                          power = 9,                # <= power here
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.999,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          numericLabels = F,
                          verbose = 3)
cor <- temp_cor
mergedColors = netwk$colors

module_df <- data.frame(
  gene_name = names(netwk$colors),
  module = netwk$colors
)

saveRDS(module_df,file="WGCNA_module_df.rds")
saveRDS(netwk,file"WGCNA_netwk.rds)