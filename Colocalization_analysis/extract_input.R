dir="/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs/traits/"
args = commandArgs(trailingOnly=TRUE)
library(vroom)
num=args[1]
trait=args[2]
source(paste0(dir,trait,"/QTL_config.R"))
trait_region = vroom(file=traitFilePath, col_names=TRUE)
print("trait input file successfully loaded")
trait_region[[trait_CHRcol]] <- as.integer(gsub('[a-zA-Z]', '', trait_region[[trait_CHRcol]]))
trait_region[[trait_BPcol]] <- as.numeric(trait_region[[trait_BPcol]])
print(head(trait_region))
loci=list.files(paste0("/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs/position_loci/",num,"loci/"),pattern=".input")
loci=sub(".input","",loci)

for(i in 1:length(loci)){
    locus= loci[i]
    print(locus)
    chrom=as.integer(gsub('[a-zA-Z]', '',paste(unlist(strsplit(locus,split="\\.")[[1]][[1]]))))
    colocStart=as.numeric(paste(unlist(strsplit(locus,split="\\.")[[1]][[2]])))
    colocStop=as.numeric(paste(unlist(strsplit(locus,split="\\.")[[1]][[3]])))
    rsid = read.table(paste0("/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs/position_loci/",num,"loci/",locus,".input"),header=F)
    rsid = paste(rsid$V1)

    reg = which(trait_region[[trait_CHRcol]] == chrom & trait_region[[trait_BPcol]] >= colocStart & trait_region[[trait_BPcol]] <= colocStop)
    rs = which(trait_region[[trait_SNPcol]] %in% rsid)

    retrieve = union(reg,rs)
	trait_reg = trait_region[retrieve,which(colnames(trait_region) %in% c(trait_CHRcol,trait_BPcol,trait_A1col,trait_A2col,trait_SNPcol,trait_Pcol,trait_Ncol,trait_MAFcol,trait_BETAcol,trait_SEcol))]
    print(head(trait_reg))
    if(any(is.na(trait_reg[[trait_SNPcol]]))){
		trait_reg = trait_reg[-which(is.na(trait_reg[[trait_SNPcol]])),]
	}
	if(length(which(trait_reg[[trait_SNPcol]]==""))>0){
		trait_reg = trait_reg[-which(trait_reg[[trait_SNPcol]]==""),]
	}
	if(length(which(duplicated(trait_reg[[trait_SNPcol]])))>0){
		trait_reg = trait_reg[-which(duplicated(trait_reg[[trait_SNPcol]])),]
	}
    write.table(trait_reg,file=paste0(dir,trait,"/",locus,".txt"),col.names=T,row.names=F,sep="\t",quote=F)

}


