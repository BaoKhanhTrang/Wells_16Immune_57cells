args = commandArgs(trailingOnly=TRUE)

dir="eQTLs"
num="258"
dbname=args[1]
#dbname="DICE"
#dbname="OneK1K"
#dbname="eQTL_catalogue"

loci=list.files(paste0(dir,"/position_loci/",num,"loci/"),pattern=".input")
loci = gsub(".input","",loci)
cells=readLines(paste0(dir,"/databases/",dbname,"/cells.txt"))

v2g=readRDS(paste0(dir,"/position_loci/",num,"loci_to_4125genes.rds"))
allv2g=readRDS(paste0(dir,"/allv2g.rds"))
allv2g = allv2g[which(allv2g$system=="immune"),]

tss=read.table("gencode_v40/gencode.v40.TSS.bed")
colnames(tss)=c("gene_chom","start","gene_tss","tx_id","gene_id","strand")
tss$ensembl_gene_id = paste(tss$gene_id)
tmp = unlist(strsplit(paste(grep("ENSG|\\.",tss$gene_id,value = TRUE)),split="\\."))
tss$ensembl_gene_id[grep("ENSG|\\.",tss$gene_id)]=paste(tmp[seq(1,length(tmp),by=2)])



match=read.table(paste0(dir,"/matching_celltypes.txt"),header=T, sep="\t")
match=match[which(match$source==dbname),]
done = readRDS(paste0(dir,"/position_loci/",num,"loci.rds"))
eGenes_ld=NULL
for(i in 1:nrow(done)){
  locus = unique(allv2g[which(allv2g$index_rsid %in% unlist(done$rsid[i])),c('Celltype','ensembl_gene_id','trait')])
  locus$locus = done$id[i]
  locus = merge(locus,match[,c('Celltype','filename')],by="Celltype",all.x=TRUE)

  cells = unique(na.exclude(locus$filename))
  if(length(cells)>0){
    for(e in 1:length(cells)){
      a=NULL
      tmp = readLines(paste0(dir,"/genes_locus/",dbname,"/",cells[e],".",done$id[i],".eGenes"))
      if(length(tmp)>0){
        ensembl_gene_id = tmp
        tmp2 = unlist(strsplit(paste(grep("ENSG|\\.",tmp,value = TRUE)),split="\\."))
        ensembl_gene_id[grep("ENSG|\\.",tmp)] = paste(tmp2[seq(1,length(tmp2),by=2)])
        a = data.frame(locus=done$id[i],filename=cells[e],ensembl_gene_id=ensembl_gene_id)
        a$shared = ifelse(a$ensembl_gene_id %in% locus$ensembl_gene_id[which(locus$locus == done$id[i] & locus$filename == cells[e])],"shared",paste(cells[e]))  
        a = merge(a,locus[which(locus$locus == done$id[i] & locus$filename == cells[e]),c("locus","filename","ensembl_gene_id","trait" )],by=c("locus","filename","ensembl_gene_id" ),all=TRUE)
        a$shared[which(is.na(a$shared))] = "V2G"
      }
      if(!is.null(a)){

      tmp = readLines(paste0(dir,"/genes_locus/",dbname,"/",cells[e],".",done$id[i],".outsideLD.eGenes"))
      if(length(tmp)>0){
        tmp = read.table(paste0(dir,"/genes_locus/",dbname,"/",cells[e],".",done$id[i],".outsideLD.eGenes"),sep="\t")
        ensembl_gene_id = tmp$V1
        tmp2 = unlist(strsplit(paste(grep("ENSG|\\.",tmp$V1,value = TRUE)),split="\\."))
        ensembl_gene_id[grep("ENSG|\\.",tmp$V1)] = paste(tmp2[seq(1,length(tmp2),by=2)])
        b = data.frame(locus=done$id[i],filename=cells[e],ensembl_gene_id=ensembl_gene_id)
        a = merge(a,b,by=c("locus","filename","ensembl_gene_id" ),all=TRUE) 
        a$shared[which(is.na(a$shared))] = "outsideLD"
      }
      }

      if(!is.null(a)){
        gene=tss[which(tss$ensembl_gene_id %in% a$ensembl_gene_id),]
        a$cis="trans"
        transgene=unique(gene$ensembl_gene_id[which(gene$gene_chom != done$seqnames[i])])
        a$cis[which(a$ensembl_gene_id %in% transgene)] = "diff_chrom"

        start = done$start[i] - 5000000
        end = done$end[i] + 5000000
        transgene=unique(gene$ensembl_gene_id[which(gene$gene_chom == done$seqnames[i]  & gene$gene_tss >= start & gene$gene_tss<=end)])
        a$cis[which(a$ensembl_gene_id %in% transgene)] = "cis"

        eGenes_ld = rbind(eGenes_ld,unique(a))

      }
      
    }
  }
}
saveRDS(eGenes_ld,file=paste0(dir,"/genes_locus/",dbname,".eGenes_ld.rds"))
