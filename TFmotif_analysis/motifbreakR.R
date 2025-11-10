library(motifbreakR)
library(BSgenome)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(Biostrings)
library(XML)
allv2g=readRDS("allv2g.rds")
allstat=readRDS("allstat.rds")

# create snp.mb object for motifbreakR #######
all = read.delim("all_proxies.bed",header=F,sep="\t")
colnames(all)=c("chr","start","end","rsid")
snps.mb = makeGRangesFromDataFrame(all,starts.in.df.are.0based = TRUE,seqnames.field = "chr",
                                   seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))
change.to.search.genome <- function(granges.object, search.genome) {
  sequence <- seqlevels(granges.object)
  ## sequence is in UCSC format and we want NCBI style
  newStyle <- mapSeqlevels(sequence,seqlevelsStyle(search.genome))
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  ## rename the seqlevels
  granges.object <- renameSeqlevels(granges.object,newStyle)
  seqlevels(granges.object) <- seqlevelsInUse(granges.object)
  seqinfo(granges.object) <- keepSeqlevels(seqinfo(search.genome),
                                           value = seqlevelsInUse(granges.object))
  return(granges.object)
}
attributes(snps.mb)$genome.package <- attributes(BSgenome.Hsapiens.UCSC.hg38)$pkgname
snps.mb = change.to.search.genome(snps.mb,BSgenome.Hsapiens.UCSC.hg38)

seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "NCBI"
snps.mb = change.to.search.genome(snps.mb,BSgenome.Hsapiens.UCSC.hg38)

genome(BSgenome.Hsapiens.UCSC.hg38) <- genome(SNPlocs.Hsapiens.dbSNP155.GRCh38)
snps.mb = change.to.search.genome(snps.mb,BSgenome.Hsapiens.UCSC.hg38)
a=snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38,snps.mb,genome=BSgenome.Hsapiens.UCSC.hg38)
b=as.data.frame(a)
b$start=b$pos-1
colnames(b)[7]="REF"
colnames(b)[8]="ALT"
b=b[,-c(5,6)]
colnames(b)[4]="SNP_id"
b$n=lengths(b$ALT)
a=b[which(b$n==1),]
b=b[which(b$n>1),]
c=NULL
for(i in 1:nrow(b)){
  alt=unlist(b$ALT[i])
  tmp=b[rep(i,b$n[i]),]
  tmp$ALT = paste(alt)
  c=rbind(c,tmp)
}
a=rbind(a,c)
rm(SNPlocs.Hsapiens.dbSNP155.GRCh38,all,alt,i,tmp,b,c,BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg38)
a$seqnames = paste0("chr",a$seqnames)
snps.mb = makeGRangesFromDataFrame(a,keep.extra.columns = TRUE,end.field = "pos",starts.in.df.are.0based = TRUE,seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))
attributes(snps.mb)$genome.package <- attributes(BSgenome.Hsapiens.UCSC.hg38)$pkgname
snps.mb = change.to.search.genome(snps.mb,BSgenome.Hsapiens.UCSC.hg38)
snps.mb$REF=DNAStringSet(paste(snps.mb$REF))
snps.mb$ALT=DNAStringSet(paste(snps.mb$ALT))

# run motifbreakR on only V2G SNPs ####
jaspar = subset(MotifDb,dataSource == 'JASPAR_CORE')
data(motifbreakR_motif)
snps.mb = snps.mb[which(snps.mb$SNP_id %in% allv2g$proxy_rsid)]
mb <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                  pwmList = motifbreakR_motif,
                  threshold = 0.005,
                  method = "ic",
                  bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                  BPPARAM = BiocParallel::SerialParam())
saveRDS(mb,file="allSNPs_motifBreakR_mb_v2g.rds")
mb <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                  pwmList = jaspar,
                  threshold = 0.005,
                  method = "ic",
                  bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                  BPPARAM = BiocParallel::SerialParam())
saveRDS(mb,file="allSNPs_motifBreakR_jaspar_v2g.rds")