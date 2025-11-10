options(scipen=100)

for( trait in c("ECZ", "ALG", "UC" , "IBD", "CRO", "AS",  "HT",  "MS" , "T1D", "CEL", "JIA", "PSO", "RA" , "SLE" ,"GD",  "VIT")){
        alg=read.delim(paste0("download/",trait,".tsv"),header=T,sep="\t")
        write.table(paste0("rs",alg$SNP_ID_CURRENT),file=paste0(trait,"_index.txt"),col.names=F,row.names=F,quote=F)
        # Run expansion by topLD
        # topld_api -inFile ${trait}_index.txt -outputInfo topLD/$trait.index -outputLD topLD/$trait.proxy -pop "EUR"  -thres 0.8
        ld=read.delim(paste0("topLD/",trait,".proxy"),header=T,sep="\t")
        ld$geneName=apply(ld,1,function(x) paste(unique(unlist(strsplit(paste(x[19]),split="\\|"))),collapse=",") )
        proxy=data.frame(chr=paste0("chr",ld$CHROM),start=ld$POS2.gb38. -1, end =ld$POS2.gb38.,index=paste0("locus:",ld$geneName,"|index_rsid:",ld$rsID1,"|proxy_rsid:",ld$rsID2))
        alg$geneName=gsub(" ","",alg$REPORTED.GENE.S.)
        tmp=data.frame(chr=paste0("chr",alg$CHR_ID),start=alg$CHR_POS - 1,end=alg$CHR_POS,index=paste0("locus:",alg$geneName,"|index_rsid:rs",alg$SNP_ID_CURRENT,"|proxy_rsid:rs",alg$SNP_ID_CURRENT))
        proxy=rbind(proxy,tmp)
        snps=setdiff(paste0("rs",alg$SNP_ID_CURRENT),unique(ld$rsID1))
        
        library(LDlinkR)
        setwd("./LDlinkRtmp/")
        LDproxy_batch(snps,token = Sys.getenv("LDLINK_TOKEN"),genome_build = "grch38_high_coverage",pop="EUR")
        setwd("../")
        for(i in 1:length(snps)){
                file=list.files("./LDlinkRtmp/",pattern=snps[i])
                if(length(file)==1){
                        tmp=read.table(paste0("./LDlinkRtmp/",file),header=T)
                        tmp=tmp[which(round(tmp$R2,1) >= 0.8),]
                        pos=unlist(strsplit(tmp$Coord,split=":"))
                        tmp$chr=paste(pos[seq(1,by=2,length(pos))])
                        tmp$end=as.numeric(paste(pos[seq(2,by=2,length(pos))]))
                        tmp$start = tmp$end -1
                        tmp$index=paste0("index_rsid:",snps[i],"|proxy_rsid:",tmp$RS_Number)
                        proxy=rbind(proxy,tmp[,c('chr','start','end','index')])
                }
        }
        proxy=proxy[order(as.numeric(gsub("chr","",proxy$chr)),as.numeric(proxy$start)),]
        write.table(proxy,file=paste0(trait,"_proxies.bed"),col.names=F,row.names=F,sep="\t",quote=F)
}

