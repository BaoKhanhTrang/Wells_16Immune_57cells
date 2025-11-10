.libPaths( c("/home/trangk/R/x86_64-pc-linux-gnu-library/4.4" , .libPaths() ) )
args = commandArgs(trailingOnly=TRUE)
library(foreach)
library(doParallel)
library(vroom)
cores <- detectCores()
registerDoParallel(cores)
options(timeout=1000000000)

num = args[1]
locus = args[2]
main_dir="/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs/"

chrom=paste(unlist(strsplit(locus,split="\\.")[[1]][[1]]))

MAF <- paste0(main_dir,"MAF/",chrom,".maf")


traits = c("ALG" ,"AS", "CEL", "CRO", "ECZ" ,"GD" ,"HT", "IBD", "JIA", "MS" ,"PSO", "RA", "SLE", "T1D", "UC" ,"VIT")
result <- foreach(i= 1:length(traits)) %dopar% {
    log_file <- paste0(main_dir,"log/extractMAFinput.",locus,".job", traits[i], ".out")
    sink(log_file)
    trait = traits[i]
    print(trait)
    source(paste0(main_dir,"traits/",trait,"/QTL_config.R"))

    trait_reg = vroom(file=paste0(main_dir,"traits/",trait,"/",locus,".txt"),col_names=T)
    if(nrow(trait_reg)>0){
        tmpfile=paste0("/scr1/users/trangk/",trait,".",locus,".txt")
        write.table(trait_reg[[trait_SNPcol]],file=tmpfile,col.names=F,row.names=F,quote=F)
        output=paste0(main_dir,"traits/",trait,"/",locus,".maf")
        system(sprintf("grep -wf %s %s > %s",tmpfile,MAF,output))
    }
}
registerDoSEQ()
doParallel::stopImplicitCluster()
                  