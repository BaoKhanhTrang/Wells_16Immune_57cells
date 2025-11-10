args = commandArgs(trailingOnly=TRUE)
library(vroom)
library(tidyr)
library(dplyr)
library(coloc)
library(foreach)
library(doParallel)
cores <- detectCores()
registerDoParallel(cores)
options(timeout=1000000000)


locus = args[1]
num = args[2]
eqtl = args[3]
N = as.numeric(paste(args[4]))
main_dir="eQTLs/"

chrom=as.integer(gsub('[a-zA-Z]', '',paste(unlist(strsplit(locus,split="\\.")[[1]][[1]]))))
colocStart=as.numeric(paste(unlist(strsplit(locus,split="\\.")[[1]][[2]])))
colocStop=as.numeric(paste(unlist(strsplit(locus,split="\\.")[[1]][[3]])))

rsid = read.table(paste0(main_dir,"position_loci/",num,"loci/",locus,".input"),header=F)
rsid = rsid$V1

genes = readRDS(paste0(main_dir,"position_loci/",num,"loci_to_4125genes.rds"))
genes = genes[genes$locus == locus,]
genes = genes$ensembl_gene_id

cells = read.table(paste0(main_dir,"databases/",eqtl,"/cells.txt"),header=F)
cells = cells$V1
header = read.table(paste0(main_dir,"databases/",eqtl,"/header"),header=F)
header = header$V1
traits = c("ALG" ,"AS", "CEL", "CRO", "ECZ" ,"GD" ,"HT", "IBD", "JIA", "MS" ,"PSO", "RA", "SLE", "T1D", "UC" ,"VIT")

# read traits input 
result <- foreach(i= 1:length(traits)) %dopar% {
    tryCatch({
        trait = traits[i]
        log_file <- paste0(main_dir,"log/runeQTL.",locus,".",eqtl,".", trait, ".out")
        sink(log_file)

        print(trait)
        source(paste0(main_dir,"traits/",trait,"/QTL_config.R"))
        trait_reg = vroom(file=paste0(main_dir,"traits/",trait,"/",locus,".txt"),col_names=T)
        
        if(nrow(trait_reg)>0){
            if(file.exists(paste0(main_dir,"traits/",trait,"/",locus,".maf")) & length(readLines(paste0(main_dir,"traits/",trait,"/",locus,".maf")))>0){
                #rectify MAF
                MAF = vroom(paste0(main_dir,"traits/",trait,"/",locus,".maf"),col_names=F)
                colnames(MAF)=c("source","allele","maf","rsid")
                MAF <- MAF %>%
                    filter(maf > 0, maf < 1) %>%
                    drop_na(maf)

                # Perform vectorized matching
                complement_map <- c(A = "T", T = "A", G = "C", C = "G")
                trait_reg$minor_allele_freq=NA
                for(e in 1:nrow(trait_reg)){
                    maf = MAF[which(MAF$rsid == trait_reg[[trait_SNPcol]][e] & MAF$allele == trait_reg[[trait_A1col]][e]  ),]
                    if(nrow(maf)==0){
                        maf = MAF[which(MAF$rsid == trait_reg[[trait_SNPcol]][e] & MAF$allele == trait_reg[[trait_A2col]][e]  ),]
                    }
                    if(nrow(maf)==0){
                        complement <- paste(complement_map[paste(trait_reg[[trait_A1col]][e])])
                        maf = MAF[which(MAF$rsid == trait_reg[[trait_SNPcol]][e] & MAF$allele == complement  ),]
                    }
                    if(nrow(maf)==0){
                        complement <- paste(complement_map[paste(trait_reg[[trait_A2col]][e])])
                        maf = MAF[which(MAF$rsid == trait_reg[[trait_SNPcol]][e] & MAF$allele == complement  ),]
                    }
                    if(nrow(maf)>0){
                        trait_reg$minor_allele_freq[e] = maf$maf[1]
                    }
                }
            }
            if(trait_MAFcol!=""){
                if(any(trait_reg[[trait_MAFcol]] <=0 | trait_reg[[trait_MAFcol]] >= 1 )){
                    is_outside_range <- which(trait_reg[[trait_MAFcol]] <=0 | trait_reg[[trait_MAFcol]] >= 1)
                    trait_reg[[trait_MAFcol]][is_outside_range] = trait_reg$minor_allele_freq[is_outside_range]
                }
                trait_reg$maf = trait_reg[[trait_MAFcol]]
            }else{
                if("minor_allele_freq" %in% colnames(trait_reg)){trait_reg$maf = trait_reg$minor_allele_freq}
            }
            if(any(is.na(trait_reg$maf))){
                trait_reg = trait_reg[-which(is.na(trait_reg$maf)),]
            }
            
            print(summary(trait_reg$maf))

            if(!trait_Ncol %in% colnames(trait_reg)){
                trait_Ncol=as.numeric(paste(trait_Ncol))
            }
            #remove any alphabetical characters from the chromosome column
            trait_reg[[trait_CHRcol]] <- as.character(gsub('[a-zA-Z]', '', trait_reg[[trait_CHRcol]]))
            trait_reg[[trait_BPcol]] <- as.numeric(trait_reg[[trait_BPcol]])
            trait_reg = unique(trait_reg)

            dataset1=list(snp=trait_reg[[trait_SNPcol]],type=traitType,s=traitProp)
            if(trait_Ncol %in% colnames(trait_reg)){
                dataset1$N=trait_reg[[trait_Ncol]]
            }else{
                if (is.numeric(trait_Ncol)) {
                    dataset1$N = trait_Ncol
                } else {
                    stop("trait_Ncol must be a numeric value")
                }
            }
            if(trait_Pcol!=""){
                dataset1$pvalues=trait_reg[[trait_Pcol]]
                if(any(is.na(dataset1$pvalues))){dataset1$pvalues[is.na(dataset1$pvalues)]=1}
            }
            if(trait_BETAcol!=""){
                if(any( is.na(trait_reg[[trait_BETAcol]]) | trait_reg[[trait_BETAcol]]==0 | is.na(trait_reg[[trait_SEcol]]) | trait_reg[[trait_SEcol]]==0) ){
                    minbeta = min(abs(trait_reg[[trait_BETAcol]][which(trait_reg[[trait_BETAcol]]!=0)]),na.rm=T)
                    minse = min(abs(trait_reg[[trait_SEcol]][which(trait_reg[[trait_SEcol]]!=0)]),na.rm=T)
                    trait_reg[[trait_SEcol]][which(trait_reg[[trait_SEcol]]==0 | is.na(trait_reg[[trait_SEcol]]))] = minse
                    trait_reg[[trait_BETAcol]][which(trait_reg[[trait_BETAcol]]==0 | is.na(trait_reg[[trait_BETAcol]]))] = minbeta
                }
                dataset1$beta=trait_reg[[trait_BETAcol]]
            }
            if(trait_SEcol!=""){
                dataset1$varbeta=trait_reg[[trait_SEcol]]^2
            }
            dataset1$MAF=trait_reg$maf
            check_dataset(dataset1)

            print("Finish dataset1")

            for(e in 1:length(genes)){
                print(paste("Start with ",genes[e]))

                done = list.files(path=paste0(main_dir,"output/",trait,"/",eqtl,"/"),pattern=paste0(locus,".",genes[e],"*.coloc_results_summary.txt"))
                if(length(done)>0){
                    done = gsub(genes[e],"",done)
                    done = gsub(".coloc_results_summary.txt","",done)
                    tmp = unlist(strsplit(done,split="\\."))
                    done = tmp[seq(2,by=2,length(tmp))]
                    left = setdiff(cells,done)
                }else{
                    left = cells
                }
                print(paste(length(left),"cells left"))

                for(cell in left){

                    colocInputFile = readLines(paste0(main_dir,"databases/",eqtl,"/input/",cell,".",genes[e],".txt"))
                    if(length(colocInputFile)>0){
                        colocInputFile = vroom(paste0(main_dir,"databases/",eqtl,"/input/",cell,".",genes[e],".txt"),col_names= header)

                        if("GeneSymbol" %in% header){
                            geneSymbol = unique(colocInputFile$GeneSymbol)
                        }else{
                            geneSymbol = tryCatch({clusterProfiler::bitr(genes[e], fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
				                                }, error = function(err){
						                            print("ENSEMBL ID could not be converted to HGNC Symbol")
							                        print(paste("geneSymbol will be set to the ENSEMBL ID",genes[e]))
							                        return(genes[e])
				                                })
                        }
                    
                        out_prefix = paste(locus,genes[e],geneSymbol,cell,sep=".")
                        print(paste("Input from ",cell))
                    
                        if(!file.exists(paste0(main_dir,"output/",trait,"/",eqtl,"/",out_prefix,".coloc_results_pp4_cond.txt"))){
                            if("Beta" %in% header & any(colocInputFile$Beta==0)){
                                colocInputFile$Beta[which(colocInputFile$Beta==0)] = 0.00000000000000000000000000000001
    		                }
                            if(! "varbeta" %in% header & "Beta" %in% header){
                                if("se" %in% header){
                                    colocInputFile$varbeta=colocInputFile$se^2
                                }else{
                                    Z = abs(qnorm(colocInputFile$Pvalue/2))
    		                        colocInputFile$varbeta=(colocInputFile$Beta/Z)^2
                                }
                            }
    		            
                            if(! "maf" %in% header){
                                colocInputFile = merge(colocInputFile,trait_reg[,c("maf",trait_SNPcol)],by.x="rsid",by.y=trait_SNPcol)
                                colocInputFile = unique(colocInputFile)
                                colnames(colocInputFile)[ncol(colocInputFile)]="maf"
                            }
    		            
    		                if(any(is.na(colocInputFile$maf))){
    			                colocInputFile = colocInputFile[-which(is.na(colocInputFile$maf)),]
    		                }
    		                if(any(c(0,1) %in% colocInputFile$maf)){
    			                colocInputFile = colocInputFile[-which(colocInputFile$maf %in% c(0,1)),]
    		                }
        		        
                            colocInputFile=colocInputFile[order(colocInputFile$rsid,colocInputFile$Pvalue,decreasing = T),]
        		            if(0 %in% colocInputFile$Pvalue){
                                colocInputFile$Pvalue[which(colocInputFile$Pvalue==0)] = min(colocInputFile$Pvalue[-which(colocInputFile$Pvalue==0)])
                            }
        		            if(length(which(duplicated(colocInputFile$rsid)))>0){
        		                colocInputFile = colocInputFile[-which(duplicated(colocInputFile$rsid)),]
          		            }
                        
                            if(any(is.na(colocInputFile$rsid))){
    			                colocInputFile = colocInputFile[-which(is.na(colocInputFile$rsid)),]
    		                }
                            if(any(colocInputFile$rsid =="")){
    			                colocInputFile = colocInputFile[-which(colocInputFile$rsid ==""),]
    		                }

                            if(nrow(colocInputFile)>0 &&  length(intersect(dataset1$snp,colocInputFile$rsid))>0){

                                snp_set = intersect(dataset1$snp,colocInputFile$rsid)
                                #subset_indices <- dataset1$snp %in% snp_set
                                #dataset_subset <- lapply(dataset1, function(x) if (length(x) == length(dataset1$snp)) x[subset_indices] else x)
                                colocInputFile = colocInputFile[which(colocInputFile$rsid %in% snp_set),]
    		                    if("N" %in% header){  N = colocInputFile$N }
                                if("an" %in% header){ N = colocInputFile$an/2 }
                                dataset2=list(snp=colocInputFile$rsid,pvalues=colocInputFile$Pvalue, N=N,type="quant",MAF=colocInputFile$maf)
                                if("Beta" %in% header){
                                    dataset2$beta=colocInputFile$Beta
                                    dataset2$varbeta = colocInputFile$varbeta
                                }
                                check_dataset(dataset2)
        		                print("Finish dataset2")
                            
                                coloc_results = coloc.abf(dataset1=dataset1,dataset2=dataset2)
    
                            #prepare useful outputs
                                coloc_results_summary = coloc_results$summary
                                coloc_results_full = coloc_results$results

                            #calculate pp4 / pp3 + pp4
                                PP3andPP4 = coloc_results_summary[5] + coloc_results_summary[6]
                                pp4_conditional = coloc_results_summary[6] / PP3andPP4

        		            #prep coloc output strings
        		                coloc_results_summary_outputStr = paste(out_prefix,"coloc_results_summary.txt",sep=".")
        		                coloc_results_full_outputStr = paste(out_prefix,"coloc_results_full.txt",sep=".")
        		                coloc_results_pp4_cond_outputStr = paste(out_prefix,"coloc_results_pp4_cond.txt",sep=".")
        
        		                write.table(coloc_results_summary, file=paste0(main_dir,"output/",trait,"/",eqtl,"/",coloc_results_summary_outputStr), sep="\t", row.names=TRUE, quote=FALSE)
        		                write.table(coloc_results_full, file=paste0(main_dir,"output/",trait,"/",eqtl,"/",coloc_results_full_outputStr), sep="\t", row.names=FALSE, quote=FALSE)
        		                write.table(pp4_conditional, file=paste0(main_dir,"output/",trait,"/",eqtl,"/",coloc_results_pp4_cond_outputStr), sep="\t", row.names=FALSE, quote=FALSE)
                        
                                rm(coloc_results_summary,coloc_results_full,coloc_results_summary_outputStr,coloc_results_full_outputStr,coloc_results_pp4_cond_outputStr)
                                rm(PP3andPP4,pp4_conditional,coloc_results,dataset2,colocInputFile,Z,out_prefix,geneSymbol)
                                gc()
                            }
                        }
                    }else{ print("no line in eqtl file")  }
                }
                print(paste("Done with ",genes[e]))
            }
            rm(trait_reg,dataset1)
            gc()
        }else{
            print(paste("no line in input ",trait ))
        }
        sink()
        return(NULL)
        gc()
    }, error = function(e) {
        message(paste("Error in trait:", traits[i], ":", e$message))
        NULL
    })
}
