args = commandArgs(trailingOnly=TRUE)
cell=args[1]
trait=args[2]

library(pathfindR)

dir="pathway/"
allv2g = readRDS("allv2g.rds")

file=list.files(path=paste0(dir,"DEanalysis/"),pattern=paste0(cell,"_overAll.apeglm.rds"))
a=readRDS(paste0(dir,"DEanalysis/",file))

# make sure the fold change is cell type/other, it not, reverse the log
check=strsplit(sub("Wald test p-value: Celltype ","",a@elementMetadata$description[4]),split=" vs ")[[1]][[1]]
if(check=="other"){
        a$tmp = 1/(2^(a$log2FoldChange))
        a$logFC = log2(a$tmp)
}else{
        colnames(a)[2]="logFC"
}
a=as.data.frame(a)
a$ensembl_gene_id=rownames(a)

# remove NA p-values
if(any(is.na(a$pvalue))){a=a[-which(is.na(a$pvalue)),]}

# run pathfindR 
genes = unique(allv2g[which(allv2g$Celltype==paste(cell) & allv2g$trait==trait_ord[t]),c('ensembl_gene_id','gene_name')])
if(nrow(genes)>=2){
    dir.create(paste0("findpathR/",cell,"_in_",trait_ord[t]),recursive=TRUE)
    input = merge(a[,c('ensembl_gene_id','logFC','pvalue')],genes,by="ensembl_gene_id")
    input = unique(input[,c('gene_name','logFC','pvalue')])
    
    # save processed-input for visualization
    input_processed=NULL
    res <- tryCatch({inp1 = input_processing(input,p_val_threshold = 1,pin_name_path = "Biogrid")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp1")){input_processed=rbind(input_processed,inp1)}

    res <- tryCatch({inp2 = input_processing(input,p_val_threshold = 1,pin_name_path = "KEGG")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp2")){input_processed=rbind(input_processed,inp2)}

    res <- tryCatch({inp3 = input_processing(input,p_val_threshold = 1,pin_name_path = "STRING")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp3")){input_processed=rbind(input_processed,inp3)}


     res <- tryCatch({inp4 = input_processing(input,p_val_threshold = 1,pin_name_path = "GeneMania")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp4")){input_processed=rbind(input_processed,inp4)}


    res <- tryCatch({inp5 = input_processing(input,p_val_threshold = 1,pin_name_path = "IntAct")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp5")){input_processed=rbind(input_processed,inp5)}
    saveRDS(input_processed,"input_processed.rds")

    # run with different PIP db

    setwd(paste0("findpathR/",cell,"_in_",trait_ord[t]))
      tryCatch({combined2 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,gene_sets = "KEGG",output_dir="./Biogrid",pin_name_path="Biogrid")},
               error = function(e) {return(NULL)},
               warning = function(w) {return(NULL)})
      if(!exists("combined2")){combined2=NULL}

      tryCatch({res22 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,gene_sets = "KEGG",output_dir="./KEGG",pin_name_path="KEGG")},
               error = function(e) {return(NULL)},
               warning = function(w) {return(NULL)})
      if(exists("res22")){
        if(!is.null(combined2)){
          combined2 = combine_pathfindR_results(combined2,res22,plot_common = FALSE)
        }else{
          combined2 = res22
        }
      }
    
      tryCatch({res32 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,gene_sets = "KEGG",output_dir="./STRING",pin_name_path="STRING")},
               error = function(e) {return(NULL)},
               warning = function(w) {return(NULL)})
      if(exists("res32")){
        if(!is.null(combined2)){
          combined2 = combine_pathfindR_results(combined2,res32,plot_common = FALSE)
        }else{
          combined = res32
        }
      }
    
      tryCatch({res42 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,gene_sets = "KEGG",output_dir="./GeneMania",pin_name_path="GeneMania")},
               error = function(e) {return(NULL)},
               warning = function(w) {return(NULL)})
      if(exists("res42")){
        if(!is.null(combined2)){
          combined2 = combine_pathfindR_results(combined2,res42,plot_common = FALSE)
        }else{
          combined2 = res42
        }
      }
    
      tryCatch({res52 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,gene_sets = "KEGG",output_dir="./IntAct",pin_name_path="IntAct")},
               error = function(e) {return(NULL)},
               warning = function(w) {return(NULL)})
      if(exists("res52")){
        if(!is.null(combined2)){
          combined2 = combine_pathfindR_results(combined2,res52,plot_common = FALSE)
        }else{
          combined2 = res52
        }
      }
    
    if(exists("combined2")){
        if(!is.null(combined2)){
          final = combined2[,1:2]
          final$Fold_Enrichment=0
          final$lowest_p=0
          final$highest_p=0
          final$Up_regulated=""
          final$Down_regulated=""
          for(r in 1:nrow(combined2)){
            final$Fold_Enrichment[r] = max(unlist(combined2[r,grep("Fold_Enrichment",colnames(combined2))]),na.rm = T)
            final$lowest_p[r] = min(unlist(combined[r,grep("lowest_p",colnames(combined))]),na.rm = T)
            final$highest_p[r] = min(unlist(combined[r,grep("highest_p",colnames(combined))]),na.rm = T)
            gene = unique(unlist(strsplit(na.exclude(unlist(combined[r,grep("Up_regulated",colnames(combined))])),split=", ")))
            if(length(gene)!=0){
              final$Up_regulated[r] = paste(gene,collapse = ", ")
            }
            gene = unique(unlist(strsplit(na.exclude(unlist(combined[r,grep("Down_regulated",colnames(combined))])),split=", ")))
            if(length(gene)!=0){
              final$Down_regulated[r] = paste(gene,collapse = ", ")
            }
          }
          out=which(final$Up_regulated=="" & final$Down_regulated=="")
          if(length(out)>0){
            final = final[-out,]
          }
          if(nrow(final)>0){
            saveRDS(final,file=paste0("findpathR/",cell,"_in_",trait_ord[t],"_KEGG_pathfindR.rds"))
            pdf(paste0(cell,"_in_",trait_ord[t],"_pathways.pdf"),height = 11,width = 12)
            enrichment_chart(final,top_terms = NULL)
            UpSet_plot(final, use_description = TRUE,num_terms = NULL)
            term_gene_graph(final, use_description = TRUE,num_terms = NULL)
            dev.off()
            final$Celltype = cells[i]
            final$trait=trait_ord[t]
            pathf = rbind(pathf,final)
          }
        }
      }
}