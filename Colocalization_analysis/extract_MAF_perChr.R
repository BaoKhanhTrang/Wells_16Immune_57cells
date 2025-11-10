args = commandArgs(trailingOnly=TRUE)
library(vroom)
library(foreach)
library(doParallel)
cores <- detectCores()
registerDoParallel(cores)
options(timeout=1000000000)

chr = args[1]
splits=list.files(paste0("/scr1/users/trangk/chr",chr,"/"),pattern = "split_")


res <- foreach(s=1:length(splits)) %dopar% {
  filename = paste0("/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs/MAF/chr",chr,".",splits[s])
  vcf=vroom(paste0("/scr1/users/trangk/chr",chr,"/",splits[s]),delim="\t" ,col_select = c(3,4,5,8),col_names=F,comment="#")
  colnames(vcf)=c("ID","REF","ALT","INFO")

  if(file.exists(filename)){
    done = system(paste("gawk 'END {print}' ",filename), intern = TRUE)
    done = unlist(strsplit(done,split="\t"))[4]
    start = which(vcf$ID == paste(done)) + 1
  }else{
    start = 1
  }
  

  for(i in start:nrow(vcf)){
    info_split = unlist(strsplit(vcf$INFO[i], split = ";"))
    freq_data = info_split[grep("FREQ", info_split)]

    # Initialize REF and ALT frequencies as NA
    num_ref = unlist(strsplit(vcf$REF[i], split = ","))
    num_alt = unlist(strsplit(vcf$ALT[i], split = ","))
    
    sumAllel = c(num_ref,num_alt)

    # Process frequency data if available
    if (length(freq_data) > 0) {
      freq_values = gsub("FREQ=", "", freq_data)
      freq_split = unlist(strsplit(freq_values, split = "\\|"))

      # Extract REF and ALT frequencies
      if (length(freq_split) > 0) {
        allele_split = unlist(strsplit(freq_split, split = ":"))
        freq_split_values = unlist(strsplit(allele_split[seq(2, by = 2, length(allele_split))], split = ","))
        sourc = allele_split[seq(1, by = 2, length(allele_split))]
        
        # Assign values if they exist
        if (length(freq_split_values) >= 2) {
          tab=data.frame(matrix(nrow=length(sourc),ncol=length(sumAllel)+1))
          colnames(tab)=c("source",sumAllel)
          tab$source = sourc
          for(r in 1:length(sumAllel)){
            tab[,r+1]=as.numeric(freq_split_values[seq(r,by=length(sumAllel),length(freq_split_values))])
          }
          a = reshape2::melt(tab, id.vars="source",variable.name="Allele",value.name ="freq")
          a$ID=vcf$ID[i]
          write.table(a,file=filename,col.names=F,row.names=F,quote=F,sep="\t",append=T)
        }
      }
    }
    # Return a data frame with the necessary columns
    
  }
  return(NULL)
}
registerDoSEQ()
doParallel::stopImplicitCluster()
