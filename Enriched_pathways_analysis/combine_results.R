metaData=readRDS("metaData.rds")
allstat=readRDS("allstat.rds")

cells=unique(metaData$Celltype)
pathf=NULL

for(i in 1:length(cells)){
  for(t in 1:length(trait_ord)){
    if(file.exists(paste0("findpathR/",cells[i],"_in_",trait_ord[t],"_KEGG_pathfindR.rds"))){
      tmp=readRDS(paste0("findpathR/",cells[i],"_in_",trait_ord[t],"_KEGG_pathfindR.rds"))
      if(nrow(tmp)>0){
        tmp$Celltype = cells[i]
        tmp$trait=trait_ord[t]
        pathf = rbind(pathf,tmp)
      }
    }
  }
}
pathf=merge(pathf,allstat[,c('Celltype','Name','system')],by="Celltype")
pathf$numGene= apply(pathf,1,function(x) ifelse(grepl(",",x[['Up_regulated']]),
                                          stringr::str_count(x[['Up_regulated']],",") + 1,ifelse(x[['Up_regulated']]=="",0,1))+
                    ifelse(grepl(",",x[['Down_regulated']]),stringr::str_count(x[['Down_regulated']],",")+1,ifelse(x[['Down_regulated']]=="",0,1)))
saveRDS(pathf,file="pathfindRresult.rds")
