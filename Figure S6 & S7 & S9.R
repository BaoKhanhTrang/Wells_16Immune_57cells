library(ggplot2)
library(ggrepel)
library(ggiraph)


pathf=readRDS("pathfindRresult.rds")

# ======= Figure S6 A =======
unip=pathf[which(pathf$Term_Description == "Cytokine-cytokine receptor interaction"),] # one pathway at a time
unip$Down_regulated[which(unip$Down_regulated=="")]=NA
unip$Up_regulated[which(unip$Up_regulated=="")]=NA
genes=unip[,c('Term_Description','Name','Up_regulated','system','trait')]
colnames(genes)[3]="label"
genes$direction="Higher_expressed"
tmp=unip[,c('Term_Description','Name','Down_regulated','system','trait')]
colnames(tmp)[3]="label"
tmp$direction="Lower_expressed"
genes=rbind(genes,tmp)
genes=genes[-which(is.na(genes$label)),]
g=ggplot(unip)+
  geom_point(aes(Name,trait,size=factor(numGene)),shape=1,alpha=.6)+
  theme_bw()+
  geom_label_repel_interactive(data=genes,aes(Name,trait,label=stringr::str_wrap(label,7),color=direction),size=2,label.padding=unit(0.1, "lines"),force = 4)+
  scale_size_discrete(name="#genes")+
  facet_grid(cols = vars(system),drop=T,space = "free",scale="free")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

pdf("Cytokine-cytokine receptor interaction.pdf",height = 9,width = 14)
print(g)
dev.off()

# ======= Figure S6 B =======
up=unique(unlist(strsplit(unip$Up_regulated,split=", ")))
down=unique(unlist(strsplit(unip$Down_regulated,split=", ")))
intersect(up,down)
change_vec = structure(c(rep(1,length(up)),rep(-1,length(down))),names=c(up,down))
pathview(gene.data = change_vec,pathway.id = unique(unip$ID),out.suffix = "",same.layer=F,gene.idtype = "SYMBOL",
         limit = list(gene=c(min(change_vec),max(change_vec))),high = "red",low = "#8FD6D9",mid = "white",
         kegg.dir = "Kegg/")


# ======= Figure S7 A =======
unip=pathf[which(pathf$Term_Description == "Salmonella infection"),] # one pathway at a time
unip$Down_regulated[which(unip$Down_regulated=="")]=NA
unip$Up_regulated[which(unip$Up_regulated=="")]=NA
genes=unip[,c('Term_Description','Name','Up_regulated','system','trait')]
colnames(genes)[3]="label"
genes$direction="Higher_expressed"
tmp=unip[,c('Term_Description','Name','Down_regulated','system','trait')]
colnames(tmp)[3]="label"
tmp$direction="Lower_expressed"
genes=rbind(genes,tmp)
genes=genes[-which(is.na(genes$label)),]
g=ggplot(unip)+
  geom_point(aes(Name,trait,size=factor(numGene)),shape=1,alpha=.6)+
  theme_bw()+
  geom_label_repel_interactive(data=genes,aes(Name,trait,label=stringr::str_wrap(label,7),color=direction),size=2,label.padding=unit(0.1, "lines"),force = 4)+
  scale_size_discrete(name="#genes")+
  facet_grid(cols = vars(system),drop=T,space = "free",scale="free")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

pdf("Salmonella infection.pdf",height = 9,width = 14)
print(g)
dev.off()

# ======= Figure S7 B =======
up=unique(unlist(strsplit(unip$Up_regulated,split=", ")))
down=unique(unlist(strsplit(unip$Down_regulated,split=", ")))
intersect(up,down)
change_vec = structure(c(rep(1,length(up)),rep(-1,length(down))),names=c(up,down))
pathview(gene.data = change_vec,pathway.id = unique(unip$ID),out.suffix = "",same.layer=F,gene.idtype = "SYMBOL",
         limit = list(gene=c(min(change_vec),max(change_vec))),high = "red",low = "#8FD6D9",mid = "white",
         kegg.dir = "Kegg/")


# ======= Figure S9 =======
path = reshape2::dcast(Term_Description~trait,data=unique(pathf[,c('trait','Term_Description')]),fun.aggregate = length)
path$numtrait = apply(path[,-1],1,function(x) length(which(x!=0)))
path$Term_Description[which(path$numtrait==1)]

unip=pathf[which(pathf$Term_Description %in% path$Term_Description[which(path$numtrait==1)]),] # one pathway at a time
unip$Down_regulated[which(unip$Down_regulated=="")]=NA
unip$Up_regulated[which(unip$Up_regulated=="")]=NA
genes=unip[,c('Term_Description','Name','Up_regulated','system','trait')]
colnames(genes)[3]="label"
genes$direction="Higher_expressed"
tmp=unip[,c('Term_Description','Name','Down_regulated','system','trait')]
colnames(tmp)[3]="label"
tmp$direction="Lower_expressed"
genes=rbind(genes,tmp)
genes=genes[-which(is.na(genes$label)),]
g=ggplot(unip)+
  geom_point(aes(Term_Description,Name,size=factor(numGene)),shape=1,alpha=.6)+
  theme_bw()+
  geom_label_repel_interactive(data=genes,aes(Term_Description,Name,label=stringr::str_wrap(label,7),color=direction),size=2,label.padding=unit(0.1, "lines"),force = 4)+
  scale_size_discrete(name="#genes")+
  facet_grid(rows = vars(system),cols=vars(trait),drop=T,space = "free",scale="free")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))


pdf("31 trait-unique KEGG pathsways.pdf",height = 9,width = 14)
print(g)
dev.off()