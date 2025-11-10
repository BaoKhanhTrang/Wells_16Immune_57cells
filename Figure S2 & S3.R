trait_ord = readRDS("trait_ord.rds")
allstat=readRDS("allstat.rds")
trait=read.delim("trait_gwas.txt",header=T)
ALL=NULL
for(i in 1:nrow(trait)){
  tmp = read.delim(paste0("~/Documents/analysis/ldsc/immune_panel/narrow/summary.lia_",trait$Name[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(trait$Name[i])
  tmp$cate = "cREs"
  tmp$N = tmp$Prop._SNPs * trait$RemainSNPs[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0("~/Documents/analysis/ldsc/immune_panel/atac/summary.lia_",trait$Name[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(trait$Name[i])
  tmp$cate = "Total OCRs"
  tmp$N = tmp$Prop._SNPs * trait$RemainSNPs[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0("~/Documents/analysis/ldsc/immune_panel/Prom/summary.lia_",trait$Name[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(trait$Name[i])
  tmp$cate = "Promoter OCRs"
  tmp$N = tmp$Prop._SNPs * trait$RemainSNPs[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0("~/Documents/analysis/ldsc/immune_panel/expand/summary.lia_",trait$Name[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(trait$Name[i])
  tmp$cate = "cREs ±500bps"
  tmp$N = tmp$Prop._SNPs * trait$RemainSNPs[i]
  ALL=rbind(ALL,tmp)
  
  tmp = read.delim(paste0("~/Documents/analysis/ldsc/immune_panel/notcRE_notProm/summary.lia_",trait$Name[i],".txt"),header=T,sep="\t",fill = TRUE)
  tmp$trait = paste(trait$Name[i])
  tmp$cate = "not_cREProm_OCRs"
  tmp$N = tmp$Prop._SNPs * trait$RemainSNPs[i]
  ALL=rbind(ALL,tmp)
}
rm(i,tmp)
ALL = merge(ALL,allstat[,c('Celltype','system','Platform','Name')],by='Celltype')
ALL$material="Bulk ATAC-seq"
ALL$material[grep("klaus_pancreas_cells",ALL$Celltype)]="Single-cell ATAC-seq"
ALL$Name = factor(ALL$Name,levels = levels(allstat$Name))
ALL$trait = factor(ALL$trait,levels = trait_ord)
ALL$cate = factor(ALL$cate,levels = c("Total OCRs",
                                      "Promoter OCRs",
                                      "cREs",
                                      "cREs ±500bps",
                                      "not_cREProm_OCRs"))
ALL$sig = round(ALL$Enrichment_p,digits = 2) <= 0.05

sd=read.table("~/Documents/analysis/ldsc/immune_panel/narrow.sd.txt",header=F,sep=":",strip.white = T)
colnames(sd)=c("Celltype","sdc")
sd$cate="cREs"
tmp=read.table("~/Documents/analysis/ldsc/immune_panel/atac.sd.txt",header=F,sep=":",strip.white = T)
colnames(tmp)=c("Celltype","sdc")
tmp$cate="Total OCRs"
sd = rbind(sd,tmp)
tmp=read.table("~/Documents/analysis/ldsc/immune_panel/Prom.sd.txt",header=F,sep=":",strip.white = T)
colnames(tmp)=c("Celltype","sdc")
tmp$cate="Promoter OCRs"
sd = rbind(sd,tmp)
tmp=read.table("~/Documents/analysis/ldsc/immune_panel/expand.sd.txt",header=F,sep=":",strip.white = T)
colnames(tmp)=c("Celltype","sdc")
tmp$cate="cREs ±500bps"
sd = rbind(sd,tmp)
tmp=read.table("~/Documents/analysis/ldsc/immune_panel/notcRE_notProm.sd.txt",header=F,sep=":",strip.white = T)
colnames(tmp)=c("Celltype","sdc")
tmp$cate="not_cREProm_OCRs"
sd = rbind(sd,tmp)
rm(tmp)

ALL=merge(ALL,sd,by=c("Celltype","cate"))
ALL=merge(ALL,trait[,c("Name","RemainSNPs")],by.x="trait",by.y="Name")

h2=read.table("~/Documents/analysis/ldsc/immune_panel/traits.h2",header=F)
colnames(h2)=c("trait",'h2g','h2gSE')
ALL=merge(ALL,h2,by="trait")

ALL$tau = ALL$RemainSNPs*ALL$sdc*ALL$Coefficient/ALL$h2g
ALL$tau_se = ALL$RemainSNPs*ALL$sdc*ALL$Coefficient_std_error/ALL$h2g
ALL$tau_p <- 2 * (1 - pnorm(abs(ALL$Coefficient_z.score)))

ALL$trait = factor(ALL$trait,levels = trait_ord)
all=NULL
for(i in 1:nrow(trait)){
  tmp=ALL[ALL$trait==trait$Name[i],]
  tmp$FDR = p.adjust(tmp$tau_p,"fdr")
  all=rbind(all,tmp)
}
ALL = all

all=NULL
for(i in 1:nrow(allstat)){
  a = ALL[which(ALL$Celltype==allstat$Celltype[i]),]
  
  for(e in 1:length(trait_ord)){
    b = a[which(a$trait==trait_ord[e]),]
    cre=b[which(b$cate=="cREs"),]
    b = b[-which(b$cate=="cREs"),]
    b$pval=NA
    for(x in 1:4){
      b$pval[x] = BSDA::tsum.test(mean.x = cre$Enrichment,s.x = cre$Enrichment_std_error,n.x = round(cre$N),
                            mean.y = b$Enrichment[x],s.y = b$Enrichment_std_error[x],n.y = round(b$N[x]),
                            alternative = "greater")$p.value
      b$pval_tau[x] = BSDA::tsum.test(mean.x = cre$tau,s.x = cre$tau_se,n.x = round(cre$N),
                                  mean.y = b$tau[x],s.y = b$tau_se[x],n.y = round(b$N[x]),
                                  alternative = "greater")$p.value
    }
    cre$pval=0
    cre$pval_tau=0
    all=rbind(all,cre,b)
  }
} 
all$tsig = ifelse(all$pval < 0.05 & all$cate!="cREs","p-value < 0.05",NA)
all$tsig_tau = ifelse(all$pval_tau < 0.05 & all$cate!="cREs","p-value < 0.05",NA)
all$FDRsig = ifelse(all$FDR < 0.05 ,"perTrait FDR < 0.05",NA)
all$sig = ifelse(all$Enrichment_p < 0.05,"p-value < 0.05",NA)
all$sig_tau = ifelse(all$tau_p < 0.05,"p-value < 0.05",NA)
all$fdr = p.adjust(all$tau_p,"fdr")
all$fdrsig = ifelse(all$fdr < 0.05 ,"FDR < 0.05",NA)

# ======= Figure S2 A ======= 
tmp=reshape2::dcast(trait+cate~Name,data=all[which(all$Enrichment_p < 0.05 & all$Enrichment>-1),],value.var = "Enrichment")
tmp$good = apply(tmp[,-c(1,2)],1,function(x) length(na.exclude(x)>1))
tmp$trait = factor(tmp$trait,levels = paste(rev(tmp2$trait)))
pdf("trait_vs_cell.pdf",width = 5,height = 4)
ggplot(tmp)+geom_bar(aes(y=trait,x=good),stat = 'identity')+
  facet_grid(cols=vars(cate),space = "free_y",scales = "free")+
  theme_bw()+
  scale_x_continuous(breaks = pretty_breaks())
dev.off()

# ======= Figure S2 B ======= 
tmp=reshape2::dcast(Name+cate~trait,data=all[which(all$Enrichment_p < 0.05 & all$Enrichment>-1),],value.var = "Enrichment")
tmp$good = apply(tmp[,-c(1,2)],1,function(x) length(na.exclude(x)>1))
tmp = merge(tmp,allstat[,c('Name','system')],by="Name")
tmp$Name = factor(tmp$Name,levels = paste(rev(tmp2$Name)))
pdf("cell_vs_trait.pdf",width = 5,height = 7.5)
ggplot(tmp)+geom_bar(aes(y=Name,x=good),stat = 'identity')+
  facet_grid(cols=vars(cate),rows = vars(system),space = "free_y",scales = "free")+
  theme_bw()+
  scale_x_continuous(breaks = pretty_breaks())
dev.off()

# ======= Figure S3 A ======= 
tmp=reshape2::dcast(trait+cate~Name,data=all[which(all$tau_p < 0.05 ),],value.var = "tau")
tmp$good = apply(tmp[,-c(1,2)],1,function(x) length(na.exclude(x)>1))
tmp2=dcast(trait~cate,data=tmp,value.var = "good")
tmp2=tmp2[order(tmp2$`Total OCRs`,
                tmp2$`Promoter OCRs`,
                tmp2$`cREs`,
                tmp2$`cREs ±500bps`,
                tmp2$`not_cREProm_OCRs`,decreasing = T),]
ord = paste(tmp2$trait)
tmp2=reshape2::dcast(trait+cate~Name,data=all[which(all$FDR < 0.05 ),],value.var = "tau")
tmp2$FDRgood = apply(tmp2[,-c(1,2)],1,function(x) length(na.exclude(x)>1))
tmp=merge(tmp[,c('trait','cate','good')],tmp2[,c('trait','cate','FDRgood')],by=c('trait','cate'),all=TRUE)
tmp$pgood=tmp$good-tmp$FDRgood
tmp$pgood[which(is.na(tmp$FDRgood))] = tmp$good[which(is.na(tmp$FDRgood))]
tmp2=melt(tmp,id.vars = c('trait','cate'),measure.vars = c('FDRgood','pgood'),value.name = "num")
tmp2$variable=factor(tmp2$variable,levels = c("pgood","FDRgood"))
tmp2$trait = factor(tmp2$trait,levels = rev(ord))
pdf("trait_vs_cell_tau.pdf",width = 5,height = 4)
ggplot(tmp2)+geom_bar(aes(y=trait,x=num,fill=variable),stat = 'identity',position = "stack")+
  facet_grid(cols=vars(cate))+
  theme_bw()+
  scale_x_continuous(breaks = pretty_breaks())+theme(legend.position = "none")
dev.off()

# ======= Figure S3 B ======= 
tmp=reshape2::dcast(Name+cate~trait,data=all[which(all$tau_p < 0.05 ),],value.var = "tau")
tmp$good = apply(tmp[,-c(1,2)],1,function(x) length(na.exclude(x)>1))
tmp = merge(tmp,allstat[,c('Name','system')],by="Name")
tmp2=dcast(Name+system~cate,data=tmp,value.var = "good")
tmp2=tmp2[order(tmp2$system,
                -tmp2$`Total OCRs`,
                -tmp2$`Promoter OCRs`,
                -tmp2$`cREs`,
                -tmp2$`cREs ±500bps`,
                -tmp2$`not_cREProm_OCRs`,decreasing = F),]
ord=rev(paste(tmp2$Name))
tmp2=reshape2::dcast(Name+cate~trait,data=all[which(all$FDR < 0.05 ),],value.var = "tau")
tmp2$FDRgood = apply(tmp2[,-c(1,2)],1,function(x) length(na.exclude(x)>1))
tmp=merge(tmp[,c('Name','cate','good')],tmp2[,c('Name','cate','FDRgood')],by=c('Name','cate'),all=TRUE)
tmp$pgood=tmp$good-tmp$FDRgood
tmp$pgood[which(is.na(tmp$FDRgood))] = tmp$good[which(is.na(tmp$FDRgood))]
tmp2=melt(tmp,id.vars = c('Name','cate'),measure.vars = c('FDRgood','pgood'),value.name = "num")
tmp = merge(tmp2,allstat[,c('Name','system')],by="Name")
tmp$variable=factor(tmp$variable,levels = c("pgood","FDRgood"))
tmp$Name = factor(tmp$Name,levels = ord)
pdf("cell_vs_trait_tau.pdf",width = 5,height = 7.5)
ggplot(tmp)+geom_bar(aes(y=Name,x=num,fill=variable),stat = 'identity',position = "stack")+
  facet_grid(cols=vars(cate),rows = vars(system),space = "free_y",scales = "free")+
  theme_bw()+
  scale_x_continuous(breaks = pretty_breaks())+theme(legend.position = "none")
dev.off()

# ======= Figure S3 C ======= 
all$trait = factor(all$trait,levels = c("GD","HT",
                                        "AS",'JIA','MS',"PSO",
                                        'ALG','CEL','ECZ','T1D',"RA","SLE",
                                        'IBD','UC','CRO','VIT'))
cairo_pdf("~/Documents/analyses/ldsc/immune_panel/MAplot_trait_tau.pdf",height = 8,width = 6)
ggplot(all) + geom_point(aes(tau,-log10(tau_p)))+
  geom_hline(yintercept = -log10(0.05),color="red",linetype=2)+
  xlab(paste0("Conditional annotation effect size (", "\U1D749","*)"))+
  ylab("-log10(p-value)")+
  facet_grid(cols=vars(cate),rows = vars(trait),scale = "free")+theme_bw()
dev.off()

# ======= Figure S3 D ======= 
cairo_pdf("~/Documents/analysis/ldsc/immune_panel/MAplot_tau.pdf",height = 4,width = 6)
ggplot(all) + geom_point(aes(tau,-log10(tau_p)))+
  geom_hline(yintercept = -log10(0.05),color="red",linetype=2)+
  xlab(paste0("Conditional annotation effect size (", "\U1D749","*)"))+
  ylab("-log10(p-value)")+
  facet_grid(cols=vars(cate),rows = vars(system),scale = "free")+theme_bw()
dev.off()