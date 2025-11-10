DIR=v2g

module load R/4.2.3
cat $DIR/traits.list | while read trait ; do
	proxy=$DIR/SNPset/${trait}_proxies.bed
	prefix=$trait".topLD"
	mkdir -p $DIR/${prefix}_v2g
    for cell in $(cat $DIR/celltypes.txt) ; do
        cd $DIR/celltypes/$cell
        ANNO_DIR=$DIR/celltypes/$cell
        OUT_PREFIX=${prefix}_$cell.cis
        OUT=$DIR/$trait.topLD_v2g
        Rscript v2g_proxybed_cisanno.R $PROXIES $ANNO_DIR OUT_PREFIX $OUT 
        done
done

# external data
module load BEDTools
mkdir $DIR/external_data
for cell in "GSM4118993_T0" "GSM4118994_T20" "GSM4118995_T1H" "GSM4118996_T2H" "GSM4118997_T4H" "GSM4118998_T24H" ; do
	time=$(echo $cell | cut -d'_' -f2)
	zcat $DIR/external_data/Yang2020_promoterCaptureC/${cell}_ATAC.narrowPeak.gz | cut -f 1-3 | sortBed -i - | uniq | awk -v c="$cell" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > Yang2020_$time.atac.bed
done

for cell in "T0" "T20" "T1H" "T4H" "T24H" ; do
	zcat $DIR/external_data/Yang2020_promoterCaptureC/*_${cell}_CHiC_pool*_interactions.txt.gz | awk '{if($1==$4 && $7>5) print}' | sort -k1,1n -k2,2 | uniq > Yang2020_$cell.cis.bedpe
done

for cell in "T0" "T20" "T1H" "T4H" "T24H" ; do
	Rscript annotate_bedpe2geneOCR_hg38.R Yang2020_$cell.cis.bedpe Yang2020_$cell.atac.bed Yang2020_$cell.anno
done
for cell in "T0" "T20" "T1H" "T4H" "T24H"; do
        sbatch -c 8 --mem-per-cpu 16G -o log/$cell.bedpe2geneOCR.out -J $cell.bedpe2geneOCR.out $DIR/run_R.sh  annotate_bedpe2geneOCR_hg38.R $DIR/external_data/Yang2020_$cell.cis.bedpe $DIR/external_data/Yang2020_$cell.atac.bed $DIR/external_data/Yang2020_$cell.anno
done

sed '1d' Javierre_pchic.txt | awk '{if($1==$6) print "chr" $1 "\t" $2 "\t" $3 "\t" "chr" $6 "\t" $7 "\t" $8}' > jav.cis.hg19.bedpe
python liftOverBedpe.py --lift liftOver --chain hg19ToHg38.over.chain.gz --i jav.cis.hg19.bedpe --o jav.cis.bedpe --h F
sed '1d' Jarrier.ActivePromoterEnhancerLinks.tsv | cut -f1-3 > tmp
sed '1d' Jarrier.ActivePromoterEnhancerLinks.tsv | cut -f5-7 >> tmp
liftOver <(sortBed -i tmp | uniq ) hg19ToHg38.over.chain.gz jav.atac.bed fail 
mv jav.atac.bed tmp
awk -v c="JAV" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' tmp > jav.atac.bed
Rscript annotate_bedpe2geneOCR_hg38.R jav.cis.hg19.bedpe jav.atac.bed jav.anno

hc=read.table("external_data/Burren_2017_GB_PCHiC_ActivatedTCells.txt")
hc$baitEnd = hc$baitStart+ hc$baitLength
hc$oeEnd = hc$oeStart + hc$oeLength
write.table(hc[which(hc$baitChr==hc$oeChr & hc$b2b==FALSE),c('baitChr' ,'baitStart','baitEnd','oeChr','oeStart','oeEnd')],file="Burren_phic.hg19.bedpe",col.names=F,row.names=F,sep = "\t",quote=F)
atac=hc[,4:6]
colnames(atac)=colnames(hc)[1:3]
atac=unique(rbind(atac,hc[,1:3]))
atac=atac[order(atac$baitChr,atac$baitStart),]
atac$name = paste0("Burren.atac_",seq(1,nrow(atac)))
write.table(atac,file="Burren.atac.hg19.bed",col.names=F,row.names=F,sep = "\t",quote=F)
python liftOverBedpe.py --lift liftOver --chain hg19ToHg38.over.chain.gz --i Burren_phic.hg19.bedpe --o Burren_phic.cis.bedpe --h F
sort -u Burren_phic.cis.bedpe > tmp
sort -k1,1n -k2,2n tmp | awk '{if($1==$4) print}' > Burren_phic.cis.bedpe
liftOver Burren.atac.hg19.bed hg19ToHg38.over.chain.gz Burren.atac.bed fail


cell="Burren"
sbatch -c 8 --mem 264G -o log/$cell.bedpe2geneOCR.out -J $cell.bedpe2geneOCR.out $DIR/run_R.sh annotate_bedpe2geneOCR_hg38.R $DIR/external_data/Burren_phic.cis.bedpe $DIR/external_data/Burren.atac.bed $DIR/external_data/Burren.anno
Rscript annotate_bedpe2geneOCR_hg38.R Burren_phic.cis.bedpe Burren.atac.bed Burren.anno

sed '1d' external_data/Gate_2018_NatGen/HiC_Arrowhead.txt | awk '{if($1==$4) print "chr" $1 "\t" $2 "\t" $3 "\t" "chr" $4 "\t" $5 "\t" $6}' | sort -k1,1 -k2,2n -u   > Gate.cis.hg19.bedpe
python liftOverBedpe.py --lift liftOver --chain hg19ToHg38.over.chain.gz --i Gate.cis.hg19.bedpe --o Gate.cis.bedpe --h F
mv Gate.cis.bedpe tmp
uniq tmp > Gate.cis.bedpe
sed '1d' external_data/Gate_2018_NatGen/ATAC_coaccessiblity.txt | cut -f1 | sed "s/-/\t/g;s/:/\t/g" > atac.tmp
sed '1d' external_data/Gate_2018_NatGen/ATAC_coaccessiblity.txt | cut -f2 | sed "s/-/\t/g;s/:/\t/g" >> atac.tmp
liftOver <(sortBed -i atac.tmp | uniq) hg19ToHg38.over.chain.gz Gate.atac.bed fail
mv Gate.atac.bed atac.tmp
awk -v c="Gate" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' atac.tmp > Gate.atac.bed
Rscript annotate_bedpe2geneOCR_hg38.R Gate.cis.bedpe Gate.atac.bed Gate.anno


proxy=$DIR/SNPset/all_proxies.bed
ls *.gene2OCR.txt | while read ANNO_FILE ; do
	cell=$(basename $ANNO_FILE ".anno.gene2OCR.txt")
	Rscript v2g_proxybed_cisanno.R $proxy $ANNO_FILE $cell.v2g $DIR/external_data
done