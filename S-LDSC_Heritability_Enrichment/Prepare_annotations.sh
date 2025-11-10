#!/bin/bash

scriptdir=~/tools
cdir=ldsc
DIRIN=v2g/celltypes

# create cREs set from ATAC-seq peaks and hi-C/Capture-C chromatin loops
# ATAC-seq peaks in BED format atac.bed
# hi-C/Capture-C chromatin loops in BEDPE format ${cell}_cis.bedpe (only consider cis-contact on same chromosome)
# gene promoters annotation from Gencode V40
for cell in $(cat cells.txt ); do
	Rscript scripts/annotate_bedpe2geneOCR_hg38.R $DIRIN/$cell/${cell}_cis.bedpe $DIRIN/$cell/atac.bed $DIRIN/$cell/$cell.cis.anno
done

# Prepare BED files #########################
DIROUT=$cdir/BEDfiles

for cell in $(cat cells.txt ); do
	# all OCRs from atac.bed
        ln -s $DIRIN/$cell/atac.bed $DIROUT/atac/$cell.bed 
	# promoter OCRs
        intersectBed -wa -u -a $DIROUT/atac/$cell.bed -b gencode.v40.promoter.bed > $DIROUT/Prom/$cell.bed
	# cREs
		grep -v "close" $DIRIN/$cell/$cell.cis.anno.gene2OCR.txt | cut -f1 | sort -u | sed "s/:/\t/g" > $DIROUT/narrow/$cell.bed
	# extend cREs ± 500 bp
		grep -v "close" $DIRIN/$cell/$cell.cis.anno.gene2OCR.txt | cut -f1 | sort -u | sed "s/:/\t/g"  | awk '{OFS = "\t" ;$2=$2-500; $3=$3+500; print}' > $DIROUT/expand/$cell.bed
done


# Make annotation #########################
conda activate ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC

for base in $(cat cells.txt ); do
	dir="$base""_ldscore"
	for cate in "narrow" "atac" "Prom" "expand" ; do
        annodir=$cdir/annotationfile_$cate
        outdir=$cdir/out_$cate
        file=$cdir/BEDfiles/$cate/$base.bed
		
		mkdir -p $cdir/$annodir/$dir
		# For each chromosome
		for j in {1..22}; do
			## Generate annotation files for each BED files
   			python $scriptdir/ldsc/make_annot.py  \
				--bed-file $file \
				--bimfile $plink.$j.bim  \
				--annot-file $annodir/$dir/$base.$j.annot.gz 

			## Then Calculate ld for each annotation
			python $scriptdir/ldsc/ldsc.py \
				--l2 --bfile $plink.$j \
				--ld-wind-cm 1 \
				--annot $annodir/$dir/$base.$j.annot.gz \
				--thin-annot  \
				--print-snps $cdir/list.txt \
				--out $annodir/$dir/$base.$j
		done
	done
done

# calculate annotation SD for tau-star
for base in $(cat cells.txt ); do
    	dir="$base""_ldscore"
        echo $base ":" $(zcat $cdir/annotationfile_narrow/$dir/$base.*.l2.ldscore.gz | cut -f4 | grep -v "L2" | Rscript $cdir/sd_calc.R) >> $cdir/annotationfile_narrow/sd.txt
        echo $base ":" $(zcat $cdir/annotationfile_atac/$dir/$base.*.l2.ldscore.gz | cut -f4 | grep -v "L2" | Rscript $cdir/sd_calc.R) >> $cdir/annotationfile_atac/sd.txt
        echo $base ":" $(zcat $cdir/annotationfile_Prom/$dir/$base.*.l2.ldscore.gz | cut -f4 | grep -v "L2" | Rscript $cdir/sd_calc.R) >> $cdir/annotationfile_Prom/sd.txt
        echo $base ":" $(zcat $cdir/annotationfile_expand/$dir/$base.*.l2.ldscore.gz | cut -f4 | grep -v "L2" | Rscript $cdir/sd_calc.R) >> $cdir/annotationfile_expand/sd.txt
        echo $base ":" $(zcat $cdir/annotationfile_notcRE_notProm/$dir/$base.*.l2.ldscore.gz | cut -f4 | grep -v "L2" | Rscript $cdir/sd_calc.R) >> $cdir/annotationfile_notcRE_notProm/sd.txt
done

