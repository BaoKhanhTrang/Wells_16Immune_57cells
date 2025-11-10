#!/bin/bash

scriptdir=~/tools
cdir=ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC

# Run partitioned LD score regression #########################
gwas_sumstats=($(echo $cdir/snps/non_merged/*.sumstats.gz))

for gwas_sumstat in gwas_sumstats[@]}; do
    gwas_base=($(basename $gwas_sumstat | cut -d . -f 1))
    sample_prev=$(grep "$gwas_base" traits_specs.txt | cut -f2 )
    pop_prev=$(grep "$gwas_base" traits_specs.txt | cut -f3 )
    cat cells.txt | while read base ; do
        dir="$base""_ldscore"
	    basesrc="pldsr_""$base"
        python ~/tools/ldsc/ldsc.py --h2 $gwas_sumstat --print-coefficients --ref-ld-chr $cdir/annotationfile_narrow/$dir/$base.,$baseline. --samp-prev $sample_prev --pop-prev $pop_prev --w-ld-chr $weight.  --overlap-annot  --frqfile-chr $frq.  --out $cdir/out_narrow/$base_$gwas_base
        python ~/tools/ldsc/ldsc.py --h2 $gwas_sumstat --print-coefficients --ref-ld-chr $cdir/annotationfile_atac/$dir/$base.,$baseline. --samp-prev $sample_prev  --pop-prev $pop_prev  --w-ld-chr $weight.  --overlap-annot  --frqfile-chr $frq.  --out $cdir/out_atac/$base_$gwas_base
        python ~/tools/ldsc/ldsc.py --h2 $gwas_sumstat --print-coefficients --ref-ld-chr $cdir/annotationfile_Prom/$dir/$base.,$baseline. --samp-prev $sample_prev  --pop-prev $pop_prev  --w-ld-chr $weight.  --overlap-annot  --frqfile-chr $frq.  --out $cdir/out_Prom/$base_$gwas_base
        python ~/tools/ldsc/ldsc.py --h2 $gwas_sumstat --print-coefficients --ref-ld-chr $cdir/annotationfile_expand/$dir/$base.,$baseline. --samp-prev $sample_prev  --pop-prev $pop_prev  --w-ld-chr $weight.  --overlap-annot  --frqfile-chr $frq.  --out $cdir/out_expand/$base_$gwas_base
        python ~/tools/ldsc/ldsc.py --h2 $gwas_sumstat --print-coefficients --ref-ld-chr $cdir/annotationfile_notcRE_notProm/$dir/$base.,$baseline. --samp-prev $sample_prev  --pop-prev $pop_prev  --w-ld-chr $weight.  --overlap-annot  --frqfile-chr $frq.  --out $cdir/out_notcRE_notProm/$base_$gwas_base

    done
done

# Summarize results
for gwas_sumstat in gwas_sumstats[@]}; do
    gwas_base=($(basename $gwas_sumstat | cut -d . -f 1))
    for cate in "narrow" "atac" "Prom" "expand" "notcRE" "notcRE_notProm" ; do
        outdir=$cdir/out_$cate
       	echo "Celltype        Prop._SNPs      Prop._h2        Prop._h2_std_error      Enrichment      Enrichment_std_error    Enrichment_p    Coefficient     Coefficient_std_error   Coefficient_z-score" > $cdir/$outdir/summary.lia_$gwas_base.txt
		cat cells.txt | while read base ; do
        	dir="$base""_ldscore"
        	basesrc="pldsr_""$base"
		    res=$(grep "L2_0" $cdir/$outdir/${base}_${gwas_base}.results | cut -f2- )
		    echo $base $res | sed "s/ /\t/g" >> $cdir/$outdir/summary.lia_$gwas_base.txt
	    done
    done
done

# extract h2
for gwas_sumstat in gwas_sumstats[@]}; do
	gwas_base=($(basename $gwas_sumstat | cut -d . -f 1))
    for cate in "narrow" "atac" "Prom" "expand" "notcRE" "notcRE_notProm" ; do
        outdir=$cdir/out_$cate
	    cat cells.txt | while read base ; do
		    echo $base ":" $(grep "Total Liability scale h2" $cdir/$outdir/base}_gwas_base}_lia.log | cut -d' ' -f 5) >> $cdir/$outdir/h2.lia_$gwas_base.txt
	done
done
