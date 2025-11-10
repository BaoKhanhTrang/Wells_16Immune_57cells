#!/bin/bash

cdir=ldsc
cd $cdir/snps

# Download sumstat and transform

zcat 34594039-GCST90018855-EFO_0003779.h.tsv.gz | cut -f 2,5,6,7,11,20 | sed 's/hm_//g' | gzip - > Hashimoto.sumstat.gz
zcat 34594039-GCST90018847-EFO_0004237.h.tsv.gz | cut -f 2,5,6,7,11,20 | sed 's/hm_//g' | gzip - > Grave.sumstat.gzz
zcat 29083406-GCST005038-EFO_0003785.h.tsv.gz | cut -f 2,5,6,7,11,19 | sed 's/hm_//g' | gzip - > ALG.sumstat.gz
zcat 23749187-GCST005529-EFO_0003898.h.tsv.gz | cut -f 2,5,6,7,11,18 | sed 's/hm_//g' | gzip - > AS.sumstat.gz
zcat 20190752-GCST000612-EFO_0001060.h.tsv.gz | cut -f 2,5,6,7,11,18 | sed 's/hm_//g' | gzip - > CEL.sumstat.gz
zcat 28067908-GCST004131-EFO_0003767.h.tsv.gz | cut -f 2,5,6,7,11,17 | sed 's/hm_//g' | gzip - > IBD.sumstat.gz
zcat 28067908-GCST004132-EFO_0000384.h.tsv.gz | cut -f 2,5,6,7,11,17 | sed 's/hm_//g' | gzip - > CRO.sumstat.gz
zcat 28067908-GCST004133-EFO_0000729.h.tsv.gz | cut -f 2,5,6,7,11,17 | sed 's/hm_//g' | gzip - > UC.sumstat.gz
cut -f1,4,5,6,8,10 EAGLE_AD_no23andme_results_29072015.txt | sed 's/eaf/freq/g' | gzip - > ECZ.sumstat.gz
zcat 23603761-GCST005528-EFO_1001999.h.tsv.gz | cut -f 2,5,6,7,11,18 | sed 's/hm_//g' | gzip - > JIA.sumstat.gz
zcat 24076602-GCST005531-EFO_0003885.h.tsv.gz | cut -f 2,5,6,7,11,18 | sed 's/hm_//g' | gzip - > MS.sumstat.gz
zcat 23143594-GCST005527-EFO_0000676.h.tsv.gz | cut -f 2,5,6,7,11,18 | sed 's/hm_//g' | gzip - > PSO.sumstat.gz
zcat 24390342-GCST002318-EFO_0000685.h.tsv.gz | cut -f 2,16,17,18,11,19 | awk '{if($1!="NA")   print $0}' | sed 's/hm_//g' | gzip - > RA.sumstat.gz
zcat 34012112-GCST90014023-EFO_0001359.h.tsv.gz | cut -f 2,5,6,7,11,14 | sed 's/hm_//g' | awk 'NR==1 || $6< 1e300 && $1!="NA" ' | gzip - > T1D.sumstat.gz # not completely remove small p, have to go in with R and cap p by min of the res
sed 1d /mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome/comparison/immunePanel/snp/sumStat/downloads/SLE/SLE_finalSumStat.tsv | sed "1iSNP\talt_id\tchromosome\tposition\tA2\tA1\tBETA\tSE\tP" > SLE.sumstat
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/disease/sumstat/downloads/VIT/VIT.oriSumstat.txt

# munge the sumstat 
conda activate ldsc
python ~/ldsc/munge_sumstats.py --sumstats ALG.sumstat.gz --N-cas 180129 --N-con 180709 --maf-min 0.01 --out non_merged/ALG > ALG_munged.log
python ~/ldsc/munge_sumstats.py --sumstats AS.sumstat.gz --N-cas 9069 --N-con 13578 --out non_merged/AS > AS_munged.log
python ~/ldsc/munge_sumstats.py --sumstats CEL.sumstat.gz --N-cas 4533 --N-con 10750 --out non_merged/CEL > CEL_munged.log
python ~/ldsc/munge_sumstats.py --sumstats IBD.sumstat.gz --N-cas 25042 --N-con 34915 --out non_merged/IBD > IBD_munged.log
python ~/ldsc/munge_sumstats.py --sumstats CRO.sumstat.gz --N-cas 12194 --N-con 28072 --out non_merged/CRO > CRO_munged.log
python ~/ldsc/munge_sumstats.py --sumstats UC.sumstat.gz --N-cas 12366 --N-con 33609 --out non_merged/UC > UC_munged.log
python ~/ldsc/munge_sumstats.py --sumstats ECZ.sumstat.gz --N-cas 18900 --N-con 84166 --maf-min 0.0003 --out non_merged/ECZ > ECZ_munged.log
python ~/ldsc/munge_sumstats.py --sumstats JIA.sumstat.gz --N-cas 2816 --N-con 13056 --out non_merged/JIA > JIA_munged.log
python ~/ldsc/munge_sumstats.py --sumstats MS.sumstat.gz --N-cas 14498 --N-con 24091 --out non_merged/MS > MS_munged.log
python ~/ldsc/munge_sumstats.py --sumstats PSO.sumstat.gz --N-cas 10588 --N-con 22806 --out non_merged/PSO > PSO_munged.log
python ~/ldsc/munge_sumstats.py --sumstats RA.sumstat.gz --N-cas 14361 --N-con 42923 --out non_merged/RA > RA_munged.log
python ~/ldsc/munge_sumstats.py --sumstats T1D.sumstat.gz --N-cas 18942 --N-con 501638 --maf-min 0.0001 --out non_merged/T1D > T1D_munged.log
python ~/ldsc/munge_sumstats.py --sumstats Grave.sumstat.gz --N-cas 1678 --N-con 456942 --maf-min 0.0001 --out non_merged/GD > GD_munged.log
python ~/ldsc/munge_sumstats.py --sumstats Hashimoto.sumstat.gz --N-cas 15654 --N-con 379986 --maf-min 0.0001 --out non_merged/HT > HT_munged.log
python ~/ldsc/munge_sumstats.py --sumstats VIT.oriSumstat.txt --N-cas 2853 --N-con 37405 --maf-min 0.0001 --signed-sumstats ORX,1 --out non_merged/VIT > VIT_munged.log
python ~/ldsc/munge_sumstats.py --sumstats SLE.sumstat.gz --N-cas 6748 --N-con 11516 --out non_merged/SLE > SLE_munged.log

# merge with w_hm3.snplist
python ~/ldsc/munge_sumstats.py --sumstats ALG.sumstat --N-cas 180129 --N-con 180709 --merge-alleles $cdir/w_hm3.snplist --out merged/ALG_merged > ALG_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats AS.sumstat --N-cas 9069 --N-con 13578 --merge-alleles $cdir/w_hm3.snplist --out merged/AS_merged > AS_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats CEL.sumstat --N-cas 4533 --N-con 10750 --merge-alleles $cdir/w_hm3.snplist --out merged/CEL_merged > CEL_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats IBD.sumstat --N-cas 25042 --N-con 34915 --merge-alleles $cdir/w_hm3.snplist --out merged/IBD_merged > IBD_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats CRO.sumstat --N-cas 12194 --N-con 28072 --merge-alleles $cdir/w_hm3.snplist --out merged/CRO_merged > CRO_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats UC.sumstat --N-cas 12366 --N-con 33609 --merge-alleles $cdir/w_hm3.snplist  --out merged/UC_merged > UC_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats ECZ.sumstat --N-cas 18900 --N-con 84166 --merge-alleles $cdir/w_hm3.snplist --maf-min 0.0001 --out merged/ECZ_merged > ECZ_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats JIA.sumstat --N-cas 2816 --N-con 13056 --merge-alleles $cdir/w_hm3.snplist  --out merged/JIA_merged > JIA_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats MS.sumstat --N-cas 14498 --N-con 24091 --merge-alleles $cdir/w_hm3.snplist  --out merged/MS_merged > MS_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats PSO.sumstat --N-cas 10588 --N-con 22806 --merge-alleles $cdir/w_hm3.snplist  --out merged/PSO_merged > PSO_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats RA.sumstat --N-cas 14361 --N-con 42923  --merge-alleles $cdir/w_hm3.snplist  --out merged/RA_merged > RA_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats T1D.sumstat --N-cas 18942 --N-con 501638 --merge-alleles $cdir/w_hm3.snplist --maf-min 0.0001 --out merged/T1D_merged > T1D_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats Grave.sumstat --N-cas 1678 --N-con 456942 --merge-alleles $cdir/w_hm3.snplist --maf-min 0.0001 --out merged/GD_merged > GD_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats Hashimoto.sumstat --N-cas 15654 --N-con 379986 --merge-alleles $cdir/w_hm3.snplist --maf-min 0.0001 --out merged/HT_merged > HT_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats /mnt/isilon/sfgi/suc1/analyses/grant/disease/sumstat/downloads/VIT/VIT.oriSumstat.txt --N-cas 2853 --N-con 37405  --merge-alleles $cdir/w_hm3.snplist --maf-min 0.0001 --signed-sumstats ORX,1 --out merged/VIT_merged > VIT_merged_munged.log
python ~/ldsc/munge_sumstats.py --sumstats SLE.sumstat --N-cas 6748 --N-con 11516 --merge-alleles $cdir/w_hm3.snplist --out merged/SLE_merged > SLE_merged_munged.log

