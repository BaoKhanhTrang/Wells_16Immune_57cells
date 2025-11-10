#!/bin/bash

#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=trangk@chop.edu

dir=/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs
trait=$1
db=$2

printf "locus\tgene_name\tensembl_gene_id\tCell\tPPID\tPP\n" > $dir/output/$trait/$db.coloc_results_all_summary.txt
cd $dir/output/$trait/$db

ls $dir/position_loci/258loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")
    cat $dir/databases/$db/cells.txt | while read cell ; do
        for colocOut in ./$loc.*.$cell.coloc_results_summary.txt ; do
            colocfilename=`basename $colocOut`
            geneID=$(echo $colocfilename | cut -d'.' -f4)
            genename=$(echo $colocfilename | cut -d'.' -f5)
            sed "s/^/$cell /"  $colocOut | sed "s/^/$geneID /" | sed "s/^/$genename /" | sed "s/^/$loc /" | tr " " "\t"  | tail -n+2 >> $dir/output/$trait/$db.coloc_results_all_summary.txt
        done
    done
done

module load R/4.3.3
cd $dir/output/$trait
Rscript $dir/summarize_qtl_coloc_PP3_PP4_results.R $dir/output/$trait/$db.coloc_results_all_summary.txt $db