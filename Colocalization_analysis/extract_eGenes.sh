#!/bin/bash

#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=trangk@chop.edu

dir=/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs
num=$1
eGene_db=$2
dbname=$3
name=$4
tail=$5
gene_col=$6
pos_col=$7
snp_col=$8
chrom_col=$9


ls $dir/position_loci/${num}loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")
    chrom=$(echo $loc | cut -d'.' -f1)
    chr=$(echo $chrom | sed "s/chr//g")
    start=$(echo $loc | cut -d'.' -f2)
    stop=$(echo $loc | cut -d'.' -f3)
    if [[ ! -s $dir/genes_locus/$dbname/$name.$loc.eGenes ]] ; then
        if [[ $tail =~ ".gz" ]]; then
            zcat $eGene_db/${name}${tail} | grep -wf $locus - | cut -f"$gene_col" | sort -u | cut -d'.' -f1 | sort -u > $dir/genes_locus/$dbname/$name.$loc.eGenes
        else
            grep -wf $locus $eGene_db/${name}${tail} | cut -f"$gene_col" | sort -u | cut -d'.' -f1 | sort -u > $dir/genes_locus/$dbname/$name.$loc.eGenes
        fi
    fi
    if [[ ! -s $dir/genes_locus/$dbname/$name.$loc.outsideLD.eGenes ]] ; then
        if [[ $tail =~ ".gz" ]]; then
            zcat $eGene_db/${name}${tail} |  awk -F'\t' -v chrom="$chrom" -v chr="$chr" -v chrom_col="$chrom_col" -v start="$start" -v stop="$stop" -v col="$pos_col" -v gene="$gene_col" -v snp="$snp_col" '($chrom_col == chr || $chrom_col == chrom) && ($col >= start && $col <= stop) { print $gene "\t" $snp}'  > $dir/genes_locus/$dbname/$name.$loc.outsideLD.eGenes
        else
            awk -F'\t' -v chrom="$chrom" -v chr="$chr" -v chrom_col="$chrom_col" -v start="$start" -v stop="$stop" -v col="$pos_col" -v gene="$gene_col" -v snp="$snp_col" '($chrom_col == chr || $chrom_col == chrom) && ($col >= start && $col <= stop) { print $gene "\t" $snp}' $eGene_db/${name}${tail} > $dir/genes_locus/$dbname/$name.$loc.outsideLD.eGenes
        fi
    fi
done
