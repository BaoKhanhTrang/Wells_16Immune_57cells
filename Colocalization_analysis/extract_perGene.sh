#!/bin/bash

#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=trangk@chop.edu
#SBATCH --open-mode=append

name=$1
db=$2
tail=$3
dirin=$4
dir=/mnt/isilon/sfgi/trangk/analyses/wells/eQTLs

lastfile=$(ls -l -t $dir/databases/$db/input/$name.ENSG*.txt | head -n 1 | awk '{print $NF}' | xargs basename | cut -d'.' -f2)
lin=$(grep -n "$lastfile" $dir/genes_locus/V2G_ensembl_gene_id.txt | cut -d':' -f1)
tail -n +"$lin" $dir/genes_locus/V2G_ensembl_gene_id.txt | while read gene ; do
	#if [ ! -s $dir/databases/$db/input/$name.$gene.txt ] ; then
		if [[ $tail =~ ".gz" ]]; then
			zcat $dirin/${name}${tail} | grep "$gene" > $dir/databases/$db/input/$name.$gene.txt
		else
			grep "$gene" $dirin/${name}${tail} > $dir/databases/$db/input/$name.$gene.txt
		fi
	#fi
	echo $(wc -l $dir/databases/$db/input/$name.$gene.txt)
	#if [ ! -s $dir/databases/$db/input/$name.$gene.txt ] ; then
	#	echo "no " $gene " in " $name
	#fi 
done
