dir=eQTLs


cd $dir
for trait in "ALG" "AS" "CEL" "CRO" "ECZ" "GD" "HT" "IBD" "JIA" "MS" "PSO" "RA" "SLE" "T1D" "UC" "VIT" ; do
    mkdir -p traits/$trait
    for db in "GTEx" "OneK1K" "DICE" "eQTL_catalogue" "ImmuNexUT" ; do
        mkdir -p databases/$db
        mkdir -p output/$trait/$db
    done
done

# expand proxies for each locus
cd $dir/position_loci
Rscript createProxies.R 

# extract MAF
for i in {1..22} "X" "Y" ; do
    echo -e "ID\tREF\tALT\tREF_FREQ\tALT_FREQ" > $dir/MAF/chr$i.maf
    sbatch -t 24:00:00 -c 4 --mem 100G -o $dir/log/extractMAF.chr$i -J extractMAF.chr$i $dir/run_R.sh extract_MAF_perChr.R $i
done

# extract eGenes from input ###############################################################################
num="258"

dbname="DICE"
tail=".hg38.gz"
eGene_db="DICE/eGENE"
gene_col=8
pos_col=2
snp_col=3
chrom_col=1
for name in "B_CELL_NAIVE" "CD4_NAIVE" "CD4_STIM" "M2" "MONOCYTES" "NK" "TFH" "TH1" "TH2" "TH17" "TREG_NAIVE" "TREG_MEM" ; do
    sbatch -c 8 -t 12:00:00 -o $dir/log/ext.eGen.$dbname.$name.out -J ext.eGen.$dbname.$name $dir/extract_eGenes.sh $num $eGene_db $dbname $name $tail $gene_col $pos_col $snp_col $chrom_col
done

dbname="OneK1K"
tail="_esnp_table.tsv.gz"
eGene_db="OneK1K/eSNP"
gene_col=6
pos_col=8
snp_col=3
chrom_col=7
for name in "bin" "bmem" "plasma" "cd4nc" "cd4et" "monoc" "mononc" "nk" "dc" ; do
    sbatch -c 8 -t 12:00:00 -o $dir/log/ext.eGen.$dbname.$name.out -J ext.eGen.$dbname.$name $dir/extract_eGenes.sh $num $eGene_db $dbname $name $tail $gene_col $pos_col $snp_col $chrom_col
done

dbname="eQTL_catalogue"
tail=".permuted.tsv.gz"
eGene_db="eQTL_catalogue"
gene_col=1
pos_col=7
snp_col=11
chrom_col=6
for name in "QTD000021" "QTD000031" "QTD000036" "QTD000105" "QTD000409"  "QTD000419" "QTD000424" "QTD000439" "QTD000444" "QTD000449" "QTD000454" "QTD000459" "QTD000464" "QTD000469" "QTD000474" "QTD000479" "QTD000484" "QTD000499" "QTD000504" "QTD000509" ; do
    sbatch -c 8 -t 12:00:00 -o $dir/log/ext.eGen.$dbname.$name.out -J ext.eGen.$dbname.$name $dir/extract_eGenes.sh $num $eGene_db $dbname $name $tail $gene_col $pos_col $snp_col $chrom_col
done

dbname="GTEx"
tail=".v8.egenes.txt.gz"
eGene_db="GTEx/GTEx_Analysis_v8_eQTL/"
gene_col=1
pos_col=15
snp_col=19
chrom_col=14
for name in "Whole_Blood" ; do
    sbatch -c 8 -t 12:00:00 -o $dir/log/ext.eGen.$dbname.$name.out -J ext.eGen.$dbname.$name $dir/extract_eGenes.sh $num $eGene_db $dbname $name $tail $gene_col $pos_col $snp_col $chrom_col
done

# count eGenes
for db in "ImmuNexUT" "DICE" "eQTL_catalogue" "OneK1K" "GTEx" ; do
    sbatch -c 12 --mem 50G -o $dir/log/count.eGenes.$db.out -J count.eGenes.$db $dir/run_R.sh genes_locus/compare_eGene_within_locus.R $db 
done


# extract input from traits per locus
num="258"
for trait in "ALG" "AS" "CEL" "CRO" "ECZ" "GD" "HT" "IBD" "JIA" "MS" "PSO" "RA" "SLE" "T1D" "UC" "VIT" ; do
    sbatch -c 8 -o $dir/log/extract_input.$trait.out -J extract_input.$trait $dir/run_R.sh extract_input.R $num $trait
done

num="258"
ls $dir/position_loci/${num}loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")
    sbatch -c 8 -o $dir/log/extractMAFinput.$loc.out -J extractMAF_input.$loc $dir/run_R.sh extract_MAF_forInput.R $num $loc
done

# check completion
ls $dir/position_loci/${num}loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")
    co=$(wc -l < traits/ALG/$loc.maf)
    if [[ $co -eq 0 ]]; then 
        echo $loc
    fi
done


# extract input from eQTLs per genes
db="DICE"
tail=".hg38.gz"
dirin="DICE/unfiltered"
for name in "B_CELL_NAIVE" "CD4_NAIVE" "CD4_STIM" "M2" "MONOCYTES" "NK" "TFH" "TH1" "TH2" "TH17" "TREG_NAIVE" "TREG_MEM" ; do
    sbatch -c 8 -t 60:00:00 -J extract.$name.$db -o $dir/log/extract.$name.$db.out $dir/genes_locus/extract_perGene.sh $name $db $tail $dirin
done


db="OneK1K"
tail="_eqtl_table.tsv.gz"
dirin="OneK1K/eQTL"
for name in "bin" "bmem" "plasma" "cd4nc" "cd4et" "monoc" "mononc" "nk" "dc" ; do
    sbatch -c 8 -t 60:00:00 -J extract.$name.$db -o $dir/log/extract.$name.$db.out $dir/genes_locus/extract_perGene.sh $name $db $tail $dirin
done

db="eQTL_catalogue"
dirin="eQTL_catalogue"
tail=".all.tsv.gz"
for name in "QTD000021" "QTD000031" "QTD000036" "QTD000105" "QTD000409" "QTD000439" "QTD000444" "QTD000449" "QTD000454" "QTD000459" "QTD000464" "QTD000469" "QTD000474" "QTD000479" "QTD000484" "QTD000499" "QTD000504" "QTD000509" ; do
    sbatch -c 8 -t 60:00:00 -J extract.$name.$db -o $dir/log/extract.$name.$db.out $dir/genes_locus/extract_perGene.sh $name $db $tail $dirin
done
#for name in "B-cell_naive" "CD4_T-cell_naive" "T-cell" "CD4_T-cell_anti-CD3-CD28" "monocyte_naive" "monocyte" "monocyte_CD16_naive" "NK-cell_naive" "Tfh_memory" "Th1_memory" "Th2_memory" "Th17_memory" "Th1-17_memory" "Treg_naive" "Treg_memory"

# run eQTL

db="DICE"
N=91
ls $dir/position_loci/${num}loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")  
    sbatch -c 10 --mem 50G -o $dir/log/runeQTL.$loc.$db.out -J $db.$loc.runeQTL $dir/run_R.sh eQTL_perLocus.R $loc $num $db $N
done

db="OneK1K"
N=982
ls $dir/position_loci/${num}loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")  
    sbatch -c 10 --mem 50G -o $dir/log/runeQTL.$loc.$db.out -J $db.$loc.runeQTL $dir/run_R.sh eQTL_perLocus.R $loc $num $db $N
done

db="eQTL_catalogue"
N="NULL"
ls $dir/position_loci/${num}loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")  
    sbatch -c 10 --mem 50G -o $dir/log/runeQTL.$loc.$db.out -J $db.$loc.runeQTL $dir/run_R.sh eQTL_perLocus.R $loc $num $db $N
done

# for FDFT1 only
for trait in "ALG" "AS" "CEL" "CRO" "ECZ" "GD" "HT" "IBD" "JIA" "MS" "PSO" "RA" "SLE" "T1D" "UC" "VIT" ; do
  sbatch -c 14 --mem 80G -o $dir/log/runeQTL.ENSG00000079459.$db.$trait.out -J $db.$trait.runeQTL.ENSG00000079459 $dir/run_R.sh eQTL_perLocus.oneGene.R "chr8.9921161.11972641" "258" $db "ENSG00000079459" $N $trait
done

# summary
for trait in "ALG" "AS"  "CEL" "CRO" "ECZ" "GD"  "HT"  "IBD" "JIA" "MS"  "PSO" "RA"  "SLE" "T1D" "UC"  "VIT" ; do
    for db in "ImmuNexUT" "DICE" "eQTL_catalogue" "OneK1K" ; do
        sbatch -c 8 -o $dir/log/summarize.$trait.$db -J summarize.$trait.$db $dir/summarize.sh $trait $db
    done
done

# summary for FDFT1 only
printf "ensembl_gene_id\tgene_name\tfilename\tTrait\tsource\tPPID\tPP\n" > $dir/output/ENSG00000079459.coloc_results_all_summary.txt
for trait in "ALG" "AS"  "CEL" "CRO" "ECZ" "GD"  "HT"  "IBD" "JIA" "MS"  "PSO" "RA"  "SLE" "T1D" "UC"  "VIT" ; do
    for db in "ImmuNexUT" "DICE" "eQTL_catalogue" "OneK1K" ; do
        cd $dir/output/$trait/$db
        
        for colocOut in ./ENSG00000079459.FDFT1.*coloc_results_summary.txt ; do
            colocfilename=$(basename $colocOut |  sed "s/\.coloc_results_summary.txt//;s/\./\t/g")
            sed "s/^/$db /" $colocOut | sed "s/^/$trait /" | sed "s/^/$colocfilename /"  | tr " " "\t"  | tail -n+2 >> $dir/output/ENSG00000079459.coloc_results_all_summary.txt
        done
    done 
done
Rscript $dir/summarize_qtl_coloc_PP3_PP4_results.R $dir/output/ENSG00000079459.coloc_results_all_summary.txt

# RUN Susie for fine-mapping
ls $dir/position_loci/258loci/*.input | while read locus ; do
    loc=$(basename $locus ".input")
    sbatch -t 48:00:00 -c 8 --mem 50G -o $dir/log/susie.$loc.out -J susie.$loc $dir/run_R.sh add_Susie.R $loc
done