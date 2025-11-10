dir="pathway"

# run DE analysis consider each cell type versus all other
for cell in $(cat cells.txt ); do
    Rscript DE_analysis.R $cell
done

# Use DE results as input for pathfindR
# run for each cell type
for cell in $(cat cells.txt ); do
    for trait in "ECZ" "ALG" "UC"  "IBD" "CRO" "AS"  "HT"  "MS"  "T1D" "CEL" "JIA" "PSO" "RA"  "SLE" "GD"  "VIT"; do
        Rscript pathfindR_analysis.R $cell $trait
    done
done

#combine the results
Rscrript combine_results.R

# create WGCNA modules
Rscript WGCNA_analysis.R
