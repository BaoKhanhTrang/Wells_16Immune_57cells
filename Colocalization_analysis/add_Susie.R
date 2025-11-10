library(vroom)
library(tidyr)
library(dplyr)
library(coloc)
library(Rfast)
library(foreach)
library(doParallel)
library(data.table)
library(susieR)
options(timeout=1000000000)
options(scipen=100)

# Set up parallel backend
cores <- detectCores()
registerDoParallel(cores)

# Main directory and input arguments
main_dir <- "eQTLs/"
num=258
args <- commandArgs(trailingOnly = TRUE)
locus <- args[1]
chrom=gsub('chr', '',paste(unlist(strsplit(locus,split="\\.")[[1]][[1]])))
ld_ref=paste0("1000G_v3/1000G.EUR.hg38.",chrom)
rsid = readLines(paste0(main_dir,"position_loci/",num,"loci/",locus,".input"))

Ns=list(OneK1K=982,DICE=91,eQTL_catalogue=NULL)
traits = c("ALG" ,"AS", "CEL", "CRO", "ECZ" ,"GD" ,"HT", "IBD", "JIA", "MS" ,"PSO", "RA", "SLE", "T1D", "UC" ,"VIT")
genelist = readRDS(paste0(main_dir,"position_loci/",num,"loci_to_4125genes.rds"))
genelist = genelist[genelist$locus == locus,]
genelist = genelist$ensembl_gene_id
db=c("OneK1K","DICE","eQTL_catalogue")

try_susie <- function(dataset2) {
  result <- NULL

  # Calculate prior_variance
  prior_variance <- ifelse(dataset2$type == "cc", 0.2^2,
                             (0.15 / sqrt(sum(dataset2$MAF * (1 - dataset2$MAF) * dataset2$varbeta)))^2)

  # Try option A: runsusie with max_iter and no convergence repeat
  try_result_A <- tryCatch({
    runsusie(dataset2, max_iter = 1000, repeat_until_convergence = FALSE)
  }, error = function(e) {
    return(e)
  })

  # Check for "did not converge" error in A
  if (inherits(try_result_A, "error") && grepl("did not converge", try_result_A$message)) {
    # Try option B: susie_rss with specified parameters
    try_result_B <- tryCatch({
      susie_rss(bhat = dataset2$beta, shat = sqrt(dataset2$varbeta), R = dataset2$LD,n=dataset2$N, L = 5, tol = 1e-2)
    }, error = function(e) {
      return(e)
    })
    if (!inherits(try_result_B, "error")) {
      result <- try_result_B
    } else {
      warning("SuSiE run failed in both options A and B.")
      return(NULL)
    }
  } else {
    # Try option C: runsusie with fixed prior variance
    try_result_C <- tryCatch({
      runsusie(dataset2, prior_variance = prior_variance, estimate_prior_variance = FALSE, max_iter = 1000)
    }, error = function(e) {
      return(e)
    })

    # Check for "estimated prior variance is unreasonably large" error in C
    if (inherits(try_result_C, "error") && grepl("estimated prior variance is unreasonably large", try_result_C$message)) {
      # Try option D: susie_rss with fixed prior variance
      try_result_D <- tryCatch({
        susie_rss(bhat = dataset2$beta, shat = sqrt(dataset2$varbeta), R = dataset2$LD, n=dataset2$N,prior_variance = prior_variance, estimate_prior_variance = FALSE, L = 5, tol = 1e-2)
      }, error = function(e) {
        return(e)
      })
      if (!inherits(try_result_D, "error")) {
        result <- try_result_D
      } else {
        warning("SuSiE run failed in options C and D.")
        return(NULL)
      }
    } else if (!inherits(try_result_C, "error")) {
      result <- try_result_C
    } else if (!inherits(try_result_A, "error")) {
      # If C failed for a different reason and A succeeded
      result <- try_result_A
    } else {
      warning("SuSiE run failed in options A and C, and subsequent options were not attempted or failed.")
      return(NULL)
    }
  }

  return(result)
}

# read traits input 
foreach(trait = traits, .packages = c("vroom", "dplyr", "susieR", "coloc", "data.table")) %dopar% {
    tryCatch({
        log_file <- paste0(main_dir,"log/susie",locus,".", trait, ".out")
        sink(log_file)

        source(paste0(main_dir,"traits/",trait,"/QTL_config.R"))
        if(file.exists(paste0(main_dir,"traits/",trait,"/",locus,".rds"))){
            trait_reg = readRDS(paste0(main_dir,"traits/",trait,"/",locus,".rds"))
        }else{
            trait_reg = vroom(file=paste0(main_dir,"traits/",trait,"/",locus,".txt"),col_names=T)

            # rectify MAF
            MAF = vroom(paste0(main_dir,"traits/",trait,"/",locus,".maf"),col_names=F)
            colnames(MAF)=c("source","allele","maf","rsid")
            # Ensure minor allele frequency is between 0 and 1
            MAF <- MAF %>%
                filter(maf > 0, maf < 1) %>%
                drop_na(maf)

            # Perform vectorized matching
            complement_map <- c(A = "T", T = "A", G = "C", C = "G")
            MAF <- as_tibble(MAF)  # Ensure MAF is a tibble
            trait_reg <- trait_reg %>%
                mutate(
                    A1 = as.character(!!sym(trait_A1col)),
                    A2 = as.character(!!sym(trait_A2col)),
                    SNP_ID = as.character(!!sym(trait_SNPcol)),
                    A1_comp = ifelse(A1 %in% names(complement_map), complement_map[A1], NA),
                    A2_comp = ifelse(A2 %in% names(complement_map), complement_map[A2], NA)
                ) %>%
                left_join(
                    MAF %>%
                        rename(MAF_SNP = rsid, MAF_allele = allele, MAF_freq = maf),
                    by = c("SNP_ID" = "MAF_SNP")
                ) %>%
                filter(MAF_allele %in% c(A1, A2, A1_comp, A2_comp)) %>%
                    group_by(SNP_ID) %>%
                    slice(1) %>%  # Ensure only the first matching row is taken
                    ungroup() %>%
                    select(-A1, -A2, -A1_comp, -A2_comp, -MAF_allele) %>%
                    rename(minor_allele_freq = MAF_freq)


            # Handle missing or invalid MAF values
            if (trait_MAFcol != "") {
                trait_reg <- trait_reg %>%
                    mutate(
                    maf = ifelse(
                        is.na(!!sym(trait_MAFcol)) | !!sym(trait_MAFcol) <= 0 | !!sym(trait_MAFcol) >= 1,
                        minor_allele_freq,
                        !!sym(trait_MAFcol)
                    )
                    )
            } else {
                trait_reg <- trait_reg %>%
                    mutate(maf = minor_allele_freq)
            }

            # Remove rows with invalid MAF
            if ("maf" %in% colnames(trait_reg) && any(!is.na(trait_reg$maf) & trait_reg$maf > 0 & trait_reg$maf < 1)) {
                trait_reg <- trait_reg %>%
                    filter(!is.na(maf), maf > 0, maf < 1)
            } else {
                warning("No valid MAF values found in trait_reg. Skipping MAF filtering.")
            }

            # Ensure numeric and unique values for required columns
            if (!is.numeric(trait_reg[[trait_BPcol]])) {
                trait_reg <- trait_reg %>%
                    filter(!is.na(as.numeric(!!trait_BPcol))) %>%  # Remove non-numeric values
                    mutate(
                        !!trait_BPcol := as.numeric(!!trait_BPcol)
                    )
            }
            trait_reg <- trait_reg %>%
                distinct()

            # Handle missing or zero Beta/SE values
            if (trait_BETAcol != "") {
                min_beta <- min(abs(trait_reg[[trait_BETAcol]][trait_reg[[trait_BETAcol]] != 0]), na.rm = TRUE)
                min_se <- min(abs(trait_reg[[trait_SEcol]][trait_reg[[trait_SEcol]] != 0]), na.rm = TRUE)
                trait_reg <- trait_reg %>%
                    mutate(
                        !!trait_BETAcol := ifelse(is.na(!!sym(trait_BETAcol)) | !!sym(trait_BETAcol) == 0, min_beta, !!sym(trait_BETAcol)),
                        !!trait_SEcol := ifelse(is.na(!!sym(trait_SEcol)) | !!sym(trait_SEcol) == 0, min_se, !!sym(trait_SEcol))
                    )
            }

            # Save processed trait_reg
            rownames(trait_reg)=trait_reg[[trait_SNPcol]]
            saveRDS(trait_reg,file=paste0(main_dir,"traits/",trait,"/",locus,".rds"))   
        }

        foreach(eqtl = db, .packages = c("vroom", "dplyr", "susieR", "coloc", "data.table")) %do% {
            tryCatch({
                cells = readLines(paste0(main_dir,"databases/",eqtl,"/cells.txt"))
                N = Ns[[eqtl]]
                header = readLines(paste0(main_dir,"databases/",eqtl,"/header"))
                
                foreach(cell = cells, .packages = c("vroom", "dplyr", "susieR", "coloc", "data.table")) %do% {
                    tryCatch({
                        foreach(gene = genelist, .packages = c("vroom", "dplyr", "susieR", "coloc", "data.table")) %do% {
                            tryCatch({
                                print(paste0("Processing: ",locus," ",trait," ",eqtl," ",cell," ",gene))
                                colocInputFile = vroom(paste0(main_dir,"databases/",eqtl,"/input/",cell,".",gene,".txt"),col_names= header)
                                if("GeneSymbol" %in% header){
                                    geneSymbol = unique(colocInputFile$GeneSymbol)
                                }else{
                                    geneSymbol = tryCatch({clusterProfiler::bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
                                        }, error = function(err){
                                            print("ENSEMBL ID could not be converted to HGNC Symbol")
                                            print(paste("geneSymbol will be set to the ENSEMBL ID",gene))
                                            return(gene)
                                        })
                                }

                                susie_output_path <- paste0(main_dir, "output/susie_results/", locus, ".", trait, ".", eqtl, ".", cell, ".", geneSymbol, ".susie.rds")
                                if (!file.exists(susie_output_path)) {
                                    if (!"maf" %in% colnames(colocInputFile)) {
                                        colocInputFile <- merge(
                                            colocInputFile,
                                            trait_reg[, c("maf", trait_SNPcol)],
                                            by.x = "rsid",
                                            by.y = trait_SNPcol
                                        )
                                        colocInputFile <- unique(colocInputFile)
                                        colnames(colocInputFile)[ncol(colocInputFile)] <- "maf"
                                    }

                                    colocInputFile <- colocInputFile %>%
                                    mutate(Beta = ifelse(Beta == 0, 1e-30, Beta)) %>%
                                    mutate(varbeta = ifelse(!"varbeta" %in% colnames(.), ifelse("se" %in% colnames(.), se^2, (Beta / abs(qnorm(Pvalue / 2)))^2), varbeta)) 
                                    if (any(!is.na(colocInputFile$maf) & colocInputFile$maf > 0 & colocInputFile$maf < 1)) {
                                        colocInputFile <- colocInputFile %>%
                                        filter(!is.na(maf), maf > 0, maf < 1)
                                    }

                                    # Ensure no duplicate or invalid SNPs
                                    colocInputFile <- colocInputFile %>% distinct(rsid, .keep_all = TRUE)
                                    colocInputFile <- colocInputFile %>% filter(!duplicated(rsid), !is.na(rsid), rsid != "")
                                    rownames(colocInputFile) <- colocInputFile$rsid

                                    # Create dataset
                                    dataset1 <- list(
                                        snp = trait_reg[[trait_SNPcol]],
                                        type = traitType,
                                        s = traitProp,
                                        N = if (trait_Ncol %in% colnames(trait_reg)) trait_reg[[trait_Ncol]] else as.numeric(paste(trait_Ncol)),
                                        pvalues = if (trait_Pcol != "") {
                                            pvals <- trait_reg[[trait_Pcol]]
                                            if (any(is.na(pvals))) pvals[is.na(pvals)] <- 1
                                            pvals
                                        } else NULL,
                                        beta = if (trait_BETAcol != "") trait_reg[[trait_BETAcol]] else NULL,
                                        varbeta = if (trait_SEcol != "") trait_reg[[trait_SEcol]]^2 else NULL,
                                        MAF = trait_reg$maf
                                    )
                                    check_dataset(dataset1)

                                    dataset2 <- list(
                                        snp = colocInputFile$rsid,
                                        pvalues = colocInputFile$Pvalue,
                                        N = N,
                                        type = "quant",
                                        MAF = colocInputFile$maf,
                                        beta = if ("Beta" %in% header) colocInputFile$Beta else NULL,
                                        varbeta = if ("Beta" %in% header) colocInputFile$varbeta else NULL
                                    )
                                    check_dataset(dataset2)
                                    print("Finish dataset2")
                                    coloc_results = coloc.abf(dataset1=dataset1,dataset2=dataset2)
                                    saveRDS(coloc_results,paste0(main_dir,"output/",trait,"/",eqtl,"/",locus,".",trait,".",eqtl,".",cell,".",geneSymbol,".coloc.rds"))

                                    # Generate LD matrix using PLINK
                                    snp_set = intersect(trait_reg[[trait_SNPcol]],colocInputFile$rsid)
                                    if(length(snp_set)>=2){
                                        plink_tmp <- paste0(main_dir, "position_loci/", locus, ".", trait, ".", eqtl, ".", geneSymbol, ".tmp")
                                        write.table(snp_set, file = plink_tmp, col.names = FALSE, row.names = FALSE, quote = FALSE)
                                        plink_command <- paste0(
                                            "plink --bfile ", ld_ref, " --maf 5e-16 --nonfounders --allow-no-sex --extract ", plink_tmp,
                                            " --r2 square --write-snplist --out ", main_dir, "position_loci/plink.", locus, ".", trait, ".", eqtl, ".", geneSymbol
                                        )
                                        system(plink_command)   
                                        
                                        # Read LD matrix and SNP list
                                        ld <- fread(paste0(main_dir, "position_loci/plink.", locus, ".", trait, ".", eqtl, ".", geneSymbol, ".ld"))
                                        ld <- as.matrix(ld)
                                        snps <- readLines(paste0(main_dir, "position_loci/plink.", locus, ".", trait, ".", eqtl, ".", geneSymbol, ".snplist"))
                                        colnames(ld) <- rownames(ld) <- snps
                                        snp_set <- intersect(snp_set, snps)
                                    
                                        dataset_subset <- trait_reg[snp_set, ]
                                        dataset1 <- list(
                                            snp = dataset_subset[[trait_SNPcol]],
                                            type = traitType,
                                            s = traitProp,
                                            N = if (trait_Ncol %in% colnames(trait_reg)) dataset_subset[[trait_Ncol]] else as.numeric(paste(trait_Ncol)),
                                            pvalues = if (trait_Pcol != "") {
                                            pvals <- dataset_subset[[trait_Pcol]]
                                            if (any(is.na(pvals))) pvals[is.na(pvals)] <- 1
                                            pvals
                                            } else NULL,
                                            beta = if (trait_BETAcol != "") dataset_subset[[trait_BETAcol]] else NULL,
                                            varbeta = if (trait_SEcol != "") dataset_subset[[trait_SEcol]]^2 else NULL,
                                            MAF = dataset_subset$maf
                                        )
                                        check_dataset(dataset1)
                                    
                                        colocInputFile <- colocInputFile[snp_set, ]
                                        # Determine sample size (N)
                                        N <- if ("N" %in% header) {
                                            max(colocInputFile$N, na.rm = TRUE)
                                        } else if ("an" %in% header) {
                                            max(colocInputFile$an / 2, na.rm = TRUE)
                                        } else {
                                            N
                                        }
                                        # Create dataset2
                                        dataset2 <- list(
                                            snp = colocInputFile$rsid,
                                            pvalues = colocInputFile$Pvalue,
                                            N = N,
                                            type = "quant",
                                            MAF = colocInputFile$maf,
                                            beta = if ("Beta" %in% header) colocInputFile$Beta else NULL,
                                            varbeta = if ("Beta" %in% header) colocInputFile$varbeta else NULL
                                        )
                                        check_dataset(dataset2)
                                        print("Finish dataset2")
                                        dataset1$LD = dataset2$LD = as.matrix(ld[snp_set,snp_set])
                                        dataset1$position = dataset2$position = dataset_subset[[trait_BPcol]]        
                                    
                                        if (!("beta" %in% names(dataset1)) || !("varbeta" %in% names(dataset1))) {
                                            message("Computing beta and varbeta from p-values...")
                                            library(stats)  # For qnorm()
                                            # Compute Z-scores
                                            Z <- sign(log(dataset1$pvalues)) * qnorm(1 - dataset1$pvalues / 2)
                                            
                                            # Compute beta and varbeta and Add beta and varbeta to dataset1
                                            dataset1$beta <- Z / sqrt(2 * dataset1$MAF * (1 - dataset1$MAF) * dataset1$N * dataset1$s * (1 - dataset1$s))
                                            dataset1$varbeta <- 1 / (2 * dataset1$MAF * (1 - dataset1$MAF) * dataset1$N * dataset1$s * (1 - dataset1$s))
                                        }
                                        dataset1$N = max(dataset1$N)
                                        
                                        S3 <- try_susie(dataset1)
                                        S4 <- try_susie(dataset2)
                                    
                                        susie.res=coloc.susie(S3,S4)
                                        saveRDS(susie.res,susie_output_path)
                                    
                                    }else{
                                        message(paste("Not enough SNPs for SuSiE in cell:", cell, "for gene:", gene))
                                        # Save empty result
                                        empty_result <- list()
                                        saveRDS(empty_result, susie_output_path)
                                    } 

                                }else{
                                    message(paste("Susie output already exists for:", susie_output_path))
                                }
                            },    error = function(e) {
                                message(paste("Error processing gene:", gene, "in cell:", cell, "of",eqtl,"for",trait,"within",locus,":", e$message))
                            })
                        }
                    }, error = function(e) {
                        message(paste("Error in processing cell:", cell, "for", eqtl, "in", trait, "within", locus, ":", e$message))
                    })
                }
                rm(trait_reg,dataset1)
                gc()
            }, error = function(e) {
                message(paste("Error in processing database:", eqtl, "for", trait, "within", locus, ":", e$message))
            })
        }
        sink()
        return(NULL)
        gc()
    }, error = function(e) {
        message(paste("Error in processing trait:", trait, "for locus:", locus, ":", e$message))
        sink()  # Ensure the sink is closed even if there's an error
    })
}
# Stop parallel backend
stopImplicitCluster()