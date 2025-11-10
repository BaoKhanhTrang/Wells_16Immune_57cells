### Required Packages 
require(ggnewscale)
require(patchwork)
require(dplyr)
require(Gviz)
require(GenomicRanges)
require(biomaRt)
require(ggpubr)
require(LDheatmap)
require(ggplotify)


pairGWASplot <- function(GWAS.df, eQTL.df, LD.df = TRUE, trait1,trait2,
                         chromosome, startpos, stoppos,
                    sigpvalue_GWAS = 0.05, sigpvalue_eQTL = 0.05,
                    NESeQTLRange = c(NA,NA), 
                    congruence = FALSE, R2min = 0.2, LDmin = 10, leadSNP = TRUE,
                    LDcolor = "color", ylima = NA, ylimd = NA, xlimd = NA,
                    genometrackheight = 2, gbuild = "hg38",
                    res = 300, wi = "wi", 
                    getplot = TRUE, saveplot = TRUE){
  
  ######## eQTpLot Main Function ######################################

  
  ### Check Data ########################
  print("Checking input data...")
  
  
  if(all(c("CHR", "BP", "SNP", "BETA", "P") %in% colnames(GWAS.df))==FALSE) {
    stop("The data supplied to GWAS.df must contain columns 'CHR', 'BP', 'SNP', 'BETA', and 'P'")
  }
  
  if("PHE" %in% colnames(GWAS.df)){
    print(paste(sep = "", 'PHE column found in GWAS.df. Analyzing data for phenotype ', trait1))
  } else {
    print(paste(sep = "", 'PHE column not found in GWAS.df. Assuming all data in GWAS.df is for phenotype ', trait1))
    GWAS.df$PHE <- trait1
  }
  
  if((trait1 %in% GWAS.df$PHE) == FALSE) {
    stop('Sorry, the  phenotype ', paste(trait1), ' does not exist in the PHE column of the GWAS.df dataframe. Phenotypes included in the data supplied to GWAS.df are:\n',
         paste("'",as.character(unique(GWAS.df$PHE)),"'",collapse=", ",sep=""))
  }
  
  if(all(c("SNP.Id", "P.Value", "NES") %in% colnames(eQTL.df))==FALSE) {
    stop("The data supplied to eQTL.df must contain columns 'SNP.Id', 'P.Value', 'NES'")
  }
  
  
  if(is.numeric(eQTL.df$P.Value) == FALSE | is.numeric(eQTL.df$NES) == FALSE) {
    stop('Sorry, the eQTL.df dataframe must contain only numeric data for P.Value and NES')
  }
  
  if(all("N" %in% colnames(eQTL.df))==TRUE){
    if(is.numeric(eQTL.df$N) == FALSE & is.integer(eQTL.df$N) == FALSE){
      stop('Sorry, the  column N in eQTL.df must contain only numeric values')
    }}
  
  if(is.numeric(GWAS.df$P) == FALSE | is.numeric(GWAS.df$BETA) == FALSE | (is.integer(GWAS.df$BP) == FALSE & is.numeric(GWAS.df$BP) == FALSE) | (is.integer(GWAS.df$CHR) == FALSE & is.numeric(GWAS.df$CHR) == FALSE)) {
    stop('Sorry, the GWAS.df dataframe must contain only numeric data for CHR, BP, P, and BETA (Note: chromosomes must be coded numerically)')
  }
 
  
  if(LDcolor != "color" & LDcolor != "black"){
    stop('Sorry, the argument LDcolor must be set to either "color" or "black"')
  }
  
  
  if(isTRUE(LD.df) == FALSE){
    if(all(c("BP_A", "SNP_A", "BP_B", "SNP_B", "R2") %in% colnames(LD.df))==FALSE) {
      stop("The data supplied to LD.df must contain columns 'BP_A', 'SNP_A', 'BP_B', 'SNP_B', and 'R2'")
    }
    
    if((is.integer(LD.df$BP_A) == FALSE | is.integer(LD.df$BP_B) == FALSE)) {
      stop('Sorry, the LD.df dataframe must contain only integer values for BP_A and BP_B')
    }
    
    if(is.numeric(LD.df$R2) == FALSE){
      stop('Sorry, the LD.df dataframe must contain only numeric values for R2')
    }
    
    if(isTRUE(leadSNP) == FALSE & any(leadSNP %in% LD.df$SNP_A) == FALSE & any(leadSNP %in% LD.df$SNP_B) == FALSE){
      stop('Sorry, the specified leadSNP is not present in your LD.df')
    }
    
    if(isTRUE(leadSNP) == FALSE & any(leadSNP %in% GWAS.df$SNP) == FALSE){
      stop('Sorry, the specified leadSNP is not present in your GWAS.df')
    }
    
    if(isTRUE(leadSNP) == FALSE & any(leadSNP %in% eQTL.df$SNP.Id) == FALSE){
      stop('Sorry, the specified leadSNP is not present in your eQTL.df')
    }
  }
  
  
  
  ### Compile Data ########################
  print("Compiling 2 GWAS data...")
  
  ### Subset GWAS.df for gene of interest, check to make sure data is present
  gwas.data <- GWAS.df[which(GWAS.df$CHR == chromosome & GWAS.df$BP >= startpos & GWAS.df$BP <= stoppos & !(is.na(GWAS.df$P)) & !(is.na(GWAS.df$BETA))), ]
  if(dim(gwas.data)[1] == 0) stop('Sorry, GWAS.df 1 does not contain data for any SNPs in the range ',chromosome, ':', startpos,"-",stoppos,  ' for the trait ', paste(trait1))
  if(dim(gwas.data[which(gwas.data$P <= sigpvalue_GWAS), ])[1] == 0) {
    NoFisher<-TRUE; print(paste(sep = "", 'WARNING: GWAS.df does not contain any SNPs with p-value < sigpvalue_GWAS within the range',
                                chromosome, ':', startpos,"-",stoppos, ' for the trait ',
                                trait1, '. Enrcihment Plot statistics will not be calculated'))
  } else {NoFisher <- FALSE}
  
  
  ### Subset eQTL.df for tissue and gene of interest, check to make sure data is present
 
  eqtl.data <- eQTL.df[which( !(is.na(eQTL.df$NES)) & !(is.na(eQTL.df$P.Value))), ]
  if(dim(eqtl.data)[1] == 0) stop('Sorry, GWAS.df 2 does not have any SNPs with p-value < sigpvalue_GWAS within the range',
                                  chromosome, ':', startpos,"-",stoppos)

  ### Complete eQTL meta-analysis, if requested
  eqtl.data$P.Value <- ifelse(eqtl.data$P.Value  <= 1e-300, 1e-300,eqtl.data$P.Value)
  
  if(dim(eqtl.data)[1] == 0) stop('Sorry, there are no eQTLs for the tissue', paste(tissue), ' with a p-value < sigeQTL')
  eqtl.data <- dplyr::ungroup(eqtl.data)
  
  
  ### Join GWAS and eQTL data, check to make sure there is at least some overlap in SNPs between the two datasets
  gwas.data$SNP <- as.factor(gwas.data$SNP)
  eqtl.data$SNP.Id <- as.factor(eqtl.data$SNP.Id)
  combinedSNPS <- sort(union(levels(gwas.data$SNP), levels(eqtl.data$SNP.Id)))
  Combined.eQTL.GWAS.Data <- dplyr::left_join(dplyr::mutate(gwas.data, SNP=factor(SNP, levels=combinedSNPS)),
                                              dplyr::mutate(eqtl.data, SNP.Id=factor(SNP.Id, levels=combinedSNPS)) %>% dplyr::rename(SNP = SNP.Id), by = "SNP")
  if(dim(Combined.eQTL.GWAS.Data)[1] == 0) {
    stop('Sorry, for the range ', chromosome, ':', startpos,"-",stoppos, ' there is no overlap between the SNPs in your GWAS.df 1 and 2')
  }
  
  
  ### Determine directions of effect and congruence
  Combined.eQTL.GWAS.Data$DirectionOfEffect_GWAS <- ifelse(Combined.eQTL.GWAS.Data$BETA < 0, "Negative", ifelse(Combined.eQTL.GWAS.Data$BETA > 0, "Positive", NA))
  Combined.eQTL.GWAS.Data$DirectionOfEffect_eQTL <- ifelse(Combined.eQTL.GWAS.Data$NES < 0, "Negative", ifelse(Combined.eQTL.GWAS.Data$NES > 0, "Positive", NA))
  Combined.eQTL.GWAS.Data$Congruence <- (Combined.eQTL.GWAS.Data$BETA*Combined.eQTL.GWAS.Data$NES)
  Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence < 0, "Incongruent", ifelse(Combined.eQTL.GWAS.Data$Congruence > 0, "Congruent", NA))
  
  if(congruence == FALSE) {
    Combined.eQTL.GWAS.Data$Congruence <- ifelse(Combined.eQTL.GWAS.Data$Congruence == "Incongruent", "Congruent", ifelse(Combined.eQTL.GWAS.Data$Congruence == "Congruent", "Congruent", NA))
  }
  
  ### Build final dataframe
  Combined.eQTL.GWAS.Data$NeglogeQTLpValue <- -(log10(Combined.eQTL.GWAS.Data$P.Value))
  Combined.eQTL.GWAS.Data$Neglog10pvalue_GWAS <- -(log10(Combined.eQTL.GWAS.Data$P))
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data[which(!(is.na(Combined.eQTL.GWAS.Data$P))), ]
  Combined.eQTL.GWAS.Data$significance <- ifelse(Combined.eQTL.GWAS.Data$P >= sigpvalue_GWAS, "Non-significant", "Significant")
  Combined.eQTL.GWAS.Data$Congruence2 = Combined.eQTL.GWAS.Data$Congruence
  Combined.eQTL.GWAS.Data$Congruence[which(Combined.eQTL.GWAS.Data$P.Value >= sigpvalue_eQTL)] <- "Non-Significant for eQTL"
  Combined.eQTL.GWAS.Data$Congruence2[which(Combined.eQTL.GWAS.Data$P >= sigpvalue_GWAS)] <- "Non-Significant for GWAS"
  Combined.eQTL.GWAS.Data <- Combined.eQTL.GWAS.Data %>% dplyr::mutate(Congruence=factor(Congruence, levels=c("Non-Significant for eQTL", "Congruent", "Incongruent"), ordered=TRUE))
  if(dim(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ])[1] == 0){
    Congruentdata <- FALSE} else {Congruentdata <- TRUE
    }
  if(dim(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ])[1] == 0){
    Incongruentdata <- FALSE} else {Incongruentdata <- TRUE
    }
  
  
  ### LD Calculations ########################
  if(isTRUE(LD.df) == FALSE){
    print("Compiling LD information...")
    Combined.eQTL.GWAS.Data$pvaluemult <- Combined.eQTL.GWAS.Data$P.Value*Combined.eQTL.GWAS.Data$P
    LD.df <- LD.df %>% dplyr::filter(SNP_A %in% levels(Combined.eQTL.GWAS.Data$SNP))
    LD.df <- LD.df %>% dplyr::filter(SNP_B %in% levels(Combined.eQTL.GWAS.Data$SNP))
    
    ### Select lead SNP for congruent data
    if(Congruentdata == TRUE){
      #mostsigsnp.cong <-  as.character(Combined.eQTL.GWAS.Data %>% dplyr::filter(!is.na(pvaluemult)) %>% dplyr::filter(Congruence == "Congruent") %>% dplyr::filter(pvaluemult == min(pvaluemult)) %>% dplyr::pull(SNP))
      #sample(mostsigsnp.cong, 1) -> mostsigsnp.cong
      #if(isTRUE(leadSNP) == FALSE){
      #  if((as.character(Combined.eQTL.GWAS.Data%>% dplyr::filter(SNP %in% leadSNP) %>% dplyr::pull(Congruence)) == "Congruent") == TRUE){
          mostsigsnp.cong <- leadSNP
      #  }}
      Combined.eQTL.GWAS.Data <- dplyr::left_join(Combined.eQTL.GWAS.Data, LD.df %>% dplyr::filter(SNP_A %in% mostsigsnp.cong) %>% dplyr::select(c("SNP_B","R2")), by = c("SNP" = "SNP_B"))
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "R2"] <- "R2cong"
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "SNP_B"] <- "SNP_Bcong"
      Combined.eQTL.GWAS.Data.cong <- subset(Combined.eQTL.GWAS.Data, SNP %in% mostsigsnp.cong)
      Combined.eQTL.GWAS.Data.cong=Combined.eQTL.GWAS.Data.cong[Combined.eQTL.GWAS.Data.cong$DirectionOfEffect_GWAS==Combined.eQTL.GWAS.Data.cong$DirectionOfEffect_eQTL,]
      Combined.eQTL.GWAS.Data.cong = unique(Combined.eQTL.GWAS.Data.cong)
      if(any(duplicated(Combined.eQTL.GWAS.Data.cong$SNP))){
        Combined.eQTL.GWAS.Data.cong = Combined.eQTL.GWAS.Data.cong[-which(duplicated(Combined.eQTL.GWAS.Data.cong$SNP)),]
      }
      
    }
    
    ### Select lead SNP for incongruent data    
    if(Incongruentdata == TRUE){
      #mostsigsnp.incong <-  as.character(Combined.eQTL.GWAS.Data %>% dplyr::filter(!is.na(pvaluemult)) %>% dplyr::filter(Congruence == "Incongruent") %>% dplyr::filter(pvaluemult == min(pvaluemult)) %>% dplyr::pull(SNP))
      #sample(mostsigsnp.incong, 1) -> mostsigsnp.incong
      #if(isTRUE(leadSNP) == FALSE){
      #  if((as.character(Combined.eQTL.GWAS.Data%>% dplyr::filter(SNP %in% leadSNP) %>% dplyr::pull(Congruence)) == "Incongruent") == TRUE){
          mostsigsnp.incong <- leadSNP
      #  }}
      Combined.eQTL.GWAS.Data <- dplyr::left_join(Combined.eQTL.GWAS.Data, LD.df %>% dplyr::filter(SNP_A %in% mostsigsnp.incong) %>% dplyr::select(c("SNP_B","R2")), by = c("SNP" = "SNP_B"))
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "R2"] <- "R2incong"
      names(Combined.eQTL.GWAS.Data)[names(Combined.eQTL.GWAS.Data) == "SNP_B"] <- "SNP_Bincong"
      Combined.eQTL.GWAS.Data.incong <- subset(Combined.eQTL.GWAS.Data, SNP %in% mostsigsnp.incong)
      Combined.eQTL.GWAS.Data.incong=Combined.eQTL.GWAS.Data.incong[Combined.eQTL.GWAS.Data.incong$DirectionOfEffect_GWAS!=Combined.eQTL.GWAS.Data.incong$DirectionOfEffect_eQTL,]
      Combined.eQTL.GWAS.Data.incong=unique(Combined.eQTL.GWAS.Data.incong)
      if(any(is.na(Combined.eQTL.GWAS.Data.incong$R2incong))){
        Combined.eQTL.GWAS.Data.incong=Combined.eQTL.GWAS.Data.incong[-which(is.na(Combined.eQTL.GWAS.Data.incong$R2incong)),]
      }
      if(any(duplicated(Combined.eQTL.GWAS.Data.incong$SNP))){
        Combined.eQTL.GWAS.Data.incong = Combined.eQTL.GWAS.Data.incong[-which(duplicated(Combined.eQTL.GWAS.Data.incong$SNP)),]
      }
      
    }
    
    
    ### Filter and combine LD data, generate square LD matrix
    LD.df$R2[LD.df$R2 <= R2min] = NA 
    LD.df.matrix <- as.data.frame(tidyr::spread(LD.df[(!duplicated(LD.df[,c("SNP_B","SNP_A")])),c("SNP_A","SNP_B","R2")], SNP_A, R2))
    rownames(LD.df.matrix) <- LD.df.matrix$SNP_B
    LD.df.matrix$SNP_B <- NULL
    LD.df.matrix$startpos <-  NA
    LD.df.matrix$stoppos <- NA
    dat2 <- data.frame(matrix(nrow = 2, ncol = ncol(LD.df.matrix)))
    rownames(dat2) <- c("startpos", "stoppos")
    colnames(dat2) <- colnames(LD.df.matrix)
    dplyr::bind_rows(LD.df.matrix, dat2) -> LD.df.matrix
    LD.df.matrix[,c("startpos", "stoppos")] <- NA
    matrix <- as.matrix(LD.df.matrix)
    un1 <- unique(sort(c(colnames(LD.df.matrix), rownames(LD.df.matrix))))
    matrix2 <- matrix(NA, length(un1), length(un1), dimnames = list(un1, un1))
    matrix2[row.names(matrix), colnames(matrix)] <- matrix
    matrix <- t(matrix2)
    LD.df.matrix <- dplyr::coalesce(as.data.frame(matrix), as.data.frame(matrix2))
    rownames(LD.df.matrix[rowSums(!is.na(LD.df.matrix)) >= LDmin, ]) -> SNPsWithLDData1
    c(SNPsWithLDData1,"startpos", "stoppos") -> SNPsWithLDData
    if(length(SNPsWithLDData) < 4){stop('Sorry, after filtering the LD.df data by the supplied R2min and LDmin thresholds, fewer than 2 SNPs remain that are also present in GWAS.df')}
    LD.df.matrix[rownames(LD.df.matrix) %in% SNPsWithLDData, colnames(LD.df.matrix) %in% SNPsWithLDData] -> LD.df.matrix
    LD.df.matrix[is.na(LD.df.matrix)] = 0
    SNPPositions <- unique(dplyr::bind_rows(unique(LD.df[which(LD.df$SNP_A %in% colnames(LD.df.matrix)), c('SNP_A', 'BP_A')]), unique(LD.df[which(LD.df$SNP_B %in% colnames(LD.df.matrix)), c('SNP_B', 'BP_B')] %>% dplyr::rename(SNP_A = 1, BP_A = 2))))
    SNPPositions2 <- data.frame(matrix(nrow = 2, ncol = 2))
    colnames(SNPPositions2) <- colnames(SNPPositions)
    SNPPositions2$SNP_A <- c("startpos", "stoppos")
    SNPPositions2$BP_A <- c(startpos, stoppos)
    dplyr::bind_rows(SNPPositions, SNPPositions2) -> SNPPositions
    SNPPositions[order(SNPPositions$BP_A),] -> SNPPositions
    rownames(SNPPositions) <- SNPPositions$SNP_A
    SNPPositions[order(SNPPositions$BP_A),]$SNP_A -> SNPorder
    SNPPositions[order(SNPPositions$BP_A),]$BP_A -> positions
    LD.df.matrix[SNPorder, SNPorder] -> LD.df.matrix
  }
  
  
  
  
  ### Generate main plot ########################
  print("Generating main plot...")
  
  ### Set plot limits
  if(is.na(ylima) == TRUE){ylima <- (max(Combined.eQTL.GWAS.Data %>% dplyr::select (Neglog10pvalue_GWAS), na.rm = TRUE) + 1)}
  minpos <- min(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE)
  maxpos <- max(Combined.eQTL.GWAS.Data$BP, na.rm = TRUE)
  
  ### Generate plot 1 ####
  if(is.na(ylima) == TRUE){ylima <- (max(Combined.eQTL.GWAS.Data %>% dplyr::select (Neglog10pvalue_GWAS), na.rm = TRUE) + 1)}
  
  p1 <-
    ggplot2::ggplot (data=Combined.eQTL.GWAS.Data, aes(stroke = 0)) +
    ggplot2::coord_cartesian(xlim = c(minpos, maxpos), expand = FALSE) +
    ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Non-Significant for eQTL"), shape = 15, color = "black", alpha = 0.2,
                        aes(x=BP, y=Neglog10pvalue_GWAS)) +
    ggplot2::xlab("") +
    ggplot2::ylab(bquote(-log10(P[.(trait1)]))) +
    ggplot2::scale_y_continuous(limits=c(0,ylima)) +
    ggplot2::ggtitle(paste("GWAS of ", trait1, ", colored by p-value data for ", trait2,sep = "")) +
    ggplot2::scale_shape_manual("GWAS Direction\nof Effect", values=c("Negative" = 25, "Positive" = 24), na.value = 22) +
    ggplot2::guides(alpha = FALSE,
                    size = guide_legend("Beta",
                                        override.aes = list(shape = 24, color = "black", fill ="grey"),
                                        title.position = "top", order = 2,direction = "vertical"),
                    shape = guide_legend(title.position = "top",
                                         direction = "vertical",
                                         order = 1,
                                         override.aes = list(size= 3,
                                                             fill = "grey"))) +
    ggplot2::theme(legend.direction = "horizontal", legend.key = element_rect(fill = NA, colour = NA, size = 0.25)) +
    ggplot2::geom_hline(yintercept=-log10(sigpvalue_GWAS), linetype="solid", color = "red", linewidth=0.5) +
    ggplot2::theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggplot2::theme(plot.margin = unit(c(0,1,-0.8,0), "cm"))
  
  
  ### If congruence is TRUE, add congruent and incongruent data
  if(congruence == TRUE){
    if(Congruentdata == TRUE){
      p1 <- p1 + 
        ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Congruent"),
                            alpha = 1,
                            aes(x=BP, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1,
                                shape = DirectionOfEffect_GWAS, size = abs(BETA))) +
        ggplot2::scale_fill_gradient(bquote(atop(-log10(P),paste("Congruous SNPs"))),
                                     low="#000099", high="#33FFFF",
                                     guide = guide_colorbar(title.position = "top"),
                                     limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                    dplyr::select (NeglogeQTLpValue,Neglog10pvalue_GWAS), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (NeglogeQTLpValue,Neglog10pvalue_GWAS), na.rm = TRUE)))
    }
    if(Congruentdata == TRUE & Incongruentdata == TRUE){
      p1 <- p1 + ggnewscale::new_scale_fill()
    }
    if(Incongruentdata == TRUE){
      p1 <- p1 + 
        ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Incongruent"),
                            alpha = 1,
                            aes(x=BP, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1,
                                shape = DirectionOfEffect_GWAS, size = abs(BETA))) +
        ggplot2::scale_fill_gradient(bquote(atop(-log10(P),paste("Incongruous SNPs"))),
                                     low="#990000", high="#FFCC33",
                                     guide = guide_colorbar(title.position = "top"),
                                     limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                    dplyr::select (NeglogeQTLpValue,Neglog10pvalue_GWAS), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (NeglogeQTLpValue,Neglog10pvalue_GWAS), na.rm = TRUE)))
    }
  }
  
  ### If congruence is FALSE, add all data
  if(Congruentdata == TRUE & Incongruentdata == FALSE & congruence != TRUE) {
    p1 <- p1 + 
      ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Congruent"),
                          alpha = 1,
                          aes(x=BP, y=Neglog10pvalue_GWAS, fill=NeglogeQTLpValue, alpha = 1, shape = DirectionOfEffect_GWAS, size = abs(BETA))) +
      ggplot2::scale_fill_viridis_c((bquote(-log10(P))),
                                    option = "C",
                                    guide = guide_colorbar(title.position = "top"),
                                    limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                   dplyr::select (NeglogeQTLpValue), na.rm = TRUE), max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (NeglogeQTLpValue), na.rm = TRUE)))
  }
  
  ### Add size scale
  p1 <- p1 + ggplot2::scale_size_continuous(limits = c(min(abs(c(Combined.eQTL.GWAS.Data$BETA,Combined.eQTL.GWAS.Data$NES))),
                                                       max(abs(c(Combined.eQTL.GWAS.Data$BETA,Combined.eQTL.GWAS.Data$NES)))))

  ### Add lead SNP labels, if LD.df is supplied
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE){
    p1 <- p1 + ggrepel::geom_label_repel(aes(x=BP, y=Neglog10pvalue_GWAS, label = ifelse(SNP %in% mostsigsnp.cong, paste(SNP), ""), fontface = "bold"),
                                         color = ifelse(congruence == T, "#000099", "black"), size =4, data = Combined.eQTL.GWAS.Data.cong, max.overlaps = Inf, force = 1, box.padding = 3, min.segment.length = unit(0, 'lines'))
  }
  
  if(isTRUE(LD.df) == FALSE & Incongruentdata == TRUE){
    p1 <- p1 + ggrepel::geom_label_repel(aes(x=BP, y=Neglog10pvalue_GWAS, label = ifelse(SNP %in% mostsigsnp.incong, paste(SNP), ""), fontface = "bold"),
                                         color = "#990000", size =4, data = Combined.eQTL.GWAS.Data.incong, max.overlaps = Inf, force = 1, box.padding = 3, min.segment.length = unit(0, 'lines'))
  }
  
  p1 <- p1 + theme_minimal()
  
  
  
  
  
  
  ### Generate plot 1b ####
  ylima <- (max(Combined.eQTL.GWAS.Data %>% dplyr::select (NeglogeQTLpValue), na.rm = TRUE) + 1)
  
  p1b <-
    ggplot2::ggplot (data=Combined.eQTL.GWAS.Data, aes(stroke = 0)) +
    ggplot2::coord_cartesian(xlim = c(minpos, maxpos), expand = FALSE) +
    ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence2 == "Non-Significant for GWAS1"), shape = 15, color = "black", alpha = 0.2,
                        aes(x=BP, y=NeglogeQTLpValue)) +
    ggplot2::xlab("") +
    ggplot2::ylab(bquote(-log10(P[.(trait2)]))) +
    ggplot2::scale_y_continuous(limits=c(0,ylima)) +
    ggplot2::ggtitle(paste("GWAS of ", trait2, ", colored by p-value data for ", trait1,sep = "")) +
    ggplot2::scale_shape_manual("GWAS Direction\nof Effect", values=c("Negative" = 25, "Positive" = 24), na.value = 22) +
    ggplot2::guides(alpha = FALSE,
                    size = guide_legend("Beta",
                                        override.aes = list(shape = 24, color = "black", fill ="grey"),
                                        title.position = "top", order = 2),
                    shape = guide_legend(title.position = "top",
                                         direction = "vertical",
                                         order = 1,
                                         override.aes = list(size= 3,
                                                             fill = "grey"))) +
    ggplot2::theme(legend.direction = "horizontal", legend.key = element_rect(fill = NA, colour = NA, size = 0.25)) +
    ggplot2::geom_hline(yintercept=-log10(sigpvalue_eQTL), linetype="solid", color = "red", linewidth=0.5) +
    ggplot2::theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggplot2::theme(plot.margin = unit(c(0,1,-0.8,0), "cm"))
  
  
  ### If congruence is TRUE, add congruent and incongruent data
  if(congruence == TRUE){
    if(Congruentdata == TRUE){
      p1b <- p1b + 
        ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence2 == "Congruent"),
                            alpha = 1,
                            aes(x=BP, y=NeglogeQTLpValue , fill=Neglog10pvalue_GWAS, alpha = 1,
                                shape = DirectionOfEffect_eQTL, size = abs(NES))) +
        ggplot2::scale_fill_gradient(bquote(atop(-log10(P),paste("Congruous SNPs"))),
                                     low="#000099", high="#33FFFF",
                                     guide = guide_colorbar(title.position = "top"),
                                     limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                    dplyr::select (Neglog10pvalue_GWAS,NeglogeQTLpValue), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (Neglog10pvalue_GWAS,NeglogeQTLpValue), na.rm = TRUE)))
    }
    if(Congruentdata == TRUE & Incongruentdata == TRUE){
      p1b <- p1b + ggnewscale::new_scale_fill()
    }
    if(Incongruentdata == TRUE){
      p1b <- p1b + 
        ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence2 == "Incongruent"),
                            alpha = 1,
                            aes(x=BP, y=NeglogeQTLpValue , fill=Neglog10pvalue_GWAS, alpha = 1,
                                shape = DirectionOfEffect_eQTL, size = abs(NES))) +
        ggplot2::scale_fill_gradient(bquote(atop(-log10(P),paste("Incongruous SNPs"))),
                                     low="#990000", high="#FFCC33",
                                     guide = guide_colorbar(title.position = "top"),
                                     limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                    dplyr::select (Neglog10pvalue_GWAS,NeglogeQTLpValue), na.rm = TRUE),max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (Neglog10pvalue_GWAS,NeglogeQTLpValue), na.rm = TRUE)))
    }
  }
  
  ### If congruence is FALSE, add all data
  if(Congruentdata == TRUE & Incongruentdata == FALSE & congruence != TRUE) {
    p1b <- p1b + 
      ggplot2::geom_point(data = subset(Combined.eQTL.GWAS.Data, Congruence == "Congruent"),
                          alpha = 1,
                          aes(x=BP, y=NeglogeQTLpValue, fill=Neglog10pvalue_GWAS , alpha = 1, shape = DirectionOfEffect_eQTL, size = abs(NES))) +
      ggplot2::scale_fill_viridis_c((bquote(-log10(P[.(trait1)]))),
                                    option = "C",
                                    guide = guide_colorbar(title.position = "top"),
                                    limits=c(min(Combined.eQTL.GWAS.Data %>%
                                                   dplyr::select (Neglog10pvalue_GWAS,NeglogeQTLpValue), na.rm = TRUE), max(Combined.eQTL.GWAS.Data %>%
                                                                                                          dplyr::select (Neglog10pvalue_GWAS,NeglogeQTLpValue), na.rm = TRUE)))
  }
  
  ### Add size scale
  p1b <- p1b + ggplot2::scale_size_continuous(limits = c(min(abs(c(Combined.eQTL.GWAS.Data$BETA,Combined.eQTL.GWAS.Data$NES))),
                                                         max(abs(c(Combined.eQTL.GWAS.Data$BETA,Combined.eQTL.GWAS.Data$NES)))))

  ### Add lead SNP labels, if LD.df is supplied
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE){
    p1b <- p1b + ggrepel::geom_label_repel(aes(x=BP, y=NeglogeQTLpValue, label = ifelse(SNP %in% mostsigsnp.cong, paste(SNP), ""), fontface = "bold"),
                                         color = ifelse(congruence == T, "#000099", "black"), size =4, data = Combined.eQTL.GWAS.Data.cong , max.overlaps = Inf, force = 1, box.padding = 3, min.segment.length = unit(0, 'lines'))
  }
  
  if(isTRUE(LD.df) == FALSE & Incongruentdata == TRUE){
    p1b <- p1b + ggrepel::geom_label_repel(aes(x=BP, y=NeglogeQTLpValue, label = ifelse(SNP %in% mostsigsnp.incong, paste(SNP), ""), fontface = "bold"),
                                         color = "#990000", size =4, data = Combined.eQTL.GWAS.Data.incong, max.overlaps = Inf, force = 1, box.padding = 3, min.segment.length = unit(0, 'lines'))
  }
  
  p1b <- p1b + theme_minimal() + theme(legend.position = "none")
  
  
  ### Generate Gene Tracks ########################
  print("Generating gene tracks...")
  
  ### Generate Gene Track Plot
  if(gbuild == "hg19"){hostname <- "https://grch37.ensembl.org"}
  if(gbuild == "hg38"){hostname <- "https://apr2020.archive.ensembl.org"}
  
  bm <- biomaRt::useMart(host = hostname,
                         biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl")
  
  biomTrack <- Gviz::BiomartGeneRegionTrack(genome = gbuild,
                                            chromosome = chromosome,
                                            start = minpos, 
                                            end = maxpos,
                                            filter = list("with_refseq_mrna"=TRUE),
                                            name = "ENSEMBL",
                                            background.panel="gray95",
                                            biomart = bm,
                                            margin = c(-3,-3))
  
  #gtrack <- Gviz::GenomeAxisTrack(fontcolor="#000000", fontsize=14, margin = c(-3,-3)) 
  
  genetracks <- patchwork::wrap_elements(panel = (grid::grid.grabExpr(Gviz::plotTracks(list(biomTrack),
                                                                                       collapseTranscripts = "meta",
                                                                                       transcriptAnnotation = "symbol",
                                                                                       chromosome = chromosome,
                                                                                       from = minpos,
                                                                                       to= maxpos,
                                                                                       showTitle = FALSE,
                                                                                       distFromAxis = 10,
                                                                                       innermargin = 0,
                                                                                       maxHeight = (genometrackheight*10),
                                                                                       minHeight = (genometrackheight*10),
                                                                                       sizes=genometrackheight,
                                                                                       margin = c(-3,-3)))))
  genetracks <- genetracks + theme_minimal()
  
  ### Generate LDHeatMap ########################  
  if(isTRUE(LD.df) == FALSE){
    if(length(SNPsWithLDData) > 1000 & interactive()){
      notrun <- askYesNo(default = TRUE,
                         msg = paste(sep = "", length(SNPsWithLDData), ' variants being used to generate LDHeatMap.\nUsing more than 1000 variants can take a long time to run.\nYou can increase the values for LDmin and R2min to use fewer variants.\nDo you want to continue with ', length(SNPsWithLDData),  ' variants?'))
    } else {
      notrun <- TRUE }
    
    if(notrun == FALSE){
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      print("Stopping analysis")
      stop()
    } else {
      if(notrun == TRUE){
        print(paste(sep = "", 'Generating LDHeatMap with ', length(SNPsWithLDData), ' variants'))
      }
      
      if(LDcolor == "color"){
        colorscale <-c(viridisLite::viridis(30, option = "C", direction = -1), "white")
      } else {if(LDcolor == "black"){
        colorscale <- c("grey10", "grey20", "grey30", "grey40","grey50", "grey60", "grey70", "grey80", "grey90","grey100", "white")
      }
      }
      
      
      LDheatmap::LDheatmap(as.matrix(LD.df.matrix), genetic.distances = positions, color = colorscale, flip = TRUE, add.map= TRUE, title = "", geneMapLabelX = NA, geneMapLabelY = NA, newpage = FALSE) -> LDmap
      dev.off()
    }
    
    ###Update heatmap graphics
    
    if(length(SNPsWithLDData) >= 50){
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.69, "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.69, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.69, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.69, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.69, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.69, "snpc")}
    
    if(length(SNPsWithLDData) < 50 & length(SNPsWithLDData) >= 25){
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.7, "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.7, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.7, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.7, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.7, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.7, "snpc")}
      
    if(length(SNPsWithLDData) < 25 & length(SNPsWithLDData) >= 15){
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.725, "snpc")
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.725, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.725, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.725, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.725, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.725, "snpc")}
    
    if(length(SNPsWithLDData) < 15 & length(SNPsWithLDData) >= 9){
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.75, "snpc")
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.75, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.75, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.75, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.75, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.75, "snpc")}
    
    if(length(SNPsWithLDData) < 9 & length(SNPsWithLDData) >= 5){
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.8, "snpc")
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.8, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.8, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.8, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.8, "snpc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.8, "snpc")}
    
    if(length(SNPsWithLDData) < 5){
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$width <- unit(0.875, "snpc")
      LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$height <- unit(0.875, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$width <- unit(0.875, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$height <- unit(0.875, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$width <- unit(0.875, "snpc")
      LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$height <- unit(0.875, "snpc")}
    
    LDmap$LDheatmapGrob$children$heatMap$children$heatmap$vp$y <- unit(0.875, "npc")
    LDmap$LDheatmapGrob$children$geneMap$children$diagonal$vp$y <- unit(0.875, "npc")
    LDmap$LDheatmapGrob$children$geneMap$children$segments$vp$y <- unit(0.875, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$y <- unit(0.65, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$x <- unit(0.925, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$width <- unit(0.15, "npc")
    LDmap$LDheatmapGrob$children$Key$vp$height <- unit(0.05, "npc")
    LDmap$flipVP$justification <- c('center', 'top')
    #LDmap$LDheatmapGrob$children$geneMap$children$segments <- NULL
    g.ld <- ggplotGrob(ggplotify::as.ggplot(LDmap$LDheatmapGrob) +theme(plot.margin=unit(c(-0.1,0,-0.8,0),"npc"))) 
    g.p1 <- ggplotGrob(p1)
    g.p1$widths -> g.ld$widths
  }

  
  
  ### Generate eQTL Enrichment Plot ########################
  print("Generating eQTL enrichment plot...")
  
  ###Fisher's exact test of eQTL enrichment
  Combined.eQTL.GWAS.Data$iseQTL <- Combined.eQTL.GWAS.Data$Congruence
  Combined.eQTL.GWAS.Data$iseQTL <- ifelse(Combined.eQTL.GWAS.Data$iseQTL=="Non-Significant for eQTL", "Non-Significant for eQTL", "Significant for eQTL")
  if(nrow(as.table(table(Combined.eQTL.GWAS.Data$iseQTL, Combined.eQTL.GWAS.Data$significance))) < 2 | ncol(as.table(table(Combined.eQTL.GWAS.Data$iseQTL, Combined.eQTL.GWAS.Data$significance))) < 2){
    NoFisher <- TRUE; print ("Not enough data to compute enrichment significance for Plot C")
  }
  if(NoFisher==FALSE){fisher <- fisher.test(table(Combined.eQTL.GWAS.Data$iseQTL,
                                                  Combined.eQTL.GWAS.Data$significance))}
  if(NoFisher==FALSE){fpvalue <- fisher$p.value}
  
  ### Generate plot 2
  p2 <- ggplot2::ggplot(Combined.eQTL.GWAS.Data) +
    ggplot2::aes(x=significance, y = 1, fill = (Congruence)) +
    ggplot2::geom_bar(stat = 'identity', position = "fill") +
    ggplot2::ggtitle(paste("Enrichment of ",trait2," among\n",trait1,"GWAS-significant SNPs")) +
    ggplot2::ylab(paste("Proportion of", trait1 ,"SNPs\nthat are significant for",trait2)) +
    ggplot2::xlab(paste(trait1,"GWAS significance\n(threshold p <", sigpvalue_GWAS,")"))+
    ggplot2::ylim(0,1.2) +
    if(NoFisher==FALSE){ggpubr::geom_signif(y_position=c(1.1,1.1),
                                            xmin=c("Non-significant"),
                                            xmax=c("Significant"),
                                            annotation=(paste("p =", formatC(fpvalue, format = "e", digits = 2))),
                                            tip_length=0.05)}
  
  if(Congruentdata == TRUE & Incongruentdata == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c("Congruent" = "Congruent eQTL", "Incongruent" = "Incongruent eQTL", "Non-Significant for eQTL" = "Non-Significant for eQTL"), values = c("Congruent" = "#000099", "Incongruent" = "#990000", "Non-Significant for eQTL" = "#C0C0C0")) + 
      ggplot2::guides(fill = guide_legend(title = NULL))}
  
  if(Congruentdata == TRUE & Incongruentdata == FALSE & congruence == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c("Congruent" = "Congruent eQTL", "Non-Significant for eQTL" = "Non-Significant for eQTL", ), values = c("Congruent" = "#000099", "Non-Significant for eQTL" = "#C0C0C0")) +
      ggplot2::guides(fill = guide_legend(title = NULL))}
  
  if(Congruentdata == TRUE & Incongruentdata == FALSE & congruence == FALSE){
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c("Congruent" = "Congruent eQTL", "Non-Significant for eQTL" = "Non-Significant for eQTL"), values = c("Congruent" = "#ffee00", "Non-Significant for eQTL" = "#360f70")) + 
      ggplot2::guides(fill = guide_legend(""))}
  
  if(Congruentdata == FALSE & Incongruentdata == TRUE){
    p2 <- p2 + ggplot2::scale_fill_manual(labels = c("Incongruent" = "Incongruent eQTL", "Non-Significant for eQTL" = "Non-Significant for eQTL"), values = c("Incongruent" = "#990000", "Non-Significant for eQTL" = "#C0C0C0")) +
      ggplot2::guides(fill = guide_legend(title = NULL))}
  
  p2 <- p2 + theme_minimal()
  
  
  ### Generate P-P plot ########################
  print("Generating P-P plot...")
  
  p3 <- ggplot2::ggplot(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES)) , ]) +
    ggplot2::guides(color = guide_legend("Direction of Effect")) +
    ggplot2::xlab((bquote(-log10(P[.(trait2)])))) +
    ggplot2::ylab((bquote(-log10(P[.(trait1)])))) +
    ggplot2::scale_y_continuous(limits=c(NA,ylimd))+
    ggplot2::scale_x_continuous(limits=c(NA,xlimd))+
    ggplot2::theme(legend.position = "right") +
    theme(legend.spacing.y = unit(0.1, "cm")) + 
    theme(legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25))
  
  ### Add data and annotations to plot 3 if LD info is supplied
  ### For congruent data
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]) >=2){
    pearson.congruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$NeglogeQTLpValue,
                                  Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$Neglog10pvalue_GWAS,
                                  method = "pearson")
    
    p3.1 <- p3 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Congruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == TRUE), ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2cong), shape = 21, size = 3)  + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Congruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == FALSE) , ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2cong), shape = 21, size = 3)  +
      ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                           aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), color = "black", method = "lm", formula = (y ~ x)) 
    
    p3.1 <- p3.1 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data.cong,
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, shape = SNP), size = 3, fill = "#33FF33") +  guides(shape = guide_legend(order = 2, title = NULL, override.aes = list(size = 3, shape = 23, fill = "#33FF33")))
    p3.1 <- p3.1 + ggplot2::geom_text(x= -Inf, y= Inf, 
                                      label = paste(sep = "", "r = ", round(pearson.congruent$estimate, 2), "\np = ", 
                                                    formatC(pearson.congruent$p.value, format = "e", digits = 2)), hjust = -0.05, vjust= 1.1, color = "black")
    p3.1 <- p3.1 + theme_minimal()
    
  } else {if(isTRUE(LD.df) == FALSE){
    print("Not enough data to generate P-P plot for Congruent eQTLs")
  }}
  
  ### For incongruent data
  if(isTRUE(LD.df) == FALSE & Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]) >=2){
    pearson.incongruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$NeglogeQTLpValue,
                                    Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$Neglog10pvalue_GWAS,
                                    method = "pearson")
    
    p3.2 <- p3 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Incongruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == TRUE) , ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2incong), shape = 21, size = 3)  + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Incongruent" & is.na(Combined.eQTL.GWAS.Data$R2cong) == FALSE) , ],
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, fill=R2incong), shape = 21, size = 3)  +
      ggplot2::scale_fill_gradient(bquote(R^2 ~ "with" ~{.(mostsigsnp.incong)}), limits = c(0,1), breaks = c(0.2,0.4,0.6,0.8), low="#990000", high="#FFCC33", na.value = "grey80", guide = guide_colorbar(order =1, direction = "horizontal",title.position = "top", label.position = "bottom" )) +
      ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which(Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                           aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), color = "black", method = "lm", formula = (y ~ x)) 
    
    p3.2 <- p3.2 + 
      ggplot2::geom_point(data=Combined.eQTL.GWAS.Data.incong,
                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, shape = SNP), size = 3, fill = "#33FF33") + guides(shape = guide_legend(order = 2, title = NULL, override.aes = list(size = 3, shape = 23, fill = "#33FF33"))) +
      ggplot2::geom_text(x= -Inf, y= Inf, 
                         label = paste(sep = "", "r = ", round(pearson.incongruent$estimate, 2), "\np = ", 
                                       formatC(pearson.incongruent$p.value, format = "e", digits = 2)), hjust = -0.05, vjust = 1.1, color = "black") +
      ggplot2::ggtitle(paste("P-P plot, incongruent SNPs"))
    p3.2 <- p3.2 + theme_minimal()
    
  } else {if(isTRUE(LD.df) == FALSE & congruence == TRUE){
    print("Not enough data to generate P-P plot for Incongruent eQTLs")
  }}
  
  ### Add color scale to plot 3.1 if LD info is supplied
  if(isTRUE(LD.df) == FALSE & congruence == FALSE){p3.1 <- p3.1 +ggplot2::scale_fill_viridis_c(bquote(R^2 ~ "with" ~{.(mostsigsnp.cong)}), limits = c(0,1), breaks = c(0.2,0.4,0.6,0.8), option = "C", na.value = "grey80", guide = guide_colorbar(order = 1, direction = "horizontal",title.position = "top", label.position = "bottom", title.vjust = 0 )) +
    ggplot2::ggtitle(paste("P-P plot"))
  p3.1 <- p3.1 + theme_minimal()
  
  } else {
    if(isTRUE(LD.df) == FALSE & congruence == TRUE){p3.1 <- p3.1 + ggplot2::scale_fill_gradient(bquote(R^2 ~ "with" ~{.(mostsigsnp.cong)}), limits = c(0,1), breaks = c(0.2,0.4,0.6,0.8), low="#000099", high="#33FFFF", na.value = "grey80", guide = guide_colorbar(order =1, direction = "horizontal",title.position = "top", label.position = "bottom", title.vjust = 0)) + 
      ggplot2::ggtitle(paste("P-P plot, congruent SNPs"))
    }}
  
  ### Add data and annotations to plot 3 if LD data is not supplied
  if(isTRUE(LD.df) == TRUE){
    
    ### Annotations for congruent data points
    if(Congruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]) >=2){
      pearson.congruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$NeglogeQTLpValue,
                                    Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ]$Neglog10pvalue_GWAS,
                                    method = "pearson")
      
      ### Regression line for congruent data points  
      if(congruence == TRUE){
        p3 <- p3 + ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                                        aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence), method = "lm", formula = (y ~ x))
      } else {
        if(congruence == FALSE)     
          p3 <- p3 + ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                                          aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence), method = "lm", formula = (y ~ x), show.legend = FALSE, color = "#ffee00")
      }
      
      ### Data for congruent data points
      if(congruence == TRUE){
        p3 <- p3 +
          ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                              aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence)) 
      } else {
        p3 <- p3 + 
          ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Congruent") , ],
                              aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue), color = "#360f70") 
        p3 <- p3 + 
          ggplot2::geom_text(x= -Inf, y= Inf,
                             label = paste(sep = "", "r = ", 
                                           round(pearson.congruent$estimate, 3),
                                           "\np = ", 
                                           formatC(pearson.congruent$p.value, format = "e", digits = 2)),
                             color="#000099", hjust = -0.05, vjust= 1.1)
      }
    } else {
      print("Not enough data to generate P-P plot for Congruent eQTLs")
    }
    
    ### Annotations and data for incongruent data points
    if(Incongruentdata == TRUE & nrow(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]) >=2){
      pearson.incongruent <- cor.test(Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$NeglogeQTLpValue,
                                      Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ]$Neglog10pvalue_GWAS,
                                      method = "pearson")
      p3 <- p3 + ggplot2::geom_point(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                                     aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence)) +
        ggplot2::geom_smooth(data=Combined.eQTL.GWAS.Data[ which( !is.na(Combined.eQTL.GWAS.Data$NES) & Combined.eQTL.GWAS.Data$Congruence == "Incongruent") , ],
                             aes(y=Neglog10pvalue_GWAS, x=NeglogeQTLpValue, color=Congruence), method = "lm", formula = (y ~ x)) + 
        ggplot2::geom_text(x= -Inf, y= Inf,
                           label = paste(sep = "", "\n\nr = ",
                                         round(pearson.incongruent$estimate, 3),
                                         "\np = ",
                                         formatC(pearson.incongruent$p.value, format = "e", digits = 2)),
                           color="#990000", hjust = -0.05, vjust= 1.1)
      p3 <- p3 + 
        ggplot2::geom_text(x= -Inf, y= Inf,
                           label = paste(sep = "", "r = ", 
                                         round(pearson.congruent$estimate, 3),
                                         "\np = ", 
                                         formatC(pearson.congruent$p.value, format = "e", digits = 2)),
                           color="#000099", hjust = -0.05, vjust= 1.1)
    } else {if(congruence == "TRUE"){
      print("Not enough data to generate P-P Plot for Incongruent eQTLs")
    }}
    
    ### Add correct color scales 
    if(Congruentdata == TRUE & Incongruentdata == TRUE){p3 <- p3 + scale_color_manual(values = c("Congruent" = "#000099", "Incongruent" = "#990000")) }
    
    if(Congruentdata == TRUE & Incongruentdata == FALSE){p3 <- p3 + scale_color_manual(values = c("Congruent" = "#000099"))}
    
    if(Congruentdata == FALSE & Incongruentdata == TRUE){p3 <- p3 + scale_color_manual(values = c("Incongruent" = "#990000"))} 
    
    p3 <- p3 + ggplot2::ggtitle(paste("P-P plot"))
  }
  
  
  ### Combine plots ########################
  print("Merging and plotting...")

  
  if(isTRUE(LD.df) == FALSE & Incongruentdata == FALSE){
    p4 <- ((p1+p1b + genetracks + g.ld + patchwork::plot_layout(nrow = 4, byrow = FALSE, heights = c(2,2,(genometrackheight/5),1))) |
             (p2 / patchwork::plot_spacer()  / p3.1 / patchwork::plot_spacer() / patchwork::plot_spacer() + patchwork::plot_layout(heights = c(5,0.5,5,0.5,5)))) +
      patchwork::plot_layout(ncol = 2, widths = c(2.5,1)) +
      patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12), plot.tag.position = c(0, 0.995))
  }
  
  if(isTRUE(LD.df) == FALSE & Congruentdata == FALSE){
    p4 <- ((p1 + p1b + genetracks + g.ld + patchwork::plot_layout(nrow = 4, byrow = FALSE, heights = c(2,2,(genometrackheight/5),1))) |
             (p2 / patchwork::plot_spacer() / p3.2 / patchwork::plot_spacer() / patchwork::plot_spacer() + patchwork::plot_layout(heights = c(5,0.5,5,0.5,5)))) +
      patchwork::plot_layout(ncol = 2, widths = c(2.5,1)) +
      patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12), plot.tag.position = c(0, 0.995))
  }
  
  if(isTRUE(LD.df) == FALSE & Congruentdata == TRUE & Incongruentdata == TRUE){
    p4 <- ((p1  +p1b+genetracks + g.ld + patchwork::plot_layout(nrow = 4, byrow = FALSE, heights = c(3,3,(genometrackheight/5),4),guides = "collect") &
              theme(legend.direction = "horizontal",legend.justification = "top", legend.key = element_rect(fill = NA, colour = NA, size = 0.25) ) ) |
             (p2 / patchwork::plot_spacer() / p3.1 / patchwork::plot_spacer() / p3.2) + patchwork::plot_layout(heights = c(5,0.5,5,0.5,5))) +
      patchwork::plot_layout(ncol = 2, widths = c(2,1)) +
      patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12), plot.tag.position = c(0, 0.995))
    
    }
  
  
  if(isTRUE(LD.df) == TRUE) {
    p4 <- (p1 + p1b+ genetracks + patchwork::plot_spacer() + (p2 + p3 + patchwork::plot_layout(ncol = 2, widths = c(2,3))) + 
             patchwork::plot_layout(ncol=1, height = c(4,4,genometrackheight, 0.5, 3)) + 
             patchwork::plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 18), text = element_text(size = 12)))
  }
  
  pfinal <- p4 + patchwork::plot_annotation(title = paste("pairwise GWAS analysis for ", trait1," and ", trait2, "\nin locus ", chromosome,":", minpos,"-", maxpos, "\n", sep = ""), theme=theme(plot.title = element_text(size = 19, face = "bold")))
  plots <- list(p1 = p1, p1b = p1b, ld = g.ld, p2 = p2, p3 = p3, genetracks = genetracks)
  
  
  ### Export final plot ########################  
  if(isTRUE(LD.df) == FALSE & wi == "wi"){
    wi <- 16
  }
  if(isTRUE(LD.df) == TRUE & wi == "wi"){
    wi <-14
  }
  if(isTRUE(LD.df) == FALSE){
    hgt <- ((1.3*(wi) - 0.01*(wi)^2 - 5.5))
  }
  if(isTRUE(LD.df) == TRUE){
    hgt <- wi*0.85
  }
  if(congruence == TRUE){congruence <- "WithCongruenceData"} else {congruence <- "WithoutCongruenceData"}
  if(isTRUE(LD.df) == FALSE){LDinfo <- "WithLinkageData"} else {LDinfo <- "WithoutLinkageData"}
  if(saveplot == TRUE){
    ggsave(pfinal, filename=paste(chromosome,startpos, stoppos, trait1, trait2, congruence, LDinfo,"GWASpair_plot" ,"pdf", sep="."), dpi=res, units="in", height=hgt, width=wi)
    }
  if(getplot == TRUE){
    return(plots)}
  
}
