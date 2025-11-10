#This script will intersect a proxylist given in 4 column bed format (chr, start, end, id) with the output from HiC ATAC-seq interection Chun's pipeline

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Three arguments supplied. Proxy bed file (chr, start, end, id; tab-deliminated, noheader;0base), annotation_directory (HiC pipeline output, 1based), output_prefix", call.=FALSE)
}

#install library dependencies if needed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("GenomicRanges", quietly = TRUE))
BiocManager::install("GenomicRanges")

if (!require("tidyverse", quietly = TRUE))
install.packages("tidyverse")


#Assumes bed file in 0base and .cis.anno.gene2OCR.txt in 1base
library(tidyverse)
library(GenomicRanges)

proxies_file = args[[1]]
anno_dir = args[[2]]
output_prefix = args[[3]]
outdir=args[4]

#move to annotation directory
setwd(anno_dir)

#Get files from annotation directory
anno_files = list.files(anno_dir)
anno_files= anno_files[grepl(".cis.anno.gene2OCR.txt", anno_files)]

#Get proxy file (assumes 4 column bedfile)
proxies = read.table(
	proxies_file, 
	header=F,sep="\t",
	col.names= c("chr", "start", "end", "variant_id")) %>% 
	mutate(start = start +1)

#Loop through annotation files, find SNPs located in gene connected OCRs for each file
out = lapply(anno_files, function(anno_file){
	
	anno_path = paste0(anno_dir, "/", anno_file)

	anno = read.delim(anno_path, 
		header=F,
		col.names= c("ocr","ocrgene" ,"gene_name","tx_id", "gene_id","type"))

	anno = anno %>% separate(ocr, into = c("chr", "start", "end", "ocr_id") , sep=":") 
	
	index = as.data.frame(GenomicRanges::findOverlaps(GRanges(anno), GRanges(proxies)))

	out_intersection = data.frame(anno[index[,1],], proxies[index[,2],]) %>% 
 		select(-chr.1, -start.1) %>% 
 		dplyr::rename(variant_pos = end.1, ocr_start = start, ocr_end = end) %>% 
 		dplyr::relocate(variant_id, chr, variant_pos,  variant_pos, gene_name, gene_id) %>%
 		mutate(source = gsub(".cis.anno.gene2OCR.txt", "", anno_file))
})

#Convert list to dataframe
out_df = do.call("rbind", out)

#make output file name 
out_name = paste0(output_prefix, "_v2g.txt")

#add a v2g directory and write table

write.table(out_df , file = paste0(outdir,"/", out_name), quote=F, row.names=F, sep= "\t")
