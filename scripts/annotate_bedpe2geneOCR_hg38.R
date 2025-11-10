#R/4.0.2
library(tidyverse)
if (! require("parseIbed")){
        devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "parseIbed")
}
if (! require("DFbedtools")){
	devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "DFbedtools")
}
if (! require("vroom")){
        install.packages("vroom")
}
library(parseIbed)
library(vroom)
options(scipen = 1000)
args=commandArgs(trailingOnly=T)

bedpe_file=args[1]
atac_bed_file=args[2]
outfile_prefix=args[3]

bedpe = vroom(
        bedpe_file,
        col_names=c("chr_a","start_a","end_a","chr_b","start_b","end_b","val1","val2"),
        comment="#"
)

ocr_bed=vroom(
        atac_bed_file,
        col_names=c("chr","start","end","id"),
        comment="#"
)

if (all(is.na(ocr_bed$id))) {
        ocr_bed = ocr_bed %>%
        mutate(id=paste0("peak_",row_number()))
}

prom_bed = read.delim("gencode.v40.promoter.txt")
bedpe_geneOCR = annotate_bedpe2geneOCR(bedpe, ocr_bed, prom_bed)

bedpe_geneOCR = bedpe_geneOCR %>%
mutate_at(vars(contains("start|end")), function(x){gsub(" ", "", format(x,scientific=F))}) %>%
select(-contains("X"))



write.table(
        bedpe_geneOCR,
        file=paste0(outfile_prefix,".bedpe"),
        sep="\t", row.names=F, col.names=F, quote=F
)


savefile=bedpe_geneOCR

bedpe_geneOCR = bind_rows(savefile %>% filter(!is.na(anno_a) & !is.na(ocr_a) & is.na(anno_b) & !is.na(ocr_b)),
			  savefile %>% filter(is.na(anno_a) & !is.na(ocr_a) & !is.na(anno_b) & !is.na(ocr_b)))
gene2ocr1 = bind_rows( 
        bedpe_geneOCR %>% 
        select(ocr=ocr_a,ocrgene=anno_a, anno=anno_b),
        bedpe_geneOCR %>% 
        select(ocr=ocr_b, ocrgene=anno_b,anno=anno_a)
) %>%
filter(!is.na(ocr), !is.na(anno)) %>%
left_join(
        prom_bed %>%
        select(anno=pro_anno, gene_id) %>%
        distinct()
) %>%
separate(anno,c("gene_name","transcript"), sep="\\+") %>%
distinct()
gene2ocr1$type="ocr2openGene"


bedpe_geneOCR = savefile %>% filter(!is.na(anno_a) & !is.na(ocr_a) & !is.na(anno_b) & !is.na(ocr_b))
gene2ocr2 = bind_rows(
        bedpe_geneOCR %>%
        select(ocr=ocr_a, ocrgene=anno_a,anno=anno_b),
        bedpe_geneOCR %>%
        select(ocr=ocr_b, ocrgene=anno_b ,anno=anno_a)
) %>%
filter(!is.na(ocr), !is.na(anno)) %>%
left_join(
        prom_bed %>%
        select(anno=pro_anno, gene_id) %>%
        distinct()
) %>%
separate(anno,c("gene_name","transcript"), sep="\\+") %>%
distinct()
openProm = overlap_df(ocr_bed,prom_bed,df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end", df1_chr_col="chr", df1_start_col="start", df1_end_col="end",df1_0base=T, df2_0base=T, minoverlap=1L)
openProm = as.data.frame(openProm)
openProm$ocr=paste0(openProm$overlap_df1.chr,":",openProm$overlap_df1.start,":",openProm$overlap_df1.end,":",openProm$overlap_df1.id)
gene2ocr2$type="ocrnearopenGene2openGene"
gene2ocr2$type[which(gene2ocr2$ocr %in% openProm$ocr)]="openGene2openGene"


bedpe_geneOCR = bind_rows(savefile %>% filter(is.na(anno_a) & !is.na(ocr_a) & !is.na(anno_b) & is.na(ocr_b)),
                          savefile %>% filter(!is.na(anno_a) & is.na(ocr_a) & is.na(anno_b) & !is.na(ocr_b)))
gene2ocr3= bind_rows(
        bedpe_geneOCR %>%
        select(ocr=ocr_a,ocrgene=anno_a, anno=anno_b),
        bedpe_geneOCR %>%
        select(ocr=ocr_b, ocrgene=anno_b,anno=anno_a)
) %>%
filter(!is.na(ocr), !is.na(anno)) %>%
left_join(
        prom_bed %>%
        select(anno=pro_anno, gene_id) %>%
        distinct()
) %>%
separate(anno,c("gene_name","transcript"), sep="\\+") %>%
distinct()
gene2ocr3$type="ocr2closeGene"

bedpe_geneOCR = bind_rows(savefile %>% filter(!is.na(anno_a) & !is.na(ocr_a) & !is.na(anno_b) & is.na(ocr_b)),
                          savefile %>% filter(!is.na(anno_a) & is.na(ocr_a) & !is.na(anno_b) & !is.na(ocr_b)))
gene2ocr4= bind_rows(
        bedpe_geneOCR %>%
        select(ocr=ocr_a, ocrgene=anno_a,anno=anno_b),
        bedpe_geneOCR %>%
        select(ocr=ocr_b, ocrgene=anno_b,anno=anno_a)
) %>%
filter(!is.na(ocr), !is.na(anno)) %>%
left_join(
        prom_bed %>%
        select(anno=pro_anno, gene_id) %>%
        distinct()
) %>%
separate(anno,c("gene_name","transcript"), sep="\\+") %>%
distinct()
gene2ocr4$type="openGene2closeGene"


gene2ocr=bind_rows(gene2ocr1,gene2ocr2,gene2ocr3,gene2ocr4)

write.table(
        gene2ocr,
        file=paste0(outfile_prefix,".gene2OCR.txt"),
        sep="\t", row.names=F, col.names=F, quote=F
)

### summarize
summary_out = tibble(
        file=bedpe_file,
        total_LoopN=bedpe %>%
        select(contains("chr"),contains("start"), contains("end")) %>%
        distinct() %>%
        count() %>% pull(n),
        anno_LoopN = bedpe_geneOCR %>%
        select(contains("chr"),contains("start"), contains("end")) %>%
        distinct() %>%
        count() %>% pull(n),
        geneOCR_pairN = nrow(gene2ocr),
        gene_N = gene2ocr %>% distinct(gene_id) %>% count() %>% pull(n),
        ocr_N = gene2ocr %>% distinct(ocr) %>% count() %>% pull(n)
)

write.table(
        summary_out,
        file=paste0(outfile_prefix,".summary.txt"),
        sep="\t", row.names=F, quote=F
)


