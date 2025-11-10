fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script_name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))
script_dir <- dirname(script_name)

suppressMessages(library(dplyr))
library(stringr)

# constants
.VERSION = "1.0.3"
CIRCOS_CONF = file.path(script_dir, "config",      "circos.conf")
CIRCOS_PATH = "circos-0.69-6-kanai/bin/circos"
# CIRCOS_DEBUG_GROUP = "text,textplace"
CIRCOS_DEBUG_GROUP = "summary"
OUTPUT_BARPLOT = TRUE
SCATTER_BACKGROUND_COLOR_ALPHA = 0.3
LARGE_POINT_SIZE = 16
SMALL_POINT_SIZE = 8

# intermediate files
COLOR_CONF =              file.path(script_dir, "config",      "color.conf")
SCATTER_BACKGROUND_CONF = file.path(script_dir, "config",      "scatter_background.conf")
CHR_DATA =                file.path(script_dir, "config", "chromosome_left.conf")
IDEO_DATA =               file.path(script_dir, "config", "ideogram_left.conf")
HIGHLIGHT_DATA =          file.path(script_dir, "data_tracks", "highlights.txt")
SCATTER_DATA =            file.path(script_dir, "data_tracks", "scatter.txt")
STACKED_DATA =            file.path(script_dir, "data_tracks", "stacked.txt")
LABEL_DATA =              file.path(script_dir, "data_tracks", "label.txt")

################################################################################
writeLines(c("*********************************************************************",
             "* Fuji plot -- a circos representation of multiple GWAS results",
     sprintf("* Version %s", .VERSION),
             "* Masahiro Kanai (mkanai@g.harvard.edu)",
             "* Harvard Medical School / RIKEN IMS / Osaka Univerisity",
             "* GNU General Public License v3",
             "*********************************************************************"
           ))

if (identical(args, character(0))) {
  args = file.path(script_dir, "input_example", c("input.txt", "traitlist.txt"))
}
input_fname = normalizePath(args[1])
traitlist_fname = normalizePath(args[2])
if (length(args) > 2){
  output_dir = normalizePath(args[3])
} else {
  output_dir = file.path(script_dir, 'output_example')
}

################################################################################
# helper func
most_common = function(x) {tail(names(sort(table(x))), 1)}

################################################################################
# load data
message("Loading input files...")

df = read.table(input_fname, T, sep = '\t', as.is = T, quote = '', comment.char = '')
traitlist = read.table(traitlist_fname, T, sep = '\t', as.is = T, quote = '', comment.char = '', fileEncoding='utf-8')

n_loci = length(unique(df$LOCUS_ID))
writeLines(c(
     sprintf("* Input data: %s", input_fname),
     sprintf("* Number of significant SNPs: %d", nrow(df)),
     sprintf("* Number of unique loci: %d", n_loci),
             "",
     sprintf("* Trait list: %s", traitlist_fname),
     sprintf("* Number of traits: %d (%s)", nrow(traitlist), str_c(traitlist$TRAIT, collapse = ',')),
     sprintf("* Number of categories: %d (%s)", length(unique(traitlist$CATEGORY)), str_c(unique(traitlist$CATEGORY), collapse = ',')),
             "",
     sprintf("* Output dir: %s", output_dir)
           ))

if ( ! dir.exists(output_dir) ){
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

input_traits = traitlist$TRAIT
if (!all(df$TRAIT %in% input_traits)) {
    missing_traits = setdiff(df$TRAIT, input_traits)
    n_missing = length(missing_traits)
    stop(sprintf("TRAIT columns mismatch.\n%d trait%s in %s %s missing from %s (%s).",
                 n_missing, ifelse(n_missing > 1, "s", ""), input_fname,
                 ifelse(n_missing > 1, "are", "is"), traitlist_fname,
                 str_c(missing_traits, collapse=",")))
}

traitlist = traitlist %>% filter(TRAIT %in% df$TRAIT) %>%
                           mutate(idx = 1:n(),  # Keep original index for ordering
                                  category_lower = str_replace_all(str_to_lower(CATEGORY), '[^a-z0-9_]', '_'), # Keep for background
                                  trait_lower = str_replace_all(str_to_lower(TRAIT), '[^a-z0-9_]', '_')) %>% # Add unique ID for trait color
                    mutate(parameters = str_c('fill_color=', trait_lower)) # Use trait_lower for point color
excluded_traits = setdiff(input_traits, traitlist$TRAIT)
writeLines(c("",
     sprintf("Excluded %d traits because of no significant SNPs (%s).", length(excluded_traits), str_c(excluded_traits, collapse = ',')),
             ""
           ))

################################################################################
message("Generating configuration and data files for circos...")

# output chr confi
chroms <- unique(df$CHR) # Get unique chromosomes
chroms <- sort(chroms)  # Optional: sort them
hs_chroms <- paste0("hs", chroms) # Add "hs" prefix
chromosomes_line <- paste("chromosomes =", paste(hs_chroms, collapse = ";"))
chromosomes_radius <- paste0("chromosomes_radius = ", paste(paste0(hs_chroms, ":1.0r"), collapse = ";"))
chromosomes_scale <- paste0("chromosomes_scale = ", paste(paste0(hs_chroms, "=1.0"), collapse = ","))
# Combine all lines
output_lines <- c(
  "chromosomes_display_default = no",
  chromosomes_line,
  "chromosomes_reverse = ",
  chromosomes_radius,
  chromosomes_scale,
  "<<include ./chromosome_order.conf>>"
)
# Write to file
writeLines(output_lines, CHR_DATA)

# output ideogram config
conf <- readLines(file.path(script_dir, "config", "ideogram.conf")) # Load file content
new_pairwise_line <- paste0("\t\t<pairwise ", hs_chroms[length(hs_chroms)], " ", hs_chroms[1], ">") # Construct new pairwise line
conf <- gsub("<pairwise hsX hs1>", new_pairwise_line, conf, fixed = TRUE) # Replace existing pairwise line
writeLines(conf, IDEO_DATA)


# output color config
str_c_comma = function(x){str_c(x, collapse = ",")}
# Define TRAIT colors (used for points/bars and now background alpha)
trait_cols = traitlist %>% select(trait_lower, COLOR) %>%
                           unique() %>%
                           mutate(rgb = apply(t(col2rgb(COLOR)), 1, str_c_comma),
                                  # Calculate alpha version for each trait color
                                  rgba = apply(floor(t((1 - SCATTER_BACKGROUND_COLOR_ALPHA) * 255 + SCATTER_BACKGROUND_COLOR_ALPHA * col2rgb(COLOR))), 1, str_c_comma))

# Generate color definitions lines
color_lines = c(
    # Add solid color definition for each trait
    sprintf("\t%s = %s", trait_cols$trait_lower, trait_cols$rgb),
    # Add alpha color definition for each trait
    sprintf("\talpha_%s = %s", trait_cols$trait_lower, trait_cols$rgba),
    # Add predefined gene label colors
    "\tgene_label_v2g  = 190,190,190",
    "\tgene_label_eqtl = 130,106,237",
    "\tgene_label_both = 246,205,48"
)

writeLines(c("<colors>",
             color_lines,
             "</colors>"), COLOR_CONF)
message(sprintf("* Color configuration (trait-specific + alpha): %s", COLOR_CONF))

################################################################################
# output scatterplot background config (Use category colors for background)
# Ensure traitlist is ordered by idx as determined earlier (usually category-grouped)
traitlist_ordered = traitlist %>% arrange(idx)

# Generate background lines for each trait
background_lines = sapply(1:nrow(traitlist_ordered), function(i) {
  trait_info = traitlist_ordered[i, ]
  # Calculate y0 and y1 based on the trait's index (idx)
  # y-axis runs from 0 (bottom) to nrow(traitlist) (top)
  # Plot value is y = nrow(traitlist) - idx
  # Background band is y +/- 0.5
  y0 = nrow(traitlist_ordered) - trait_info$idx - 0.5
  y1 = nrow(traitlist_ordered) - trait_info$idx + 0.5
  sprintf("<background>\n\tcolor = alpha_%s\n\ty0 = %.1f\n\ty1 = %.1f\n</background>",
          trait_info$trait_lower, y0, y1)
})

writeLines(c("<backgrounds>",
             background_lines,
             "</backgrounds>"), SCATTER_BACKGROUND_CONF)
message(sprintf("* Scatter background configuration (trait-specific): %s", SCATTER_BACKGROUND_CONF))


################################################################################
# output pleiotropy highlight data
nsnps_per_locus = df %>% group_by(LOCUS_ID) %>% summarize(n = n())
df = df %>% mutate(CHR = str_c("hs", CHR),
                   nsnps = nsnps_per_locus$n[match(LOCUS_ID, nsnps_per_locus$LOCUS_ID)])

inter_categorical = df %>% group_by(LOCUS_ID) %>% summarize(CHR = most_common(CHR),
                                                            BP = most_common(BP),
                                                            n = length(unique(CATEGORY))) %>%
                                                  filter(n > 1)
write.table(inter_categorical[c("CHR", "BP", "BP")], HIGHLIGHT_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Highlights data (inter-categorical pleiotropic loci): %s", HIGHLIGHT_DATA))


################################################################################
# output outer scatter plot data
scatter = merge(df, traitlist, by = "TRAIT", all.x = T)
scatter$value = nrow(traitlist) - scatter$idx
scatter$parameters = str_c(scatter$parameters, str_c('z=', scatter$nsnps), str_c('glyph_size=', ifelse(scatter$nsnps > 1, LARGE_POINT_SIZE, SMALL_POINT_SIZE)), sep = ",")
scatter = scatter[order(scatter$nsnps, decreasing=T),]
write.table(scatter[c("CHR", "BP", "BP", "value", "parameters")], SCATTER_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Scatter plot data (significant loci): %s", SCATTER_DATA))


################################################################################
# output inner stacked scatter plot data
stacked = list()
stacked_y = rep(0, n_loci)
names(stacked_y) = sort(unique(scatter$LOCUS_ID))

for (i in 1:nrow(traitlist)) {
  x = subset(scatter, idx == i)
  x$value = stacked_y[x$LOCUS_ID]
  stacked_y[x$LOCUS_ID] = stacked_y[x$LOCUS_ID] + 1
  stacked[[i]] = x
}

stacked = do.call(rbind, stacked)
stacked$parameters = str_c(stacked$parameters, str_c('z=', stacked$nsnps), sep = ",")
stacked = stacked[order(stacked$nsnps, decreasing=T),]
write.table(stacked[c("CHR", "BP", "BP", "value", "parameters")], STACKED_DATA, sep = "\t", row.names = F, col.names = F, quote = F)
message(sprintf("* Stacked bar plot data (# significant SNPs per locus): %s", STACKED_DATA))

################################################################################
# output label data for ALL genes, colored by Type (embedding RGB values)
message("Generating label data for all genes with embedded colors...")

# Ensure CHR has "hs" prefix
if (!all(startsWith(df$CHR, "hs"))) {
    df <- df %>% mutate(CHR = ifelse(startsWith(as.character(CHR), "hs"), CHR, str_c("hs", CHR)))
}

# Determine type status for each LOCUS_ID and GENE combination
gene_labels <- df %>%
  select(LOCUS_ID, GENE, Type, CHR, BP) %>%
  filter(!is.na(GENE) & GENE != "") %>%
  distinct() %>%
  group_by(LOCUS_ID, GENE) %>%
  summarize(
    CHR = most_common(CHR),
    BP = most_common(BP),
    has_v2g = "V2G" %in% Type,
    has_eqtl = "eQTL" %in% Type,
    .groups = 'drop'
  ) %>%
  mutate(
  gene_type_status = case_when(
    has_v2g & has_eqtl ~ "Both",
    has_v2g            ~ "V2G",
    has_eqtl           ~ "eQTL",
    TRUE               ~ "Other"
  ),
  # --- *** Create parameter string with NAMED colors *** ---
  label_color_param = case_when(
     gene_type_status == "Both" ~ "color=gene_label_both", # Use names like gene_label_v2g etc.
     gene_type_status == "V2G"  ~ "color=gene_label_v2g",
     gene_type_status == "eQTL" ~ "color=gene_label_eqtl",
     TRUE                      ~ "black" # Default if needed
  )
) %>%
filter(gene_type_status != "Other") # Keep filter if desired

# Check if successful
if (nrow(gene_labels) == 0) {
  warning("No gene labels found to generate label data.")
  file.create(LABEL_DATA)
} else {
  # Select columns for output: CHR, BP_start, BP_end, Label, Parameters
  label_output_prep <- gene_labels %>%
    mutate(
      BP_START = BP, # Explicitly name the start column
      BP_END = BP    # Explicitly name the end column (same as start for point label)
    )

  # --- Select ONLY the required columns in the SPECIFIC order for Circos text plots ---
  label_output_final <- label_output_prep %>%
     select(
         CHR,        # 1st column: Chromosome
         BP_START,   # 2nd column: Start position
         BP_END,     # 3rd column: End position
         GENE,       # 4th column: The label text itself
         label_color_param # 5th column: Parameters (like color=R,G,B)
     )

  write.table(label_output_final, LABEL_DATA, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message(sprintf("* Label data (all genes, with embedded RGB colors): %s", LABEL_DATA))
}

################################################################################
# call circos
cmd = sprintf("%s -conf %s %s", CIRCOS_PATH, CIRCOS_CONF, ifelse(CIRCOS_DEBUG_GROUP == "", "", sprintf("-debug_group %s", CIRCOS_DEBUG_GROUP)))
writeLines(c("",
             "Calling circos to plot...",
     sprintf("* Call: %s", cmd),
             ""
           ))
setwd(script_dir) # circos config files are specified with the relative paths
system(cmd)

# move the output file from circos to the specified location
if (output_dir != file.path(script_dir, 'output')){
  for(ext in c('png', 'svg')){
    system(sprintf("mv %s %s", file.path(script_dir, 'output', sprintf('circos.%s', ext)), file.path(output_dir, sprintf('circos.%s', ext))))
  }
}

# clean-up the intermediate files
#for(f in c(
#  COLOR_CONF,
#  SCATTER_BACKGROUND_CONF,
#  HIGHLIGHT_DATA,
#  SCATTER_DATA,
#  STACKED_DATA,
#  LABEL_DATA
#)){
#  if (file.exists(f)) file.remove(f)
#}

################################################################################
# output bar plot
if (OUTPUT_BARPLOT) {
  bar = df %>% group_by(TRAIT) %>%
               summarize(total = length(MARKER),
                         pleiotropic = length(MARKER[nsnps > 1]),
                         inter_categorical = length(MARKER[LOCUS_ID %in% inter_categorical$LOCUS_ID])) %>%
               mutate(single = total - pleiotropic,
                      intra_categorical = pleiotropic - inter_categorical)
  bar = merge(bar, traitlist, by = "TRAIT")
  bar = bar[order(bar$idx, decreasing=T),]
  rownames(bar) = bar$TRAIT

  cairo_pdf(file.path(output_dir, "barplot.pdf"), width = 8, height = 8, family = "Helvetica")
    barplot(t(bar[,c("inter_categorical", "intra_categorical", "single")]), ylim = c(0, 100), space = 0, col = c("black", "grey50", "white"))
  . = dev.off()
}

writeLines(c("", "",
             sprintf("* Final circos outputs: %s.{png,svg}.", file.path(output_dir, 'circos')),
             "",
     sprintf("Finished at %s.", Sys.time())
           ))
