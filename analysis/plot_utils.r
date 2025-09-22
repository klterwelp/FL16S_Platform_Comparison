# data for all plots
all_ps <- readRDS(all_ps_path)

# set up analysis variable
all_ps <- all_ps %>% ps_mutate(analysis = factor(analysis, levels = c("MGX", "FL-16S", "V1-V3")))
analysis_labs <- c("V1-V3", "FL-16S", "Metagenomics")
names(analysis_labs) <- c("V1-V3", "FL-16S", "MGX")
# useful plotting functions
# all_ps is relative but we can transform it back into counts using the total_counts
f_relative_to_counts <- function(ps_rel) {
  total_reads <- sample_data(ps_rel)$total_reads
  otu_rel <- as(otu_table(ps_rel), "matrix") # extract matrix otu table
  otu_counts <- sweep(otu_rel, 2, total_reads, `*`) # multiply rel abx by read counts
  otu_counts <- otu_counts %>% round() # round to keep whole numbers (@ mgx data)
  otu_table(ps_rel) <- otu_table(otu_counts, taxa_are_rows = taxa_are_rows(ps_rel))
  return(ps_rel)
}

# function to set up phyloseq object for plots
f_plot_prep_phyloseq <- function(ps, min_prev=2, counts=FALSE, tax_select=1:7, rank="Species",
                                 min_abundance=0) {
  # filter out rare taxa by prevalence
  new_ps <- ps %>%
    tax_filter(min_prevalence = min_prev, use_counts = counts, undetected = 0,
               min_sample_abundance = min_abundance)
  # remove taxonomy columns not used by plot
  new_ps@tax_table <- new_ps@tax_table[,tax_select]
  # use tax_fix from microViz to set up unknowns as Last Taxa Rank (ex: Streptococcus Genus)
  new_ps <- new_ps %>%
    tax_fix(
      min_length = 1,
      unknowns = NA,
      sep = " ",
      anon_unique = TRUE,
      suffix_rank = "classified") %>%
    tax_agg(rank = rank, force = TRUE) # ignores convergent taxa (ie; conflicts)
}
# export plot
f_export_plot <- function(plot, out_file, width = 8, height = 4, dpi = 300) {
  # make figures folder if doesn't exist
  if (!dir.exists("../figures")) {
    dir.create("../figures")
  }
  # create full output path
  out_path <- file.path("../figures/", out_file)
  ggsave(out_path, plot = plot, width = width, height = height, dpi = dpi)
}
