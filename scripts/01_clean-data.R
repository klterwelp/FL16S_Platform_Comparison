# This script prepares the phyloseq objects for downstream analyses
# Inputs:
# data/16S_V1-V3/
#   - table.qza: feature table from QIIME2
#   - taxonomy_HOMD_V15.22.qza: taxonomy info from QIIME2
#   - V1V3_meta.tsv: sample metadata
# data/FL-16S/
#   - table.qza: feature table from QIIME2
#   - taxonomy_HOMD_V15.22.qza: taxonomy info from QIIME2
#   - FL16S_meta.tsv: sample metadata
# data/
#   - metagenomics-phy.relative.rds: relative abundance phyloseq object
# intermediates/NCBI-16S-blast-results/
#   - FL16S_NCBI-16S-blast-results.tsv: BLAST results for FL-16S
#   - V1V3_NCBI-16S-blast-results.tsv: BLAST results for V1-V3
# Outputs:
# intermediates/RDS/
#   - mgx.RDS: cleaned rel. abundance phyloseq object for metagenomics
#   - fl.RDS: cleaned abs. abundance phyloseq object for full-length 16S
#   - v1v3.RDS: cleaned abs. abundance phyloseq object for 16S V1-V3
#   These RDS objects are used in script 03
# intermediates/TSV/
#   - mgx_last_tax.tsv: last identified taxa and its rank per ASV for mgx
#   - fl_last_tax.tsv: last identified taxa and its rank per ASV for FL-16S
#   - v1v3_last_tax.tsv: last identified taxa and its rank per ASV for 16S V1-V3
#   These tables will be used in script 02


# Set up ------------------------------------------------------------------

# load libraries
library(qiime2R) # import qza files
library(tidyverse) # work with tibbles etc
library(phyloseq) # work with phyloseq objects
library(stringr) # manipulate strings, like sample names
library(skimr) # create summaries of dataframes
library(janitor) # great for cleaning up tibbles and col names
library(microbiome) # creates useful summaries of phyloseq data
library(microViz) # useful functions ex: ps_mutate
# setting input and output path variables:
v1v3_path <- "data/16S_V1-V3"
fl_path <- "data/FL-16S"
mgx_path <- "data"
int_fol <- "intermediates/NCBI-16S-blast-results"
v1v3_blast_path <- file.path(int_fol, "V1V3_NCBI-16S-blast-results.tsv")
fl_blast_path <- file.path(int_fol, "FL16S_NCBI-16S-blast-results.tsv")
out_path <- "intermediates/RDS"
out_tsv <- "intermediates/TSV"
out_mgx <- file.path(out_path, "mgx.RDS")
out_fl <- file.path(out_path, "fl.RDS")
out_v1v3 <- file.path(out_path, "v1v3.RDS")

# Functions ---------------------------------------------------------------

# Extract sample DNA from sample name
#' Create sample DNA metadata from sample name
#'
#' @param sample_names the sample names of the phyloseq object
#'
#' @returns `sample_dna` the sample dna metadata column that is the
#' unified sample names across all three datasets
f_sample_dna <- function(sample_names) {
  sample_names %>%
    str_replace_all("-", ".") %>%
    # removes V1-V3 primer info from sample name
    str_replace_all(".CH", "") %>%
    str_replace_all(".CM", "") %>%
    # remove platform from sample name
    str_replace_all(".V1V3", "") %>%
    str_replace_all(".MGX", "") %>%
    str_replace_all(".FL", "")
}

#' Shorten sample names and potentially add analysis suffix
#'
#' @param sample_names the sample names of the phyloseq object
#' @param suffix the suffix to indicate the analysis type (mgx, v1v3, fl)
#'
#' @returns `sample_names` updated to be shorter and to optionally contain the
#' analysis suffix.
f_shrt_names <- function(sample_names, suffix) {
  sample_names %>%
    str_replace_all("\\-", "\\.") %>% # replaces "-" with "."
    str_replace_all("\\_", "\\.") %>%  # replaces "_" with "."
    #str_replace_all(".CH", "") %>% # replaces "CH" with nothing
    str_replace_all(".INF.NPS", "") %>% # replaces ".INF.NPS" with nothing
    str_replace_all("L.", "") %>% # replaces L. with nothing
    str_c(suffix) # appends a suffix to the end of the names
}

#' Add genus to species column for HOMD data
#'
#' @param tax HOMD taxonomy with only "species" name in the Species column.
#'
#' @returns `tax` HOMD taxonomy with the Species column updated to
#' "Genus species"
f_add_genus <- function(tax) {
  # only add Genus to assigned species
  species_tax <- !is.na(tax[, "Genus"]) & !is.na(tax[, "Species"])
  tax[species_tax][, "Species"] <- paste(tax[species_tax][, "Genus"],
                                         tax[species_tax][, "Species"],
                                         sep = " ")
  return(tax)
}

#' Replace unknowns with NA across all datasets
#'
#' @param tax taxonomy table from phyloseq object
#'
#' @returns `tax` with the unknown list replaced with NA
f_replace_unknwn <- function(tax) {
  # set of unknowns to replace
  unknwns <- c(";[pcofgs]", " sp\\.", " bacterium",
               "HMT", "_\\[C", "_\\[O", "_\\[F", "_\\[G", "\\]\\[")
  unknwn_grep <- paste(unknwns, collapse = "|")
  tax_df <- as.data.frame(tax)
  # Replace values with NA for specified patterns
  tax_df <- apply(tax_df,
                  2,
                  function(col) ifelse(grepl(unknwn_grep, col), NA, col))
  rownames(tax_df) <- rownames(tax) # re-add sample names
  tax <- tax_table(tax_df) # return as tax table
}

#' Filter to top-scoring hits per ASV and adds consensus species classification
#'
#' @param blast_tbl the imported blast results from 00_blast-species.sh
#'
#' @returns `blast_tbl` that's been filtered to only contain matches per ASV
#' that have the maximum bitscore for that ASV. Adds the `consensus_species`
#' column which is NA if there's more than one species with a top bitscore for
#' that ASV or is the scientific name of the species with the top bitscore.
f_shrtn_hits <- function(blast_tbl) {
  blast_tbl %>%
    dplyr::group_by(query_asv) %>%
    dplyr::mutate(max_bitscore = max(bitscore),
                  sci_species = word(sciname, 1, 2)) %>%
    filter(bitscore == max_bitscore) %>%
    mutate(num_species_per_asv = length(unique(sci_species)),
           uniq_species_per_asv = toString(unique(sci_species))) %>%
    mutate(consensus_species = case_when(
      num_species_per_asv == 1 ~ sci_species,
      .default = NA
    )) %>%
    dplyr::ungroup() %>%
    dplyr::rename(ASV = query_asv)
}

#' Replace HOMD classified species with BLAST consensus species
#'
#' @param tax_table phyloseq taxonomy table
#'
#' @returns `tax_table` with Species replaced with BLAST `consensus_species`
#' from the `f_shrtn_hits` function. Adds a `HOMD_Species` column that contains
#' the HOMD classifier species. Also adds a `ASV` column which contains the ASV
#' name from the rownames of the phyloseq object.
f_replace_blast_species <- function(tax_table) {
  tax_table_ranks <- colnames(tax_table) # save table ranks
  tax_table <- tax_table %>% tax_mutate(HOMD_Species = Species) # add HOMD col
  tax_table <- tax_names2rank(tax_table, colname = "ASV") # add ASV col
  tax_table_asvs <- rownames(tax_table) # save ASV from table
  if (any(tax_table_asvs %in% blast_all_asvs)) {
    # replace matching phyloseq ASV results with blast table consensus Species
    tax_df <- as.data.frame(tax_table, stringsAsFactors = FALSE) %>%
      dplyr::select(-Species)
    tax_df <- dplyr::left_join(tax_df, blast_all) %>%
      dplyr::select(dplyr::any_of(tax_table_ranks), HOMD_Species, ASV)
    rownames(tax_df) <- tax_table_asvs
    tax_table <- as.matrix.data.frame(tax_df) %>% phyloseq::tax_table()
  }
  return(tax_table)
}

#' Get last taxon identified and last taxonomic rank from a row
#'
#' @param taxon_row taxonomy table row from phyloseq taxonomy table
#' @param ranks ranks of the current phyloseq taxonomy table
#'
#' @returns A list of `last_value` and `last_col`, where `last_value` is the
#' taxon name of the last non-NA cell and `last_col` is the rank of the last
#' identified taxon.
get_last_ided <- function(taxon_row, ranks) {
  tax_levels <- taxon_row
  non_empty_idx <- which(tax_levels != "" & !is.na(tax_levels))
  if (length(non_empty_idx) == 0) {
    return(list(last_value = NA, last_col = NA))
  } else {
    last_idx <- tail(non_empty_idx, 1) # last identified taxa
    last_value <- tax_levels[last_idx] # extract name
    last_col <- ranks[last_idx] # extract rank
    return(list(last_value = last_value, last_col = last_col))
  }
}

#' Create TSV of the last identified taxa for each ASV
#'
#' @param tax_table taxonomy table from phyloseq object
#' @param out_file TSV path to write final table
#'
#' @returns `out_file` TSV file used in 02_update_taxonomy.sh with the
#' following columns: `ASV`, `last_taxonomic_level`, and
#' `last_taxon_identified`. `last_taxonomic_level` is the last identified rank
#' for a given `ASV` while `last_taxon_identified` is the scientific name of
#' the last identified (not NA) taxon for that `ASV`.
f_write_last_tax <- function(tax_table, out_file) {
  ranks <- rank_names(tax_table)[1:7] # extract taxonomic ranks King-Species
  tax_df <- as.data.frame(tax_table, stringsAsFactors = FALSE)
  # removes anything after _ for taxons
  tax_df <- tax_df %>%
    select(-HOMD_Species) %>%
    mutate(across(where(is.character), ~ sub("_.*", "", .))) %>%
    rowwise() %>%
    mutate(last_info = list(get_last_ided(c_across(all_of(ranks)), ranks))) %>%
    transmute(
      ASV = ASV,
      last_taxonomic_level = last_info$last_col,
      last_taxon_identified = last_info$last_value
    ) %>%
    ungroup()
  write_tsv(tax_df, file = out_file, na = "")
}

# Step One: Import data ---------------------------------------------------

# importing V1-V3 data as phyloseq object
v1v3_ps <- qiime2R::qza_to_phyloseq(
  features = file.path(v1v3_path, "table.qza"),
  taxonomy = file.path(v1v3_path, "taxonomy-HOMD_V15.22.qza"),
  metadata = file.path(v1v3_path, "V1V3_meta.tsv")
)
# importing FL-16S data as phyloseq object
fl_ps <- qiime2R::qza_to_phyloseq(
  features = file.path(fl_path, "table.qza"),
  taxonomy = file.path(fl_path, "taxonomy-HOMD_V15.22.qza"),
  metadata = file.path(fl_path, "FL16S_meta.tsv")
)
# importing metagenomics data
mgx_ps <- readRDS(file.path(mgx_path, "metagenomics-phy.relative.rds"))

# Step Two: Keep samples in all three datasets ----------------------------

# set sample names for full-length to match Name col in metadata
fl_ps@sam_data$old_name <- sample_names(fl_ps)
sample_names(fl_ps) <- fl_ps@sam_data$Name
# extract names in the format expected of mgx names
names_mgx <- sample_names(mgx_ps)
names_fl <- sample_names(fl_ps) %>%
  str_replace_all("-", ".")
names_v1v3 <- sample_names(v1v3_ps) %>%
  str_replace_all("-CH|-CM", "") %>%
  str_replace_all("-", ".")
keep_names_mgx <- intersect(names_fl, names_v1v3)
# prune mgx samples and taxa that are not included in FL and V1-V3
mgx_ps <- prune_samples(samples = keep_names_mgx, x = mgx_ps)
mgx_ps <- prune_taxa(taxa_sums(mgx_ps) > 0, mgx_ps)
# prune v1v3 samples to only those that are included in FL and metagenomics
keep_names_v1v3 <- sample_names(mgx_ps) %>%
  str_replace_all("\\.", "-")
# keep both CH and CM primers
keep_names_v1v3_ch <- keep_names_v1v3 %>% str_c("-CH")
keep_names_v1v3_cm <- keep_names_v1v3 %>% str_c("-CM")
keep_names_v1v3 <- c(keep_names_v1v3_ch, keep_names_v1v3_cm)
v1v3_ps <- prune_samples(samples = keep_names_v1v3, x = v1v3_ps)
v1v3_ps <- prune_taxa(taxa_sums(v1v3_ps) > 0, v1v3_ps)
# prune fl samples to only those that are included in V1-V3 and metagenomics
keep_names_fl <- sample_names(mgx_ps) %>%
  str_replace_all("\\.", "-")
fl_ps <- prune_samples(samples = keep_names_fl, x = fl_ps)
fl_ps <- prune_taxa(taxa_sums(fl_ps) > 0, fl_ps)

# Step Three: Clean metadata and names ------------------------------------

# set old_name for V1-V3 and metagenomics (FL already added)
v1v3_ps@sam_data$old_name <- sample_names(v1v3_ps)
mgx_ps@sam_data$old_name <- sample_names(mgx_ps)
# create sample DNA metadata column from sample names
v1v3_ps@sam_data$sample_dna <- f_sample_dna(sample_names(v1v3_ps))
fl_ps@sam_data$sample_dna <- f_sample_dna(sample_names(fl_ps))
mgx_ps@sam_data$sample_dna <- f_sample_dna(sample_names(mgx_ps))
# set up analysis variable
v1v3_ps@sam_data$analysis <- "V1-V3"
fl_ps@sam_data$analysis <- "FL-16S"
mgx_ps@sam_data$analysis <- "MGX"
# set up shrt_sample_dna variable
v1v3_ps@sam_data$shrt_sample_dna <- f_shrt_names(v1v3_ps@sam_data$sample_dna, "")
fl_ps@sam_data$shrt_sample_dna <- f_shrt_names(fl_ps@sam_data$sample_dna, "")
mgx_ps@sam_data$shrt_sample_dna <- f_shrt_names(mgx_ps@sam_data$sample_dna, "")
# set up shrt_names variable
v1v3_ps@sam_data$shrt_names <- f_shrt_names(sample_names(v1v3_ps), ".V1V3")
fl_ps@sam_data$shrt_names <- f_shrt_names(sample_names(fl_ps), ".FL")
mgx_ps@sam_data$shrt_names <- f_shrt_names(sample_names(mgx_ps), ".MGX")
# update phyloseq sample names to be shrt names
sample_names(v1v3_ps) <- v1v3_ps@sam_data$shrt_names
sample_names(fl_ps) <- fl_ps@sam_data$shrt_names
sample_names(mgx_ps) <- mgx_ps@sam_data$shrt_names

# Step Four: Clean taxonomy tables ----------------------------------------

# HOMD database taxonomy uses "species" only, change to "Genus species"
fl_ps@tax_table <- f_add_genus(fl_ps@tax_table)
v1v3_ps@tax_table <- f_add_genus(v1v3_ps@tax_table)
# Deal with NA values consistently
fl_ps@tax_table <- f_replace_unknwn(fl_ps@tax_table)
v1v3_ps@tax_table <- f_replace_unknwn(v1v3_ps@tax_table)
mgx_ps@tax_table <- f_replace_unknwn(mgx_ps@tax_table)
# use BLAST NCBI taxonomy results for species level identification
# (HOMD results for the upper taxonomy levels)
# read in the BLAST results for V1-V3 and FL-16S
cols <- c("query_asv",
          "query_len",
          "seq_id",
          "seq_len",
          "alignment_len",
          "perc_id",
          "num_id",
          "evalue",
          "bitscore",
          "mismatch",
          "query_coverage",
          "gaps",
          "sciname")
blast_v1v3 <- read_tsv(v1v3_blast_path, col_names = cols)
blast_fl <- read_tsv(fl_blast_path, col_names = cols)
# shorten results to only the best hits per ASV
# then add consensus species if there's only one species per ASV
blast_fl <- blast_fl %>% f_shrtn_hits()
blast_v1v3 <- blast_v1v3 %>% f_shrtn_hits()
blast_all <- bind_rows(blast_fl, blast_v1v3) %>%
  select(ASV, consensus_species) %>%
  rename(Species = consensus_species) %>%
  distinct()
blast_all_asvs <- blast_all$ASV
# replace Species with BLAST consensus_Species
fl_ps@tax_table <- f_replace_blast_species(fl_ps@tax_table)
v1v3_ps@tax_table <- f_replace_blast_species(v1v3_ps@tax_table)
mgx_ps@tax_table <- f_replace_blast_species(mgx_ps@tax_table)

# Step Five: Export phyloseq objects --------------------------------------

# save absolute abundance list
if (!dir.exists(out_path)) {
  dir.create(out_path)
}
saveRDS(fl_ps, out_fl)
saveRDS(v1v3_ps, out_v1v3)
saveRDS(mgx_ps, out_mgx)

# Step Six: Export last identified taxon table ----------------------------

# create output directory if it doesn't exist
if (!dir.exists(out_tsv)) {
  dir.create(out_tsv)
}
# export for use in the 02_update_taxonomy.sh
f_write_last_tax(fl_ps@tax_table, file.path(out_tsv, "fl_last_tax.tsv"))
f_write_last_tax(v1v3_ps@tax_table, file.path(out_tsv, "v1v3_last_tax.tsv"))
f_write_last_tax(mgx_ps@tax_table, file.path(out_tsv, "mgx_last_tax.tsv"))
