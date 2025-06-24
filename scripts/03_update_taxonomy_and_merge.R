# uses 02_update_taxonomy.sh results to update to current (6/2025) NCBI Taxonomy
# then creates merged versions of the phyloseq objects 
# Inputs:
# intermediates/TSV/
#   - mgx_last_tax.tsv: ASV, name of last taxon identified, and its rank
#   - up_mgx_tax.tsv: old taxon name, NCBI taxid, new taxon name and lineage
#   - same for fl (full-length 16S) and v1v3 (16S variable regions 1-3)
# Outputs:
# intermediates/TSV/
#   - not_updated_asvs.tsv: a list of ASVs that could not have updated taxonomy
# intermediates/RDS/
#   - all_ps: merged relative abundance phyloseq of all 3 data types
#   - amp_ps: merged absolute abundance phyloseq of amplicon datasets
#   - abs_ps_ls: list of amplicon phyloseq objects
#   - rel_ps_ls: list of all phyloseq objects


# Set up ------------------------------------------------------------------

# load libraries
library(tidyverse)
library(phyloseq)
library(microViz)
# set up folder paths for inputs and outputs
in_fol <- "intermediates/TSV"
mgx_tax_path <- file.path(in_fol, "mgx_last_tax.tsv")
mgx_lin_path <- file.path(in_fol, "up_mgx_tax.tsv")
v1v3_tax_path <- file.path(in_fol, "v1v3_last_tax.tsv")
v1v3_lin_path <- file.path(in_fol, "up_v1v3_tax.tsv")
fl_tax_path <- file.path(in_fol, "fl_last_tax.tsv")
fl_lin_path <- file.path(in_fol, "up_fl_tax.tsv")
out_no_updates_path <- file.path(in_fol, "not_updated_asvs.tsv")
out_path <- "intermediates/RDS"
out_abs_path <- file.path(out_path, "abs_ps_ls.RDS")
out_rel_path <- file.path(out_path, "rel_ps_ls.RDS")
out_all_path <- file.path(out_path, "all_ps.RDS")
out_amp_path <- file.path(out_path, "amp_ps.RDS")
in_rds <- "intermediates/RDS"


# Functions ---------------------------------------------------------------

#' Generate updated taxonomy table merged with last identified taxa
#'
#' @param last_tax_path the path to the TSV which has three columns: 
#' `ASV`, `last_taxonomic_level`, `last_taxon_identified`, where `ASV` is the 
#' ASV from the taxonomy table, `last_taxonomic_level` is the last non-NA rank
#' for that ASV, and `last_taxon_identified` is the scientific name in that
#' rank. This table was the input for the shell script that runs taxonkit to
#' update the taxonomies. 
#' @param lineage_path the path to the TSV output of the taxonkit shell script, 
#' which has the following columns: `last_taxon_identified`, `ncbi_taxid`, 
#' `last_taxonomic_level`, `full_lineage`, `ncbi_sciname`, `lineage`. 
#' `last_taxon_identified` and `last_taxonomic_level` are used to merge this
#' TSV with the `last_tax` table. `lineage` is the updated lineage that only
#' contains Kingdom, Phylum, Class, Order, Family, Genus, and Species while 
#' `full_lineage` includes the taxonomic ranks outside of these seven. 
#'
#' @returns `joined` table that combines the last identified taxa with the
#' updated lineage information. It includes all of their columns, along with 
#' `num_taxids_per_asv` which is the number of unique NCBI taxids associated 
#' with each ASV and `missing_taxid` which is TRUE/FALSE depending on whether
#' there's an NA for the `ncbi_taxid` column from the taxonkit script. 
#'
f_update_lineage <- function(last_tax_path, lineage_path) {
  # column format for updated taxonomy TSVs
  lineage_cols <- c("last_taxon_identified",
                    "ncbi_taxid",
                    "last_taxonomic_level",
                    "full_lineage",
                    "ncbi_sciname",
                    "lineage")
  lineage <- read_tsv(lineage_path, col_names = lineage_cols) # from 02 script
  last_tax <- read_tsv(last_tax_path) # from 01 script
  # set up last tax and lineage tables to be merged
  last_tax <- last_tax %>%
    mutate(last_taxonomic_level = tolower(last_taxonomic_level))
  lineage <- lineage %>%
    mutate(different = case_when(last_taxon_identified != ncbi_sciname ~ TRUE,
                                 .default = FALSE))  %>%
    separate(lineage, into = c("Kingdom",	"Phylum",	"Class",	"Order",
                               "Family",	"Genus",	"Species"), sep = ";")
  joined <- left_join(last_tax, lineage) %>%
    group_by(ASV) %>%
    mutate(num_taxids_per_asv = length(unique(ncbi_taxid))) %>%
    mutate(missing_taxid = is.na(ncbi_taxid)) %>%
    ungroup()
  return(joined)
}

#' Create list of ASVs that could not be updated to current taxonomy
#'
#' @param joined_tbl the table generated from the `f_update_lineage` function. 
#' Used to identify the ASVs that taxonkit could not find a matching NCBI 
#' taxid. 
#'
#' @returns a list of ASVs that taxonkit could not find a matching NCBI taxid 
#' for the last identified taxon. 
f_no_update <- function(joined_tbl) {
  no_update_str <- joined_tbl %>%
    filter(is.na(ncbi_taxid)) %>%
    pull(ASV) %>%
    unique()
  return(no_update_str)
}

#' Update phyloseq table taxonomy to newest NCBI taxonomy 
#'
#' @param old_table the taxonomy table from one of the phyloseq objects
#' @param new_table the joined table generated from using taxonkit to
#' identify the NCBI taxid (and its current taxonomy) associated with the last 
#' identified taxon. Created with the `f_update_lineage` function. 
#'
#' @returns `old_table` the updated taxonomy table. All ASVs with non-NA NCBI
#' taxids were updated with the current taxonomy. For the ASVs that had no 
#' taxids (and therefore NA kingdom), the old taxonomy will be used. 
f_update_taxonomy <- function(old_table, new_table) {
  new_table[new_table == "unknown"] <- NA # replace unknown with na
  # keep only relevant cols, removing duplicates
  new_table <- new_table %>% 
    select(c("Kingdom",	"Phylum",	"Class",	"Order",
             "Family",	"Genus",	"Species", "ASV")) %>%
    distinct()
  old_table <- as.data.frame(as.matrix(old_table))
  old_table <- old_table %>%
    left_join(new_table, by = "ASV", suffix = c("_old", "_new")) %>%
    mutate(
      # If Kingdom from new_table is NOT NA, use new_table values
      Kingdom = if_else(!is.na(Kingdom_new), Kingdom_new, Kingdom_old),
      Phylum  = if_else(!is.na(Kingdom_new), Phylum_new,  Phylum_old),
      Class   = if_else(!is.na(Kingdom_new), Class_new,   Class_old),
      Order   = if_else(!is.na(Kingdom_new), Order_new,   Order_old),
      Family  = if_else(!is.na(Kingdom_new), Family_new,  Family_old),
      Genus   = if_else(!is.na(Kingdom_new), Genus_new,   Genus_old),
      Species = if_else(!is.na(Kingdom_new), Species_new, Species_old),
      Updated = if_else(!is.na(Kingdom_new), TRUE, FALSE)
    ) %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species,
           HOMD_Species, ASV, Updated)
  rownames(old_table) <- old_table$ASV
  as.matrix.data.frame(old_table) %>% tax_table() #convert to tax table
}
#' Update phyloseq object with new taxonomy table
#'
#' @param old_ps old phyloseq object used in the `f_update_taxonomy` function. 
#' @param new_tax new taxonomy table generated with the `f_update_taxonomy` 
#' function. 
#' @returns `new_ps` the newly merged phyloseq object with `old_ps` otu table, 
#' `old_ps` sample data, and the `new_tax` taxonomy table. 
f_update_ps <- function(old_ps, new_tax) {
  new_ps <- merge_phyloseq(otu_table(old_ps),
                           sample_data(old_ps),
                           new_tax)
}

#' Convert feature to relative abundance
#'
#' @param x feature
#'
#' @returns `x` divided by the sum of all values of `x`
f_rel_transform <- function(x) {
  x / sum(x)
}

#' Merge metadata from phyloseq objects
#'
#' @param ls list of phyloseqs whose metadata will be combined. 
#'
#' @returns `sd_all` all metadata combined using `dplyr::bind_rows`. 
f_add_metadata <- function(ps_ls) {
  keep_cols <- c("shrt_names", "old_name", "sample_dna", "shrt_sample_dna",
                 "database", "analysis", "total_reads")
  # extract sample data as tibbles
  tib_ls <- lapply(ps_ls, function(ps) {
    sd <- sample_data(ps)
    tib <- as_tibble(sd)
    tib <- tib %>% select(any_of(keep_cols))
    return(tib)})
  # bind all tibbles and convert back to sample data
  tib_all <- bind_rows(tib_ls)
  sd_all <- sample_data(tib_all)
  sample_names(sd_all) <- tib_all$shrt_names
  sd_all <- sample_data(sd_all)
  return(sd_all)
}

#' Merge otu tables or taxonomy tables
#'
#' @param ps_ls list of phyloseq objects whose otu or taxonomy tables will be
#' combined. 
#' @param otu_or_tax whether to combine the otu or taxonomy table. Accepts
#' "otu" if combining otu tables or "tax" if combining taxonomy tables. 
#'
#' @returns `all_tables` a combined phyloseq table, either taxonomy or otu. 
f_merge_otu_or_tax <- function(ps_ls, otu_or_tax = "otu") {
  if (otu_or_tax == "otu") {
    f_extract_table <- function(ps) {otu_table(ps)}
  } else {
    f_extract_table <- function(ps) {tax_table(ps)}
  }
  table_ls <- lapply(ps_ls, f_extract_table)
  all_tables <- ""
  for (table in table_ls) {
    all_tables <- merge_phyloseq(all_tables, table)
  }
  return(all_tables)
}

# Step One: Importing data ------------------------------------------------

mgx_ps <- readRDS(file.path(in_rds, "mgx.RDS"))
v1v3_ps <- readRDS(file.path(in_rds, "v1v3.RDS"))
fl_ps <- readRDS(file.path(in_rds, "fl.RDS"))

# Step Two: Generate updated lineage tables -------------------------------

mgx_joined <- f_update_lineage(mgx_tax_path, mgx_lin_path)
v1v3_joined <- f_update_lineage(v1v3_tax_path, v1v3_lin_path)
fl_joined <- f_update_lineage(fl_tax_path, fl_lin_path)

# Step Three: Record ASVs with no associated NCBI taxid -------------------

# several taxa had no taxids according to taxonkit
# we'll keep a record of ASVs that cannot have updated taxonomy
no_updated_taxonomy_asvs <- c(f_no_update(mgx_joined),
                    f_no_update(v1v3_joined),
                    f_no_update(fl_joined))

write_tsv(x = as.data.frame(no_updated_taxonomy_asvs),
          file = out_no_updates_path)


# Step Four: Trim ( off of metagenomics Species ---------------------------

# one ASV in metagenomics has two NCBI taxids but they are the same species:
# ex: Genus species (lit 2021) = Genus species (lit 2022)
# removing the () will allow them to be merged
mgx_joined <- mgx_joined %>% 
  mutate(Species = str_trim(str_remove(Species, "\\s*\\(.*"))) %>%
  group_by(ASV) %>% 
  mutate(num_species_per_ASV = length(unique(Species)))

# Step Five: Update taxonomy tables ---------------------------------------

# extract old taxonomy tables
old_tax_v1v3 <- v1v3_ps@tax_table
old_tax_fl <- fl_ps@tax_table
old_tax_mgx <- mgx_ps@tax_table
# update taxonomy if ncbi taxid was found, otherwise use old taxonomy
new_tax_v1v3 <- f_update_taxonomy(old_tax_v1v3, v1v3_joined)
new_tax_fl <- f_update_taxonomy(old_tax_fl, fl_joined)
new_tax_mgx <- f_update_taxonomy(old_tax_mgx, mgx_joined)
# add new taxonomy table to phyloseq object
mgx_ps <- f_update_ps(mgx_ps, new_tax_mgx)
fl_ps <- f_update_ps(fl_ps, new_tax_fl)
v1v3_ps <- f_update_ps(v1v3_ps, new_tax_v1v3)

# Step Six: Export relative and absolute phyloseq lists -------------------

# add sample sums for full-length and V1-V3
fl_ps <- fl_ps %>% microViz::ps_mutate(total_reads = sample_sums(fl_ps))
v1v3_ps <- v1v3_ps %>% microViz::ps_mutate(total_reads = sample_sums(v1v3_ps))
# make out directory if it doesn't exist
if (!dir.exists(out_path)) {
  dir.create(out_path)
}
# save absolute abundance list
abs_ps_ls <- list(fl_ps = fl_ps, v1v3_ps = v1v3_ps)
saveRDS(abs_ps_ls, file = out_abs_path)
# convert to relative abundance
fl_ps_rel <- transform_sample_counts(fl_ps, fun = f_rel_transform)
v1v3_ps_rel <- transform_sample_counts(v1v3_ps, fun = f_rel_transform)
rel_ps_ls <- list(fl_ps = fl_ps_rel, v1v3_ps = v1v3_ps_rel, mgx_ps = mgx_ps)
# save relative abundance phyloseq list
saveRDS(rel_ps_ls, file = out_rel_path)

# Step Seven: Merge phyloseq objects and export ---------------------------

# merge metadata
all_sdata <- f_add_metadata(rel_ps_ls)
amp_sdata <- f_add_metadata(abs_ps_ls)
# merge tax and otu tables
all_otu <- f_merge_otu_or_tax(rel_ps_ls, otu_or_tax = "otu")
amp_otu <- f_merge_otu_or_tax(abs_ps_ls, otu_or_tax = "otu")
all_tax <- f_merge_otu_or_tax(rel_ps_ls, otu_or_tax = "tax")
amp_tax <- f_merge_otu_or_tax(abs_ps_ls, otu_or_tax = "tax")
# merge phyloseq objects
all_ps <- merge_phyloseq(all_tax, all_otu, all_sdata)
amp_ps <- merge_phyloseq(amp_tax, amp_otu, amp_sdata)
# save as RDS
saveRDS(all_ps, file = out_all_path)
saveRDS(amp_ps, file = out_amp_path)