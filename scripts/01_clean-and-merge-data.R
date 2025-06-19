# This script prepares combined phyloseq objects for downstream analyses
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
# Outputs:
# RDS/
#   - all.ps: merged relative abundance phyloseq of all 3 data types
#   - amp.ps: merged absolute abundance phyloseq of amplicon datasets
#   - abs.ps.ls: list of amplicon phyloseq objects
#   - rel.ps.ls: list of all phyloseq objects

# load libraries
library(qiime2R) # import qza files
library(tidyverse) # work with tibbles etc
library(phyloseq) # work with phyloseq objects
library(stringr) # manipulate strings, like sample names
library(skimr) # create summaries of dataframes
library(janitor) # great for cleaning up tibbles and col names
library(microbiome) # creates useful summaries of phyloseq data
library(microViz) # useful functions ex: ps_mutate

# STEP ONE: IMPORTING DATA FROM ALL THREE DATASETS

# setting path variables:
v1v3_path <- "data/16S_V1-V3"
fl_path <- "data/FL-16S"
mgx_path <- "data"
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
# remove paths from environment
rm(list = ls(pattern = ".path"))

# STEP TWO: FILTER OUT SAMPLES NOT IN ALL THREE DATASETS

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
  str_replace_all("\\.", "-") %>%
  str_c("-CH") #only keeping CH primer results
v1v3_ps <- prune_samples(samples = keep_names_v1v3, x = v1v3_ps)
v1v3_ps <- prune_taxa(taxa_sums(v1v3_ps) > 0, v1v3_ps)
# prune fl samples to only those that are included in V1-V3 and metagenomics
keep_names_fl <- sample_names(mgx_ps) %>%
  str_replace_all("\\.", "-")
fl_ps <- prune_samples(samples = keep_names_fl, x = fl_ps)
fl_ps <- prune_taxa(taxa_sums(fl_ps) > 0, fl_ps)

# STEP THREE: CLEAN METADATA

# set old_name for V1-V3 and metagenomics (FL already added)
v1v3_ps@sam_data$old_name <- sample_names(v1v3_ps)
mgx_ps@sam_data$old_name <- sample_names(mgx_ps)
# extract sample DNA from sample name
f_sample_dna <- function(sample_names) {
  sample_names %>%
    str_replace_all("-", ".") %>%
    str_replace_all(".CH", "")
}
v1v3_ps@sam_data$sample_dna <- f_sample_dna(sample_names(v1v3_ps))
fl_ps@sam_data$sample_dna <- f_sample_dna(sample_names(fl_ps))
mgx_ps@sam_data$sample_dna <- f_sample_dna(sample_names(mgx_ps))
# set up analysis variable
v1v3_ps@sam_data$analysis <- "V1-V3"
fl_ps@sam_data$analysis <- "FL-16S"
mgx_ps@sam_data$analysis <- "MGX"
# shrt_name function
f_shrt_names <- function(sample_names, suffix) {
  sample_names %>%
    str_replace_all("\\-", "\\.") %>% # replaces "-" with "."
    str_replace_all("\\_", "\\.") %>%  # replaces "_" with "."
    str_replace_all(".CH", "") %>% # replaces "CH" with nothing
    str_replace_all(".INF.NPS", "") %>% # replaces ".INF.NPS" with nothing
    str_replace_all("L.", "") %>% # replaces L. with nothing
    str_c(suffix) # appends a suffix to the end of the names
}
# shrt_sample_dna
v1v3_ps@sam_data$shrt_sample_dna <- f_shrt_names(sample_names(v1v3_ps), "")
fl_ps@sam_data$shrt_sample_dna <- f_shrt_names(sample_names(fl_ps), "")
mgx_ps@sam_data$shrt_sample_dna <- f_shrt_names(sample_names(mgx_ps), "")
# shrt_sampleNames
v1v3_ps@sam_data$shrt_names <- f_shrt_names(sample_names(v1v3_ps), ".V1V3")
fl_ps@sam_data$shrt_names <- f_shrt_names(sample_names(fl_ps), ".FL")
mgx_ps@sam_data$shrt_names <- f_shrt_names(sample_names(mgx_ps), ".MGX")
# update sample names to be shrt names
sample_names(v1v3_ps) <- v1v3_ps@sam_data$shrt_names
sample_names(fl_ps) <- fl_ps@sam_data$shrt_names
sample_names(mgx_ps) <- mgx_ps@sam_data$shrt_names