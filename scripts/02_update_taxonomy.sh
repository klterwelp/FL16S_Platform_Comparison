#!/bin/bash
# This script updates the taxonomy tables from 01 to 2025-06-01 NCBI taxonomy.
# Inputs:
# intermediates/TSV/
#   - mgx_last_tax.tsv: last identified taxa and its rank per ASV for mgx
#   - fl_last_tax.tsv: last identified taxa and its rank per ASV for FL-16S
#   - v1v3_last_tax.tsv: last identified taxa and its rank per ASV for 16S V1-V3
# references/NCBI_blast/
#   - names.dmp: NCBI taxonomy names file
#   - nodes.dmp: NCBI taxonomy nodes file
#   - from the NCBI taxonomy archive for taxdmp_2025-06-01.zip
# Environment:
#  - taxonkit conda environment from environments/taxonkit.yml
#  - requires wget installed
# Outputs:
# intermediates/TSV/
#   - up_fl_tax.tsv: updated taxonomy for FL-16S
#   - up_v1v3_tax.tsv: updated taxonomy for 16S V1-V3
#   - up_mgx_tax.tsv: updated taxonomy for MGX
#   - Format of TSV Output Files:
#       old_name (from 01), ncbi_taxid, ncbi taxonomic level,
#       full lineage, scientific NCBI name, short lineage (Kingdom;Phylum;Class;Order;Family;Genus;Species)

# Download the latest NCBI taxonomy files
NCBI_FOL="../references/NCBI_blast"
# if names.dmp already exist, skip download
if [ -f "${NCBI_FOL}/names.dmp" ]; then
    echo "NCBI taxonomy files already exist in $NCBI_FOL. Skipping download."
else
    echo "Downloading NCBI taxonomy files to $NCBI_FOL"
    wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2025-06-01.zip -P $NCBI_FOL
    unzip $NCBI_FOL/taxdmp_2025-06-01.zip -d $NCBI_FOL
    rm $NCBI_FOL/taxdmp_2025-06-01.zip
fi

# Path to input files
INT_TSV="../intermediates/TSV"
FL_TAX="${INT_TSV}/fl_last_tax.tsv"
V1V3_TAX="${INT_TSV}/v1v3_last_tax.tsv"
MGX_TAX="${INT_TSV}/mgx_last_tax.tsv"
# Path to output files
FL_OUT="${INT_TSV}/up_fl_tax.tsv"
V1V3_OUT="${INT_TSV}/up_v1v3_tax.tsv"
MGX_OUT="${INT_TSV}/up_mgx_tax.tsv"

# Activate the conda environment
# Initialize conda for non-interactive shells and activate taxonkit environment
conda activate taxonkit
# set TAXONKIT_DB to the path to the NCBI taxonomy files
export TAXONKIT_DB="$NCBI_FOL"
taxonkit version

NAMES="tmp_names.txt"
NAME2TAXID="tmp_name2taxid.tsv"
LINEAGE="tmp_lineage.tsv"

update_taxonomy() {
    local INFILE="$1"
    local OUTFILE="$2"
    tail -n +2 "$INFILE" | cut -f 3 | sort -u > "$NAMES"
    taxonkit name2taxid -n 10 -r "$NAMES" > "$NAME2TAXID"
    taxonkit lineage -i 2 -n "$NAME2TAXID" | grep -v "Eukaryota;" > "$LINEAGE"
    taxonkit reformat2 -I 2 "$LINEAGE" -r "unknown" --format "{domain};{phylum};{class};{order};{family};{genus};{species}" > "$OUTFILE"
}
# Update the taxonomy for each file
update_taxonomy "$FL_TAX" "$FL_OUT"
update_taxonomy "$V1V3_TAX" "$V1V3_OUT"
update_taxonomy "$MGX_TAX" "$MGX_OUT"

# Clean up temporary files
rm "$NAMES" "$NAME2TAXID" "$LINEAGE"
# Print completion message
echo "Taxonomy update completed."
