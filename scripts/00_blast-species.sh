#!/bin/bash
# This script performs BLAST searches against the NCBI 16S ribosomal RNA database
# Used to assign species-level taxonomy to the FL-16S and V1-V3 16S reference sequences.
# Inputs:
#   - data/FL-16S/rep-seqs.qza
#   - data/16S_V1-V3/rep-seqs.qza
#   - references/NCBI_blast/16S_ribosomal_RNA.tar.gz
#       - version May 24, 2025 via:
#       - wget "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz"
# Environment:
#   - BLAST 2.16 conda environment from environments/blast.yml
# Outputs:
#   - V1V3_NCBI-16S-blast-results.tsv: the BLAST hit results for the V1-V3 region
#   - FL16S_NCBI-16S-blast-results.tsv: the BLAST hit results for the FL-16S region
# all of these TSVs will have the following BLAST columns:
#   - qseqid
#   - qlen
#   - sseqid
#   - slen
#   - length
#   - pident
#   - nident
#   - evalue
#   - bitscore
#   - mismatch
#   - qcovs
#   - gaps
#   - sscinames
# Number of threads
NUM_THREADS=10
# Locations of the input files
V1V3_REP_SEQS="../data/16S_V1-V3/rep-seqs.qza"
FL16S_REP_SEQS="../data/FL-16S/rep-seqs.qza"
BLAST_TAR="../references/NCBI_blast/16S_ribosomal_RNA.tar.gz"
# Location of the output files
V1V3_FASTA="../data/16S_V1-V3/rep-seqs.fasta"
FL16S_FASTA="../data/FL-16S/rep-seqs.fasta"
BLAST_DB="../references/NCBI_blast/16S_ribosomal_RNA"
OUT_FOL="../intermediates/NCBI-16S-blast-results"
V1V3_BLAST_OUTPUT="${OUT_FOL}/V1V3_NCBI-16S-blast-results.tsv"
FL_BLAST_OUTPUT="${OUT_FOL}/FL16S_NCBI-16S-blast-results.tsv"
# Check that the input files exist
if [ ! -f $V1V3_REP_SEQS ]; then
    echo "Error: $V1V3_REP_SEQS does not exist."
    exit 1
fi
if [ ! -f $FL16S_REP_SEQS ]; then
    echo "Error: $FL16S_REP_SEQS does not exist."
    exit 1
fi
if [ ! -f $BLAST_TAR ]; then
    echo "Error: $BLAST_TAR does not exist."
    exit 1
fi
# Extract fasta sequences from the QZA files (if not already done)
TMPDIR="tmp_unzip"
if [ ! -f $V1V3_FASTA ]; then
    echo "Extracting V1-V3 sequences from $V1V3_REP_SEQS"
    mkdir -p $TMPDIR
    unzip -q $V1V3_REP_SEQS -d $TMPDIR
    find $TMPDIR -name "dna-sequences.fasta" -exec cp {} $V1V3_FASTA \;
    rm -rf $TMPDIR
fi
if [ ! -f $FL16S_FASTA ]; then
    echo "Extracting FL-16S sequences from $FL16S_REP_SEQS"
    mkdir -p $TMPDIR
    unzip -q $FL16S_REP_SEQS -d $TMPDIR
    find $TMPDIR -name "dna-sequences.fasta" -exec cp {} $FL16S_FASTA \;
    rm -rf $TMPDIR
fi
# Untar the BLAST database (if not already done)
if [ ! -f "../references/NCBI_blast/16S_ribosomal_RNA.nto" ]; then
    echo "Extracting BLAST database from $BLAST_TAR"
    tar -xzf $BLAST_TAR -C ../references/NCBI_blast/
fi
# Activate the BLAST conda environment
conda activate blast
export BLASTDB="../references/NCBI_blast/" # location of NCBI taxdb database
# make the output directory if it does not exist
mkdir -p $OUT_FOL
if [ ! -f $V1V3_BLAST_OUTPUT ]; then
    echo "Running BLAST for V1-V3 region"
    blastn -query "$V1V3_FASTA" \
        -db "$BLAST_DB" \
        -outfmt "6 qseqid qlen sseqid slen length pident nident evalue bitscore mismatch qcovs gaps sscinames" \
        -perc_identity "99" \
        -evalue 1e-50 \
        -num_threads $NUM_THREADS \
    > "$V1V3_BLAST_OUTPUT"
    echo "BLAST results for V1-V3 region saved to $V1V3_BLAST_OUTPUT"
else
    echo "V1-V3 BLAST results already exist at $V1V3_BLAST_OUTPUT"
fi
if [ ! -f $FL_BLAST_OUTPUT ]; then
    echo "Running BLAST for FL-16S region"
    blastn -query "$FL16S_FASTA" \
        -db "$BLAST_DB" \
        -outfmt "6 qseqid qlen sseqid slen length pident nident evalue bitscore mismatch qcovs gaps sscinames" \
        -perc_identity "99" \
        -evalue 1e-50 \
        -num_threads $NUM_THREADS \
    > "$FL_BLAST_OUTPUT"
    echo "BLAST results for FL-16S region saved to $FL_BLAST_OUTPUT"
else
    echo "FL-16S BLAST results already exist at $FL_BLAST_OUTPUT"
fi
