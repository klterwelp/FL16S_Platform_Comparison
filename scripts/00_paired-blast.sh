#!/bin/bash
# This script performs paired BLAST searches for each amplicon dataset against the other.
# Used to determine which ASVs across V1-V3 and FL-16S are paired,
# where pairs share sequence identity.
# Inputs:
#   - data/FL-16S/dada2_seq.fasta.qza: fasta sequences after DADA2 denoising
#   - data/16S_V1-V3/dada2_seq.fasta.qza: fasta sequences after DADA2 denoising
# Environment:
#   - BLAST 2.16 conda environment from environments/blast.yml
# Outputs:
#   - data/16S_V1-V3/dada2_rep-seqs.fasta: fasta sequences for the V1-V3 region
#   - data/FL-16S/dada2_rep-seqs.fasta: fasta sequences for the FL-16S region
#   - intermediates/paired-blast-results/
#       - V1V3_blast_db/: BLAST database for V1-V3 region built with the V1-V3 DADA2 fasta sequences
#       - FL_blast_db/: BLAST database for the FL-16S region built with the FL-16S DADA2 fasta sequences
#       - V1V3_paired-blast-results.tsv: the paired BLAST hit results for the V1-V3 region against FL-16S
#       - FL16S_paired-blast-results.tsv: the paired BLAST hit results for the FL-16S region against V1-V3
#       - Format of TSV Output Files:
#           - qseqid: query sequence ID
#           - qlen: query sequence length
#           - sseqid: subject sequence ID
#           - slen: subject sequence length
#           - length: alignment length
#           - pident: percentage of identical matches
#           - nident: number of identical matches
#           - evalue: expect value
#           - bitscore: bit score of the alignment
#           - mismatch: number of mismatches
#           - qcovs: query coverage of the subject sequence
#           - gaps: number of gaps in the alignment
# Number of threads
NUM_THREADS=10
# Locations of the input files
V1V3_REP_SEQS="../data/16S_V1-V3/dada2_seq.fasta.qza"
FL16S_REP_SEQS="../data/FL-16S/dada2_seq.fasta.qza"
# Location of the output files
V1V3_FASTA="../data/16S_V1-V3/dada2_rep-seqs.fasta"
FL16S_FASTA="../data/FL-16S/dada2_rep-seqs.fasta"
OUT_FOL="../intermediates/paired-blast-results"
V1V3_BLAST_DB="${OUT_FOL}/V1V3_blast_db"
FL_BLAST_DB="${OUT_FOL}/FL_blast_db"
V1V3_BLAST_OUTPUT="${OUT_FOL}/V1V3_paired-blast-results.tsv"
FL_BLAST_OUTPUT="${OUT_FOL}/FL16S_paired-blast-results.tsv"
# Check that the input files exist
if [ ! -f $V1V3_REP_SEQS ]; then
    echo "Error: $V1V3_REP_SEQS does not exist."
    exit 1
fi
if [ ! -f $FL16S_REP_SEQS ]; then
    echo "Error: $FL16S_REP_SEQS does not exist."
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
# Activate the BLAST conda environment
conda activate blast
# make the output directory if it does not exist
mkdir -p $OUT_FOL
# create the BLAST databases if they do not exist
if [ ! -d $V1V3_BLAST_DB ]; then
    echo "Creating BLAST database for V1-V3 region"
    makeblastdb -dbtype nucl -in "$V1V3_FASTA" -out "$V1V3_BLAST_DB"
fi
if [ ! -d $FL_BLAST_DB ]; then
    echo "Creating BLAST database for FL-16S region"
    makeblastdb -dbtype nucl -in "$FL16S_FASTA" -out "$FL_BLAST_DB"
fi
# run paired BLAST searches
if [ ! -f $V1V3_BLAST_OUTPUT ]; then
    echo "Running BLAST for V1-V3 region against FL-16S data"
    blastn -query "$V1V3_FASTA" \
        -db "$FL_BLAST_DB " \
        -outfmt "6 qseqid qlen sseqid slen length pident nident evalue bitscore mismatch qcovs gaps" \
        -perc_identity "99" \
        -evalue 1e-50 \
        -num_threads $NUM_THREADS \
    > "$V1V3_BLAST_OUTPUT"
    echo "Paired BLAST results for V1-V3 region saved to $V1V3_BLAST_OUTPUT"
else
    echo "V1-V3 paired BLAST results already exist at $V1V3_BLAST_OUTPUT"
fi
if [ ! -f $FL_BLAST_OUTPUT ]; then
    echo "Running BLAST for FL-16S region against V1-V3 data"
    blastn -query "$FL16S_FASTA" \
        -db "$V1V3_BLAST_DB" \
        -outfmt "6 qseqid qlen sseqid slen length pident nident evalue bitscore mismatch qcovs gaps" \
        -perc_identity "99" \
        -evalue 1e-50 \
        -num_threads $NUM_THREADS \
    > "$FL_BLAST_OUTPUT"
    echo "Paired BLAST results for FL-16S region saved to $FL_BLAST_OUTPUT"
else
    echo "FL-16S paired BLAST results already exist at $FL_BLAST_OUTPUT"
fi
