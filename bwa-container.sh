#!/bin/bash
# Script name: bwa-container.sh

# Setting BWA location
export PATH=/bwa:$PATH

echo "Indexing E. coli genome"
bwa index ecoli_rel606.fasta.gz

echo "Starting bwa alignment for SRR2584863"
bwa mem ecoli_rel606.fasta.gz SRR2584863_1.trim.sub.fastq SRR2584863_2.trim.sub.fastq > SRR2584863.aligned.sam

echo "Done with bwa alignment for SRR2584863!"

echo "Cleaning up files generated from genome indexing"
rm ecoli_rel606.fasta.gz.amb
rm ecoli_rel606.fasta.gz.ann
rm ecoli_rel606.fasta.gz.bwt
rm ecoli_rel606.fasta.gz.pac
rm ecoli_rel606.fasta.gz.sa
