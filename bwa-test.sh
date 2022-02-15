#!/bin/bash
# Script name: bwa-test.sh

echo "Unpacking software"

echo "Setting PATH for bwa" 

echo "Indexing E. coli genome"

echo "Starting bwa alignment for SRR2584863"

echo "Done with bwa alignment for SRR2584863!"

echo "Cleaning up files generated from genome indexing"
rm ecoli_rel606.fasta.gz.amb
rm ecoli_rel606.fasta.gz.ann
rm ecoli_rel606.fasta.gz.bwt
rm ecoli_rel606.fasta.gz.pac
rm ecoli_rel606.fasta.gz.sa
